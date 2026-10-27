"""Solve for (n_units, torsion_deg) that make a homopolymer chain approximate a target
shape (straight line, regular polygon, or circle) of a user-given size, instead of the
user picking those two numbers by hand.

The geometric basis: chem.builder.build_polymer places every unit with the *same*
rigid transform (see _place_incoming_unit), so the angle between consecutive
tail-to-tail step vectors (the "turn angle") and the step length (the "chord") are
each roughly constant along the chain and are both pure functions of torsion_deg
alone, for a given monomer/head/tail choice. Setting turn angle = 0 gives a fully
extended, effectively straight chain; setting turn angle = 360/n matches the exterior
angle of a regular n-gon *at every junction*.

IMPORTANT geometric caveat, found by testing against a real monomer rather than
assuming it away: matching the per-junction turn angle to 360/n does NOT generally
make the chain actually close into a flat n-sided polygon. That would additionally
require the chain to have zero "pitch" (no net advance along the rotation axis each
step) at that same torsion value, which is a second, independent condition -- and with
only one free parameter (torsion) to satisfy two conditions, exact flat closure is
only possible for whichever specific angle a given monomer's own fixed bond geometry
happens to allow, not for an arbitrarily chosen shape. Decided with the user: this
module deliberately produces an OPEN, approximate shape -- the right per-junction turn
angle and edge count, but not a guaranteed flat/closed loop -- rather than rejecting
every shape a monomer can't exactly close into. Report/UI text must say "open
approximation", not "closed polygon" or "ring", so this isn't mistaken for an exact
closed macrocycle.

Every solve_* function measures turn_angle(torsion) and chord(torsion) directly with a
3-unit probe build and root-finds on the measurement (never assumes the relationship
analytically), and `_check_constant_step` verifies the constant-repeated-transform
premise itself (comparing the probe's two step lengths) before trusting anything the
scan finds, so an unusual monomer where that premise breaks down fails loudly instead
of silently producing a nonsense shape.

Edge/unit bookkeeping: a repeat unit occupies one *junction*, not one *tail atom* --
n_units repeat units have n_units tail atoms but only (n_units - 1) tail-to-tail steps
(edges). "Requested length" for a line is the distance between the first and last
unit's tail atoms (n_units - 1 edges). A polygon has n_sides - 1 "regular" sides (each
ending in a corner turn into the next) plus one distinguished "closing" side with no
turn after it -- see `_closing_run_length` for why the closing side, not every side,
is where any leftover closure error goes. The last tail atom is simply wherever that
lands, not bonded back to the first (no actual closing bond is ever added) -- but see
`closure_gap_nm` below for how close it actually gets, measured, not assumed.

Scope: only shapes reachable by one constant global turn (line, regular polygon,
circle) are supported. A cube or other true 3D polyhedron needs a different turn at
each junction, which this app's single global torsion cannot produce -- out of scope.

Terminology: a polygon's "side length" is the geometric target the user requests. One
repeat unit's own bond geometry at the required turn angle usually can't reach an
arbitrary requested side length on its own (e.g. it might naturally be 0.25 nm), so a
side is built from a "straight run" of several repeat units (each at the "most
extended" torsion from `_solve_most_extended_torsion`) ending in one "corner" junction
(at the turn-matching torsion) into the next side -- not a separate "edge object",
still one continuous chain. `chem.builder.build_polymer` accepts this as a per-junction
torsion list; see its docstring. `achieved_length_nm` is always a *pre-optimization
construction estimate* from the chosen integer units-per-side, not a promise about the
final optimized structure. The one exception is the *closing* side (see
`_closing_run_length`): its length is chosen by direct measurement to best return the
chain to its starting point, not fit to the requested side length -- `torsion_deg` is
therefore always a per-junction list for polygons (never a bare scalar, even when every
regular side is a single unit), and `closing_side_length_nm`/`closure_gap_nm` report
what that side actually came out to and how close the shape got to closing.

A circle is deliberately NOT built the same way. A polygon's sides are meant to be
long and straight with sharp corners, so straight-run padding is the right tool. A
circle is meant to look round, which needs many small, evenly spread turns -- one
repeat unit per facet, no straight runs at all -- so `solve_circle` solves directly for
how many facets (repeat units) the requested circumference needs, the same way
`solve_line` solves for how many units a requested length needs. See `solve_circle`'s
own docstring for why the straight-run model was tried and rejected for circles.

`flatness_score` (see `_edge_flatness_score`): every `solve_polygon`/`solve_circle`
result reports how consistent the turn axis stays corner-to-corner (1.0 = planar,
lower = corkscrewing in 3D). This is measured, not assumed, and can come out low for a
monomer whose only available torsion roots are inherently non-planar at that turn
angle -- that is a property of the monomer's chemistry, not a bug this module can code
around with a single scalar torsion per junction. Report it to the user rather than
silently building a shape known in advance to drift.

`closure_gap_nm` (polygon only, see `_closing_run_length`): the measured distance
between the chain's start and end after the closing side, relative to the requested
side length, is the real-world consequence of a low `flatness_score` -- a good
flatness score means the closing side's direction happens to point back near the
start, so sizing it by measurement closes the loop well; a poor one means no side
length can fully close it, and `closure_gap_nm` says by how much, instead of leaving
the user to discover it from the rendered shape.
"""

import math
from dataclasses import dataclass

from chem.builder import build_polymer
from chem.validation import MAX_REPEAT_UNITS, validate_pre_build

MIN_POLYGON_EDGES = 3
DEFAULT_MAX_RELATIVE_ERROR = 0.25  # above this, "achieved" no longer counts as satisfying the request
_STEP_CONSISTENCY_TOLERANCE = 0.03  # relative gap between the probe's two step lengths


@dataclass
class ShapeSolveResult:
    ok: bool
    message: str = ""
    n_units: int = None
    torsion_deg: float = None  # float (uniform) or list[float] (per-junction schedule)
    units_per_edge: int = None  # repeat units per polygon side/arc segment (line: not set)
    turn_angle_deg: float = None  # the per-junction turn angle this result targets
    requested_length_nm: float = None  # side length for line/polygon, radius for circle
    achieved_length_nm: float = None  # pre-optimization construction estimate, not a guarantee
    relative_error: float = None
    flatness_score: float = None  # see _edge_flatness_score; None for solve_line (always flat)
    closure_gap_nm: float = None  # polygon only; measured gap between chain start and end after the closing side
    closing_side_length_nm: float = None  # polygon only; the closing side's actual length, may differ from the rest
    corner_torsion_deg: float = None  # polygon only; the schedule's last entry is the closing side's straight
    # torsion, not the corner -- callers that need "the corner angle" (e.g. a UI slider
    # reference value) should read this field, not torsion_deg[-1].


def probe_step_geometry(monomer, head_idx, tail_idx, torsion_deg: float) -> tuple[float, float, float]:
    """Build a 3-unit probe chain and measure (turn_angle_deg, chord_angstrom,
    step_consistency) from the two tail-to-tail step vectors. `step_consistency` is
    the relative difference between the two steps' lengths -- near 0 confirms
    build_polymer really is applying one constant repeated transform for this
    monomer/head/tail choice, which is what makes measuring 2 steps enough to
    characterize the whole chain."""
    result = build_polymer(monomer, head_idx, tail_idx, n_units=3, torsion_deg=torsion_deg)
    conf = result.mol.GetConformer()
    p0, p1, p2 = (conf.GetAtomPosition(i) for i in result.unit_tail_indices)
    v1 = [p1.x - p0.x, p1.y - p0.y, p1.z - p0.z]
    v2 = [p2.x - p1.x, p2.y - p1.y, p2.z - p1.z]
    n1 = math.sqrt(sum(c * c for c in v1))
    n2 = math.sqrt(sum(c * c for c in v2))
    dot = sum(a * b for a, b in zip(v1, v2))
    cos_turn = max(-1.0, min(1.0, dot / (n1 * n2)))
    turn_angle_deg = math.degrees(math.acos(cos_turn))
    chord_angstrom = (n1 + n2) / 2.0
    step_consistency = abs(n1 - n2) / max(n1, n2, 1e-9)
    return turn_angle_deg, chord_angstrom, step_consistency


def _check_constant_step(monomer, head_idx, tail_idx):
    """One-off sanity probe (at an arbitrary torsion, 180 deg): the constant-repeated-
    transform premise this whole module relies on is a structural property of the
    monomer/head/tail choice, not something that varies with torsion, so it only
    needs checking once per solve_* call rather than on every scan sample. Returns an
    error message if the premise doesn't hold, else None."""
    _, _, step_consistency = probe_step_geometry(monomer, head_idx, tail_idx, 180.0)
    if step_consistency > _STEP_CONSISTENCY_TOLERANCE:
        return (
            f"This monomer/head/tail choice does not produce a consistent repeat "
            f"step (consecutive step lengths differ by {step_consistency:.1%}); shape "
            f"presets assume a uniformly repeated geometry and aren't reliable here."
        )
    return None


def _scan_turn_angle(monomer, head_idx, tail_idx, n_samples=36):
    """turn_angle(phi) and chord(phi) depend only on the monomer/head/tail choice, not
    on which shape is being solved for -- scanned once per solve_* call and reused for
    every candidate target turn angle (e.g. every candidate circle resolution) via
    `_find_roots_for_turn_angle`, instead of rescanning from scratch each time."""
    phis = [-180.0 + 360.0 * i / n_samples for i in range(n_samples + 1)]
    turns, chords = [], []
    for phi in phis:
        t, c, _ = probe_step_geometry(monomer, head_idx, tail_idx, phi)
        turns.append(t)
        chords.append(c)
    return phis, turns, chords


def _bisect_for_turn_angle(monomer, head_idx, tail_idx, target_turn_deg, lo, hi, iters=25):
    f_lo = probe_step_geometry(monomer, head_idx, tail_idx, lo)[0] - target_turn_deg
    for _ in range(iters):
        mid = (lo + hi) / 2.0
        f_mid = probe_step_geometry(monomer, head_idx, tail_idx, mid)[0] - target_turn_deg
        if (f_mid < 0) == (f_lo < 0):
            lo, f_lo = mid, f_mid
        else:
            hi = mid
    mid = (lo + hi) / 2.0
    _, chord, _ = probe_step_geometry(monomer, head_idx, tail_idx, mid)
    return mid, chord


def _find_roots_for_turn_angle(monomer, head_idx, tail_idx, target_turn_deg, scan):
    """Every torsion in `scan` where turn_angle(phi) == target_turn_deg, refined via
    bisection. There can be more than one (e.g. mirror-image torsions with different
    chord lengths) -- callers pick among them by whichever chord best matches size."""
    phis, turns, chords = scan
    roots = []
    for i in range(len(phis) - 1):
        d0, d1 = turns[i] - target_turn_deg, turns[i + 1] - target_turn_deg
        if d0 == 0:
            roots.append((phis[i], chords[i]))
        elif d0 * d1 < 0:
            roots.append(_bisect_for_turn_angle(monomer, head_idx, tail_idx, target_turn_deg, phis[i], phis[i + 1]))
    if not roots and turns[-1] - target_turn_deg == 0:
        roots.append((phis[-1], chords[-1]))
    return roots


def solve_torsion_for_turn_angle(monomer, head_idx, tail_idx, target_turn_deg, n_samples=36):
    """Single-target convenience wrapper (used by solve_polygon/solve_circle): scans
    once, finds roots for `target_turn_deg`, and returns (roots, achievable_range)."""
    scan = _scan_turn_angle(monomer, head_idx, tail_idx, n_samples)
    roots = _find_roots_for_turn_angle(monomer, head_idx, tail_idx, target_turn_deg, scan)
    return roots, (min(scan[1]), max(scan[1]))


def _end_to_end_per_unit(monomer, head_idx, tail_idx, torsion_deg, n_probe=9):
    result = build_polymer(monomer, head_idx, tail_idx, n_units=n_probe, torsion_deg=torsion_deg)
    conf = result.mol.GetConformer()
    p0 = conf.GetAtomPosition(result.unit_tail_indices[0])
    p_last = conf.GetAtomPosition(result.unit_tail_indices[-1])
    return p0.Distance(p_last) / (n_probe - 1)


def _solve_most_extended_torsion(monomer, head_idx, tail_idx, n_probe=9, coarse_samples=36):
    """"Straight" for a homopolymer backbone does not mean each single repeat unit
    points the same direction as the last -- a classic all-trans alkane zigzags
    atom-by-atom around a fixed (~109 deg) bond angle and is still the maximally
    extended, "straight fiber axis" conformation. So instead of targeting a
    per-junction turn angle of zero (which, per the module docstring's caveat,
    conflates "genuinely bent" with "regularly zigzagging," and for many monomers
    isn't reachable at any torsion at all), this directly maximizes end-to-end
    distance over a several-unit probe -- the standard definition of "extended" in
    polymer chemistry -- via a coarse scan refined by a finer local scan."""
    phis = [-180.0 + 360.0 * i / coarse_samples for i in range(coarse_samples)]
    best_phi = max(phis, key=lambda phi: _end_to_end_per_unit(monomer, head_idx, tail_idx, phi, n_probe))
    step = 360.0 / coarse_samples
    fine_phis = [best_phi - step + 2 * step * i / 20 for i in range(21)]
    best_phi = max(fine_phis, key=lambda phi: _end_to_end_per_unit(monomer, head_idx, tail_idx, phi, n_probe))
    chord = _end_to_end_per_unit(monomer, head_idx, tail_idx, best_phi, n_probe)
    return best_phi, chord


def _edge_flatness_score(monomer, head_idx, tail_idx, straight_torsion, corner_torsion, units_per_edge):
    """1.0 = consecutive side-to-side turns share a consistent rotation axis (the
    shape stays close to planar); lower values mean each side's turn twists off in a
    new direction (a corkscrewing, non-planar result -- see the module docstring's
    flatness caveat). There is usually more than one torsion value that reaches the
    same target turn angle (mirror-image roots); they are not equally flat, so when
    more than one gives an acceptable size, callers prefer the flattest.

    Measured directly on a short 3-side probe built with the actual candidate
    torsion values, not estimated from corner_torsion alone: whether a corner root
    precesses depends on how its rotation composes with `units_per_edge` straight-run
    steps, which is specific to this (straight_torsion, corner_torsion, units_per_edge)
    combination.
    """
    edge = [straight_torsion] * (units_per_edge - 1) + [corner_torsion]
    n_units = 3 * units_per_edge + 1
    result = build_polymer(monomer, head_idx, tail_idx, n_units=n_units, torsion_deg=edge * 3)
    conf = result.mol.GetConformer()
    vertices = [conf.GetAtomPosition(i) for i in result.unit_tail_indices[::units_per_edge]]
    if len(vertices) < 4:
        return 1.0

    side_vecs = [
        (vertices[i + 1].x - vertices[i].x, vertices[i + 1].y - vertices[i].y, vertices[i + 1].z - vertices[i].z)
        for i in range(len(vertices) - 1)
    ]
    axes = []
    for a, b in zip(side_vecs, side_vecs[1:]):
        cx, cy, cz = a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]
        norm = math.sqrt(cx * cx + cy * cy + cz * cz)
        if norm < 1e-6:
            return 1.0
        axes.append((cx / norm, cy / norm, cz / norm))
    if len(axes) < 2:
        return 1.0

    dots = [abs(sum(a * b for a, b in zip(x, y))) for x, y in zip(axes, axes[1:])]
    return min(dots)


def _pick_flattest_viable(monomer, head_idx, tail_idx, straight_torsion, candidates, is_viable):
    """`candidates` is a list of (corner_torsion, units_per_edge, ..., error) tuples
    (error last). Among whichever are viable, pick the flattest (see
    _edge_flatness_score); if none are viable, just pick the smallest error so the
    caller has something concrete to report in its rejection message. Always measures
    and returns the winning flatness score too (not just when there's a real choice
    to make) so callers can report it to the user honestly, even when only one root
    exists.

    Returns (chosen_candidate, flatness_score).
    """
    viable = [c for c in candidates if is_viable(c)]
    pool = viable if viable else candidates
    scored = [(c, _edge_flatness_score(monomer, head_idx, tail_idx, straight_torsion, c[0], c[1])) for c in pool]
    return max(scored, key=lambda item: item[1])


def _relative_error(requested: float, achieved: float) -> float:
    if requested <= 0:
        return 0.0
    return abs(achieved - requested) / requested


def _choose_units_per_edge(corner_chord, straight_chord, target_edge_length, max_units_per_edge):
    """How many repeat units should make up one side: 1 corner-turn unit plus however
    many "most extended" straight-run units close the gap to `target_edge_length`.
    `units_per_edge` is necessarily an integer, so this evaluates floor/round/ceil
    candidates (clipped to [1, max_units_per_edge]) rather than assuming round() is
    always the closest-error choice, and returns whichever minimizes
    |edge_length - target_edge_length|. Returns None if even a single corner-only unit
    doesn't fit the caller's overall unit budget.

    Returns (units_per_edge, edge_length_angstrom) or None.
    """
    if max_units_per_edge < 1:
        return None
    if straight_chord < 1e-6:
        return 1, corner_chord
    raw = 1 + (target_edge_length - corner_chord) / straight_chord
    candidates = {max(1, min(max_units_per_edge, int(f(raw)))) for f in (math.floor, round, math.ceil)}
    best_k = min(candidates, key=lambda k: abs(corner_chord + (k - 1) * straight_chord - target_edge_length))
    return best_k, corner_chord + (best_k - 1) * straight_chord


def _vsub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vdot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _vnorm(a):
    return math.sqrt(_vdot(a, a))


def _vnormalize(a):
    n = _vnorm(a)
    return (a[0] / n, a[1] / n, a[2] / n) if n > 1e-9 else (0.0, 0.0, 0.0)


def _build_positions(monomer, head_idx, tail_idx, torsion_schedule):
    """Every repeat unit's tail-atom position, actually built and measured -- not
    estimated from chord/turn formulas -- for a given per-junction torsion schedule."""
    n_units = len(torsion_schedule) + 1
    result = build_polymer(monomer, head_idx, tail_idx, n_units=n_units, torsion_deg=torsion_schedule)
    conf = result.mol.GetConformer()
    return [tuple(conf.GetAtomPosition(i)) for i in result.unit_tail_indices]


def _closing_run_length(monomer, head_idx, tail_idx, straight_torsion, corner_torsion, units_per_edge, n_sides, max_extra_units, heading_probe_units=5):
    """How many straight-run units the final side needs to best close the loop back to
    the start -- found by direct measurement, not by assuming the final side should be
    the same length as the others. Matching every corner's turn *magnitude* (see
    solve_torsion_for_turn_angle) does not guarantee the turns share a consistent
    rotation axis (see _edge_flatness_score) -- so after n_sides - 1 real turns, the
    direction the chain ends up heading is an empirical fact that has to be measured,
    not a formula. Given that heading, length is the only free parameter left, and it
    can only bring the endpoint as close to the start as that heading allows -- when
    flatness is good, the heading points back near the start and this closes the loop
    almost exactly (possibly at a slightly different side length than the rest); when
    flatness is poor, this still finds the best achievable approach, honestly reported
    via the returned gap rather than assumed away.

    The heading itself is measured over several units (`heading_probe_units`), not
    just the single last step -- a straight run zigzags atom-by-atom (see module
    docstring), so a single-step heading estimate is noisy and can point noticeably
    off the run's real net direction. That estimate is then only a starting point: a
    direct hill-climb over measured gaps (not the linear estimate) picks the actual
    minimum, since the estimate assumes perfectly straight travel that the real
    zigzagging chain doesn't exactly deliver.

    Returns (closing_units, achieved_gap_angstrom, closing_side_length_angstrom).
    `closing_units` is clipped to >= 1 (a real closing side must exist) and
    <= max_extra_units (the caller's remaining unit budget).
    """
    max_extra_units = max(1, max_extra_units)
    base_schedule = ([straight_torsion] * (units_per_edge - 1) + [corner_torsion]) * (n_sides - 1)
    positions = _build_positions(monomer, head_idx, tail_idx, base_schedule)
    origin, p_last = positions[0], positions[-1]

    probe_extra = min(heading_probe_units, max_extra_units)
    probe_positions = _build_positions(monomer, head_idx, tail_idx, base_schedule + [straight_torsion] * probe_extra)
    heading_vec = _vsub(probe_positions[-1], p_last)
    step_length = _vnorm(heading_vec) / probe_extra
    heading = _vnormalize(heading_vec)

    def measure(k):
        positions_k = _build_positions(monomer, head_idx, tail_idx, base_schedule + [straight_torsion] * k)
        gap = _vnorm(_vsub(positions_k[-1], origin))
        closing_length = _vnorm(_vsub(positions_k[-1], p_last))
        return gap, closing_length

    if step_length < 1e-6:
        k0 = 1
    else:
        projected_steps = _vdot(_vsub(origin, p_last), heading) / step_length
        k0 = max(1, min(max_extra_units, round(projected_steps)))

    best_k = k0
    best_gap, best_length = measure(k0)
    for direction in (1, -1):
        k = best_k
        while True:
            k_next = k + direction
            if k_next < 1 or k_next > max_extra_units:
                break
            gap_next, length_next = measure(k_next)
            if gap_next < best_gap - 1e-9:
                best_k, best_gap, best_length, k = k_next, gap_next, length_next, k_next
            else:
                break

    return best_k, best_gap, best_length


def solve_line(
    monomer,
    head_idx,
    tail_idx,
    length_nm: float,
    max_units: int = MAX_REPEAT_UNITS,
    max_relative_error: float = DEFAULT_MAX_RELATIVE_ERROR,
) -> ShapeSolveResult:
    if length_nm <= 0:
        return ShapeSolveResult(False, message="Length must be positive.")

    bad = _check_constant_step(monomer, head_idx, tail_idx)
    if bad:
        return ShapeSolveResult(False, message=bad)

    torsion_deg, chord = _solve_most_extended_torsion(monomer, head_idx, tail_idx)
    if chord < 1e-6:
        return ShapeSolveResult(
            False, message="This monomer has no meaningfully extended conformation to build a straight chain from."
        )

    # `length_nm` is the end-to-end distance between the first and last unit's tail
    # atoms: n_units tail atoms make (n_units - 1) steps, not n_units.
    n_steps = max(1, round(length_nm * 10.0 / chord))
    n_units = min(n_steps + 1, max_units)
    achieved_length_nm = (n_units - 1) * chord / 10.0

    rel_err = _relative_error(length_nm, achieved_length_nm)
    if rel_err > max_relative_error:
        return ShapeSolveResult(
            False,
            message=(
                f"Best achievable straight-chain length within the {max_units}-unit limit is "
                f"{achieved_length_nm:.2f} nm, too far from the requested {length_nm:.2f} nm."
            ),
        )

    pre = validate_pre_build(monomer, head_idx, tail_idx, n_units=n_units, max_units=max_units)
    if not pre.ok:
        return ShapeSolveResult(False, message=pre.message)

    return ShapeSolveResult(
        True,
        n_units=n_units,
        torsion_deg=torsion_deg,
        requested_length_nm=length_nm,
        achieved_length_nm=achieved_length_nm,
        relative_error=rel_err,
    )


def solve_polygon(
    monomer,
    head_idx,
    tail_idx,
    n_sides: int,
    side_length_nm: float,
    max_units: int = MAX_REPEAT_UNITS,
    max_relative_error: float = DEFAULT_MAX_RELATIVE_ERROR,
) -> ShapeSolveResult:
    if side_length_nm <= 0:
        return ShapeSolveResult(False, message="Side length must be positive.")
    if n_sides + 1 > max_units:
        return ShapeSolveResult(
            False,
            message=f"An open {n_sides}-sided-corner shape needs at least {n_sides + 1} repeat units, above the {max_units}-unit limit.",
        )

    bad = _check_constant_step(monomer, head_idx, tail_idx)
    if bad:
        return ShapeSolveResult(False, message=bad)

    target_turn = 360.0 / n_sides
    roots, (lo, hi) = solve_torsion_for_turn_angle(monomer, head_idx, tail_idx, target_turn)
    if not roots:
        return ShapeSolveResult(
            False,
            message=(
                f"This monomer's local geometry can only produce turn angles between "
                f"{lo:.1f}° and {hi:.1f}° at this head/tail choice; a regular {n_sides}-sided "
                f"polygon needs {target_turn:.1f}° per junction, which isn't achievable."
            ),
        )

    straight_torsion, straight_chord = _solve_most_extended_torsion(monomer, head_idx, tail_idx)
    target_edge_length = side_length_nm * 10.0
    max_units_per_edge = (max_units - 1) // n_sides

    candidates = []  # (corner_torsion, units_per_edge, edge_length, relative_error)
    for corner_torsion, corner_chord in roots:
        choice = _choose_units_per_edge(corner_chord, straight_chord, target_edge_length, max_units_per_edge)
        if choice is None:
            continue
        units_per_edge, edge_length = choice
        err = _relative_error(side_length_nm, edge_length / 10.0)
        candidates.append((corner_torsion, units_per_edge, edge_length, err))

    if not candidates:
        return ShapeSolveResult(
            False,
            message=(
                f"A regular {n_sides}-sided shape needs at least {n_sides * max_units_per_edge + 1} "
                f"repeat units to reach any side length within the {max_units}-unit limit."
            ),
        )

    # Multiple torsion values can hit the same target turn angle (mirror-image
    # roots) but aren't equally flat (see _edge_flatness_score) -- among whichever
    # give an acceptable size, prefer the flattest rather than just the best size fit.
    (corner_torsion, units_per_edge, edge_length, rel_err), flatness_score = _pick_flattest_viable(
        monomer, head_idx, tail_idx, straight_torsion, candidates, lambda c: c[3] <= max_relative_error
    )
    achieved_length_nm = edge_length / 10.0

    if rel_err > max_relative_error:
        return ShapeSolveResult(
            False,
            message=(
                f"Best achievable side length within the {max_units}-unit limit is "
                f"{achieved_length_nm:.2f} nm, too far from the requested {side_length_nm:.2f} nm."
            ),
        )

    # Matching every corner's turn magnitude (above) does not guarantee the chain
    # actually comes back around to its starting point -- that depends on whether the
    # turns share a consistent axis (see _edge_flatness_score), which is a fact about
    # this monomer's real 3D geometry, not something the target turn angle alone
    # determines. So build every side but the last with the schedule above, MEASURE
    # where that leaves the chain, and size the final side by direct measurement to
    # best close the gap back to the start -- rather than assuming it should be the
    # same nominal length as the rest.
    base_units_used = (n_sides - 1) * units_per_edge
    max_extra_units = max_units - 1 - base_units_used
    if max_extra_units < 1:
        return ShapeSolveResult(
            False,
            message=f"Closing the {n_sides}-sided shape needs at least one more repeat unit beyond the {max_units}-unit limit.",
        )

    closing_units, closure_gap_ang, closing_side_ang = _closing_run_length(
        monomer, head_idx, tail_idx, straight_torsion, corner_torsion, units_per_edge, n_sides, max_extra_units
    )
    n_units = base_units_used + closing_units + 1

    pre = validate_pre_build(monomer, head_idx, tail_idx, n_units=n_units, max_units=max_units)
    if not pre.ok:
        return ShapeSolveResult(False, message=pre.message)

    edge_pattern = [straight_torsion] * (units_per_edge - 1) + [corner_torsion]
    torsion_result = edge_pattern * (n_sides - 1) + [straight_torsion] * closing_units

    return ShapeSolveResult(
        True,
        n_units=n_units,
        torsion_deg=torsion_result,
        units_per_edge=units_per_edge,
        turn_angle_deg=target_turn,
        requested_length_nm=side_length_nm,
        achieved_length_nm=achieved_length_nm,
        relative_error=rel_err,
        flatness_score=flatness_score,
        closure_gap_nm=closure_gap_ang / 10.0,
        closing_side_length_nm=closing_side_ang / 10.0,
        corner_torsion_deg=corner_torsion,
    )


def solve_circle(
    monomer,
    head_idx,
    tail_idx,
    radius_nm: float,
    max_units: int = MAX_REPEAT_UNITS,
    min_edges: int = MIN_POLYGON_EDGES,
    max_relative_error: float = DEFAULT_MAX_RELATIVE_ERROR,
) -> ShapeSolveResult:
    """Approximates the circle as a smooth many-faceted regular polygon: ONE repeat
    unit per facet, every junction turning by the same small angle (360 / n_facets)
    -- deliberately NOT solve_polygon's straight-run-plus-corner model. A real circle
    should look round (many small, evenly spread turns); this module used to let a
    straight run pad every facet out to whatever length hit the requested radius,
    but that meant almost *any* n_facets could match the radius equally well (the
    straight run just absorbs the difference), so the search had no signal pushing
    it toward "many small facets" over "a handful of long straight stretches with
    sharp corners" -- confirmed by a live report where a requested circle rendered
    as a visible pentagon/hexagon, not a circle.

    Here radius is NOT independently tunable from n_facets: for a given monomer,
    each n_facets has exactly one turn angle (360/n_facets) and hence one fixed
    natural chord (from the turn-angle root-find) -- no free length parameter -- so
    the search is over which integer n_facets' resulting achieved_radius lands
    closest to the request, the same way solve_line searches over n_steps. Same
    open-approximation caveat as solve_polygon (see module docstring): matching the
    turn angle at every junction does not guarantee a flat, closed ring -- check
    the returned flatness_score."""
    if radius_nm <= 0:
        return ShapeSolveResult(False, message="Radius must be positive.")

    bad = _check_constant_step(monomer, head_idx, tail_idx)
    if bad:
        return ShapeSolveResult(False, message=bad)

    max_facets = max_units - 1
    if max_facets < min_edges:
        return ShapeSolveResult(
            False, message=f"The {max_units}-unit limit doesn't allow even a {min_edges}-facet ring."
        )

    target_radius = radius_nm * 10.0
    scan = _scan_turn_angle(monomer, head_idx, tail_idx)

    def evaluate(n_facets):
        target_turn = 360.0 / n_facets
        roots = _find_roots_for_turn_angle(monomer, head_idx, tail_idx, target_turn, scan)
        if not roots:
            return None
        half_angle = math.radians(target_turn) / 2.0

        candidates = []  # (corner_torsion, units_per_edge=1, achieved_radius, radius_error)
        for corner_torsion, corner_chord in roots:
            achieved_radius = corner_chord / (2 * math.sin(half_angle))
            candidates.append((corner_torsion, 1, achieved_radius, abs(achieved_radius - target_radius)))

        # No straight-run padding to fall back on here (unlike solve_polygon) --
        # among the (usually 1-2) mirror-image roots, prefer the flatter one, same
        # tie-break as solve_polygon.
        return _pick_flattest_viable(
            monomer, head_idx, tail_idx, 0.0, candidates,
            lambda c: c[3] <= max_relative_error * target_radius,
        )

    evaluated = {}
    coarse = list(range(min_edges, max_facets + 1, 5))
    if max_facets not in coarse:
        coarse.append(max_facets)
    for n in coarse:
        result = evaluate(n)
        if result is not None:
            evaluated[n] = result

    if not evaluated:
        return ShapeSolveResult(False, message="No repeat count in the allowed range produces a usable ring approximation for this monomer.")

    best_n = min(evaluated, key=lambda n: abs(evaluated[n][0][2] - target_radius))
    for n in range(max(min_edges, best_n - 5), min(max_facets, best_n + 5) + 1):
        if n not in evaluated:
            result = evaluate(n)
            if result is not None:
                evaluated[n] = result

    best_n = min(evaluated, key=lambda n: abs(evaluated[n][0][2] - target_radius))
    (corner_torsion, _, achieved_radius_ang, _), flatness_score = evaluated[best_n]
    achieved_length_nm = achieved_radius_ang / 10.0
    n_units = best_n + 1

    rel_err = _relative_error(radius_nm, achieved_length_nm)
    if rel_err > max_relative_error:
        return ShapeSolveResult(
            False,
            message=(
                f"Closest achievable ring is a {best_n}-facet approximation with radius "
                f"{achieved_length_nm:.2f} nm, too far from the requested {radius_nm:.2f} nm."
            ),
        )

    pre = validate_pre_build(monomer, head_idx, tail_idx, n_units=n_units, max_units=max_units)
    if not pre.ok:
        return ShapeSolveResult(False, message=pre.message)

    return ShapeSolveResult(
        True,
        n_units=n_units,
        torsion_deg=corner_torsion,
        units_per_edge=1,
        turn_angle_deg=360.0 / best_n,
        requested_length_nm=radius_nm,
        achieved_length_nm=achieved_length_nm,
        relative_error=rel_err,
        flatness_score=flatness_score,
    )
