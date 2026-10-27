import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from chem.builder import build_polymer
from chem.shape_solver import (
    _edge_flatness_score,
    _solve_most_extended_torsion,
    probe_step_geometry,
    solve_circle,
    solve_line,
    solve_polygon,
    solve_torsion_for_turn_angle,
)


def _ethane():
    """Pure geometry/algorithm test case, not a claim about real polymer chemistry --
    same rationale as tests/test_builder.py's _ethane() helper."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _hexane():
    """A longer, more flexible test case than _ethane(): ethane's step-vector turn
    angle can never go below ~54.6 deg (confirmed by direct scan -- see
    chem/shape_solver.py's docstring), so it structurally cannot form a many-facet
    smooth circle no matter the algorithm. Hexane's extra backbone length reaches
    turn angles near 0, which is what a large-radius circle actually needs."""
    mol = Chem.MolFromSmiles("CCCCCC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def _tail_positions(build_result):
    conf = build_result.mol.GetConformer()
    return np.array([list(conf.GetAtomPosition(i)) for i in build_result.unit_tail_indices])


def test_probe_step_geometry_is_self_consistent_for_ethane():
    monomer = _ethane()
    turn, chord, step_consistency = probe_step_geometry(monomer, head_idx=0, tail_idx=1, torsion_deg=180.0)
    assert step_consistency < 0.01, "the two probe steps should match closely for a genuinely repeated transform"
    assert chord > 0.5


def test_solve_line_hits_requested_length_within_one_chord():
    monomer = _ethane()
    result = solve_line(monomer, head_idx=0, tail_idx=1, length_nm=1.0)
    assert result.ok, result.message
    assert result.n_units >= 2
    chord_nm = result.achieved_length_nm / (result.n_units - 1)
    assert abs(result.achieved_length_nm - 1.0) <= chord_nm + 1e-6


def test_solve_line_result_actually_builds_an_extended_chain():
    """"Straight" doesn't mean every atom is literally collinear -- a classic
    all-trans alkane zigzags atom-by-atom around its fixed bond angle and is still
    the maximally extended conformation (see chem/shape_solver.py's docstring). The
    real invariant to check is that end-to-end distance accounts for nearly all of
    the contour length, i.e. the chain doesn't fold back on itself."""
    monomer = _ethane()
    result = solve_line(monomer, head_idx=0, tail_idx=1, length_nm=1.5)
    assert result.ok, result.message

    build = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=result.n_units, torsion_deg=result.torsion_deg)
    positions = _tail_positions(build)

    end_to_end = np.linalg.norm(positions[-1] - positions[0])
    contour_length = np.linalg.norm(np.diff(positions, axis=0), axis=1).sum()
    # A tetrahedral (~109.5 deg) all-trans zigzag's theoretical end-to-end/contour
    # ratio is ~0.82-0.87, not 1.0 -- collinearity was never the right bar (see the
    # docstring above); this threshold rules out a coiled/folded conformation, not
    # literal straightness.
    assert end_to_end / contour_length > 0.8


def test_most_extended_torsion_is_near_maximal_over_a_short_run_too():
    """_solve_most_extended_torsion fits its torsion/chord over a 9-unit probe (see
    its docstring); confirm the same torsion is also close to maximally extended over
    a much shorter 2-unit run, i.e. it's genuinely a near-locally-maximal torsion, not
    an artifact of the specific probe length it happened to be fit on. This matters
    because solve_polygon/solve_circle reuse it to build multi-unit straight runs of
    varying length."""
    monomer = _ethane()
    torsion_deg, chord = _solve_most_extended_torsion(monomer, head_idx=0, tail_idx=1)
    build = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=2, torsion_deg=torsion_deg)
    positions = _tail_positions(build)
    two_unit_extension = np.linalg.norm(positions[-1] - positions[0])
    assert two_unit_extension > 0.9 * chord, (
        f"expected the straight-run torsion to also be near-maximal over 2 units, "
        f"got {two_unit_extension:.3f} vs the 9-unit-probe chord {chord:.3f}"
    )


def _bend_angles_deg(positions):
    steps = np.diff(positions, axis=0)
    norms = np.linalg.norm(steps, axis=1)
    cos_bend = np.clip(np.sum(steps[:-1] * steps[1:], axis=1) / (norms[:-1] * norms[1:]), -1.0, 1.0)
    return np.degrees(np.arccos(cos_bend))


def _regular_side_vertices(positions, units_per_edge, n_regular_sides):
    """The n_regular_sides + 1 vertices of solve_polygon's n_sides - 1 "regular"
    (equal-length) sides -- NOT including the final closing side, whose own unit
    count can differ from units_per_edge (see chem/shape_solver.py's
    _closing_run_length), so it can't be swept up by a uniform stride like the rest."""
    return positions[: n_regular_sides * units_per_edge + 1 : units_per_edge]


def _corner_bend_angles_deg(positions, units_per_edge, n_sides):
    """The turn angle actually solved for lives at each *corner junction* -- NOT at
    the angle between aggregated side-to-side vectors, which is a different quantity
    once a straight run wobbles (a zigzagging run's average direction isn't exactly
    its final step's direction, so collapsing to side vertices and measuring the
    angle between those undercounts the real per-corner turn; confirmed empirically
    while writing this test). `bends[i]` (from `_bend_angles_deg`) reflects the
    junction that placed step i+1, i.e. torsion-schedule index i+1; corners sit at
    schedule indices `units_per_edge - 1, 2*units_per_edge - 1, ...`, so their bends
    are one index lower (skipping the very first corner, which has no preceding step
    to measure a bend against)."""
    bends = _bend_angles_deg(positions)
    corner_schedule_indices = [k * units_per_edge - 1 for k in range(1, n_sides + 1)]
    bend_indices = [c - 1 for c in corner_schedule_indices if c - 1 >= 0]
    return bends[bend_indices]


def _regular_polygon_check(monomer, n_sides, side_length_nm):
    """Verifies what solve_polygon actually guarantees: the right side count, equal
    length among the n_sides - 1 "regular" sides, the correct constant per-corner turn
    angle on those sides, and a closing side whose reported length/gap match direct
    measurement of the actual built structure. The closing side is deliberately NOT
    required to match the other sides' length -- see chem/shape_solver.py's
    _closing_run_length docstring for why forcing that would be assuming away exactly
    the non-planarity this module measures instead."""
    result = solve_polygon(monomer, head_idx=0, tail_idx=1, n_sides=n_sides, side_length_nm=side_length_nm)
    assert result.ok, result.message
    assert isinstance(result.torsion_deg, list)
    assert len(result.torsion_deg) == result.n_units - 1

    build = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=result.n_units, torsion_deg=result.torsion_deg)
    positions = _tail_positions(build)
    assert len(positions) == result.n_units

    n_regular_sides = n_sides - 1
    vertices = _regular_side_vertices(positions, result.units_per_edge, n_regular_sides)
    assert len(vertices) == n_regular_sides + 1

    side_lengths = np.linalg.norm(np.diff(vertices, axis=0), axis=1)
    assert side_lengths.std() < 0.05 * side_lengths.mean(), "the n_sides - 1 regular sides should be nearly equal length"

    target_turn = 360.0 / n_sides
    corner_bends = _corner_bend_angles_deg(positions, result.units_per_edge, n_regular_sides)
    assert np.allclose(corner_bends, target_turn, atol=2.0), (
        f"expected every regular corner to turn by {target_turn:.1f} deg, got {corner_bends}"
    )

    # The closing side's reported length/gap must match what's actually in the built
    # conformer, not just the internal estimate used to choose it.
    measured_closing_length = np.linalg.norm(positions[-1] - vertices[-1])
    measured_closure_gap = np.linalg.norm(positions[-1] - positions[0])
    assert abs(measured_closing_length / 10.0 - result.closing_side_length_nm) < 1e-3
    assert abs(measured_closure_gap / 10.0 - result.closure_gap_nm) < 1e-3

    return result, positions


def test_solve_square_has_correct_edge_count_and_turn_angle():
    monomer = _ethane()
    # Probe the side length this monomer naturally achieves for a 90 deg turn rather
    # than guessing an achievable size in advance -- the test only needs to confirm
    # the *turn/edge geometry* is correct, not fight the solver's size-error rejection.
    roots, (lo, hi) = solve_torsion_for_turn_angle(monomer, head_idx=0, tail_idx=1, target_turn_deg=90.0)
    assert roots, f"ethane can't reach a 90 degree turn (range [{lo:.1f}, {hi:.1f}]); pick a different target_turn_deg"
    natural_side_nm = max(r[1] for r in roots) / 10.0
    _regular_polygon_check(monomer, n_sides=4, side_length_nm=natural_side_nm)


def test_solve_polygon_prefers_the_flattest_of_equally_sized_roots():
    """A target turn angle is usually reachable via more than one torsion (mirror-
    image roots with the same chord, hence the same size fit) -- they are not equally
    flat (see _edge_flatness_score's module docstring), so solve_polygon should pick
    whichever one keeps consecutive side-to-side turns closest to a consistent axis,
    not just the first/arbitrary one."""
    monomer = _ethane()
    straight_torsion, _ = _solve_most_extended_torsion(monomer, head_idx=0, tail_idx=1)
    roots, _ = solve_torsion_for_turn_angle(monomer, head_idx=0, tail_idx=1, target_turn_deg=90.0)
    assert len(roots) >= 2, "this test needs at least two equally-sized roots to be meaningful"

    result = solve_polygon(monomer, head_idx=0, tail_idx=1, n_sides=4, side_length_nm=roots[0][1] / 10.0)
    assert result.ok, result.message

    scores = [
        _edge_flatness_score(monomer, 0, 1, straight_torsion, corner_torsion, result.units_per_edge)
        for corner_torsion, _ in roots
    ]
    chosen_score = _edge_flatness_score(monomer, 0, 1, straight_torsion, result.corner_torsion_deg, result.units_per_edge)
    assert chosen_score >= max(scores) - 1e-6, (
        f"solve_polygon picked a root scoring {chosen_score:.3f}, but {max(scores):.3f} was available"
    )


def test_solve_hexagon_has_correct_edge_count_and_turn_angle():
    monomer = _ethane()
    roots, (lo, hi) = solve_torsion_for_turn_angle(monomer, head_idx=0, tail_idx=1, target_turn_deg=60.0)
    assert roots, f"ethane can't reach a 60 degree turn (range [{lo:.1f}, {hi:.1f}]); pick a different target_turn_deg"
    natural_side_nm = max(r[1] for r in roots) / 10.0
    _regular_polygon_check(monomer, n_sides=6, side_length_nm=natural_side_nm)


def test_solve_triangle_closing_side_actually_closes_the_loop():
    """Regression test for the live bug report: a 5 nm ethane triangle's first two
    sides connected fine but the third never came back around to meet the first --
    matching every corner's turn *magnitude* doesn't guarantee the turns share a
    consistent axis (see _edge_flatness_score), so after 2 real turns the chain's
    heading isn't guaranteed to point back at the start. The closing side is now
    sized by direct measurement (see _closing_run_length) to close that gap as well
    as this monomer's geometry allows, rather than assuming it should equal the
    nominal side length."""
    monomer = _ethane()
    result = solve_polygon(monomer, head_idx=0, tail_idx=1, n_sides=3, side_length_nm=5.0)
    assert result.ok, result.message
    assert result.closure_gap_nm is not None
    assert result.closing_side_length_nm is not None

    build = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=result.n_units, torsion_deg=result.torsion_deg)
    positions = _tail_positions(build)
    measured_gap = np.linalg.norm(positions[-1] - positions[0])
    assert abs(measured_gap / 10.0 - result.closure_gap_nm) < 1e-3

    # This monomer's best available corner root (flatness ~0.94, i.e. the two real
    # corners' turn axes are already ~20 deg apart) has a real, measured geometric
    # floor on how well ANY closing-side length can close the loop -- the residual
    # gap is dominated by the *perpendicular* component of that axis mismatch, which
    # length alone cannot correct (confirmed via a direct scratch measurement: the
    # closing side converges to very nearly the nominal 5 nm side length itself,
    # with ~1.8 nm of unclosable gap left over). So this is a bounded-improvement
    # check, not a "closes almost perfectly" claim -- it must be well under a full
    # side length (which is what an unconnected two-arm shape -- the original bug --
    # would show), not near-zero.
    assert result.closure_gap_nm < 0.5 * result.requested_length_nm, (
        f"closure gap {result.closure_gap_nm:.2f} nm is too large relative to the "
        f"{result.requested_length_nm:.2f} nm side length -- the triangle isn't closing"
    )


def test_solve_polygon_rejects_unachievable_size():
    monomer = _ethane()
    result = solve_polygon(monomer, head_idx=0, tail_idx=1, n_sides=4, side_length_nm=1000.0)
    assert not result.ok
    assert result.message


def test_solve_polygon_supports_side_lengths_beyond_the_natural_single_unit_chord():
    """Regression test for the live bug report: a triangle side length far beyond
    what one repeat unit's own bond geometry can naturally reach (ethane's natural
    chord is a few tenths of a nm) used to be rejected outright ("side length isn't
    independently adjustable once the shape and monomer are fixed"). A straight
    multi-unit run per side should now reach it."""
    monomer = _ethane()
    result = solve_polygon(monomer, head_idx=0, tail_idx=1, n_sides=3, side_length_nm=5.0)
    assert result.ok, result.message
    assert isinstance(result.torsion_deg, list)
    assert len(result.torsion_deg) == result.n_units - 1
    assert result.units_per_edge > 1
    assert result.relative_error < 0.05


def test_solve_polygon_reports_best_achievable_when_still_infeasible():
    """Even with multi-unit edges, a request far beyond what's reachable within
    max_units must still be cleanly rejected (not silently truncated), with the
    message stating what's actually achievable at that ceiling."""
    monomer = _ethane()
    result = solve_polygon(monomer, head_idx=0, tail_idx=1, n_sides=3, side_length_nm=500.0, max_units=10)
    assert not result.ok
    assert "nm" in result.message


def test_solve_circle_matches_requested_radius_within_tolerance():
    monomer = _ethane()
    # Probe a plausible radius first (rather than guessing) and then ask the solver
    # for exactly that -- isolates "does the solver hit its own target" from "did the
    # test happen to pick an achievable radius."
    probe = solve_circle(monomer, head_idx=0, tail_idx=1, radius_nm=0.3)
    assert probe.ok, probe.message
    target_radius_nm = probe.achieved_length_nm

    result = solve_circle(monomer, head_idx=0, tail_idx=1, radius_nm=target_radius_nm)
    assert result.ok, result.message
    assert result.relative_error < 0.1

    # What solve_circle actually guarantees (see module docstring): a constant
    # per-junction turn angle matching result.turn_angle_deg (every junction is a
    # "facet", one repeat unit each -- no straight runs), not a guaranteed flat/
    # closed ring.
    build = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=result.n_units, torsion_deg=result.torsion_deg)
    positions = _tail_positions(build)
    corner_bends = _bend_angles_deg(positions)
    assert np.allclose(corner_bends, result.turn_angle_deg, atol=2.0), (
        f"expected every junction to turn by {result.turn_angle_deg:.1f} deg, got {corner_bends}"
    )


def test_solve_circle_uses_one_unit_per_facet_not_straight_runs():
    """Regression test for the live bug report: circles used to reuse solve_polygon's
    straight-run-plus-corner model, which let a straight run pad any n_edges out to
    the requested radius -- so the search had no reason to prefer many small facets
    over a handful of big ones, and a requested circle rendered as a visible
    pentagon/hexagon. solve_circle now always uses exactly one repeat unit per facet
    (see its docstring) and a single uniform torsion, never a per-junction schedule.
    Needs a flexible monomer (see _hexane()'s docstring) -- ethane structurally can't
    reach the small turn angles a many-facet circle needs, which is its own separate,
    correctly-rejected case below."""
    monomer = _hexane()
    result = solve_circle(monomer, head_idx=0, tail_idx=5, radius_nm=3.0)
    assert result.ok, result.message
    assert result.units_per_edge == 1
    assert isinstance(result.torsion_deg, float)
    # Should need many more facets than any regular polygon in this suite (3-6
    # sides) -- otherwise it's indistinguishable from a crude polygon, the bug.
    assert result.n_units - 1 > 10, (
        f"expected a smooth many-facet circle, got only {result.n_units - 1} facets"
    )
    assert result.relative_error < 0.05


def test_solve_circle_rejects_cleanly_when_monomer_too_rigid_for_requested_radius():
    """Ethane's step-vector turn angle can never go below ~54.6 deg (direct scan,
    see chem/shape_solver.py's docstring), so at most 6 facets are geometrically
    reachable with no straight-run padding -- capping the achievable radius far
    below 3 nm. The new no-padding model must say so cleanly, not silently return
    a small/wrong ring or (as the old straight-run model did) fake a good-looking
    size fit that was actually a crude few-sided polygon in disguise."""
    monomer = _ethane()
    result = solve_circle(monomer, head_idx=0, tail_idx=1, radius_nm=3.0)
    assert not result.ok
    assert "nm" in result.message


def test_solve_polygon_and_solve_circle_report_flatness_score():
    ethane = _ethane()
    roots, (lo, hi) = solve_torsion_for_turn_angle(ethane, head_idx=0, tail_idx=1, target_turn_deg=90.0)
    assert roots, f"ethane can't reach a 90 degree turn (range [{lo:.1f}, {hi:.1f}])"
    natural_side_nm = max(r[1] for r in roots) / 10.0

    polygon_result = solve_polygon(ethane, head_idx=0, tail_idx=1, n_sides=4, side_length_nm=natural_side_nm)
    assert polygon_result.ok, polygon_result.message
    assert polygon_result.flatness_score is not None
    assert 0.0 <= polygon_result.flatness_score <= 1.0

    circle_result = solve_circle(_hexane(), head_idx=0, tail_idx=5, radius_nm=3.0)
    assert circle_result.ok, circle_result.message
    assert circle_result.flatness_score is not None
    assert 0.0 <= circle_result.flatness_score <= 1.0
