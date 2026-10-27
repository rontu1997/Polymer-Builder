"""Feasibility checks for polymer construction: can two atoms actually be bonded
together, and once built, is the resulting geometry sane?

This module intentionally does not try to infer polymerization chemistry -- see
the project README's scope statement. It only checks whether the two
user-designated atoms can structurally accept the new bond, and whether the
built/optimized geometry that results is plausible.
"""

from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

MAX_ATOMS_HARD_CAP = 100_000
MAX_REPEAT_UNITS = 1_000

# Common backbone bond lengths (Angstrom), keyed by a frozenset of atomic symbols.
# Anything not in this table falls back to summed covalent radii.
_BOND_LENGTH_TABLE = {
    frozenset(("C", "C")): 1.54,
    frozenset(("C", "N")): 1.47,
    frozenset(("C", "O")): 1.43,
    frozenset(("C", "S")): 1.82,
    frozenset(("N", "N")): 1.45,
    frozenset(("Si", "O")): 1.63,
    frozenset(("Si", "Si")): 2.33,
}

_PT = Chem.GetPeriodicTable()


@dataclass
class ValenceCheck:
    ok: bool
    message: str = ""
    # How the new bond's valence should be freed up on this atom:
    #  - "explicit_h": delete the real H atom at explicit_h_idx, then add the bond
    #  - "implicit":   just add the bond; RDKit's implicit-H accounting absorbs it
    slot_kind: str = None
    explicit_h_idx: int = None


@dataclass
class ValidationResult:
    ok: bool
    message: str = ""
    warnings: list = field(default_factory=list)
    torsion_report: list = field(default_factory=list)


def _max_valence(atom) -> int | None:
    """Default allowed valence for this atom, adjusted for formal charge. Uses
    RDKit's *typical* valence for the element (not the full extended-valence
    table), so it can be conservative for hypervalent elements like S or P --
    e.g. a sulfone sulfur may be flagged as "full" even though sulfur can extend
    past its default valence. Acceptable for V1; revisit with an element-specific
    table if that turns out to matter in practice.
    """
    default = _PT.GetDefaultValence(atom.GetAtomicNum())
    if default < 0:
        return None
    return default + atom.GetFormalCharge()


def validate_linkage_atom(mol: Chem.Mol, atom_idx: int) -> ValenceCheck:
    """Can this atom accept one more bond?

    A saturated atom with a removable hydrogen (explicit or implicit) normally
    already has `current_valence == max_valence` -- that H *is* what's using up
    its last valence slot. So the feasibility test is not "is there valence
    headroom" but "is there a hydrogen to give up, or is there headroom already
    (an under-coordinated atom with no H at all)". This also means a bare
    aromatic ring C-H is a legitimate linkage point (ordinary aromatic
    substitution chemistry) -- only atoms with neither headroom nor a removable H
    (e.g. a quaternary carbon, or a fully-substituted ring position) are rejected.
    """
    if atom_idx is None or atom_idx < 0 or atom_idx >= mol.GetNumAtoms():
        return ValenceCheck(False, f"Atom index {atom_idx} does not exist on this molecule.")

    atom = mol.GetAtomWithIdx(atom_idx)
    max_v = _max_valence(atom)
    if max_v is None:
        return ValenceCheck(
            False,
            f"Cannot determine a reliable valence limit for element "
            f"{atom.GetSymbol()} (atom {atom_idx}); please pick a different atom.",
        )

    current_v = atom.GetTotalValence()
    if current_v > max_v:
        return ValenceCheck(
            False,
            f"Atom {atom_idx} ({atom.GetSymbol()}) has valence {current_v}, above "
            f"its expected maximum of {max_v} -- this structure may be malformed.",
        )

    explicit_h_idx = next(
        (nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1),
        None,
    )
    has_implicit_h = atom.GetNumImplicitHs() > 0

    if current_v < max_v:
        # Genuine open valence -- no hydrogen needs to be sacrificed.
        return ValenceCheck(True, "OK", slot_kind="none")
    if explicit_h_idx is not None:
        return ValenceCheck(True, "OK", slot_kind="explicit_h", explicit_h_idx=explicit_h_idx)
    if has_implicit_h:
        return ValenceCheck(True, "OK", slot_kind="implicit")

    return ValenceCheck(
        False,
        f"Atom {atom_idx} ({atom.GetSymbol()}) is already at its maximum valence "
        f"({current_v}/{max_v}) and has no hydrogen to give up for a new bond.",
    )


def validate_pre_build(
    monomer: Chem.Mol,
    head_idx: int,
    tail_idx: int,
    n_units: int | None = None,
    max_atoms: int = MAX_ATOMS_HARD_CAP,
    max_units: int = MAX_REPEAT_UNITS,
) -> ValidationResult:
    if head_idx == tail_idx:
        return ValidationResult(False, "Head and tail atom must be different atoms.")

    for role, idx in (("Head", head_idx), ("Tail", tail_idx)):
        check = validate_linkage_atom(monomer, idx)
        if not check.ok:
            return ValidationResult(False, f"{role} atom invalid: {check.message}")

    if monomer.GetNumConformers() == 0:
        return ValidationResult(False, "Monomer has no 3D conformer to build from.")

    frags = Chem.GetMolFrags(monomer)
    if len(frags) > 1:
        return ValidationResult(
            False,
            f"Structure contains {len(frags)} disconnected components -- please "
            f"upload a single monomer (no solvent/counter-ions).",
        )

    if n_units is not None:
        if n_units < 1:
            return ValidationResult(False, "Repeat count must be at least 1.")
        if n_units > max_units:
            return ValidationResult(
                False, f"Repeat count {n_units} exceeds the maximum of {max_units}."
            )
        estimated_atoms = monomer.GetNumAtoms() * n_units
        if estimated_atoms > max_atoms:
            return ValidationResult(
                False,
                f"This polymer would contain approximately {estimated_atoms} atoms, "
                f"above the {max_atoms}-atom limit. Please reduce the repeat count.",
            )

    return ValidationResult(True, "OK")


def expected_bond_length(symbol_a: str, symbol_b: str) -> float:
    table_hit = _BOND_LENGTH_TABLE.get(frozenset((symbol_a, symbol_b)))
    if table_hit is not None:
        return table_hit
    ra = _PT.GetRcovalent(_PT.GetAtomicNumber(symbol_a))
    rb = _PT.GetRcovalent(_PT.GetAtomicNumber(symbol_b))
    return ra + rb


def _normalize_angle(deg: float) -> float:
    return ((deg + 180.0) % 360.0) - 180.0


def _angle_delta(a: float, b: float) -> float:
    """Smallest signed difference between two angles, periodic-aware."""
    return _normalize_angle(a - b)


def validate_post_build(
    mol: Chem.Mol,
    new_bonds: list[tuple[int, int]],
    torsion_specs: list[tuple[int, int, int, int, float]],
    energy_before: float,
    energy_after: float,
    bond_length_range=(0.6, 1.6),
    clash_factor: float = 0.6,
) -> ValidationResult:
    warnings = []

    mol_copy = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(mol_copy)
    except Exception as exc:
        return ValidationResult(False, f"Built structure failed sanitization: {exc}")

    frags = Chem.GetMolFrags(mol_copy)
    if len(frags) > 1:
        return ValidationResult(
            False, f"Built structure has {len(frags)} disconnected fragments (expected 1)."
        )

    conf = mol_copy.GetConformer()

    for a1, a2 in new_bonds:
        pos1, pos2 = conf.GetAtomPosition(a1), conf.GetAtomPosition(a2)
        dist = pos1.Distance(pos2)
        sym1 = mol_copy.GetAtomWithIdx(a1).GetSymbol()
        sym2 = mol_copy.GetAtomWithIdx(a2).GetSymbol()
        expected = expected_bond_length(sym1, sym2)
        low, high = bond_length_range[0] * expected, bond_length_range[1] * expected
        if not (low <= dist <= high):
            warnings.append(
                f"New bond {a1}-{a2} ({sym1}-{sym2}) has length {dist:.2f} Å, "
                f"outside the expected range [{low:.2f}, {high:.2f}] Å."
            )

    bonded_pairs = {
        frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in mol_copy.GetBonds()
    }
    one_three_pairs = set()
    for atom in mol_copy.GetAtoms():
        nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                one_three_pairs.add(frozenset((nbrs[i], nbrs[j])))

    n_atoms = mol_copy.GetNumAtoms()
    positions = np.array([list(conf.GetAtomPosition(i)) for i in range(n_atoms)])
    radii = np.array(
        [_PT.GetRcovalent(mol_copy.GetAtomWithIdx(i).GetAtomicNum()) for i in range(n_atoms)]
    )
    diffs = positions[:, None, :] - positions[None, :, :]
    dist_matrix = np.sqrt((diffs**2).sum(axis=-1))
    min_allowed = clash_factor * (radii[:, None] + radii[None, :])

    clash_count = 0
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            pair = frozenset((i, j))
            if pair in bonded_pairs or pair in one_three_pairs:
                continue
            if dist_matrix[i, j] < min_allowed[i, j]:
                clash_count += 1
    if clash_count:
        warnings.append(f"{clash_count} pair(s) of non-bonded atoms are clashing (too close).")

    torsion_report = []
    for a1, a2, a3, a4, target in torsion_specs:
        actual = rdMolTransforms.GetDihedralDeg(conf, a1, a2, a3, a4)
        delta = _angle_delta(actual, target)
        torsion_report.append(
            {"atoms": (a1, a2, a3, a4), "target": target, "actual": actual, "delta": delta}
        )
        if abs(delta) > 15.0:
            warnings.append(
                f"Junction {a2}-{a3} torsion is {actual:.1f}° vs target {target:.1f}° "
                f"(off by {delta:.1f}°)."
            )

    if energy_after > energy_before + 1.0:
        warnings.append(
            f"Energy increased during optimization ({energy_before:.1f} -> "
            f"{energy_after:.1f} kcal/mol) -- the torsion restraint may be fighting "
            f"the force field."
        )

    return ValidationResult(True, "OK", warnings=warnings, torsion_report=torsion_report)
