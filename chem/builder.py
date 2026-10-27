"""Build a linear homopolymer chain from a monomer by repeatedly bonding a
designated tail atom of the growing chain to the designated head atom of a fresh
monomer copy, at a user-chosen backbone torsional angle.

Placement uses an explicit rigid-body transform (not a bare two-atom alignment):
for each new unit, the incoming monomer copy is (1) rotated so its head atom's
"leaving direction" aligns with the growing chain's tail "leaving direction", then
(2) translated so the head atom lands at the correct bond length, then (3) twisted
about the new bond axis so the T_ref-T-H-H_ref dihedral matches the target torsion.
That's the full six rigid-body degrees of freedom, so the torsion setting is
well-defined rather than left to whatever the optimizer happens to settle on.
"""

import bisect
from dataclasses import dataclass, replace

import numpy as np
from rdkit import Chem

from chem.validation import (
    ValenceCheck,
    expected_bond_length,
    validate_linkage_atom,
    validate_pre_build,
)


class PolymerBuildError(Exception):
    pass


@dataclass
class BuildResult:
    mol: Chem.Mol
    new_bonds: list  # list[(atom_idx_a, atom_idx_b)] in the final mol's numbering
    torsion_specs: list  # list[(a1, a2, a3, a4, target_deg)] in the final mol's numbering
    unit_tail_indices: list  # final index of each unit's tail atom, in unit order


def _normalize_angle(deg: float) -> float:
    return ((deg + 180.0) % 360.0) - 180.0


def _vec(conf, idx) -> np.ndarray:
    p = conf.GetAtomPosition(idx)
    return np.array([p.x, p.y, p.z])


def _pick_reference_neighbor(mol: Chem.Mol, atom_idx: int) -> int:
    """A bonded neighbor to use as the 4th atom of a backbone dihedral. Prefers a
    heavy atom (deterministically, lowest index) since that's almost always the
    chemically meaningful choice; falls back to any neighbor (including H) only
    if the atom has no heavy-atom neighbor at all.
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    neighbors = sorted(atom.GetNeighbors(), key=lambda a: a.GetIdx())
    heavy = [n.GetIdx() for n in neighbors if n.GetAtomicNum() > 1]
    if heavy:
        return heavy[0]
    if neighbors:
        return neighbors[0].GetIdx()
    raise PolymerBuildError(
        f"Atom {atom_idx} has no bonded neighbors; cannot define a backbone dihedral."
    )


def _leave_direction(mol: Chem.Mol, atom_idx: int, slot: ValenceCheck) -> np.ndarray:
    """Unit vector pointing from `atom_idx` towards where its new bond partner
    should go. When the atom has a real leaving hydrogen, that hydrogen's position
    already marks the chemically correct direction (it's the exact spot the new
    substituent replaces). Otherwise, fall back to "away from existing heavy
    neighbors" as a generic open-valence heuristic.
    """
    conf = mol.GetConformer()
    origin = _vec(conf, atom_idx)

    if slot.slot_kind == "explicit_h" and slot.explicit_h_idx is not None:
        vec = _vec(conf, slot.explicit_h_idx) - origin
    else:
        atom = mol.GetAtomWithIdx(atom_idx)
        nbr_dirs = []
        for nbr in atom.GetNeighbors():
            d = _vec(conf, nbr.GetIdx()) - origin
            n = np.linalg.norm(d)
            if n > 1e-6:
                nbr_dirs.append(d / n)
        vec = -np.mean(nbr_dirs, axis=0) if nbr_dirs else np.array([1.0, 0.0, 0.0])

    norm = np.linalg.norm(vec)
    if norm < 1e-6:
        return np.array([1.0, 0.0, 0.0])
    return vec / norm


def _rotation_aligning(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation matrix R such that R @ a == b, for unit vectors a, b."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = np.dot(a, b)
    if s < 1e-8:
        if c > 0:
            return np.eye(3)
        # a and b are anti-parallel: rotate 180 deg about any axis perpendicular to a
        perp = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, perp)
        axis = axis / np.linalg.norm(axis)
        return _axis_angle_rotation(axis, 180.0)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s**2))


def _axis_angle_rotation(axis: np.ndarray, theta_deg: float) -> np.ndarray:
    """Rodrigues' rotation formula: rotate by theta_deg (right-hand rule) about
    `axis` (unit vector)."""
    axis = axis / np.linalg.norm(axis)
    theta = np.radians(theta_deg)
    k = axis
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(theta) * kx + (1 - np.cos(theta)) * (kx @ kx)


def _dihedral_deg(p1, p2, p3, p4) -> float:
    """Matches rdkit.Chem.rdMolTransforms.GetDihedralDeg's sign convention
    (verified numerically against it) -- negated relative to the textbook
    Praxeolitic formula, whose raw sign runs the other way."""
    b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x, y = np.dot(n1, n2), np.dot(m1, n2)
    return -np.degrees(np.arctan2(y, x))


def _place_incoming_unit(
    growing: Chem.Mol,
    tail_idx: int,
    tail_ref_idx: int,
    tail_slot: ValenceCheck,
    fresh: Chem.Mol,
    head_idx: int,
    head_ref_idx: int,
    head_slot: ValenceCheck,
    torsion_deg: float,
) -> Chem.Mol:
    """Return a copy of `fresh` with its conformer rigidly transformed into world
    coordinates so its head atom bonds correctly onto the growing chain's tail."""
    growing_conf = growing.GetConformer()
    tail_pos = _vec(growing_conf, tail_idx)
    tail_ref_pos = _vec(growing_conf, tail_ref_idx)
    tail_leave_dir = _leave_direction(growing, tail_idx, tail_slot)

    bond_length = expected_bond_length(
        growing.GetAtomWithIdx(tail_idx).GetSymbol(),
        fresh.GetAtomWithIdx(head_idx).GetSymbol(),
    )
    target_head_pos = tail_pos + bond_length * tail_leave_dir

    fresh_conf = fresh.GetConformer()
    head_pos0 = _vec(fresh_conf, head_idx)
    head_leave_dir0 = _leave_direction(fresh, head_idx, head_slot)

    n_atoms = fresh.GetNumAtoms()
    coords0 = np.array([_vec(fresh_conf, i) for i in range(n_atoms)])
    shifted = coords0 - head_pos0

    r1 = _rotation_aligning(head_leave_dir0, -tail_leave_dir)
    rotated = shifted @ r1.T
    placed = rotated + target_head_pos  # head atom now sits exactly at target_head_pos

    head_ref_pos_current = placed[head_ref_idx]
    current_dihedral = _dihedral_deg(tail_ref_pos, tail_pos, target_head_pos, head_ref_pos_current)
    delta = _normalize_angle(torsion_deg - current_dihedral)

    axis = target_head_pos - tail_pos
    r2 = _axis_angle_rotation(axis, delta)
    final = (placed - target_head_pos) @ r2.T + target_head_pos

    result = Chem.Mol(fresh)
    result_conf = result.GetConformer()
    for i in range(n_atoms):
        result_conf.SetAtomPosition(i, final[i].tolist())
    return result


def _get_base_residue_number(mol: Chem.Mol):
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is not None:
            return info.GetResidueNumber()
    return None


def _set_unit_residue_info(mol: Chem.Mol, atom_indices, residue_number: int) -> None:
    """Give this unit its own residue number and mark its atoms as hetero atoms.

    Without this, every unit -- being a copy of the same monomer -- keeps the
    *same* residue number, so a multi-unit chain looks like one overlapping
    residue rather than a sequence of units. Separately, many viewers trust a
    biomolecule's own residue/atom-name templates over raw CONECT records for
    plain ATOM entries; since our new inter-unit bond is whatever the user
    designated (not necessarily a real peptide/glycosidic bond a template would
    recognize), forcing every atom to HETATM makes viewers fall back to reading
    the (correct, symmetric) CONECT records instead of a template that has no
    idea this bond exists.
    """
    for idx in atom_indices:
        info = mol.GetAtomWithIdx(idx).GetPDBResidueInfo()
        if info is None:
            continue
        info.SetResidueNumber(residue_number)
        info.SetIsHeteroAtom(True)


def build_polymer(monomer: Chem.Mol, head_idx: int, tail_idx: int, n_units: int, torsion_deg) -> BuildResult:
    """`torsion_deg` is either a single angle (applied at every junction, as always) or a
    list/tuple of length `n_units - 1` giving each junction its own target angle -- e.g. a
    straight run of junctions at ~180 deg ending in one sharp corner junction, for shape
    presets that need a polygon side longer than one repeat unit can naturally reach.
    """
    pre = validate_pre_build(monomer, head_idx, tail_idx, n_units=n_units)
    if not pre.ok:
        raise PolymerBuildError(pre.message)

    if isinstance(torsion_deg, (list, tuple)):
        torsion_list = [_normalize_angle(t) for t in torsion_deg]
        if len(torsion_list) != n_units - 1:
            raise PolymerBuildError(
                f"Expected {n_units - 1} torsion values (one per junction), got {len(torsion_list)}."
            )
    else:
        torsion_list = [_normalize_angle(torsion_deg)] * (n_units - 1)

    head_slot = validate_linkage_atom(monomer, head_idx)
    tail_slot = validate_linkage_atom(monomer, tail_idx)

    combined = Chem.RWMol(monomer)
    base_resnum = _get_base_residue_number(monomer)
    if base_resnum is not None:
        _set_unit_residue_info(combined, range(monomer.GetNumAtoms()), base_resnum)

    if n_units == 1:
        mol = combined.GetMol()
        Chem.SanitizeMol(mol)
        return BuildResult(mol=mol, new_bonds=[], torsion_specs=[], unit_tail_indices=[tail_idx])

    current_tail_idx = tail_idx
    current_tail_unit_offset = 0  # offset of whichever unit currently owns the growing chain's tail atom
    current_tail_ref_idx = _pick_reference_neighbor(combined, tail_idx)

    new_bonds_raw = []
    torsion_specs_raw = []
    unit_tail_indices_raw = [tail_idx]  # unit 0's (raw) tail atom index; one appended per unit thereafter
    removal_set = set()

    for unit_index in range(1, n_units):
        current_torsion = torsion_list[unit_index - 1]
        fresh = Chem.Mol(monomer)
        fresh_head_ref_idx = _pick_reference_neighbor(fresh, head_idx)

        # tail_slot.explicit_h_idx was computed once against the original monomer
        # template (offset 0) -- for every unit past the first, the tail atom's
        # own leaving H lives at that same relative position *within its own
        # unit*, i.e. offset by current_tail_unit_offset. Without this, _leave_
        # direction() looks up unit 1's hydrogen instead of the current tail
        # unit's, which silently corrupts the incoming unit's placement direction
        # (the dihedral still gets forced to target by the later twist step, so
        # this doesn't show up as a bad dihedral -- only as bad global geometry).
        current_tail_slot = tail_slot
        if tail_slot.slot_kind == "explicit_h":
            current_tail_slot = replace(
                tail_slot, explicit_h_idx=tail_slot.explicit_h_idx + current_tail_unit_offset
            )

        placed_fresh = _place_incoming_unit(
            growing=combined.GetMol(),
            tail_idx=current_tail_idx,
            tail_ref_idx=current_tail_ref_idx,
            tail_slot=current_tail_slot,
            fresh=fresh,
            head_idx=head_idx,
            head_ref_idx=fresh_head_ref_idx,
            head_slot=head_slot,
            torsion_deg=current_torsion,
        )

        offset = combined.GetNumAtoms()
        combined.InsertMol(placed_fresh)
        if base_resnum is not None:
            _set_unit_residue_info(
                combined, range(offset, offset + placed_fresh.GetNumAtoms()), base_resnum + unit_index
            )

        new_head_idx = head_idx + offset
        new_head_ref_idx = fresh_head_ref_idx + offset
        new_tail_idx = tail_idx + offset

        combined.AddBond(current_tail_idx, new_head_idx, Chem.BondType.SINGLE)
        new_bonds_raw.append((current_tail_idx, new_head_idx))
        torsion_specs_raw.append(
            (current_tail_ref_idx, current_tail_idx, new_head_idx, new_head_ref_idx, current_torsion)
        )

        if tail_slot.slot_kind == "explicit_h":
            removal_set.add(tail_slot.explicit_h_idx + current_tail_unit_offset)
        if head_slot.slot_kind == "explicit_h":
            removal_set.add(head_slot.explicit_h_idx + offset)

        current_tail_idx = new_tail_idx
        current_tail_unit_offset = offset
        current_tail_ref_idx = _pick_reference_neighbor(combined, new_tail_idx)
        unit_tail_indices_raw.append(new_tail_idx)

    for idx in sorted(removal_set, reverse=True):
        combined.RemoveAtom(idx)

    removed_sorted = sorted(removal_set)

    def remap(idx: int) -> int:
        return idx - bisect.bisect_left(removed_sorted, idx)

    new_bonds = [(remap(a), remap(b)) for a, b in new_bonds_raw]
    torsion_specs = [
        (remap(a1), remap(a2), remap(a3), remap(a4), t) for a1, a2, a3, a4, t in torsion_specs_raw
    ]
    unit_tail_indices = [remap(idx) for idx in unit_tail_indices_raw]

    mol = combined.GetMol()
    Chem.SanitizeMol(mol)

    return BuildResult(
        mol=mol, new_bonds=new_bonds, torsion_specs=torsion_specs, unit_tail_indices=unit_tail_indices
    )
