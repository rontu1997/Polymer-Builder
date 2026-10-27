"""Join two structures at explicitly selected atoms with one single bond."""

import bisect
import math
from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np
from rdkit import Chem

from chem.builder import _leave_direction, _place_incoming_unit, _vec
from chem.validation import MAX_ATOMS_HARD_CAP, validate_linkage_atom


class MergeError(Exception):
    """An input cannot form the requested junction."""


@dataclass
class MergeResult:
    mol: Chem.Mol
    new_bonds: list
    torsion_specs: list


def _validate_input(mol, atom_idx, label):
    if mol is None or mol.GetNumAtoms() == 0:
        raise MergeError(f'{label}: upload a nonempty PDB structure.')
    if isinstance(atom_idx, bool) or not isinstance(atom_idx, Integral):
        raise MergeError(f'{label}: select a connection atom.')
    if mol.GetNumConformers() == 0:
        raise MergeError(f'{label}: the structure has no 3D coordinates.')
    if mol.GetNumConformers() != 1:
        raise MergeError(f'{label}: supply a single model/conformer per PDB.')
    if not np.isfinite(mol.GetConformer().GetPositions()).all():
        raise MergeError(f'{label}: all atom coordinates must be finite.')
    if len(Chem.GetMolFrags(mol)) != 1:
        raise MergeError(f'{label}: upload one connected component (no solvent or counter-ions).')
    slot = validate_linkage_atom(mol, int(atom_idx))
    if not slot.ok:
        raise MergeError(f'{label}: {slot.message}')
    return slot


def _reference_atom(mol, atom_idx, slot):
    """Choose a retained, noncollinear neighbor for the junction dihedral."""
    neighbors = sorted(
        mol.GetAtomWithIdx(atom_idx).GetNeighbors(),
        key=lambda atom: (atom.GetAtomicNum() == 1, atom.GetIdx()),
    )
    origin = _vec(mol.GetConformer(), atom_idx)
    direction = _leave_direction(mol, atom_idx, slot)
    for neighbor in neighbors:
        if neighbor.GetIdx() == slot.explicit_h_idx:
            continue
        reference = _vec(mol.GetConformer(), neighbor.GetIdx()) - origin
        if np.linalg.norm(np.cross(reference, direction)) > 1e-6:
            return neighbor.GetIdx()
    raise MergeError(
        f'Atom {atom_idx} has no retained neighbor with usable geometry for a junction '
        'torsion. Choose another atom or supply a PDB with explicit hydrogens and valid coordinates.'
    )


def _label_source(mol, chain):
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            info = Chem.AtomPDBResidueInfo()
            info.SetName(f'{atom.GetSymbol():>2}{atom.GetIdx() % 100:02}')
            info.SetResidueName('MOL')
            info.SetResidueNumber(1)
            atom.SetMonomerInfo(info)
            info = atom.GetPDBResidueInfo()
        info.SetChainId(chain)
        info.SetIsHeteroAtom(True)


def _complete_hydrogens(mol):
    # Assign new H metadata ourselves: RDKit's addResidueInfo path can hang on
    # nonstandard residue numbering. Existing atom metadata is copied intact.
    result = Chem.AddHs(mol, addCoords=True)
    for idx in range(mol.GetNumAtoms(), result.GetNumAtoms()):
        atom = result.GetAtomWithIdx(idx)
        parent_info = atom.GetNeighbors()[0].GetPDBResidueInfo()
        if parent_info is not None:
            atom.SetMonomerInfo(parent_info)
            atom.GetPDBResidueInfo().SetName(f'H{idx % 1000:03}')
    return result


def merge_structures(first, second, first_idx, second_idx, torsion_deg=180.0):
    """Place a copy of second at first's selected atom, then form one single bond.

    Missing hydrogens are added on working copies. Leaving H atoms are removed
    and all returned indices are remapped. Both inputs are unchanged.
    """
    if isinstance(torsion_deg, bool) or not isinstance(torsion_deg, Real) or not math.isfinite(torsion_deg):
        raise MergeError('Junction torsion must be a finite number of degrees.')
    if first is not None and second is not None:
        limit = min(MAX_ATOMS_HARD_CAP, 99_999)
        if first.GetNumAtoms() + second.GetNumAtoms() > limit:
            raise MergeError(f'The combined inputs exceed the {limit:,}-atom PDB limit.')
    _validate_input(first, first_idx, 'PDB A')
    _validate_input(second, second_idx, 'PDB B')
    first = _complete_hydrogens(first)
    second = _complete_hydrogens(second)
    if first.GetNumAtoms() + second.GetNumAtoms() > min(MAX_ATOMS_HARD_CAP, 99_999):
        raise MergeError('The combined inputs exceed the atom limit after adding hydrogens.')
    first_slot = _validate_input(first, first_idx, 'PDB A')
    second_slot = _validate_input(second, second_idx, 'PDB B')
    first_ref = _reference_atom(first, first_idx, first_slot)
    second_ref = _reference_atom(second, second_idx, second_slot)
    torsion = ((float(torsion_deg) + 180.0) % 360.0) - 180.0

    placed = _place_incoming_unit(
        first, first_idx, first_ref, first_slot,
        second, second_idx, second_ref, second_slot, torsion,
    )
    first_copy = Chem.Mol(first)
    _label_source(first_copy, 'A')
    _label_source(placed, 'B')
    offset = first_copy.GetNumAtoms()
    combined = Chem.RWMol(Chem.CombineMols(first_copy, placed))
    combined.AddBond(first_idx, offset + second_idx, Chem.BondType.SINGLE)
    removed = sorted(
        slot.explicit_h_idx + base
        for slot, base in ((first_slot, 0), (second_slot, offset))
        if slot.slot_kind == 'explicit_h'
    )
    for idx in reversed(removed):
        combined.RemoveAtom(idx)

    def remap(idx):
        return idx - bisect.bisect_left(removed, idx)

    result = combined.GetMol()
    try:
        Chem.SanitizeMol(result)
    except (ValueError, RuntimeError) as exc:
        raise MergeError(f'The selected atoms cannot form a valid single bond: {exc}') from exc
    if not np.isfinite(result.GetConformer().GetPositions()).all():
        raise MergeError('The selected atoms produced invalid junction geometry.')
    return MergeResult(
        mol=result,
        new_bonds=[(remap(first_idx), remap(second_idx + offset))],
        torsion_specs=[(
            remap(first_ref), remap(first_idx), remap(second_idx + offset),
            remap(second_ref + offset), torsion,
        )],
    )
