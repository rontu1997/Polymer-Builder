import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rdkit import Chem
from rdkit.Chem import AllChem

from chem.validation import validate_linkage_atom, validate_pre_build, MAX_REPEAT_UNITS


def _embed(mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    return mol


def test_sp3_ch3_carbon_is_valid_explicit_h_slot():
    mol = _embed(Chem.MolFromSmiles("CC"))  # ethane
    # atom 0 is a carbon with 3 H neighbors
    check = validate_linkage_atom(mol, 0)
    assert check.ok
    assert check.slot_kind == "explicit_h"
    assert mol.GetAtomWithIdx(check.explicit_h_idx).GetAtomicNum() == 1


def test_sp3_ch2_carbon_is_valid():
    mol = _embed(Chem.MolFromSmiles("CCC"))  # propane, middle CH2 carbon is atom 1
    check = validate_linkage_atom(mol, 1)
    assert check.ok
    assert check.slot_kind == "explicit_h"


def test_quaternary_carbon_is_rejected():
    mol = _embed(Chem.MolFromSmiles("CC(C)(C)C"))  # neopentane, central C is atom 1
    check = validate_linkage_atom(mol, 1)
    assert not check.ok
    assert "maximum valence" in check.message


def test_aromatic_ch_is_valid():
    mol = _embed(Chem.MolFromSmiles("c1ccccc1"))  # benzene, any ring atom has an H
    check = validate_linkage_atom(mol, 0)
    assert check.ok
    assert check.slot_kind == "explicit_h"


def test_fully_substituted_aromatic_carbon_is_rejected():
    # toluene: ipso ring carbon (attached to the methyl group) has no H
    mol = Chem.MolFromSmiles("Cc1ccccc1")
    ipso_idx = [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetIsAromatic() and any(n.GetSymbol() == "C" and not n.GetIsAromatic() for n in a.GetNeighbors())
    ][0]
    mol = _embed(mol)
    check = validate_linkage_atom(mol, ipso_idx)
    assert not check.ok


def test_implicit_h_slot_without_explicit_hydrogens():
    # No AddHs -- atom valence is satisfied by implicit Hs only.
    mol = Chem.MolFromSmiles("CC")
    check = validate_linkage_atom(mol, 0)
    assert check.ok
    assert check.slot_kind == "implicit"


def test_out_of_range_index_rejected():
    mol = _embed(Chem.MolFromSmiles("CC"))
    check = validate_linkage_atom(mol, 999)
    assert not check.ok


def test_pre_build_rejects_same_head_and_tail():
    mol = _embed(Chem.MolFromSmiles("CC"))
    result = validate_pre_build(mol, 0, 0)
    assert not result.ok


def test_pre_build_rejects_disconnected_fragments():
    mol = _embed(Chem.MolFromSmiles("CC.O"))  # ethane + a stray water
    result = validate_pre_build(mol, 0, 1)
    assert not result.ok
    assert "disconnected" in result.message


def test_pre_build_rejects_excessive_repeat_count():
    mol = _embed(Chem.MolFromSmiles("CC"))
    result = validate_pre_build(mol, 0, 1, n_units=MAX_REPEAT_UNITS + 1)
    assert not result.ok


def test_pre_build_rejects_atom_count_over_cap():
    mol = _embed(Chem.MolFromSmiles("CC"))
    result = validate_pre_build(mol, 0, 1, n_units=50, max_atoms=10)
    assert not result.ok
    assert "atoms" in result.message


def test_pre_build_accepts_valid_case():
    mol = _embed(Chem.MolFromSmiles("CC"))
    result = validate_pre_build(mol, 0, 1, n_units=5)
    assert result.ok
