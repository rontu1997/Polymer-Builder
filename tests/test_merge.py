import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from chem.merge import MergeError, merge_structures
from chem.optimize import optimize_geometry
from chem.viewer import mol_to_pdb_block, pdb_block_to_mol


def molecule(smiles, hydrogens=True):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=7) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    return mol if hydrogens else Chem.RemoveHs(mol)


@pytest.mark.parametrize('torsion', [-180.0, -60.0, 0.0, 70.0, 180.0])
def test_merge_preserves_inputs_and_builds_selected_bond_and_torsion(torsion):
    first, second = molecule('CC'), molecule('CO')
    original_a, original_b = Chem.MolToMolBlock(first), Chem.MolToMolBlock(second)
    result = merge_structures(first, second, 1, 0, torsion)
    assert Chem.MolToSmiles(Chem.RemoveHs(result.mol)) == 'CCCO'
    assert result.mol.GetNumAtoms() == 12  # C3H8O: 8 + 6 input atoms - 2 leaving H
    assert len(Chem.GetMolFrags(result.mol)) == 1
    assert len(result.new_bonds) == len(result.torsion_specs) == 1
    left, right = result.new_bonds[0]
    assert result.mol.GetBondBetweenAtoms(left, right).GetBondType() == Chem.BondType.SINGLE
    conf = result.mol.GetConformer()
    assert conf.GetAtomPosition(left).Distance(conf.GetAtomPosition(right)) == pytest.approx(1.54)
    a, b, c, d, target = result.torsion_specs[0]
    actual = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
    assert abs((actual - torsion + 180) % 360 - 180) < 1e-5
    assert np.allclose(conf.GetAtomPosition(0), first.GetConformer().GetAtomPosition(0))
    assert Chem.MolToMolBlock(first) == original_a
    assert Chem.MolToMolBlock(second) == original_b


def test_merge_implicit_hydrogens_and_optimize_roundtrip():
    result = merge_structures(molecule('CC', False), molecule('CO', False), 1, 0)
    assert result.mol.GetNumAtoms() == 12
    assert Chem.MolToSmiles(Chem.RemoveHs(result.mol)) == 'CCCO'
    optimized, report = optimize_geometry(result.mol, result.torsion_specs, result.new_bonds)
    assert report.validation.ok
    restored = pdb_block_to_mol(mol_to_pdb_block(optimized))
    assert len(Chem.GetMolFrags(restored)) == 1
    assert restored.GetBondBetweenAtoms(*result.new_bonds[0]) is not None


def test_reference_atom_is_not_a_removed_hydrogen():
    result = merge_structures(molecule('C'), molecule('C'), 0, 0, 60)
    assert result.mol.GetNumAtoms() == 8
    assert Chem.MolToSmiles(Chem.RemoveHs(result.mol)) == 'CC'
    a, b, c, d, target = result.torsion_specs[0]
    assert result.mol.GetBondBetweenAtoms(a, b) is not None
    assert result.mol.GetBondBetweenAtoms(c, d) is not None
    assert rdMolTransforms.GetDihedralDeg(result.mol.GetConformer(), a, b, c, d) == pytest.approx(60)


def test_residue_metadata_is_preserved_and_sources_have_distinct_chains():
    first, second = molecule('CC'), molecule('CO')
    for mol, resname in ((first, 'ONE'), (second, 'TWO')):
        for atom in mol.GetAtoms():
            info = Chem.AtomPDBResidueInfo()
            info.SetName(f'{atom.GetSymbol()}{atom.GetIdx():<3}'[:4])
            info.SetResidueName(resname)
            info.SetResidueNumber(10 + atom.GetIdx())
            info.SetChainId('X')
            atom.SetMonomerInfo(info)
    result = merge_structures(first, second, 1, 0)
    restored = pdb_block_to_mol(mol_to_pdb_block(result.mol))
    assert {atom.GetPDBResidueInfo().GetChainId() for atom in restored.GetAtoms()} == {'A', 'B'}
    for atom in result.mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        original = first if info.GetChainId() == 'A' else second
        original_idx = info.GetResidueNumber() - 10
        assert info.GetName() == original.GetAtomWithIdx(original_idx).GetPDBResidueInfo().GetName()
        assert info.GetResidueName() == ('ONE' if info.GetChainId() == 'A' else 'TWO')
        assert info.GetIsHeteroAtom()


@pytest.mark.parametrize('index', [None, -1, 999, 0.5, True])
def test_rejects_invalid_connection_indices(index):
    with pytest.raises(MergeError, match='[Aa]tom'):
        merge_structures(molecule('CC'), molecule('CO'), index, 0)


@pytest.mark.parametrize('angle', [float('nan'), float('inf'), None, 'bad'])
def test_rejects_invalid_torsion(angle):
    with pytest.raises(MergeError, match='torsion|Torsion'):
        merge_structures(molecule('CC'), molecule('CO'), 0, 0, angle)


def test_rejects_saturated_connection_atom():
    with pytest.raises(MergeError, match='valence'):
        merge_structures(molecule('CC(C)(C)C'), molecule('CC'), 1, 0)


def test_rejects_disconnected_or_missing_coordinates():
    first, second = molecule('CC'), molecule('CO')
    with pytest.raises(MergeError, match='connected|component'):
        merge_structures(Chem.CombineMols(first, first), second, 0, 0)
    first.RemoveAllConformers()
    with pytest.raises(MergeError, match='coordinate|conformer'):
        merge_structures(first, second, 0, 0)


def test_rejects_nonfinite_and_degenerate_geometry():
    first, second = molecule('CC'), molecule('CO')
    first.GetConformer().SetAtomPosition(0, (float('nan'), 0, 0))
    with pytest.raises(MergeError, match='coordinate|finite'):
        merge_structures(first, second, 0, 0)
    for atom in first.GetAtoms():
        first.GetConformer().SetAtomPosition(atom.GetIdx(), (0, 0, 0))
    with pytest.raises(MergeError, match='geometry|direction|torsion'):
        merge_structures(first, second, 0, 0)


def test_rejects_multiple_conformers_instead_of_exporting_unplaced_models():
    first, second = molecule('CC'), molecule('CO')
    first.AddConformer(Chem.Conformer(first.GetConformer()), assignId=True)
    with pytest.raises(MergeError, match='model|conformer'):
        merge_structures(first, second, 0, 0)
