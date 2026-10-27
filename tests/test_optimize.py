import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from chem.builder import build_polymer
from chem.optimize import TORSION_TOLERANCE_DEG, _get_force_field, optimize_geometry


def _ethane_polymer(n_units=4, torsion_deg=170.0):
    monomer = Chem.MolFromSmiles("CC")
    monomer = Chem.AddHs(monomer)
    AllChem.EmbedMolecule(monomer, randomSeed=3)
    AllChem.MMFFOptimizeMolecule(monomer)
    return build_polymer(monomer, head_idx=0, tail_idx=1, n_units=n_units, torsion_deg=torsion_deg)


def test_energy_does_not_increase():
    result = _ethane_polymer()
    _, report = optimize_geometry(result.mol, result.torsion_specs, new_bonds=result.new_bonds)
    assert report.energy_after <= report.energy_before + 1e-3


def test_final_torsion_stays_within_restraint_window():
    torsion_deg = 170.0
    result = _ethane_polymer(torsion_deg=torsion_deg)
    optimized_mol, report = optimize_geometry(result.mol, result.torsion_specs, new_bonds=result.new_bonds)
    conf = optimized_mol.GetConformer()
    for a1, a2, a3, a4, target in result.torsion_specs:
        actual = rdMolTransforms.GetDihedralDeg(conf, a1, a2, a3, a4)
        delta = ((actual - target + 180) % 360) - 180
        assert abs(delta) < TORSION_TOLERANCE_DEG + 2.0


def test_report_uses_mmff_when_available():
    result = _ethane_polymer()
    _, report = optimize_geometry(result.mol, result.torsion_specs, new_bonds=result.new_bonds, ff="MMFF94")
    assert report.force_field_used == "MMFF94"
    assert not report.fell_back_to_uff


def test_falls_back_to_uff_when_mmff_unavailable():
    mol = Chem.MolFromSmiles("B(F)(F)F")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    assert not AllChem.MMFFHasAllMoleculeParams(mol)
    assert AllChem.UFFHasAllMoleculeParams(mol)

    ff, name, fell_back = _get_force_field(mol, "MMFF94")
    assert name == "UFF"
    assert fell_back


def test_no_fallback_flag_when_uff_requested_directly():
    mol = Chem.MolFromSmiles("CCCC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    ff, name, fell_back = _get_force_field(mol, "UFF")
    assert name == "UFF"
    assert not fell_back


def test_single_unit_optimizes_without_torsion_specs():
    monomer = Chem.MolFromSmiles("CC")
    monomer = Chem.AddHs(monomer)
    AllChem.EmbedMolecule(monomer, randomSeed=3)
    result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=1, torsion_deg=180.0)
    optimized_mol, report = optimize_geometry(result.mol, result.torsion_specs, new_bonds=result.new_bonds)
    assert optimized_mol.GetNumAtoms() == monomer.GetNumAtoms()
