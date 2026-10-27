import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

from chem.builder import PolymerBuildError, build_polymer
from chem.io import load_structure

_ALA_PDB_TEXT = """\
ATOM      1  N   ALA A   1      -0.966   0.493   1.500  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.257  -0.238   1.130  1.00  0.00           C
ATOM      3  C   ALA A   1       0.161  -0.895  -0.238  1.00  0.00           C
ATOM      4  O   ALA A   1      -0.937  -1.101  -0.757  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.453   0.716   1.184  1.00  0.00           C
ATOM      6  H   ALA A   1      -1.000   1.400   1.000  1.00  0.00           H
ATOM      7  HA  ALA A   1       0.400  -1.000   1.900  1.00  0.00           H
ATOM      8  HB1 ALA A   1       1.400   1.400   2.000  1.00  0.00           H
ATOM      9  HB2 ALA A   1       2.400   0.200   1.300  1.00  0.00           H
ATOM     10  HB3 ALA A   1       1.400   1.300   0.200  1.00  0.00           H
ATOM     11  OXT ALA A   1       1.200  -1.300  -0.800  1.00  0.00           O
END
"""


def _alanine_like_monomer(tmp_path):
    """A monomer with real PDB residue info (chain/resname/resnum), the kind an
    actual protein-style upload would carry -- used to exercise the
    residue-renumbering / hetero-atom fix, not as a claim about real peptide
    chemistry (head/tail picked below are just structurally valid, not a real
    peptide bond)."""
    path = tmp_path / "ala.pdb"
    path.write_text(_ALA_PDB_TEXT)
    return load_structure(path)


def test_each_unit_gets_a_distinct_residue_number_and_is_marked_hetero(tmp_path):
    monomer = _alanine_like_monomer(tmp_path)
    result = build_polymer(monomer, head_idx=0, tail_idx=4, n_units=3, torsion_deg=180.0)

    seen_resnums = set()
    for atom in result.mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        assert info is not None
        assert info.GetIsHeteroAtom(), "every atom should be forced to HETATM for viewer compatibility"
        seen_resnums.add(info.GetResidueNumber())

    assert seen_resnums == {1, 2, 3}, f"expected 3 distinct residue numbers, got {seen_resnums}"


def _ethane():
    """A pure geometry/algorithm test case -- not a claim that ethane is a
    realistic polymer. Isolates the placement algorithm from real-monomer
    chemistry complexity, per the plan's build-order notes."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def test_2mer_atom_and_bond_counts():
    monomer = _ethane()
    result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=2, torsion_deg=180.0)
    # each of the 1 junction removes one H from each side
    assert result.mol.GetNumAtoms() == monomer.GetNumAtoms() * 2 - 2
    assert len(result.new_bonds) == 1
    assert len(result.torsion_specs) == 1


def test_5mer_atom_and_bond_counts():
    monomer = _ethane()
    result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=5, torsion_deg=180.0)
    assert result.mol.GetNumAtoms() == monomer.GetNumAtoms() * 5 - 2 * 4
    assert len(result.new_bonds) == 4
    assert len(result.torsion_specs) == 4


def test_no_global_atomic_overlap_at_scale():
    """Regression test for a real bug: a tail unit's leaving-H lookup used a
    template-relative index instead of one offset for the unit that currently
    owns the growing chain's tail atom, silently corrupting every placement past
    the first junction. Per-bond length and per-junction dihedral checks (the
    other tests here) both stayed green through that bug -- the later "twist"
    step forces the dihedral to match target regardless of the placement
    direction, and each new bond's length is correct by construction regardless
    of which direction it points. Only a *global* non-bonded distance check
    across a long-enough chain exposes it (the chain folded back and
    self-intersected once eight-plus units had accumulated)."""
    monomer = _ethane()
    result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=30, torsion_deg=180.0)
    conf = result.mol.GetConformer()
    n = result.mol.GetNumAtoms()
    positions = np.array([list(conf.GetAtomPosition(i)) for i in range(n)])

    bonded = set()
    for bond in result.mol.GetBonds():
        bonded.add(frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))

    diffs = positions[:, None, :] - positions[None, :, :]
    dist_matrix = np.sqrt((diffs**2).sum(axis=-1))
    np.fill_diagonal(dist_matrix, np.inf)

    min_dist = np.inf
    for i in range(n):
        for j in range(i + 1, n):
            if frozenset((i, j)) in bonded:
                continue
            min_dist = min(min_dist, dist_matrix[i, j])

    assert min_dist > 0.7, f"chain has severely overlapping non-bonded atoms (closest: {min_dist:.3f} A)"


def test_all_atoms_have_valid_valence():
    monomer = _ethane()
    result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=5, torsion_deg=180.0)
    pt = Chem.GetPeriodicTable()
    for atom in result.mol.GetAtoms():
        max_v = pt.GetDefaultValence(atom.GetAtomicNum()) + atom.GetFormalCharge()
        assert atom.GetTotalValence() <= max_v, (
            f"atom {atom.GetIdx()} ({atom.GetSymbol()}) has valence "
            f"{atom.GetTotalValence()} > max {max_v}"
        )


def test_new_bonds_have_reasonable_length():
    monomer = _ethane()
    result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=3, torsion_deg=180.0)
    conf = result.mol.GetConformer()
    for a1, a2 in result.new_bonds:
        dist = conf.GetAtomPosition(a1).Distance(conf.GetAtomPosition(a2))
        assert 1.2 < dist < 1.8, f"bond {a1}-{a2} length {dist:.2f} looks wrong for a C-C bond"


def test_constructed_dihedral_matches_target_before_optimization():
    monomer = _ethane()
    target = 150.0
    result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=3, torsion_deg=target)
    conf = result.mol.GetConformer()
    for a1, a2, a3, a4, spec_target in result.torsion_specs:
        assert spec_target == target
        actual = rdMolTransforms.GetDihedralDeg(conf, a1, a2, a3, a4)
        delta = ((actual - target + 180) % 360) - 180
        assert abs(delta) < 1.0, f"dihedral {a1}-{a2}-{a3}-{a4} is {actual:.2f}, expected {target}"


def test_single_unit_returns_monomer_unchanged():
    monomer = _ethane()
    result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=1, torsion_deg=180.0)
    assert result.mol.GetNumAtoms() == monomer.GetNumAtoms()
    assert result.new_bonds == []
    assert result.torsion_specs == []


def test_same_head_and_tail_rejected():
    monomer = _ethane()
    try:
        build_polymer(monomer, head_idx=0, tail_idx=0, n_units=2, torsion_deg=180.0)
        assert False, "expected PolymerBuildError"
    except PolymerBuildError:
        pass


def test_negative_torsion_normalizes_correctly():
    monomer = _ethane()
    result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=2, torsion_deg=-170.0)
    conf = result.mol.GetConformer()
    a1, a2, a3, a4, spec_target = result.torsion_specs[0]
    assert spec_target == -170.0
    actual = rdMolTransforms.GetDihedralDeg(conf, a1, a2, a3, a4)
    delta = ((actual - spec_target + 180) % 360) - 180
    assert abs(delta) < 1.0


def test_per_junction_torsion_list_matches_each_target():
    monomer = _ethane()
    targets = [150.0, -60.0, 45.0, 170.0]
    result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=5, torsion_deg=targets)
    conf = result.mol.GetConformer()
    assert [t for *_, t in result.torsion_specs] == targets
    for (a1, a2, a3, a4, target), expected in zip(result.torsion_specs, targets):
        actual = rdMolTransforms.GetDihedralDeg(conf, a1, a2, a3, a4)
        delta = ((actual - expected + 180) % 360) - 180
        assert abs(delta) < 1.0, f"junction targeting {expected} actually measured {actual:.2f}"


def test_wrong_length_torsion_list_is_rejected():
    monomer = _ethane()
    try:
        build_polymer(monomer, head_idx=0, tail_idx=1, n_units=5, torsion_deg=[180.0, 180.0])
        assert False, "expected PolymerBuildError"
    except PolymerBuildError:
        pass


def test_uniform_torsion_list_matches_equivalent_scalar():
    """A per-junction list of all-the-same-value must be a strict generalization of
    the scalar path, not a parallel code path that could silently diverge -- same
    torsion targets and same conformer geometry either way."""
    monomer = _ethane()
    n_units = 4
    scalar_result = build_polymer(monomer, head_idx=0, tail_idx=1, n_units=n_units, torsion_deg=180.0)
    list_result = build_polymer(
        monomer, head_idx=0, tail_idx=1, n_units=n_units, torsion_deg=[180.0] * (n_units - 1)
    )

    assert [t for *_, t in scalar_result.torsion_specs] == [t for *_, t in list_result.torsion_specs]

    scalar_conf = scalar_result.mol.GetConformer()
    list_conf = list_result.mol.GetConformer()
    assert scalar_result.mol.GetNumAtoms() == list_result.mol.GetNumAtoms()
    for i in range(scalar_result.mol.GetNumAtoms()):
        p1, p2 = scalar_conf.GetAtomPosition(i), list_conf.GetAtomPosition(i)
        assert p1.Distance(p2) < 1e-6, f"atom {i} differs between the scalar and uniform-list build paths"
