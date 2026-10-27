import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rdkit import Chem
from rdkit.Chem import AllChem

from chem.io import load_structure, save_structure


def _make_test_monomer():
    """Ethylene glycol-like small monomer with an embedded 3D conformer, used
    purely to exercise the file I/O round trip (not a claim about real polymer
    chemistry)."""
    mol = Chem.MolFromSmiles("OCCO")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def test_roundtrip_mol(tmp_path):
    mol = _make_test_monomer()
    path = tmp_path / "monomer.mol"
    save_structure(mol, path, ".mol")
    reloaded = load_structure(path)
    assert reloaded.GetNumAtoms() == mol.GetNumAtoms()


def test_roundtrip_pdb(tmp_path):
    mol = _make_test_monomer()
    path = tmp_path / "monomer.pdb"
    save_structure(mol, path, ".pdb")
    reloaded = load_structure(path)
    assert reloaded.GetNumAtoms() == mol.GetNumAtoms()


def test_roundtrip_mol2(tmp_path):
    mol = _make_test_monomer()
    path = tmp_path / "monomer.mol2"
    save_structure(mol, path, ".mol2")
    assert path.exists()
    text = path.read_text()
    assert "@<TRIPOS>ATOM" in text
    reloaded = load_structure(path)
    assert reloaded.GetNumAtoms() == mol.GetNumAtoms()


def test_load_structure_unsupported_extension(tmp_path):
    path = tmp_path / "monomer.xyz"
    path.write_text("not a real file")
    try:
        load_structure(path)
        assert False, "expected StructureLoadError"
    except Exception as exc:
        assert "Unsupported" in str(exc)


def test_pdb_conect_records_are_symmetric(tmp_path):
    # RDKit's PDB writer lists each bond only once, on the lower-numbered atom's
    # CONECT line -- save_structure must post-process so every atom's own record
    # lists all of its neighbors, including bonds where it's the higher index.
    mol = _make_test_monomer()
    path = tmp_path / "monomer.pdb"
    save_structure(mol, path, ".pdb")

    adjacency = {}
    for line in path.read_text().splitlines():
        if not line.startswith("CONECT"):
            continue
        serial = int(line[6:11])
        partners = [int(line[i : i + 5]) for i in range(11, len(line), 5) if line[i : i + 5].strip()]
        adjacency.setdefault(serial, []).extend(partners)

    for serial, partners in adjacency.items():
        for partner in partners:
            assert serial in adjacency.get(partner, []), (
                f"atom {partner} is listed as bonded to {serial}, but atom {partner}'s "
                f"own CONECT record does not list {serial} back"
            )


def test_load_structure_embeds_conformer_when_missing(tmp_path):
    # A MOL block with all-zero coordinates (2D-only) should still load and get a
    # real embedded 3D conformer via load_structure's fallback.
    mol_2d = Chem.MolFromSmiles("OCCO")
    molblock = Chem.MolToMolBlock(mol_2d, kekulize=True)
    path = tmp_path / "flat.mol"
    path.write_text(molblock)
    reloaded = load_structure(path)
    assert reloaded.GetNumConformers() == 1
