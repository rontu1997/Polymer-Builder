"""Load and save monomer/polymer structures across PDB, MOL, and MOL2.

RDKit handles PDB and MOL natively in both directions. RDKit can *read* MOL2 but
has no MOL2 writer, so OpenBabel bridges the MOL2 boundary in both directions
(through an RDKit MOL block, so the rest of the pipeline only ever touches RDKit
Mol objects).
"""

from pathlib import Path

from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem

from chem.viewer import PDB_WRITE_FLAVOR

SUPPORTED_EXTENSIONS = {".pdb", ".mol", ".sdf", ".mol2"}


class StructureLoadError(Exception):
    pass


class StructureSaveError(Exception):
    pass


def _mol2_file_to_molblock(path: Path) -> str:
    obmol = next(pybel.readfile("mol2", str(path)))
    return obmol.write("mol")


def _molblock_to_mol2_string(molblock: str) -> str:
    obmol = pybel.readstring("mol", molblock)
    return obmol.write("mol2")


def _ensure_conformer(mol: Chem.Mol) -> Chem.Mol:
    if mol.GetNumConformers() > 0:
        return mol
    mol = Chem.AddHs(mol)
    embed_result = AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
    if embed_result != 0:
        raise StructureLoadError(
            "Could not generate a 3D conformer for this structure; the input "
            "file has no coordinates and RDKit's embedding failed."
        )
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def load_structure(path) -> Chem.Mol:
    """Load a PDB/MOL/SDF/MOL2 file into an RDKit Mol with a 3D conformer."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext not in SUPPORTED_EXTENSIONS:
        raise StructureLoadError(
            f"Unsupported file extension '{ext}'. Supported: "
            f"{', '.join(sorted(SUPPORTED_EXTENSIONS))}."
        )

    if ext == ".pdb":
        mol = Chem.MolFromPDBFile(str(path), removeHs=False)
    elif ext in (".mol", ".sdf"):
        mol = Chem.MolFromMolFile(str(path), removeHs=False)
    else:  # .mol2
        try:
            molblock = _mol2_file_to_molblock(path)
        except StopIteration:
            raise StructureLoadError(f"OpenBabel found no molecule in '{path}'.")
        mol = Chem.MolFromMolBlock(molblock, removeHs=False)

    if mol is None:
        raise StructureLoadError(
            f"RDKit could not parse '{path}' as a valid {ext} structure."
        )

    mol = _ensure_conformer(mol)
    return mol


def save_structure(mol: Chem.Mol, path, fmt: str) -> None:
    """Save an RDKit Mol to `path`, in the format given by `fmt` (a file extension,
    e.g. '.pdb', '.mol', '.mol2'). The selected output format may differ from
    the original upload's extension.
    """
    path = Path(path)
    fmt = fmt.lower()
    if fmt not in SUPPORTED_EXTENSIONS:
        raise StructureSaveError(
            f"Unsupported save format '{fmt}'. Supported: "
            f"{', '.join(sorted(SUPPORTED_EXTENSIONS))}."
        )

    if fmt == ".pdb":
        # flavor & 4: write CONECT records in both directions. RDKit's default
        # only lists a bond once, on the lower-numbered atom's own CONECT line --
        # a consumer that looks up a specific atom's connectivity (e.g. the head
        # atom of a later polymer unit) would see it as unbonded.
        Chem.MolToPDBFile(mol, str(path), flavor=PDB_WRITE_FLAVOR)
    elif fmt in (".mol", ".sdf"):
        Chem.MolToMolFile(mol, str(path))
    else:  # .mol2
        molblock = Chem.MolToMolBlock(mol)
        mol2_text = _molblock_to_mol2_string(molblock)
        path.write_text(mol2_text)
        _verify_mol2_roundtrip(path, expected_atoms=mol.GetNumAtoms())


def _verify_mol2_roundtrip(path: Path, expected_atoms: int) -> None:
    """Re-read a just-written .mol2 file and confirm the atom count matches, so a
    silent OpenBabel conversion failure surfaces as an error rather than a corrupt
    file left on disk.
    """
    try:
        reread = next(pybel.readfile("mol2", str(path)))
    except StopIteration:
        raise StructureSaveError(
            f"Wrote '{path}' but could not re-read it back as MOL2 to verify."
        )
    if reread.OBMol.NumAtoms() != expected_atoms:
        raise StructureSaveError(
            f"MOL2 round-trip check failed for '{path}': wrote "
            f"{expected_atoms} atoms but re-read {reread.OBMol.NumAtoms()}."
        )
