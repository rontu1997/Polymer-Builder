"""Convert an RDKit Mol into the one representation the 3D viewer (and the app's
internal dcc.Store round-trips) need.

A PDB block, not a MOL block, is used as the universal internal serialization:
3Dmol.js's `addModel(block, 'pdb')` preserves file atom order, which is exactly
RDKit's 0-based atom index order -- so a click event's `atom.index` in the browser
maps directly to an RDKit atom index in Python, with no lookup table required.
Critically, unlike a MOL block, a PDB block also round-trips chain/residue/atom-
name metadata -- which matters because build_polymer() uses that metadata (via
chem.builder's residue-renumbering) to keep polymer units distinguishable and
HETATM-tagged for viewer compatibility. A MOL-block round trip would silently
discard that metadata before the builder ever saw it.
"""

from rdkit import Chem

PDB_WRITE_FLAVOR = 4  # write CONECT records in both directions


def mol_to_pdb_block(mol: Chem.Mol) -> str:
    return Chem.MolToPDBBlock(mol, flavor=PDB_WRITE_FLAVOR)


def pdb_block_to_mol(pdb_block: str) -> Chem.Mol:
    return Chem.MolFromPDBBlock(pdb_block, sanitize=True, removeHs=False, proximityBonding=False)
