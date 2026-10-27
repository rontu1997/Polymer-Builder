import base64
import tempfile
from pathlib import Path

from dash import Input, Output, State, html, no_update
from rdkit import Chem

from chem.io import StructureLoadError, load_structure
from chem.viewer import mol_to_pdb_block

MAX_FRAGMENTS = 1


def register_upload_callbacks(app):
    @app.callback(
        Output("app-state", "data"),
        Output("upload-status", "children"),
        Output("monomer-view-data", "data"),
        Input("upload-structure", "contents"),
        State("upload-structure", "filename"),
        prevent_initial_call=True,
    )
    def handle_upload(contents, filename):
        if contents is None:
            return no_update, no_update, no_update

        ext = Path(filename).suffix.lower()
        header, _, b64data = contents.partition(",")
        raw_bytes = base64.b64decode(b64data)

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_path = Path(tmpdir) / filename
                tmp_path.write_bytes(raw_bytes)
                mol = load_structure(tmp_path)
        except StructureLoadError as exc:
            return (
                {"error": str(exc)},
                html.Div(str(exc), className="error-text"),
                None,
            )

        frags = Chem.GetMolFrags(mol)
        if len(frags) > MAX_FRAGMENTS:
            message = (
                f"'{filename}' contains {len(frags)} disconnected components -- "
                f"please upload a single monomer (no solvent/counter-ions)."
            )
            return {"error": message}, html.Div(message, className="error-text"), None

        pdbblock = mol_to_pdb_block(mol)
        state = {
            "filename": filename,
            "input_format": ext,
            "monomer_pdbblock": pdbblock,
            "head_idx": None,
            "tail_idx": None,
            "error": None,
        }
        status = html.Div(f"Loaded {filename} -- {mol.GetNumAtoms()} atoms.")
        view_data = {"pdbblock": pdbblock, "head_idx": None, "tail_idx": None}
        return state, status, view_data
