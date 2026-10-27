import tempfile
from pathlib import Path

from dash import Input, Output, State, no_update

from chem.io import StructureSaveError, save_structure
from chem.viewer import pdb_block_to_mol


EXPORT_FORMATS = {'.pdb': 'chemical/x-pdb', '.mol': 'chemical/x-mdl-molfile', '.mol2': 'chemical/x-mol2'}


def export_payload(request, pdbblock, filename):
    """Convert only the requested format and match each response to its Save dialog."""
    request_id = request.get('request_id')
    fmt = request.get('format', '.pdb')
    try:
        if fmt not in EXPORT_FORMATS:
            raise StructureSaveError('Choose PDB, MOL, or MOL2.')
        mol = pdb_block_to_mol(pdbblock) if pdbblock else None
        if mol is None:
            raise StructureSaveError('Build a valid structure before saving.')
        out_filename = Path(filename).stem + fmt
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir) / out_filename
            save_structure(mol, tmp_path, fmt)
            content = tmp_path.read_text()
        return {'request_id': request_id, 'content': content, 'filename': out_filename,
                'type': EXPORT_FORMATS[fmt]}
    except (StructureSaveError, ValueError, RuntimeError, OSError) as exc:
        return {'request_id': request_id, 'error': f'Could not export structure: {exc}'}


def register_download_callbacks(app):
    @app.callback(
        Output("download-structure", "data"),
        Input("save-request", "data"),
        State("polymer-state", "data"),
        prevent_initial_call=True,
    )
    def do_save(request, polymer_state):
        if not request:
            return no_update
        state = polymer_state or {}
        stem = Path(state.get('filename', 'structure')).stem
        return export_payload(request, state.get('polymer_pdbblock'), f'{stem}_polymer.pdb')
