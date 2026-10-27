"""Independent upload/selection state and build/download callbacks for PDB merging."""

import base64
import hashlib
import math
import uuid
from pathlib import Path

import numpy as np
from dash import Input, Output, State, ctx, html, no_update
from rdkit import Chem
from callbacks.download import export_payload

from chem.io import StructureLoadError
from chem.merge import MergeError, merge_structures
from chem.optimize import OptimizationError, optimize_geometry
from chem.validation import MAX_ATOMS_HARD_CAP, validate_linkage_atom
from chem.viewer import mol_to_pdb_block, pdb_block_to_mol

MAX_UPLOAD_BYTES = 20 * 1024 * 1024
JOIN_MODES = {
    'head-head': ('head', 'head'),
    'tail-tail': ('tail', 'tail'),
    'head-tail': ('head', 'tail'),
    'tail-head': ('tail', 'head'),
}


def _upload_key(contents, filename):
    return hashlib.sha256(f'{filename}\0{contents}'.encode()).hexdigest()


def _upload_matches(state, contents, filename):
    return bool(state and state.get('upload_key') == _upload_key(contents, filename))


def _read_upload(contents, filename):
    if not filename or Path(filename).suffix.lower() != '.pdb':
        raise StructureLoadError('Please select a .pdb file.')
    if not isinstance(contents, str) or ',' not in contents:
        raise StructureLoadError('The upload is empty or malformed. Please upload the PDB again.')
    if len(contents) > MAX_UPLOAD_BYTES * 4 // 3 + 1024:
        raise StructureLoadError('PDB uploads must be 20 MB or smaller.')
    raw = base64.b64decode(contents.partition(',')[2], validate=True)
    if not raw or len(raw) > MAX_UPLOAD_BYTES:
        raise StructureLoadError('Upload a nonempty PDB file no larger than 20 MB.')
    block = raw.decode('utf-8-sig')
    lines = block.splitlines()
    for line_number, line in enumerate(lines, 1):
        if line.startswith(('ATOM  ', 'HETATM')):
            try:
                coordinates = [float(line[start:start + 8]) for start in (30, 38, 46)]
            except ValueError as exc:
                raise StructureLoadError(f'Atom coordinates are missing or invalid on line {line_number}.') from exc
            if not all(math.isfinite(value) for value in coordinates):
                raise StructureLoadError(f'Atom coordinates must be finite on line {line_number}.')
    if sum(line.startswith('MODEL ') for line in lines) > 1:
        raise StructureLoadError('Please upload a single model per PDB file.')
    # Supplied CONECT records define the graph. Infer proximity bonds only when
    # connectivity is absent, avoiding extra bonds across folded input structures.
    has_connectivity = any(line.startswith('CONECT') for line in lines)
    mol = Chem.MolFromPDBBlock(block, removeHs=False, proximityBonding=not has_connectivity)
    if mol is None or mol.GetNumAtoms() == 0:
        raise StructureLoadError('The PDB contains no atoms.')
    if mol.GetNumConformers() != 1:
        raise StructureLoadError('Please upload a single model with 3D coordinates.')
    if len(Chem.GetMolFrags(mol)) != 1:
        raise StructureLoadError('Upload one connected structure, without solvent or counter-ions.')
    if mol.GetNumAtoms() > min(MAX_ATOMS_HARD_CAP, 99_999):
        raise StructureLoadError('This PDB exceeds the supported atom limit.')
    if not np.isfinite(mol.GetConformer().GetPositions()).all():
        raise StructureLoadError('PDB atom coordinates must be finite.')
    return {
        'filename': Path(filename).name,
        'pdbblock': mol_to_pdb_block(mol),
        'atom_count': mol.GetNumAtoms(),
        'head_idx': None,
        'tail_idx': None,
        'upload_id': uuid.uuid4().hex,
        'upload_key': _upload_key(contents, filename),
    }


def _view_data(state):
    return {
        'pdbblock': state['pdbblock'], 'head_idx': state['head_idx'],
        'tail_idx': state['tail_idx'],
        'upload_id': state['upload_id'],
    }


def _has_selection(state):
    return bool(state and state.get('pdbblock') and
                state.get('head_idx') is not None and state.get('tail_idx') is not None and
                state['head_idx'] != state['tail_idx'])


def _signature(first, second, torsion, force_field, join_mode):
    return [
        [state.get('upload_id'), state.get('head_idx'), state.get('tail_idx')] if state else [None, None, None]
        for state in (first, second)
    ] + [torsion, force_field, join_mode]


def _register_input(app, side):
    prefix = f'merge-{side}'

    @app.callback(
        Output(f'{prefix}-state', 'data'),
        Output(f'{prefix}-status', 'children'),
        Output(f'{prefix}-view-data', 'data'),
        Input(f'{prefix}-upload', 'contents'),
        Input(f'{prefix}-selected-atom', 'data'),
        Input(f'{prefix}-upload', 'filename'),
        State(f'{prefix}-state', 'data'),
        State(f'{prefix}-picking-role', 'value'),
        prevent_initial_call=True,
    )
    def update_input(contents, clicked, filename, state, picking_role):
        new_upload = contents is not None and not _upload_matches(state, contents, filename)
        if ctx.triggered_id == f'{prefix}-upload' or new_upload:
            try:
                loaded = _read_upload(contents, filename)
            except (StructureLoadError, ValueError, RuntimeError, OSError) as exc:
                return None, html.Div(f'Could not load PDB: {exc}', className='error-text'), None
            message = f"Loaded {loaded['filename']} -- {loaded['atom_count']} atoms. Select its head and tail atoms."
            return loaded, html.Div(message), _view_data(loaded)

        if not clicked or not state or clicked.get('upload_id') != state.get('upload_id'):
            return no_update, no_update, no_update
        # The viewer captures the active role with the click so a later radio
        # change cannot reclassify an in-flight selection.
        role = clicked.get('role', picking_role)
        if role not in ('head', 'tail'):
            return no_update, html.Div('Choose Head or Tail before selecting an atom.', className='error-text'), no_update
        index = clicked.get('index')
        if isinstance(index, bool) or not isinstance(index, int):
            return no_update, html.Div('Please click a valid atom.', className='error-text'), no_update
        mol = pdb_block_to_mol(state['pdbblock'])
        check = validate_linkage_atom(mol, index)
        if not check.ok:
            return no_update, html.Div(check.message, className='error-text'), no_update
        other = 'tail' if role == 'head' else 'head'
        if state.get(f'{other}_idx') == index:
            return no_update, html.Div('Head and tail must be different atoms within each polymer.', className='error-text'), no_update
        selected = dict(state, **{f'{role}_idx': index})
        def describe(end):
            idx = selected.get(f'{end}_idx')
            return 'not selected' if idx is None else f'#{idx} ({mol.GetAtomWithIdx(idx).GetSymbol()})'
        message = f"Polymer {side.upper()} — Head: {describe('head')} | Tail: {describe('tail')}"
        return selected, html.Div(message), _view_data(selected)


def register_merge_callbacks(app):
    for side in ('a', 'b'):
        _register_input(app, side)

    @app.callback(
        Output('merge-button', 'disabled'), Output('merge-atom-estimate', 'children'),
        Input('merge-a-state', 'data'), Input('merge-b-state', 'data'),
        Input('merge-join-mode', 'value'),
    )
    def update_readiness(first, second, join_mode):
        count = sum(state.get('atom_count', 0) for state in (first, second) if state)
        estimate = f'Input atoms: {count}. Missing hydrogens are added before optimization.' if count else ''
        ready = _has_selection(first) and _has_selection(second) and join_mode in JOIN_MODES
        if count > min(MAX_ATOMS_HARD_CAP, 99_999):
            return True, html.Div('The combined inputs exceed the atom limit.', className='error-text')
        return not ready, estimate

    @app.callback(
        Output('merge-result', 'data'), Output('merge-result-view-data', 'data'),
        Output('merge-stats', 'children'),
        Input('merge-button', 'n_clicks'),
        Input('merge-a-state', 'data'), Input('merge-b-state', 'data'),
        Input('merge-torsion-deg', 'value'), Input('merge-force-field', 'value'),
        Input('merge-a-upload', 'contents'), Input('merge-a-upload', 'filename'),
        Input('merge-b-upload', 'contents'), Input('merge-b-upload', 'filename'),
        Input('merge-join-mode', 'value'),
        prevent_initial_call=True,
    )
    def do_merge(n_clicks, first, second, torsion, force_field, first_contents, first_filename, second_contents, second_filename, join_mode):
        # Changing an input invalidates both the old display and its download.
        if ctx.triggered_id != 'merge-button':
            return None, None, None
        try:
            if not (_has_selection(first) and _has_selection(second)):
                raise MergeError('Upload both polymers and select a distinct head and tail atom in each.')
            if join_mode not in JOIN_MODES:
                raise MergeError('Choose which head/tail ends of polymers A and B to join.')
            if not (_upload_matches(first, first_contents, first_filename) and
                    _upload_matches(second, second_contents, second_filename)):
                raise MergeError('Wait for both PDB uploads to finish, then select the connection atoms.')
            if force_field not in ('MMFF94', 'UFF'):
                raise MergeError('Choose MMFF94 or UFF for optimization.')
            first_role, second_role = JOIN_MODES[join_mode]
            first_idx, second_idx = first[f'{first_role}_idx'], second[f'{second_role}_idx']
            built = merge_structures(
                pdb_block_to_mol(first['pdbblock']), pdb_block_to_mol(second['pdbblock']),
                first_idx, second_idx, torsion,
            )
            optimized, report = optimize_geometry(
                built.mol, built.torsion_specs, new_bonds=built.new_bonds, ff=force_field,
            )
            if not report.validation.ok:
                raise MergeError(report.validation.message)
            if not (math.isfinite(report.energy_before) and math.isfinite(report.energy_after)):
                raise MergeError('Optimization produced a nonfinite energy. Check the input coordinates.')
            block = mol_to_pdb_block(optimized)
        except (MergeError, OptimizationError, ValueError, RuntimeError) as exc:
            return None, None, html.Div(str(exc), className='error-text')

        filename = f"{Path(first['filename']).stem}_{Path(second['filename']).stem}_merged.pdb"
        signature = _signature(first, second, torsion, force_field, join_mode)
        result = {'pdbblock': block, 'filename': filename, 'signature': signature}
        lines = [
            html.Div(f'Structures: 2   Atoms: {optimized.GetNumAtoms()}   New bonds: 1'),
            html.Div(f'Connection: A {first_role} atom #{first_idx} to B {second_role} atom #{second_idx}'),
            html.Div(f'Junction torsion target: {torsion:.1f} deg'),
            html.Div(f'Force field: {report.force_field_used}' + (' (MMFF94 unavailable, fell back to UFF)' if report.fell_back_to_uff else '')),
            html.Div(f'Initial energy: {report.energy_before:.1f}   Final energy: {report.energy_after:.1f}'),
        ]
        if not report.minimize_converged:
            lines.append(html.Div('Optimization reached its iteration limit.', className='warning-text'))
        lines.extend(html.Div(warning, className='warning-text') for warning in report.validation.warnings)
        return result, {'pdbblock': block}, html.Div(lines, className='stats-box')

    @app.callback(Output('merge-save-button', 'disabled'), Input('merge-result', 'data'))
    def update_save_readiness(result):
        return not bool(result and result.get('pdbblock'))

    @app.callback(
        Output('merge-download', 'data'), Input('merge-save-request', 'data'),
        State('merge-result', 'data'), State('merge-a-state', 'data'), State('merge-b-state', 'data'),
        State('merge-torsion-deg', 'value'), State('merge-force-field', 'value'),
        State('merge-a-upload', 'contents'), State('merge-a-upload', 'filename'),
        State('merge-b-upload', 'contents'), State('merge-b-upload', 'filename'),
        State('merge-join-mode', 'value'),
        prevent_initial_call=True,
    )
    def download_merge(request, result, first, second, torsion, force_field, first_contents, first_filename, second_contents, second_filename, join_mode):
        if not request:
            return no_update
        signature = _signature(first, second, torsion, force_field, join_mode)
        if (not result or result.get('signature') != signature or
                not _upload_matches(first, first_contents, first_filename) or
                not _upload_matches(second, second_contents, second_filename)):
            return {'request_id': request.get('request_id'), 'error': 'The merger inputs changed. Merge again before saving.'}
        return export_payload(request, result['pdbblock'], result['filename'])
