import base64

import pytest
from dash import Dash, html
from rdkit import Chem
from rdkit.Chem import AllChem

from callbacks.merge import register_merge_callbacks
from chem.viewer import mol_to_pdb_block, pdb_block_to_mol


@pytest.fixture
def dash_app():
    app = Dash(__name__, suppress_callback_exceptions=True)
    app.layout = html.Div()
    register_merge_callbacks(app)
    return app


def request_callback(app, output_id, values, changed):
    key, entry = next((key, entry) for key, entry in app.callback_map.items() if output_id in key)
    outputs = entry['output']
    output_spec = ([{'id': out.component_id, 'property': out.component_property} for out in outputs]
                   if isinstance(outputs, list) else {'id': outputs.component_id, 'property': outputs.component_property})
    response = app.server.test_client().post('/_dash-update-component', json={
        'output': key, 'outputs': output_spec,
        'inputs': [dict(item, value=values.get(item['id'] + '.' + item['property'])) for item in entry['inputs']],
        'state': [dict(item, value=values.get(item['id'] + '.' + item['property'])) for item in entry['state']],
        'changedPropIds': [changed],
    })
    assert response.status_code == 200, response.get_data(as_text=True)
    return response.get_json()['response']


def upload_text(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=7) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    return 'data:chemical/x-pdb;base64,' + base64.b64encode(mol_to_pdb_block(mol).encode()).decode()


def upload(app, side, contents, filename='input.pdb', state=None):
    prefix = f'merge-{side}'
    return request_callback(app, f'{prefix}-state.data', {
        f'{prefix}-upload.contents': contents, f'{prefix}-upload.filename': filename,
        f'{prefix}-state.data': state,
    }, f'{prefix}-upload.contents')


def select(app, side, state, index, upload_id=None, role='head'):
    prefix = f'merge-{side}'
    return request_callback(app, f'{prefix}-state.data', {
        f'{prefix}-selected-atom.data': {'index': index, 'upload_id': upload_id or state['upload_id'], 'role': role},
        f'{prefix}-state.data': state,
        f'{prefix}-picking-role.value': role,
    }, f'{prefix}-selected-atom.data')


def merge_inputs(app):
    first_contents, second_contents = upload_text('CC'), upload_text('CO')
    first = upload(app, 'a', first_contents, 'first.pdb')['merge-a-state']['data']
    second = upload(app, 'b', second_contents, 'second.pdb')['merge-b-state']['data']
    first = select(app, 'a', first, 0)['merge-a-state']['data']
    first = select(app, 'a', first, 1, role='tail')['merge-a-state']['data']
    second = select(app, 'b', second, 0)['merge-b-state']['data']
    second = select(app, 'b', second, 1, role='tail')['merge-b-state']['data']
    return {
        'merge-a-state.data': first, 'merge-b-state.data': second,
        'merge-torsion-deg.value': 60, 'merge-force-field.value': 'MMFF94',
        'merge-button.n_clicks': 1,
        'merge-save-request.data': {'request_id': 'merger-save-test', 'format': '.pdb'},
        'merge-join-mode.value': 'tail-head',
        'merge-a-upload.contents': first_contents, 'merge-a-upload.filename': 'first.pdb',
        'merge-b-upload.contents': second_contents, 'merge-b-upload.filename': 'second.pdb',
    }


def test_upload_select_merge_and_download_connected_pdb(dash_app):
    inputs = merge_inputs(dash_app)
    readiness = request_callback(dash_app, 'merge-button.disabled', inputs, 'merge-a-state.data')
    assert readiness['merge-button']['disabled'] is False
    merged = request_callback(dash_app, 'merge-result.data', inputs, 'merge-button.n_clicks')
    result = merged['merge-result']['data']
    mol = pdb_block_to_mol(result['pdbblock'])
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == 'CCCO'
    assert mol.GetNumAtoms() == 12
    assert len(Chem.GetMolFrags(mol)) == 1
    download = request_callback(dash_app, 'merge-download.data', {
        **inputs, 'merge-result.data': result, 'merge-save-button.n_clicks': 1,
    }, 'merge-save-button.n_clicks')['merge-download']['data']
    assert download['filename'] == 'first_second_merged.pdb'
    restored = pdb_block_to_mol(download['content'])
    assert restored.GetNumBonds() == mol.GetNumBonds()


@pytest.mark.parametrize('changed', ['merge-a-state.data', 'merge-b-state.data', 'merge-torsion-deg.value', 'merge-force-field.value', 'merge-a-upload.contents', 'merge-b-upload.filename', 'merge-join-mode.value'])
def test_input_changes_clear_previous_result_and_preview(dash_app, changed):
    inputs = merge_inputs(dash_app)
    result = request_callback(dash_app, 'merge-result.data', inputs, changed)
    assert result['merge-result']['data'] is None
    assert result['merge-result-view-data']['data'] is None
    assert result['merge-stats']['children'] is None


def test_upload_replacement_resets_selection_and_ignores_old_clicks(dash_app):
    old = upload(dash_app, 'a', upload_text('CC'))['merge-a-state']['data']
    old = select(dash_app, 'a', old, 0)['merge-a-state']['data']
    new = upload(dash_app, 'a', upload_text('CO'), state=old)['merge-a-state']['data']
    assert new['head_idx'] is None and new['tail_idx'] is None
    assert new['upload_id'] != old['upload_id']
    ignored = select(dash_app, 'a', new, 0, upload_id=old['upload_id'])
    assert 'merge-a-state' not in ignored or ignored['merge-a-state']['data']['head_idx'] is None


@pytest.mark.parametrize('contents,filename', [
    ('not base64', 'bad.pdb'), ('data:;base64,!!!!', 'bad.pdb'),
    ('data:;base64,' + base64.b64encode(b'not a pdb').decode(), 'bad.pdb'),
    (None, 'missing.pdb'),
])
def test_bad_upload_clears_input_and_view(dash_app, contents, filename):
    old = upload(dash_app, 'a', upload_text('CC'))['merge-a-state']['data']
    result = upload(dash_app, 'a', contents, filename, state=old)
    assert result['merge-a-state']['data'] is None
    assert result['merge-a-view-data']['data'] is None


def test_non_pdb_and_disconnected_uploads_are_rejected(dash_app):
    for contents, filename in [(upload_text('CC'), 'fake.mol'), (upload_text('CC.CC'), 'fragments.pdb')]:
        result = upload(dash_app, 'a', contents, filename)
        assert result['merge-a-state']['data'] is None


def test_invalid_selection_does_not_enable_merge(dash_app):
    first = upload(dash_app, 'a', upload_text('CC(C)(C)C'))['merge-a-state']['data']
    rejected = select(dash_app, 'a', first, 1)
    assert 'merge-a-state' not in rejected or rejected['merge-a-state']['data']['head_idx'] is None
    result = request_callback(dash_app, 'merge-button.disabled', {'merge-a-state.data': first}, 'merge-a-state.data')
    assert result['merge-button']['disabled'] is True


def test_stale_output_cannot_be_downloaded(dash_app):
    inputs = merge_inputs(dash_app)
    merged = request_callback(dash_app, 'merge-result.data', inputs, 'merge-button.n_clicks')
    response = request_callback(dash_app, 'merge-download.data', {
        **inputs, 'merge-result.data': merged['merge-result']['data'],
        'merge-torsion-deg.value': -30, 'merge-save-button.n_clicks': 1,
    }, 'merge-save-button.n_clicks')
    assert response['merge-download']['data']['error']
    assert 'content' not in response['merge-download']['data']


def test_old_viewer_click_during_upload_cannot_restore_old_structure(dash_app):
    old = upload(dash_app, 'a', upload_text('CC'), 'old.pdb')['merge-a-state']['data']
    result = request_callback(dash_app, 'merge-a-state.data', {
        'merge-a-upload.contents': upload_text('CO'), 'merge-a-upload.filename': 'new.pdb',
        'merge-a-state.data': old,
        'merge-a-selected-atom.data': {'index': 0, 'upload_id': old['upload_id']},
    }, 'merge-a-selected-atom.data')
    new = result['merge-a-state']['data']
    assert new['filename'] == 'new.pdb'
    assert new['head_idx'] is None and new['tail_idx'] is None
    assert Chem.MolToSmiles(Chem.RemoveHs(pdb_block_to_mol(new['pdbblock']))) == 'CO'


def test_missing_pdb_coordinates_are_rejected_before_parser_substitutes_zeros(dash_app):
    content = upload_text('CC')
    block = base64.b64decode(content.partition(',')[2]).decode()
    lines = block.splitlines()
    lines[0] = lines[0][:30] + ' ' * 24 + lines[0][54:]
    broken = 'data:;base64,' + base64.b64encode('\n'.join(lines).encode()).decode()
    result = upload(dash_app, 'a', broken)
    assert result['merge-a-state']['data'] is None


def test_multiple_model_pdb_upload_is_rejected(dash_app):
    mol = Chem.AddHs(Chem.MolFromSmiles('CC'))
    AllChem.EmbedMultipleConfs(mol, numConfs=2, randomSeed=7)
    content = 'data:;base64,' + base64.b64encode(mol_to_pdb_block(mol).encode()).decode()
    result = upload(dash_app, 'a', content)
    assert result['merge-a-state']['data'] is None


def test_explicit_pdb_connectivity_is_not_overridden_by_proximity(dash_app):
    mol = Chem.MolFromSmiles('CCC')
    conf = Chem.Conformer(3)
    for index, xyz in enumerate(((0, 0, 0), (1.5, 0, 0), (0.75, 1.3, 0))):
        conf.SetAtomPosition(index, xyz)
    mol.AddConformer(conf)
    contents = 'data:;base64,' + base64.b64encode(mol_to_pdb_block(mol).encode()).decode()
    state = upload(dash_app, 'a', contents)['merge-a-state']['data']
    assert Chem.MolToSmiles(pdb_block_to_mol(state['pdbblock'])) == 'CCC'


def test_download_and_merge_blocked_while_new_upload_is_being_parsed(dash_app):
    inputs = merge_inputs(dash_app)
    old = request_callback(dash_app, 'merge-result.data', inputs, 'merge-button.n_clicks')['merge-result']['data']
    values = {**inputs, 'merge-a-upload.contents': upload_text('CCC'), 'merge-result.data': old}
    download = request_callback(dash_app, 'merge-download.data', values, 'merge-save-button.n_clicks')
    assert download['merge-download']['data']['error']
    assert 'content' not in download['merge-download']['data']
    merge = request_callback(dash_app, 'merge-result.data', values, 'merge-button.n_clicks')
    assert merge['merge-result']['data'] is None


def test_both_endpoints_are_retained_and_highlighted(dash_app):
    state = upload(dash_app, 'a', upload_text('CC'))['merge-a-state']['data']
    state = select(dash_app, 'a', state, 0)['merge-a-state']['data']
    result = select(dash_app, 'a', state, 1, role='tail')
    assert result['merge-a-state']['data']['head_idx'] == 0
    assert result['merge-a-state']['data']['tail_idx'] == 1
    view = result['merge-a-view-data']['data']
    assert (view['head_idx'], view['tail_idx']) == (0, 1)


def test_head_and_tail_must_be_different_atoms(dash_app):
    state = upload(dash_app, 'a', upload_text('CC'))['merge-a-state']['data']
    state = select(dash_app, 'a', state, 0)['merge-a-state']['data']
    result = select(dash_app, 'a', state, 0, role='tail')
    assert 'merge-a-state' not in result
    assert 'different' in str(result['merge-a-status'])


@pytest.mark.parametrize('mode,expected', [
    ('head-head', 'NCCCO'), ('tail-tail', 'CCONC'),
    ('head-tail', 'CNCCO'), ('tail-head', 'CCOCN'),
])
def test_join_mode_uses_the_requested_a_and_b_endpoints(dash_app, mode, expected):
    inputs = merge_inputs(dash_app)
    for side, smiles, tail in [('a', 'CCO', 2), ('b', 'CN', 1)]:
        contents = upload_text(smiles)
        filename = f'{side}.pdb'
        state = upload(dash_app, side, contents, filename)[f'merge-{side}-state']['data']
        state = select(dash_app, side, state, 0)[f'merge-{side}-state']['data']
        state = select(dash_app, side, state, tail, role='tail')[f'merge-{side}-state']['data']
        inputs.update({f'merge-{side}-state.data': state, f'merge-{side}-upload.contents': contents,
                       f'merge-{side}-upload.filename': filename})
    inputs['merge-join-mode.value'] = mode
    result = request_callback(dash_app, 'merge-result.data', inputs, 'merge-button.n_clicks')['merge-result']['data']
    mol = pdb_block_to_mol(result['pdbblock'])
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == Chem.MolToSmiles(Chem.MolFromSmiles(expected))


def test_missing_tail_disables_build_and_rejects_merge(dash_app):
    inputs = merge_inputs(dash_app)
    inputs['merge-b-state.data']['tail_idx'] = None
    readiness = request_callback(dash_app, 'merge-button.disabled', inputs, 'merge-b-state.data')
    assert readiness['merge-button']['disabled'] is True
    result = request_callback(dash_app, 'merge-result.data', inputs, 'merge-button.n_clicks')
    assert result['merge-result']['data'] is None


def test_join_mode_change_blocks_stale_download(dash_app):
    inputs = merge_inputs(dash_app)
    result = request_callback(dash_app, 'merge-result.data', inputs, 'merge-button.n_clicks')['merge-result']['data']
    response = request_callback(dash_app, 'merge-download.data', {
        **inputs, 'merge-result.data': result, 'merge-join-mode.value': 'head-head',
    }, 'merge-save-button.n_clicks')
    assert response['merge-download']['data']['error']
    assert 'content' not in response['merge-download']['data']
