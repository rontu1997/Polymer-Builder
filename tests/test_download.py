import base64

import pytest
from dash import Dash, html
from rdkit import Chem

from callbacks.download import register_download_callbacks
from callbacks.merge import register_merge_callbacks
from chem.io import load_structure
from test_merge_callbacks import merge_inputs, request_callback, upload_text


@pytest.fixture
def app():
    app = Dash(__name__, suppress_callback_exceptions=True)
    app.layout = html.Div()
    register_download_callbacks(app)
    register_merge_callbacks(app)
    return app


def polymer_state():
    return {
        'polymer_pdbblock': base64.b64decode(upload_text('CCO').partition(',')[2]).decode(),
        'input_format': '.mol2', 'filename': 'polymer-input.mol2',
    }


@pytest.mark.parametrize('fmt', ['.pdb', '.mol', '.mol2'])
@pytest.mark.parametrize('source', ['homopolymer', 'merger'])
def test_selected_export_format_has_correct_extension_and_molecular_graph(app, tmp_path, fmt, source):
    request = {'request_id': 'test-save', 'format': fmt}
    if source == 'homopolymer':
        output = 'download-structure'
        values = {'save-request.data': request, 'polymer-state.data': polymer_state()}
        trigger = 'save-request.data'
        expected, filename = 'CCO', f'polymer-input_polymer{fmt}'
    else:
        output = 'merge-download'
        values = merge_inputs(app)
        result = request_callback(app, 'merge-result.data', values, 'merge-button.n_clicks')['merge-result']['data']
        values.update({'merge-result.data': result, 'merge-save-request.data': request})
        trigger = 'merge-save-request.data'
        expected, filename = 'CCCO', f'first_second_merged{fmt}'
    payload = request_callback(app, f'{output}.data', values, trigger)[output]['data']
    assert payload['request_id'] == 'test-save'
    assert payload['filename'] == filename
    path = tmp_path / filename
    path.write_text(payload['content'])
    mol = load_structure(path)
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == expected


def test_save_defaults_to_pdb_even_when_input_was_mol2(app):
    result = request_callback(app, 'download-structure.data', {
        'save-request.data': {'request_id': 'default'}, 'polymer-state.data': polymer_state(),
    }, 'save-request.data')['download-structure']['data']
    assert result['filename'] == 'polymer-input_polymer.pdb'
    assert 'HETATM' in result['content']


@pytest.mark.parametrize('state,fmt', [(None, '.pdb'), (polymer_state(), '.xyz')])
def test_export_failure_returns_error_instead_of_a_file(app, state, fmt):
    result = request_callback(app, 'download-structure.data', {
        'save-request.data': {'request_id': 'failed', 'format': fmt},
        'polymer-state.data': state,
    }, 'save-request.data')['download-structure']['data']
    assert result['request_id'] == 'failed'
    assert result['error']
    assert 'content' not in result
