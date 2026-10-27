from dash import dcc, html

from callbacks.merge import MAX_UPLOAD_BYTES
from ui.save import save_controls


def _input_panel(side):
    prefix = f'merge-{side}'
    return html.Div([
        html.H4(f'Polymer {side.upper()} (linker {side.upper()})'),
        dcc.Upload(
            id=f'{prefix}-upload', accept='.pdb', multiple=False,
            max_size=MAX_UPLOAD_BYTES,
            children=html.Div(['Drag and drop, or ', html.A('select a .pdb file')]),
            className='structure-upload',
        ),
        html.Div([
            html.Label('Picking:'),
            dcc.RadioItems(
                id=f'{prefix}-picking-role',
                options=[{'label': ' Head (blue)', 'value': 'head'},
                         {'label': ' Tail (red)', 'value': 'tail'}],
                value='head', inline=True,
            ),
        ], className='controls-row'),
        html.P('Select Head or Tail, then click its atom. Left-drag rotates; right-drag moves the view; scroll zooms.'),
        html.Div(id=f'{prefix}-status', role='status'),
        html.Div(id=f'{prefix}-viewer', className='viewer-box'),
        html.Div(id=f'{prefix}-viewer-dummy', style={'display': 'none'}),
        dcc.Store(id=f'{prefix}-state'),
        dcc.Store(id=f'{prefix}-view-data'),
        dcc.Store(id=f'{prefix}-selected-atom'),
    ], className='panel')


def merger_layout():
    return html.Div([
        html.P(
            'Upload polymers A and B, select the head and tail atoms of each, and choose which ends to join. '
            'The merger joins the chosen ends with a single bond, then optimizes the combined structure. '
            'Missing hydrogens are added; one hydrogen at each connection site is replaced when needed.'
        ),
        html.Div([_input_panel('a'), _input_panel('b')], className='panels-row'),
        html.Div([
            html.Label('Join ends:', htmlFor='merge-join-mode'),
            dcc.Dropdown(
                id='merge-join-mode',
                options=[
                    {'label': 'Head–Head: A head → B head', 'value': 'head-head'},
                    {'label': 'Tail–Tail: A tail → B tail', 'value': 'tail-tail'},
                    {'label': 'Head–Tail: A head → B tail', 'value': 'head-tail'},
                    {'label': 'Tail–Head: A tail → B head', 'value': 'tail-head'},
                ],
                value='tail-head', clearable=False, style={'width': '340px'},
            ),
        ], className='controls-row'),
        html.Div([
            html.Label('Junction torsion (deg):', htmlFor='merge-torsion-deg'),
            html.Div(dcc.Slider(
                id='merge-torsion-deg', min=-180, max=180, step=1, value=180,
                marks={-180: '-180', 0: '0', 180: '180'},
                tooltip={'placement': 'bottom'},
            ), className='merger-torsion'),
            html.Label('Force field:', htmlFor='merge-force-field'),
            dcc.Dropdown(
                id='merge-force-field', options=['MMFF94', 'UFF'], value='MMFF94',
                clearable=False, style={'width': '140px'},
            ),
        ], className='controls-row'),
        html.Div(id='merge-atom-estimate'),
        html.Button('Merge & Optimize', id='merge-button', disabled=True),
        html.H4('Merged structure'),
        dcc.Loading([
            html.Div(id='merge-result-viewer', className='viewer-box'),
            html.Div(id='merge-result-viewer-dummy', style={'display': 'none'}),
            html.Div(id='merge-stats', role='status'),
        ], custom_spinner=html.Div('Merging & optimizing...', className='loading-text')),
        save_controls('merge-save', 'merge-download'),
        dcc.Store(id='merge-result'),
        dcc.Store(id='merge-result-view-data'),
    ])
