from dash import dcc, html


def save_controls(prefix, payload_id):
    return html.Div([
        html.Div([
            html.Label('Save as:', htmlFor=f'{prefix}-format'),
            dcc.Dropdown(
                id=f'{prefix}-format', value='.pdb', clearable=False,
                options=[{'label': 'PDB (.pdb)', 'value': '.pdb'},
                         {'label': 'MOL (.mol)', 'value': '.mol'},
                         {'label': 'MOL2 (.mol2)', 'value': '.mol2'}],
                style={'width': '180px'},
            ),
            html.Button('Save…', id=f'{prefix}-button', disabled=True),
        ], className='controls-row'),
        html.Div(id=f'{prefix}-status', role='status'),
        dcc.Store(id=f'{prefix}-request'),
        dcc.Store(id=payload_id),
        html.Div(id=f'{prefix}-config', style={'display': 'none'}),
        html.Div(id=f'{prefix}-receipt', style={'display': 'none'}),
    ])
