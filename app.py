from dash import Dash, Input, Output, dcc, html

from chem.validation import MAX_REPEAT_UNITS
from callbacks.build import register_build_callbacks
from callbacks.download import register_download_callbacks
from callbacks.merge import register_merge_callbacks
from callbacks.selection import register_selection_callbacks
from callbacks.shape import register_shape_callbacks
from callbacks.upload import register_upload_callbacks
from ui.merge import merger_layout
from ui.save import save_controls

app = Dash(__name__)
app.title = "Polymer Builder"

homopolymer_layout = html.Div(
    [
        html.P(
            "Upload a monomer, pick a head and tail atom, choose a repeat count "
            "and backbone torsion, then build and optimize the chain. This tool "
            "does not infer polymerization chemistry -- you designate the two "
            "linkage atoms yourself."
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.H4("1. Monomer"),
                        dcc.Upload(
                            id="upload-structure",
                            children=html.Div(["Drag and drop, or ", html.A("select a .pdb/.mol/.mol2 file")]),
                            style={
                                "width": "100%",
                                "height": "60px",
                                "lineHeight": "60px",
                                "borderWidth": "1px",
                                "borderStyle": "dashed",
                                "borderRadius": "6px",
                                "textAlign": "center",
                            },
                        ),
                        html.Div(id="upload-status"),
                        html.Div(
                            [
                                html.Label("Picking:"),
                                dcc.RadioItems(
                                    id="picking-role",
                                    options=[
                                        {"label": " Head atom", "value": "head"},
                                        {"label": " Tail atom", "value": "tail"},
                                    ],
                                    value="head",
                                    inline=True,
                                ),
                            ],
                            className="controls-row",
                        ),
                        html.Div(id="selection-status"),
                        html.Div(id="viewer-left", className="viewer-box"),
                        html.Div(id="viewer-left-dummy", style={"display": "none"}),
                        html.Div(
                            [
                                html.Label("Shape preset:"),
                                dcc.Dropdown(
                                    id="shape-select",
                                    options=[
                                        {"label": "Linear -- enter total length (nm)", "value": "line"},
                                    ],
                                    value="line",
                                    clearable=False,
                                ),
                                dcc.Input(
                                    id="shape-length-nm",
                                    type="number",
                                    min=0,
                                    step="any",
                                    placeholder="Total length (nm)",
                                    disabled=False,
                                ),
                                html.Button("Calculate Geometry", id="calc-shape-button", disabled=True),
                            ],
                            className="controls-row",
                        ),
                        html.Div(id="shape-status"),
                        html.Div(
                            [
                                html.Label("Repeat units:"),
                                dcc.Input(id="n-units", type="number", min=1, max=MAX_REPEAT_UNITS, step=1, value=2),
                                html.Label("Torsion (deg):"),
                                dcc.Slider(
                                    id="torsion-deg",
                                    min=-180,
                                    max=180,
                                    step=1,
                                    value=180,
                                    marks={-180: "-180", 0: "0", 180: "180"},
                                    tooltip={"placement": "bottom"},
                                ),
                                html.Label("Force field:"),
                                dcc.Dropdown(
                                    id="force-field",
                                    options=[
                                        {"label": "MMFF94", "value": "MMFF94"},
                                        {"label": "UFF", "value": "UFF"},
                                    ],
                                    value="MMFF94",
                                    style={"width": "140px"},
                                    clearable=False,
                                ),
                            ],
                            className="controls-row",
                        ),
                        html.Div(id="atom-count-estimate"),
                        html.Button("Build & Optimize", id="build-button", disabled=True),
                    ],
                    className="panel",
                ),
                html.Div(
                    [
                        html.H4("2. Built polymer"),
                        dcc.Loading(
                            [
                                html.Div(id="viewer-right", className="viewer-box"),
                                html.Div(id="viewer-right-dummy", style={"display": "none"}),
                                html.Div(id="stats-panel"),
                            ],
                            type="circle",
                            overlay_style={"visibility": "visible", "opacity": 0.5},
                            custom_spinner=html.Div("Building & optimizing...", className="loading-text"),
                        ),
                        save_controls("save", "download-structure"),
                    ],
                    className="panel",
                ),
            ],
            className="panels-row",
        ),
        dcc.Store(id="app-state"),
        dcc.Store(id="monomer-view-data"),
        dcc.Store(id="polymer-view-data"),
        dcc.Store(id="polymer-state"),
        dcc.Store(id="selected-atom-store"),
        dcc.Store(id="shape-torsion-schedule"),
    ]
)

app.layout = html.Div([
    html.H2("Polymer Builder"),
    dcc.Tabs(id="builder-tabs", value="homopolymer", children=[
        dcc.Tab(label="Homopolymer Builder", value="homopolymer"),
        dcc.Tab(label="PDB Merger", value="merger"),
    ]),
    # Both panels remain mounted so uploads, selected atoms, and controls survive
    # tab changes. The viewer's ResizeObserver handles reveal/resizing.
    html.Div(homopolymer_layout, id="homopolymer-panel"),
    html.Div(merger_layout(), id="merger-panel", style={"display": "none"}),
])


@app.callback(
    Output("homopolymer-panel", "style"), Output("merger-panel", "style"),
    Input("builder-tabs", "value"),
)
def switch_builder_tab(tab):
    return ({"display": "none"}, {}) if tab == "merger" else ({}, {"display": "none"})


app.clientside_callback(
    """
    function(data) {
        return window.dash_clientside.polymer_viewer.render(
            data && data.pdbblock, 'viewer-left', 'selected-atom-store',
            data && data.head_idx, data && data.tail_idx, true
        );
    }
    """,
    Output("viewer-left-dummy", "children"),
    Input("monomer-view-data", "data"),
)

app.clientside_callback(
    """
    function(data) {
        return window.dash_clientside.polymer_viewer.render(
            data && data.pdbblock, 'viewer-right', null, null, null, false
        );
    }
    """,
    Output("viewer-right-dummy", "children"),
    Input("polymer-view-data", "data"),
)

for side in ("a", "b"):
    app.clientside_callback(
        f"""
        function(data, pickingRole) {{
            return window.dash_clientside.polymer_viewer.render(
                data && data.pdbblock, 'merge-{side}-viewer', 'merge-{side}-selected-atom',
                data && data.head_idx, data && data.tail_idx, true, data && data.upload_id, pickingRole
            );
        }}
        """,
        Output(f"merge-{side}-viewer-dummy", "children"),
        Input(f"merge-{side}-view-data", "data"),
        Input(f"merge-{side}-picking-role", "value"),
    )

app.clientside_callback(
    """
    function(data) {
        return window.dash_clientside.polymer_viewer.render(
            data && data.pdbblock, 'merge-result-viewer', null, null, null, false
        );
    }
    """,
    Output("merge-result-viewer-dummy", "children"),
    Input("merge-result-view-data", "data"),
)

for prefix, state_id, payload_id, merged in (
    ("save", "polymer-state", "download-structure", "false"),
    ("merge-save", "merge-result", "merge-download", "true"),
):
    app.clientside_callback(
        f"""function(data, format) {{
            return window.polymer_save.configure('{prefix}', data, format, {merged});
        }}""",
        Output(f"{prefix}-config", "children"),
        Input(state_id, "data"), Input(f"{prefix}-format", "value"),
    )
    app.clientside_callback(
        f"""function(payload) {{ return window.polymer_save.receive('{prefix}', payload); }}""",
        Output(f"{prefix}-receipt", "children"), Input(payload_id, "data"),
        prevent_initial_call=True,
    )

register_upload_callbacks(app)
register_selection_callbacks(app)
register_shape_callbacks(app)
register_build_callbacks(app)
register_download_callbacks(app)
register_merge_callbacks(app)

if __name__ == "__main__":
    app.run(debug=True)
