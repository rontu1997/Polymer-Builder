from dash import Input, Output, State, html, no_update

from chem.shape_solver import solve_circle, solve_line, solve_polygon
from chem.viewer import pdb_block_to_mol

_N_SIDES = {"triangle": 3, "square": 4, "pentagon": 5, "hexagon": 6}


def _flatness_note(flatness_score):
    """Report the measured flatness score honestly instead of letting the user find
    out from the rendered shape after a 20-30s build. A low score is a real property
    of this monomer's available torsion roots (see chem/shape_solver.py's
    _edge_flatness_score docstring) -- more code can't force flatness that isn't
    chemically available, so this only informs, it never blocks the build."""
    if flatness_score >= 0.9:
        verdict = "should render close to flat/closed"
    elif flatness_score >= 0.6:
        verdict = "expect noticeable 3D drift between corners (not a flat shape)"
    else:
        verdict = (
            "expect significant 3D corkscrewing -- this monomer's available torsion "
            "roots don't support a flat version of this shape at this corner count. "
            "Side lengths may also come out visibly unequal after optimization as a "
            "result, especially the last side."
        )
    return f"Corner flatness score: {flatness_score:.2f} ({verdict})"


def _closure_note(closure_gap_nm, requested_length_nm):
    """flatness_score only measures *pairwise* consistency between two adjacent
    corners' turn axes -- for shapes with more corners (pentagon, hexagon), small
    pairwise drift can still compound into a large overall gap that a good-looking
    flatness score doesn't predict (confirmed by direct measurement: a hexagon case
    scored 0.96 flatness yet closed with a gap *larger* than its own side length).
    closure_gap_nm is the actual, comprehensive measurement -- report it plainly
    rather than let a good flatness score imply a good closure."""
    ratio = closure_gap_nm / requested_length_nm if requested_length_nm else 0.0
    if ratio < 0.15:
        verdict = "closes well"
    elif ratio < 0.5:
        verdict = "leaves a visible but partial gap"
    else:
        verdict = "does not meaningfully close -- expect the shape to look open/disconnected at this corner"
    return f"Closure quality: gap is {ratio:.0%} of one side length ({verdict})"


def register_shape_callbacks(app):
    @app.callback(
        Output("shape-length-nm", "disabled"),
        Output("calc-shape-button", "disabled"),
        Output("n-units", "disabled"),
        Output("torsion-deg", "disabled"),
        Input("shape-select", "value"),
        Input("shape-length-nm", "value"),
        Input("app-state", "data"),
    )
    def update_shape_readiness(shape, length_nm, state):
        if shape == "manual":
            return True, True, False, False
        has_head_tail = bool(state) and state.get("head_idx") is not None and state.get("tail_idx") is not None
        ready = has_head_tail and length_nm is not None and length_nm > 0
        # Linear chains use one torsion value, so their repeat count can be edited
        # independently. Other presets require a fixed per-junction schedule.
        return False, not ready, shape != "line", True

    @app.callback(
        Output("n-units", "value"),
        Output("torsion-deg", "value"),
        Output("shape-status", "children"),
        Output("shape-torsion-schedule", "data"),
        Input("calc-shape-button", "n_clicks"),
        State("shape-select", "value"),
        State("shape-length-nm", "value"),
        State("app-state", "data"),
        prevent_initial_call=True,
    )
    def calculate_geometry(n_clicks, shape, length_nm, state):
        if not state or not state.get("monomer_pdbblock"):
            return no_update, no_update, no_update, no_update

        monomer = pdb_block_to_mol(state["monomer_pdbblock"])
        head_idx, tail_idx = state["head_idx"], state["tail_idx"]

        if shape == "line":
            result = solve_line(monomer, head_idx, tail_idx, length_nm)
            size_word = "total length"
        elif shape in _N_SIDES:
            result = solve_polygon(monomer, head_idx, tail_idx, _N_SIDES[shape], length_nm)
            size_word = "side length"
        elif shape == "circle":
            result = solve_circle(monomer, head_idx, tail_idx, length_nm)
            size_word = "radius"
        else:
            return no_update, no_update, no_update, no_update

        if not result.ok:
            return no_update, no_update, html.Div(result.message, className="error-text"), None

        is_schedule = isinstance(result.torsion_deg, list)
        # The slider can't represent a per-junction schedule -- when one is active it
        # shows the corner angle only, as a reference value (read from corner_torsion_deg,
        # NOT torsion_deg[-1] -- for a polygon, the schedule's last entries are the
        # closing side's straight-run torsion, not a corner), and the schedule itself
        # (the thing that's actually authoritative for the build) lives in the store.
        slider_value = round((result.corner_torsion_deg if is_schedule else result.torsion_deg), 2)

        n_sides_effective = (result.n_units - 1) // result.units_per_edge if is_schedule else None
        status_lines = [f"Repeat units: {result.n_units}" + (f" ({result.units_per_edge} per side x {n_sides_effective}, plus a closing side)" if is_schedule else "")]
        if is_schedule:
            status_lines.append(
                f"Straight-run torsion: {result.torsion_deg[0]:.1f} deg   "
                f"Corner torsion: {result.corner_torsion_deg:.1f} deg (per-junction schedule active -- "
                f"the slider above shows the corner angle only)"
            )
        else:
            status_lines.append(f"Torsion: {result.torsion_deg:.1f} deg")
        status_lines.append(
            f"Requested {size_word}: {length_nm:.2f} nm   "
            f"Predicted (pre-optimization): {result.achieved_length_nm:.2f} nm ({result.relative_error:.1%} off)"
        )
        if result.closing_side_length_nm is not None:
            status_lines.append(
                f"Closing side length: {result.closing_side_length_nm:.2f} nm (sized by measurement to close "
                f"the loop, not fit to the requested length -- it may differ from the other sides)   "
                f"Closure gap: {result.closure_gap_nm:.2f} nm"
            )
            status_lines.append(_closure_note(result.closure_gap_nm, result.requested_length_nm))
        if shape in _N_SIDES or shape == "circle":
            status_lines.append(
                "Open approximation: correct per-corner turn angle and side count, but not a "
                "guaranteed flat, closed shape -- see chem/shape_solver.py's design note."
            )
            if result.flatness_score is not None:
                status_lines.append(_flatness_note(result.flatness_score))

        stats = html.Div([html.Div(line) for line in status_lines], className="stats-box")
        schedule = result.torsion_deg if is_schedule else None
        return result.n_units, slider_value, stats, schedule
