import math

from dash import Input, Output, State, html, no_update

from chem.builder import PolymerBuildError, build_polymer
from chem.optimize import OptimizationError, optimize_geometry
from chem.validation import MAX_ATOMS_HARD_CAP, MAX_REPEAT_UNITS, validate_pre_build
from chem.viewer import mol_to_pdb_block, pdb_block_to_mol


def register_build_callbacks(app):
    @app.callback(
        Output("build-button", "disabled"),
        Output("atom-count-estimate", "children"),
        Input("app-state", "data"),
        Input("n-units", "value"),
    )
    def update_build_readiness(state, n_units):
        if not state or not state.get("monomer_pdbblock"):
            return True, ""

        mol = pdb_block_to_mol(state["monomer_pdbblock"])
        n_units = n_units or 1
        estimate_text = f"Estimated atoms: ~{mol.GetNumAtoms() * n_units} (cap {MAX_ATOMS_HARD_CAP}, max {MAX_REPEAT_UNITS} units)"

        head_idx, tail_idx = state.get("head_idx"), state.get("tail_idx")
        if head_idx is None or tail_idx is None:
            return True, estimate_text

        pre = validate_pre_build(mol, head_idx, tail_idx, n_units=n_units)
        return not pre.ok, estimate_text

    @app.callback(
        Output("polymer-view-data", "data"),
        Output("polymer-state", "data"),
        Output("stats-panel", "children"),
        Output("save-button", "disabled"),
        Input("build-button", "n_clicks"),
        State("app-state", "data"),
        State("n-units", "value"),
        State("torsion-deg", "value"),
        State("force-field", "value"),
        State("shape-torsion-schedule", "data"),
        prevent_initial_call=True,
    )
    def do_build(n_clicks, state, n_units, torsion_deg, force_field, shape_schedule):
        if not state or not state.get("monomer_pdbblock"):
            return no_update, no_update, no_update, no_update

        monomer = pdb_block_to_mol(state["monomer_pdbblock"])
        # A shape-solved per-junction schedule takes precedence over the slider's
        # scalar value, but only if its length still matches n_units -- editing
        # n-units after calculating a shape invalidates a stale schedule for free,
        # with no separate invalidation step needed.
        torsion_arg = shape_schedule if shape_schedule and len(shape_schedule) == n_units - 1 else torsion_deg

        try:
            build_result = build_polymer(
                monomer, state["head_idx"], state["tail_idx"], n_units, torsion_arg
            )
            optimized_mol, report = optimize_geometry(
                build_result.mol,
                build_result.torsion_specs,
                new_bonds=build_result.new_bonds,
                ff=force_field,
            )
        except (PolymerBuildError, OptimizationError) as exc:
            return None, {"error": str(exc)}, html.Div(str(exc), className="error-text"), True

        pdbblock = mol_to_pdb_block(optimized_mol)
        polymer_state = {
            "polymer_pdbblock": pdbblock,
            "input_format": state["input_format"],
            "filename": state["filename"],
        }
        stats = _render_stats(report, build_result, n_units)
        return {"pdbblock": pdbblock}, polymer_state, stats, False


def _render_stats(report, build_result, n_units):
    # Read the actual per-junction targets build_polymer used, rather than trusting a
    # single raw UI value -- correct whether the build used one scalar torsion or a
    # per-junction schedule (see chem/builder.py::build_polymer's torsion_deg param).
    targets = [t for *_, t in build_result.torsion_specs]
    if not targets or len(set(round(t, 3) for t in targets)) == 1:
        torsion_line = f"Target torsion: {targets[0]:.1f} deg" if targets else "Target torsion: n/a"
    else:
        counts = {}
        for t in targets:
            key = round(t, 1)
            counts[key] = counts.get(key, 0) + 1
        breakdown = ", ".join(f"{t:.1f} deg x{n}" for t, n in sorted(counts.items()))
        torsion_line = f"Torsion targets: {breakdown} (per-junction schedule)"

    lines = [
        f"Repeat units: {n_units}   Atoms: {build_result.mol.GetNumAtoms()}   New bonds: {len(build_result.new_bonds)}",
        f"Force field: {report.force_field_used}"
        + (" (MMFF94 unavailable, fell back to UFF)" if report.fell_back_to_uff else ""),
        f"Initial energy: {report.energy_before:.1f}   Final energy: {report.energy_after:.1f}",
        torsion_line,
    ]
    if report.validation.torsion_report:
        deltas = [abs(t["delta"]) for t in report.validation.torsion_report]
        avg_actual = _circular_mean_deg([t["actual"] for t in report.validation.torsion_report])
        lines.append(f"Final torsion (avg over junctions): {avg_actual:.1f} deg (max |delta| {max(deltas):.1f})")

    children = [html.Div(line) for line in lines]
    if report.validation.warnings:
        children.append(html.Div("Warnings:", className="warning-text"))
        for w in report.validation.warnings:
            children.append(html.Div(f"- {w}", className="warning-text"))

    return html.Div(children, className="stats-box")


def _circular_mean_deg(angles_deg):
    """Plain arithmetic mean breaks near the +-180 boundary -- e.g. 177 and -178
    average to -0.5, not ~179.5. Average as unit vectors instead."""
    radians = [math.radians(a) for a in angles_deg]
    mean_sin = sum(math.sin(r) for r in radians) / len(radians)
    mean_cos = sum(math.cos(r) for r in radians) / len(radians)
    return math.degrees(math.atan2(mean_sin, mean_cos))
