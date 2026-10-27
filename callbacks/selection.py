from dash import Input, Output, State, html, no_update

from chem.validation import validate_linkage_atom
from chem.viewer import pdb_block_to_mol


def register_selection_callbacks(app):
    @app.callback(
        Output("app-state", "data", allow_duplicate=True),
        Output("selection-status", "children"),
        Output("monomer-view-data", "data", allow_duplicate=True),
        Input("selected-atom-store", "data"),
        State("app-state", "data"),
        State("picking-role", "value"),
        prevent_initial_call=True,
    )
    def handle_atom_click(clicked, state, role):
        if not clicked or not state or not state.get("monomer_pdbblock"):
            return no_update, no_update, no_update

        mol = pdb_block_to_mol(state["monomer_pdbblock"])
        idx = clicked["index"]

        other_role_idx = state.get("tail_idx" if role == "head" else "head_idx")
        if other_role_idx is not None and idx == other_role_idx:
            message = "Head and tail atom must be different atoms."
            return no_update, html.Div(message, className="error-text"), no_update

        check = validate_linkage_atom(mol, idx)
        if not check.ok:
            return no_update, html.Div(check.message, className="error-text"), no_update

        new_state = dict(state)
        new_state[f"{role}_idx"] = idx

        atom = mol.GetAtomWithIdx(idx)
        status = html.Div(f"{role.capitalize()} atom set to #{idx} ({atom.GetSymbol()}).")
        view_data = {
            "pdbblock": state["monomer_pdbblock"],
            "head_idx": new_state.get("head_idx"),
            "tail_idx": new_state.get("tail_idx"),
        }
        return new_state, status, view_data
