"""Optional end-to-end check: uv run --with playwright==1.62.0 python tests/browser_merge_controls.py."""
import sys
import threading
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from playwright.sync_api import expect, sync_playwright
from rdkit import Chem
from rdkit.Chem import AllChem
from werkzeug.serving import WSGIRequestHandler, make_server

from app import app
from chem.viewer import mol_to_pdb_block, pdb_block_to_mol


class QuietHandler(WSGIRequestHandler):
    def log_request(self, *args, **kwargs):
        pass


def pdb_file(name, smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=7) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    return {'name': name, 'mimeType': 'chemical/x-pdb', 'buffer': mol_to_pdb_block(mol).encode()}


def view(page, viewer_id):
    return page.evaluate('id => document.getElementById(id)._3dmolViewer.getView()', viewer_id)


def check_drag_controls(page, viewer_id):
    page.wait_for_function('id => Boolean(document.getElementById(id)._3dmolViewer?.getModel())', arg=viewer_id)
    box = page.locator(f'#{viewer_id}').bounding_box()
    x, y = box['x'] + 50, box['y'] + 50
    before = view(page, viewer_id)
    page.mouse.move(x, y)
    page.mouse.down(button='right')
    page.mouse.move(x + 55, y + 35, steps=10)
    page.mouse.up(button='right')
    after = view(page, viewer_id)
    assert any(abs(a - b) > 1e-4 for a, b in zip(before[:3], after[:3])), 'Right-drag must translate'
    assert all(abs(a - b) < 1e-6 for a, b in zip(before[3:], after[3:])), 'Right-drag must not zoom or rotate'
    page.mouse.move(x, y)
    page.mouse.down(button='left')
    page.mouse.move(x + 35, y + 45, steps=10)
    page.mouse.up(button='left')
    rotated = view(page, viewer_id)
    assert any(abs(a - b) > 1e-4 for a, b in zip(after[4:], rotated[4:])), 'Left-drag must rotate'


def pick_atom(page, side, index, role):
    viewer_id = f'merge-{side}-viewer'
    before_role = view(page, viewer_id)
    page.locator(f'#merge-{side}-picking-role').get_by_text('Head (blue)' if role == 'head' else 'Tail (red)').click()
    page.wait_for_function("() => document.querySelectorAll('[data-dash-is-loading=\"true\"]').length === 0")
    assert all(abs(a - b) < 1e-6 for a, b in zip(before_role, view(page, viewer_id))), 'Role switching must preserve the view'
    for attempt in range(8):
        position = page.evaluate('''([id, index, attempt]) => {
            const viewer = document.getElementById(id)._3dmolViewer;
            if (attempt) { viewer.rotate(35, 'y'); viewer.render(); }
            const atom = viewer.getModel().selectedAtoms({index})[0];
            const point = viewer.modelToScreen(atom);
            return {x: point.x - window.scrollX, y: point.y - window.scrollY};
        }''', [viewer_id, index, attempt])
        before_click = view(page, viewer_id)
        page.mouse.click(position['x'], position['y'])
        try:
            expect(page.locator(f'#merge-{side}-status')).to_contain_text(f'{role.capitalize()}: #{index}', timeout=1000)
            assert all(abs(a - b) < 1e-6 for a, b in zip(before_click, view(page, viewer_id))), 'Picking must preserve translation/rotation'
            return
        except AssertionError:
            if attempt == 7:
                raise


def main():
    server = make_server('127.0.0.1', 0, app.server, threaded=True, request_handler=QuietHandler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(headless=True)
            page = browser.new_page(viewport={'width': 1500, 'height': 1400}, accept_downloads=True)
            # Exercise the ordinary download fallback; native picker integration
            # is checked separately without trying to automate an OS dialog.
            page.add_init_script('window.showSaveFilePicker = undefined;')
            errors = []
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.goto(f'http://127.0.0.1:{server.server_port}', wait_until='networkidle')
            expect(page.locator('#n-units')).to_be_enabled()
            page.locator('#n-units').fill('7')
            page.get_by_text('PDB Merger', exact=True).click()
            expect(page.locator('#merger-panel')).to_be_visible()
            expect(page.locator('#homopolymer-panel')).to_be_hidden()
            page.locator('#merge-a-upload input').set_input_files(pdb_file('first.pdb', 'CC'))
            page.locator('#merge-b-upload input').set_input_files(pdb_file('second.pdb', 'CO'))
            expect(page.locator('#merge-a-status')).to_contain_text('Loaded first.pdb')
            expect(page.locator('#merge-b-status')).to_contain_text('Loaded second.pdb')
            for side in ('a', 'b'):
                check_drag_controls(page, f'merge-{side}-viewer')
                expect(page.locator(f'#merge-{side}-status')).to_contain_text('Loaded')
                pick_atom(page, side, 0, 'head')
                expect(page.locator('#merge-button')).to_be_disabled()
                pick_atom(page, side, 1, 'tail')
            expect(page.locator('#merge-button')).to_be_enabled()
            page.locator('#merge-button').click()
            expect(page.locator('#merge-stats')).to_contain_text('A tail atom #1 to B head atom #0', timeout=30000)
            expect(page.locator('#merge-save-button')).to_be_enabled()
            with page.expect_download() as info:
                page.locator('#merge-save-button').click()
            mol = pdb_block_to_mol(Path(info.value.path()).read_text())
            assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == 'CCCO'
            page.locator('#merge-join-mode').click()
            page.get_by_text('Head–Tail: A head → B tail', exact=True).click()
            expect(page.locator('#merge-save-button')).to_be_disabled()
            expect(page.locator('#merge-stats')).to_be_empty()
            page.locator('#merge-button').click()
            expect(page.locator('#merge-stats')).to_contain_text('A head atom #0 to B tail atom #1', timeout=30000)
            with page.expect_download() as info:
                page.locator('#merge-save-button').click()
            mol = pdb_block_to_mol(Path(info.value.path()).read_text())
            assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == 'CCOC'
            page.get_by_text('Homopolymer Builder', exact=True).click()
            expect(page.locator('#homopolymer-panel')).to_be_visible()
            expect(page.locator('#n-units')).to_have_value('7')
            page.get_by_text('PDB Merger', exact=True).click()
            expect(page.locator('#merger-panel')).to_be_visible()
            expect(page.locator('#homopolymer-panel')).to_be_hidden()
            expect(page.locator('#merge-a-status')).to_contain_text('Head: #0 (C) | Tail: #1 (C)')
            screenshot = Path('artifacts/pdb-merger.png')
            screenshot.parent.mkdir(exist_ok=True)
            page.screenshot(path=str(screenshot), full_page=True)
            page.locator('#merge-a-upload input').set_input_files(pdb_file('replacement.pdb', 'CCC'))
            expect(page.locator('#merge-a-status')).to_contain_text('Loaded replacement.pdb')
            expect(page.locator('#merge-button')).to_be_disabled()
            expect(page.locator('#merge-save-button')).to_be_disabled()
            page.wait_for_function("!document.getElementById('merge-result-viewer')._3dmolViewer.getModel()")
            assert not errors, errors
            print('PASS: head/tail selection, A/B join choices, right-drag translation, left-drag rotation, stable view, PDB downloads, tab persistence, and output invalidation.')
            browser.close()
    finally:
        server.shutdown()


if __name__ == '__main__':
    main()
