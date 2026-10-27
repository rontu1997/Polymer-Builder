"""Save UI integration; only the native OS picker/file handle is substituted.

Run: uv run --with playwright==1.62.0 python tests/browser_save.py
"""
import sys
import tempfile
import threading
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from playwright.sync_api import expect, sync_playwright
from rdkit import Chem
from werkzeug.serving import make_server

from app import app
from chem.io import load_structure
from browser_merge_controls import QuietHandler, pdb_file, pick_atom
from test_download import polymer_state


PICKER = '''
window.savedFiles = [];
window.pickerCalls = [];
window.pickerMode = 'save';
window.showSaveFilePicker = async function(options) {
    window.pickerCalls.push({options, activated: navigator.userActivation.isActive});
    if (window.pickerMode === 'cancel') { throw new DOMException('Cancelled', 'AbortError'); }
    const name = window.overrideFilename || options.suggestedName;
    const target = {name, content: '', closed: false};
    return {name, createWritable: async () => {
        if (window.pickerMode === 'write-error') { throw new Error('Disk write denied'); }
        return {
            write: async text => {target.content = text;},
            close: async () => {target.closed = true; window.savedFiles.push(target);},
            abort: async () => {target.aborted = true;}
        };
    }};
};
'''


def select_format(page, prefix, fmt):
    page.locator(f'#{prefix}-format').click()
    page.get_by_text(f'{fmt[1:].upper()} ({fmt})', exact=True).filter(visible=True).last.click()


def assert_export(page, prefix, fmt, expected):
    select_format(page, prefix, fmt)
    before = page.evaluate('savedFiles.length')
    page.locator(f'#{prefix}-button').click()
    page.wait_for_function('count => savedFiles.length > count', arg=before)
    expect(page.locator(f'#{prefix}-status')).to_contain_text('Saved ', timeout=15000)
    assert page.evaluate('savedFiles.length') == before + 1
    saved = page.evaluate('savedFiles.at(-1)')
    call = page.evaluate('pickerCalls.at(-1)')
    assert call['activated'], 'Picker must open directly from a user click'
    assert list(call['options']['types'][0]['accept'].values())[0] == [fmt]
    assert saved['name'].endswith(fmt) and saved['closed']
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / saved['name']
        path.write_text(saved['content'])
        assert Chem.MolToSmiles(Chem.RemoveHs(load_structure(path))) == expected


def main():
    server = make_server('127.0.0.1', 0, app.server, threaded=True, request_handler=QuietHandler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(headless=True)
            page = browser.new_page(viewport={'width': 1500, 'height': 1400}, accept_downloads=True)
            page.add_init_script(PICKER)
            errors = []
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.goto(f'http://127.0.0.1:{server.server_port}', wait_until='networkidle')
            # Inject a real molecular result at the build/export boundary.
            page.evaluate("state => dash_clientside.set_props('polymer-state', {data: state})", polymer_state())
            page.evaluate("dash_clientside.set_props('save-button', {disabled: false})")
            expect(page.locator('#save-format')).to_contain_text('PDB (.pdb)')
            for fmt in ('.pdb', '.mol', '.mol2'):
                assert_export(page, 'save', fmt, 'CCO')

            # Choosing a different file type in the OS dialog determines conversion.
            select_format(page, 'save', '.pdb')
            page.evaluate("window.overrideFilename = 'chosen-name.mol2'")
            page.locator('#save-button').click()
            expect(page.locator('#save-status')).to_contain_text('Saved chosen-name.mol2.')
            assert '@<TRIPOS>MOLECULE' in page.evaluate('savedFiles.at(-1).content')
            page.evaluate('window.overrideFilename = null')
            count = page.evaluate('savedFiles.length')
            page.evaluate("window.pickerMode = 'cancel'")
            page.locator('#save-button').click()
            expect(page.locator('#save-status')).to_contain_text('Save cancelled.')
            assert page.evaluate('savedFiles.length') == count
            page.evaluate("window.pickerMode = 'write-error'")
            page.locator('#save-button').click()
            expect(page.locator('#save-status')).to_contain_text('Disk write denied')
            assert page.evaluate('savedFiles.length') == count
            page.evaluate("window.pickerMode = 'save'")

            page.get_by_text('PDB Merger', exact=True).click()
            expect(page.locator('#merger-panel')).to_be_visible()
            expect(page.locator('#merge-save-format')).to_contain_text('PDB (.pdb)')
            page.locator('#merge-a-upload input').set_input_files(pdb_file('first.pdb', 'CC'))
            page.locator('#merge-b-upload input').set_input_files(pdb_file('second.pdb', 'CO'))
            expect(page.locator('#merge-a-status')).to_contain_text('Loaded first.pdb')
            expect(page.locator('#merge-b-status')).to_contain_text('Loaded second.pdb')
            for side in ('a', 'b'):
                pick_atom(page, side, 0, 'head')
                pick_atom(page, side, 1, 'tail')
            page.locator('#merge-button').click()
            expect(page.locator('#merge-save-button')).to_be_enabled(timeout=30000)
            for fmt in ('.pdb', '.mol', '.mol2'):
                assert_export(page, 'merge-save', fmt, 'CCCO')

            # No native picker: the selected format still downloads correctly.
            page.evaluate('window.showSaveFilePicker = undefined; window.crypto.randomUUID = undefined;')
            select_format(page, 'merge-save', '.mol')
            with page.expect_download() as download:
                page.locator('#merge-save-button').click()
            assert download.value.suggested_filename == 'first_second_merged.mol'
            downloaded_mol = Chem.MolFromMolBlock(Path(download.value.path()).read_text(), removeHs=False)
            assert Chem.MolToSmiles(Chem.RemoveHs(downloaded_mol)) == 'CCCO'
            expect(page.locator('#merge-save-status')).to_contain_text('Ask where to save each file')
            assert not errors, errors
            print('PASS: native picker handoff with user activation, all three formats in both tabs, chosen filename, cancellation, write errors, and download fallback.')
            browser.close()
    finally:
        server.shutdown()


if __name__ == '__main__':
    main()
