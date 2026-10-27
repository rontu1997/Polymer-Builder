# PDB Merger Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan task-by-task, with a final review using superpowers:requesting-code-review.

**Goal:** Join two uploaded PDB structures at selected atoms in a new tab.

**Architecture:** Add isolated chemistry and callback modules, reuse existing bond
placement/optimization and the 3D viewer, and keep both tab panels mounted so input
state survives switching.

**Tech Stack:** Python, Dash, RDKit, NumPy, 3Dmol.js, pytest.

**Spec:** `docs/superpowers/specs/2026-10-07-pdb-merger.md`

## Global Constraints

- Python >=3.11; no new runtime dependencies.
- Preserve homopolymer behavior and current user changes (the folder has no Git repository).
- The merger creates exactly one single bond between two PDB inputs.
- Inputs must each contain one connected structure with finite coordinates.
- Existing atom cap applies; PDB output must fit five-digit serial numbers.

## Review Focus

- A removed hydrogen must not also serve as the torsion reference (Task 1).
- Invalid atom indices, degenerate geometry, and nonfinite inputs produce readable errors (Task 1).
- Upload replacement and changed parameters invalidate stale output (Task 2).
- Delayed click events from a previous upload must be ignored (Task 2).
- Tab switching retains state and resizes hidden viewers on return (Task 3).

### Task 1: Chemical merger

**Files:** Create `chem/merge.py`, `tests/test_merge.py`.

**Interfaces:** `merge_structures(first, second, first_idx, second_idx, torsion_deg=180.0) -> MergeResult`, with `mol`, `new_bonds`, and `torsion_specs`; failures raise `MergeError`.

- [x] Write and run failing graph/geometry, hydrogen, validation, metadata, and immutability tests.
- [x] Implement validation, retained reference selection, second-structure placement, bonding, hydrogen removal, and index remapping.
- [x] Run `uv run --locked python -m pytest tests/test_merge.py -q` and confirm all pass.

### Task 2: Upload, select, merge, and download

**Files:** Create `callbacks/merge.py`, `tests/test_merge_callbacks.py`.

**Interfaces:** `register_merge_callbacks(app)`; consumes Task 1 and existing `optimize_geometry`/PDB serialization. Each input store carries PDB data, selection, filename, and an upload identity.

- [x] Write and run failing tests for two real PDB uploads, validated clicks, merged topology/download, upload errors, outdated clicks, and stale-output invalidation.
- [x] Implement independent input callbacks and one result callback triggered by both input states and parameters. Merge runs only on its button; other changes clear output.
- [x] Run `uv run --locked python -m pytest tests/test_merge_callbacks.py -q` and confirm all pass.

### Task 3: Tabs, viewer integration, and documentation

**Files:** Create `ui/merge.py`, modify `app.py`, `assets/viewer.js`, `assets/style.css`, `README.md`.

- [x] Add persistent Homopolymer Builder/PDB Merger panels, controls, stores, and viewer callbacks.
- [x] Clear viewers when their data clears, include upload identity in merger click events, and resize viewers when revealed.
- [x] Verify served layout, callback wiring, and browser interaction where tooling permits.
- [x] Update usage documentation, run all tests, and obtain an independent review.

## Execution ledger

- Ruling: implement inline in the authorized workspace because there is no Git repository; no worktree or commits apply.
- Pre-flight: Task 2 consumes Task 1's molecule/bonds/torsions; Task 3 consumes Task 2's prefixed component IDs. Existing homopolymer IDs remain unchanged.
- Ruling: add missing hydrogens on working copies before placement. The initial implicit-H test exposed collinear reference geometry at terminal heavy atoms; explicit H geometry fixes torsion placement and supplies optimization parameters. Output includes the completed hydrogens.
- Task 1: complete. Chemistry tests verified graph, torsions, bond length, hydrogen handling, validation, immutability, and metadata. An RDKit addResidueInfo hang was avoided by copying new H metadata explicitly.
- Task 2: complete. Real Dash requests verified upload/select/merge/download. Reviewer findings were reproduced with failing tests and fixed: multi-model inputs, blank XYZ coordinates, proximity bonds overriding CONECT, old clicks during upload, and stale builds/downloads.
- Task 3: complete. Headless Chromium verified actual atom clicks, optimized merge, connected PDB download, both tabs' state retention, model resize, and upload invalidation with no browser errors. Screenshot: `artifacts/pdb-merger.png`.
- Final validation: `uv run --locked python -m pytest tests -q` -> 93 passed. Reviewer independently verified 12 targeted regression tests and confirmed no outstanding findings in the reviewed fixes.
- Deferred from review: inherited large-structure optimization costs; broad serializer refactoring. Intentional HETATM conversion matches the existing viewer's connectivity handling.

## Follow-up: A/B head and tail selection and mouse controls

- [x] Add independent head/tail assignments for polymers A and B, with blue/red highlights and distinct-atom validation.
- [x] Add explicit head-head, tail-tail, head-tail, and tail-head joins. Include mode and all endpoints in stale-output/download guards.
- [x] Capture picking role with click events and preserve the view during highlight updates.
- [x] Replace native right-drag zoom with view translation, preserve left-drag rotation and scroll zoom, and suppress right-click atom selection.
- [x] Verify all four resulting molecular topologies with real callback tests. Verify real browser drags, clicks, mode changes, and downloads using `tests/browser_merge_controls.py`.
- Validation: full suite passed (102 tests); browser check passed. Independent read-only review found no blocking issues.
- Minor review note: swapping already assigned head/tail labels requires reuploading or temporarily selecting another eligible atom. Join choices allow reversing the connection directly.

## Follow-up: Save As and export format choice

- [x] Shared Save controls in both tabs with PDB default and MOL/MOL2 alternatives.
- [x] Shared server export helper, using the existing format writers and returning request-matched responses/errors.
- [x] Native browser picker launched from the click before conversion; selected extension controls format. Cancellation, stale data, write errors, timeouts, and unsupported-browser download fallback handled.
- [x] Preserve merger signature/upload guards and add actual molecular round-trip tests for all formats and both tabs.
- Validation: full suite passed (111 tests); merger browser workflow and Save browser integration both passed. OS picker/file handles were substituted in the Save automation; click activation, conversions, and fallback downloads were real.
- Independent read-only review: no blocking findings. Optional delayed-write/double-response/timeout stress coverage deferred.
