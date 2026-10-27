# Polymer Builder

A local Dash web app for building polymer structures from structure files. It has
two tools:

- **Homopolymer Builder**: upload a monomer (PDB/MOL/MOL2), pick a head and a tail
  atom in an interactive 3D view, choose a repeat count and backbone torsion (or a
  target chain length), then build and force-field-optimize a linear chain.
- **PDB Merger**: upload two PDB structures (polymers or linkers), pick a head and
  tail atom in each, and join one end of A to one end of B with a single bond.

Results can be saved as PDB, MOL, or MOL2.

## Contents

- [Chemistry model / scope](#chemistry-model--scope)
- [Features](#features)
- [How to use](#how-to-use)
- [Algorithms](#algorithms)
- [Known limitations](#known-limitations)
- [Setup](#setup-windows-using-uv)
- [Project layout](#project-layout)
- [Tests](#tests)

## Chemistry model / scope

**This tool does not infer polymerization chemistry.** It does not know that
styrene polymerizes through its vinyl group, or that a monomer has a "reactive"
functional group. You explicitly designate the atoms that form each new bond (head
and tail); the tool's job is to:

1. Check those atoms can structurally accept a new bond (valence feasibility, not
   chemical reactivity).
2. Build the resulting geometry: bond the tail of unit *i* to the head of unit
   *i+1* (or A to B in the merger), at the torsion angle you specify.
3. Force-field optimize (MMFF94, falling back to UFF when MMFF94 lacks parameters
   for an atom type present) while restraining the junction torsions near your
   target.

## Features

### Homopolymer Builder

- Upload a monomer as `.pdb`, `.mol`, or `.mol2`. A file with no 3D coordinates gets
  a generated conformer (hydrogens added, RDKit embedding, MMFF94 cleanup).
- Interactive 3D viewer with click-to-select head (blue) and tail (red) atoms.
- Valence-feasibility check on every selected atom, with a specific error message
  when an atom cannot accept a bond.
- Repeat count up to 1,000 units, with a live atom-count estimate against the
  100,000-atom cap.
- Backbone torsion slider from -180° to 180°, applied at every junction.
- **Linear length preset**: enter a total chain length in nm and the app computes the
  repeat count and torsion for the most extended chain. The repeat count stays
  editable afterward.
- Force field choice: MMFF94 or UFF.
- Build report: atoms, new bonds, force field used (and whether it fell back to UFF),
  energy before and after, target vs. final torsion, and warnings.
- Built polymer is shown in a second 3D viewer. Each unit gets its own residue
  number and every atom is written as HETATM, so units stay distinguishable and
  viewers read the real `CONECT` records.

### PDB Merger

- Two independent upload/selection panels (A and B), each with head/tail picking.
- Four join modes: A head–B head, A tail–B tail, A head–B tail, A tail–B head.
- Junction torsion slider and force field choice.
- Missing hydrogens are added; one hydrogen at each connection site is replaced
  when needed.
- Source residue names, residue numbers, and atom names are kept. Chain IDs become
  `A` and `B`. Connectivity records are written.
- Your uploaded files are never modified. Changing any input, atom, join mode,
  torsion, or force field clears the previous result so a stale structure cannot be
  saved.

### Viewing and saving

- 3D viewer controls: **left-drag rotates**, **right-drag translates**, **scroll
  zooms**, **left-click selects** an atom. Switching head/tail picking or selecting
  another atom keeps your current rotation, translation, and zoom. These gestures
  never change the molecular coordinates.
- Both tabs save as **PDB (.pdb)** (default), **MOL (.mol)**, or **MOL2 (.mol2)**,
  regardless of the uploaded format.
- **Save…** opens the browser's native Save As dialog where supported
  ([`showSaveFilePicker`](https://developer.mozilla.org/en-US/docs/Web/API/Window/showSaveFilePicker),
  secure context required; the default localhost URL qualifies). The chosen
  filename's extension decides the conversion. Cancelling downloads nothing. Where
  the dialog is unavailable, the app falls back to a normal download; enable your
  browser's **Ask where to save each file** to pick a location.
- Both tabs stay mounted, so uploads, selections, and settings survive switching
  tabs.

## How to use

### Run the app

```powershell
uv run app.py
```

Open the printed local URL (default `http://127.0.0.1:8050/`).

### Build a homopolymer

1. On the **Homopolymer Builder** tab, drag a `.pdb`, `.mol`, or `.mol2` monomer
   into the upload box. It must be a single connected molecule.
2. Choose **Head atom** and click the atom in the viewer that should bond to the
   previous unit. Choose **Tail atom** and click the atom that should bond to the
   next unit. Head and tail must be different atoms.
3. Choose the size, either way:
   - **Directly**: set **Repeat units** and the **Torsion** slider. 180° is the
     default.
   - **From a length**: enter a **Total length (nm)** under *Shape preset*
     ("Linear") and click **Calculate Geometry**. The app fills in the repeat count
     and torsion and reports the predicted length. The repeat count stays editable;
     changing it changes the length from the prediction.
4. Pick **MMFF94** or **UFF**.
5. Click **Build & Optimize**. Large chains can take tens of seconds.
6. Inspect the built chain and the report. Warnings (long or short new bonds,
   atom clashes, torsion drift, energy increase) are listed under the report.
7. Choose a format and click **Save…**.

### Merge two PDB structures

1. Open the **PDB Merger** tab and upload one connected, single-model `.pdb` into
   each of Polymer A and Polymer B (20 MB limit each).
2. For each, choose **Head** and click its head atom, then **Tail** and click its
   tail atom. Select distinct head and tail atoms in both.
3. Choose the ends to join (default: **A tail–B head**), set the junction torsion,
   and choose MMFF94 or UFF.
4. Click **Merge & Optimize**, inspect the result and any warnings, then choose a
   format and click **Save…**.

Atom records must have finite XYZ coordinates. If the input has `CONECT` records,
those define its bonds; if not, bonds are inferred from atom positions. A file whose
incomplete connectivity leaves disconnected components must be corrected before
merging.

## Algorithms

All geometry work is in `chem/`. RDKit is the chemistry engine; OpenBabel is used
only to read and write MOL2, which RDKit cannot write.

### 1. Loading and internal representation (`chem/io.py`, `chem/viewer.py`)

- PDB and MOL/SDF are read by RDKit with hydrogens kept. MOL2 is read by OpenBabel,
  converted to a MOL block, then parsed by RDKit.
- A structure without coordinates gets `AddHs`, RDKit distance-geometry embedding
  (fixed seed, so results repeat), and an MMFF94 cleanup.
- Between callbacks, every structure is held as a **PDB block** written with RDKit
  flavor 4 (`CONECT` in both directions). A MOL block would drop residue and
  atom-name metadata, and the viewer's atom index maps directly to the RDKit atom
  index because PDB preserves atom order.

### 2. Linkage-atom validation (`chem/validation.py`)

A bond needs a free valence slot on the atom. Using RDKit's typical valence for the
element (adjusted for formal charge), an atom is accepted when it has:

- open valence (valence below the maximum): the bond is simply added, or
- a hydrogen to give up, explicit (that H is deleted) or implicit.

Atoms at maximum valence with no hydrogen (for example a quaternary carbon) are
rejected. Pre-build checks also reject identical head/tail atoms, a missing
conformer, more than one fragment, and repeat/atom counts above the limits.

### 3. Chain construction (`chem/builder.py`)

Each new unit is a copy of the monomer placed with a full rigid-body transform so
the new bond and torsion are exact before any optimization:

1. **Leaving direction.** The direction a new bond should point from an atom is the
   vector to its leaving hydrogen when there is one; otherwise it is the opposite of
   the mean direction to its neighbors.
2. **Align.** Rotate the incoming unit (Rodrigues' formula, with an anti-parallel
   special case) so its head's leaving direction opposes the chain tail's leaving
   direction.
3. **Translate.** Put the head atom at the tail atom's position plus the bond length
   along the leaving direction. Bond length comes from a table (C–C 1.54, C–N 1.47,
   C–O 1.43, C–S 1.82, N–N 1.45, Si–O 1.63, Si–Si 2.33 Å) and otherwise the sum of
   covalent radii.
4. **Twist.** Rotate the unit about the new bond axis so the dihedral
   (tail reference atom, tail, head, head reference atom) equals the target. The
   reference atoms are the lowest-index heavy neighbors, so the choice is
   deterministic. The dihedral sign convention matches RDKit's.
5. **Bond and clean up.** Add a single bond, record the torsion for the optimizer,
   delete the leaving hydrogens at the end, remap all indices, renumber each unit's
   residue, and sanitize.

`torsion_deg` can be one angle for every junction or a list with one angle per
junction.

### 4. Force-field optimization (`chem/optimize.py`)

- MMFF94 is used when it has parameters for every atom; otherwise UFF. If neither
  does, optimization fails with a clear message.
- Each junction torsion gets a restraint with a ±5° flat-bottom window and a force
  constant of 1000, defined **relative to the dihedral's current value**. The
  builder has already set that value, and relative restraints avoid wraparound
  problems at the ±180° default.
- Minimization runs up to 500 iterations, and energy before and after is reported.

### 5. Post-build validation (`chem/validation.py`)

The optimized structure is checked and problems are reported as warnings:

- It must sanitize and be a single fragment.
- New bonds should be between 0.6× and 1.6× their expected length.
- Non-bonded atom pairs (excluding 1-2 and 1-3 neighbors) closer than 0.6× the sum
  of their covalent radii count as clashes.
- Each junction torsion more than 15° from target is flagged.
- An energy increase of more than 1 kcal/mol during optimization is flagged.

### 6. Linear length preset (`chem/shape_solver.py`)

"Straight" for a backbone does not mean every unit points the same way; an
all-trans alkane zigzags and is still fully extended. So the solver:

1. Probes a 3-unit build and checks that consecutive step lengths agree (within 3%),
   confirming the repeated-transform model holds for this monomer and head/tail.
2. Finds the **most extended torsion**: it scans 36 torsions over a 9-unit probe
   chain to maximize end-to-end distance per unit, then refines locally around the
   best.
3. Divides the requested length by that chord length to get the repeat count
   (`n_units - 1` steps for `n_units` units), caps it at the unit limit, and
   rejects the request if the achievable length is more than 25% off.

The reported length is a pre-optimization estimate, not a guarantee about the final
optimized structure.

### 7. PDB merging (`chem/merge.py`)

1. Validate both inputs (single conformer, finite coordinates, one connected
   component, valid linkage atom) and the 99,999-atom PDB limit.
2. Add missing hydrogens to working copies, with coordinates.
3. Choose a retained, non-collinear reference neighbor on each side for the junction
   torsion.
4. Place a copy of B on A's selected atom with the same align/translate/twist
   transform used by the chain builder.
5. Combine, add one single bond, remove the leaving hydrogens, remap indices, and
   sanitize.
6. Optimize the whole merged structure with one torsion restraint at the junction,
   using the same validation as above.

### 8. Saving (`callbacks/download.py`, `assets/save.js`)

Conversion happens on the server, only for the requested format. PDB and MOL come
from RDKit. MOL2 goes through OpenBabel and is then re-read to confirm the atom
count matches, so a silent conversion failure becomes an error instead of a corrupt
file.

### Additional solvers (library only)

`chem/shape_solver.py` also contains `solve_polygon` and `solve_circle`, which pick
a per-junction torsion schedule to approximate a regular polygon or ring by
root-finding the torsion that gives a target turn angle. They produce **open
approximations** that are not guaranteed to be flat or closed, and report measured
flatness and closure gap. The UI currently offers only the **Linear** preset, so
these are not reachable from the interface.

## Known limitations

- **Valence check is a heuristic**: `chem/validation.py` uses RDKit's *typical*
  valence per element, not a full extended-valence table. It can be conservative for
  hypervalent elements (S, P) that legitimately exceed their default valence (for
  example, a sulfone sulfur may be flagged as "full").
- **MOL2 output goes through OpenBabel**: RDKit can read `.mol2` but has no MOL2
  writer. Saving to `.mol2` re-derives SYBYL atom types for the *new* structure; it
  does not patch the original file's atom typing.
- **Hard limits**: repeat count is capped at 1,000 units, and estimated atom count
  at 100,000, before any homopolymer build is attempted. The merger also respects
  PDB's 99,999-atom serial limit, including added hydrogens.
- **Single fragment only**: a multi-fragment upload (for example a monomer plus a
  counter-ion or solvent) is rejected.
- **Merger accepts PDB only**, one model per file, up to 20 MB.
- **No reaction chemistry**: the tool never decides which atoms should react, and
  never changes bond orders or removes anything other than the leaving hydrogens.

## Setup (Windows, using `uv`)

Use Python 3.11 or newer. From the project directory, run:

```powershell
uv sync
```

This creates `.venv` and installs the application and test dependencies from
`pyproject.toml`, using the versions recorded in `uv.lock`. Dependencies: Dash,
RDKit, `openbabel-wheel` (a prebuilt OpenBabel with Python bindings, so no system
OpenBabel install or C++ toolchain is needed), and NumPy. The 3D viewer uses a
vendored copy of 3Dmol.js in `assets/`.

## Project layout

```
app.py          Dash app, layout, and client-side callbacks
callbacks/      Dash callbacks: upload, selection, shape, build, download, merge
chem/           Chemistry engine: io, validation, builder, optimize, shape_solver, merge, viewer
ui/             Layout helpers for the merger tab and Save controls
assets/         3Dmol.js, viewer.js (3D view and click bridge), save.js, style.css
tests/          pytest suite and optional Playwright browser checks
```

## Tests

```powershell
uv run python -m pytest tests\ -v
```

Tests cover file I/O round-trips (`test_io.py`), valence-feasibility edge cases
(`test_validation.py`), polymer construction geometry (`test_builder.py`, using
ethane as a pure geometry/algorithm check, not a claim about realistic polymer
chemistry), force-field optimization (`test_optimize.py`), and the shape solver
(`test_shape_solver.py`). Merger tests cover bond geometry and input validation
(`test_merge.py`) and real Dash upload, selection, merge, invalidation, and download
requests (`test_merge_callbacks.py`). Save tests check the exported molecular graph
in all three formats for both tabs, the default format, and error responses
(`test_download.py`).

The optional browser checks exercise real mouse gestures, head/tail picks, join
changes, downloads, and tab persistence:

```powershell
uv run --with playwright==1.62.0 playwright install chromium
uv run --with playwright==1.62.0 python tests/browser_merge_controls.py
uv run --with playwright==1.62.0 python tests/browser_save.py
```

The Save browser check uses real click events and server conversions, with a
controlled substitute for the native OS dialog and file handle. It verifies user
activation, format choice, cancellation, write errors, and download fallback.

The homopolymer tab's click-selection bridge (`assets/viewer.js` +
`dash_clientside.set_props`) is not covered by the pytest suite; it was verified
manually end-to-end in a Playwright-driven browser session (upload, click-select
head/tail, build and optimize, save).
