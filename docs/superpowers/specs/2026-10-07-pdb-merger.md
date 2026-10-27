# PDB Merger

Add a PDB Merger tab alongside the existing Homopolymer Builder. The user confirmed
that merging means joining one selected atom in each PDB with a chemical bond.
The follow-up specifies separate head/tail assignments for both polymers and
an explicit choice of which ends to join.

The merger accepts two single-component, single-model PDB structures, displays each in its own
clickable 3D viewer, highlights head (blue) and tail (red) atoms per input, and exposes a
junction torsion (-180 to 180 degrees) and MMFF94/UFF force-field choice. Merge &
Optimize adds missing hydrogens on working copies, builds one single bond,
removes one explicit leaving hydrogen per endpoint
when needed, aligns the second structure using the existing placement code, and
optimizes with a torsion restraint. The first input's original coordinates remain
unchanged during placement; optimization may relax both structures. Inputs are not
mutated. The internal result is a PDB with connectivity and distinct source chain
IDs; users can export PDB (default), MOL, or MOL2.

Label inputs Polymer A/linker A and Polymer B/linker B. Each has a Head/Tail
picking radio and stores distinct head_idx/tail_idx values. All four assignments
are required to enable merging. Join choices explicitly map A head to B head,
A tail to B tail, A head to B tail, or A tail to B head (default). Only the chosen
pair receives a new bond. Status text identifies the chosen roles and indices;
join-mode changes invalidate old results/downloads.

Left-drag rotates, right-drag translates the view, scroll zooms, and left click
selects an atom. Selection/role changes preserve the camera view. The click event
captures the picking role. Viewer gestures never change molecular coordinates.

Keep the homopolymer controls, editable repeat count, linear-only shape, uv setup,
and existing input state when switching tabs. This is a one-to-one merge, with no
additional repetition of either uploaded structure.

Reject malformed/non-PDB uploads, disconnected inputs, invalid linkage atoms,
missing/nonfinite coordinates, nonfinite torsion, and structures above the app's
existing atom limit. Return readable errors. PDB output also respects its five-digit
atom serial limit. Failed uploads clear their selection. Any changed input or
parameter invalidates the merger preview/download. Old viewer clicks cannot select
an atom in a newer upload. Preserve residue names/numbers, names, and internal
connectivity; assign chain A to the first input and chain B to the second.
Validate XYZ fields before RDKit parsing so blank fields cannot silently become
zeros. Preserve supplied CONECT records; infer proximity bonds only when no
connectivity records are present. Compare raw upload identities as well as parsed
state to prevent stale clicks, builds, or downloads during upload replacement.

Use the existing Python >=3.11, Dash, RDKit, NumPy, and OpenBabel dependencies. No
new runtime dependencies. Chemistry feasibility checks do not infer a reaction.

Verify molecular graph, bond length and torsion, explicit and implicit hydrogen
handling, source metadata, input immutability, real upload/select/merge/download
callbacks, invalidation, and tab switching. Run the full existing test suite.

## Save As follow-up

Both tabs expose PDB (default), MOL, and MOL2 exports. Save opens the native
browser location/filename picker directly from the click to retain user activation.
Only after the user chooses a file is its requested format converted on the server.
Each response carries a request ID and stale results are not written. Cancellation
does not trigger a fallback download. Browsers without the picker use an ordinary
download and receive instructions for enabling their location-prompt setting.
Report conversion/write failures, preserve existing merger signature checks, and
reuse the existing RDKit/OpenBabel exporters without new runtime dependencies.
