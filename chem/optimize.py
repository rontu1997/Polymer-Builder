"""Force-field optimization of a built polymer chain, restraining the backbone
torsions builder.py already constructed near the user's target angle so the
optimizer relaxes bond lengths/angles without discarding that choice.
"""

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem

from chem.validation import ValidationResult, validate_post_build

TORSION_TOLERANCE_DEG = 5.0
TORSION_FORCE_CONSTANT = 1000.0
DEFAULT_MAX_ITERS = 500


class OptimizationError(Exception):
    pass


@dataclass
class OptimizationReport:
    force_field_used: str
    fell_back_to_uff: bool
    energy_before: float
    energy_after: float
    minimize_converged: bool
    validation: ValidationResult


def _get_force_field(mol: Chem.Mol, preferred_ff: str):
    """Returns (force_field, name_used, fell_back). Falls back MMFF94 -> UFF when
    the preferred force field lacks parameters for some atom type present."""
    preferred_ff = preferred_ff.upper()

    if preferred_ff == "MMFF94" and AllChem.MMFFHasAllMoleculeParams(mol):
        props = AllChem.MMFFGetMoleculeProperties(mol)
        return AllChem.MMFFGetMoleculeForceField(mol, props), "MMFF94", False

    if AllChem.UFFHasAllMoleculeParams(mol):
        fell_back = preferred_ff == "MMFF94"
        return AllChem.UFFGetMoleculeForceField(mol), "UFF", fell_back

    raise OptimizationError(
        "Neither MMFF94 nor UFF have parameters for all atoms in this structure; "
        "cannot optimize."
    )


def optimize_geometry(
    mol: Chem.Mol,
    torsion_specs,
    new_bonds=None,
    ff: str = "MMFF94",
    torsion_tolerance_deg: float = TORSION_TOLERANCE_DEG,
    force_constant: float = TORSION_FORCE_CONSTANT,
    max_iters: int = DEFAULT_MAX_ITERS,
):
    """Optimize `mol` in place (on a copy) and return (optimized_mol, report).

    Torsion constraints are added *relative* to each dihedral's current value
    (which builder.py already placed near the target), with a small tolerance
    window -- this sidesteps the +-180 degree wraparound that an absolute-angle
    constraint would need to handle near the common 180 degree default target.
    """
    mol = Chem.Mol(mol)
    force_field, ff_used, fell_back = _get_force_field(mol, ff)
    force_field.Initialize()
    energy_before = force_field.CalcEnergy()

    constraint_method = getattr(
        force_field, "MMFFAddTorsionConstraint" if ff_used == "MMFF94" else "UFFAddTorsionConstraint"
    )
    for a1, a2, a3, a4, _target in torsion_specs:
        constraint_method(
            a1, a2, a3, a4, True, -torsion_tolerance_deg, torsion_tolerance_deg, force_constant
        )

    converged = force_field.Minimize(maxIts=max_iters) == 0
    energy_after = force_field.CalcEnergy()

    validation = validate_post_build(mol, new_bonds or [], torsion_specs, energy_before, energy_after)

    report = OptimizationReport(
        force_field_used=ff_used,
        fell_back_to_uff=fell_back,
        energy_before=energy_before,
        energy_after=energy_after,
        minimize_converged=converged,
        validation=validation,
    )
    return mol, report
