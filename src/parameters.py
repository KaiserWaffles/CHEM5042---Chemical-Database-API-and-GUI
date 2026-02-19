from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, QED

# structure for RDKit
@dataclass(frozen=True)
class Parameters:
    mw: float
    logp: float
    hbd: int
    hba: int
    # polar surface area
    tpsa: float
    # rotatble bonds
    rotb: int 
    # aromatic rings
    rings: int
    # heavy atom count
    hatom: int
    qed: float
    ro5_violations: int
    ro5_pass: int
    ro3_pass: int
    leadlike_pass: int
    veber_pass: int
    bioavail_pass: int

# drug design filters
def ro5_violations(mw: float, logp: float, hbd: int, hba: int) -> int:
    """Lipinski Ro5 violations count."""
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1
    return violations


def pass_ro5(violations: int, max_violations: int = 1) -> bool:
    """
    Ro5 pass threhold <= 1.
    compound still druglike with at most one violation
    """
    return violations <= max_violations


def pass_ro3(mw: float, logp: float, hbd: int, hba: int, hatom: int, rotb: int) -> bool:
    """Rule of 3; fragment-like filter."""
    return (mw <= 300) and (logp <= 3) and (hbd <= 3) and (hba <= 3) and (hatom < 20) and (rotb <= 3)


def pass_lead_like(mw: float, logp: float, hbd: int, hba: int, rotb: int, tpsa: float) -> bool:
    """Lead-like filter"""
    return (250 <= mw <= 350) and (logp <= 3) and (hbd <= 3) and (hba <= 6) and (rotb <= 7) and (tpsa <= 120)


def pass_veber(rotb: int, tpsa: float) -> bool:
    """Veber oral bioavailability filter."""
    return (rotb <= 10) and (tpsa <= 140)


def pass_bioavailability(ro5_ok: bool, veber_ok: bool) -> bool:
    """ bioavailability filter; combines Ro5 + Veber."""
    return ro5_ok and veber_ok