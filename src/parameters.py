from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, QED, Lipinski, Crippen,  Descriptors

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

    # filter results functionality
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

def compute_properties(mol: Chem.Mol) -> Parameters:
    """
    Compute RDKit 2D descriptors and screening filters;
    All values are deterministic for a given structure.
    """
    mw = round(float(Descriptors.MolWt(mol)), 2)
    logp = round(float(Crippen.MolLogP(mol)), 2)
    hbd = int(Lipinski.NumHDonors(mol))
    hba = int(Lipinski.NumHAcceptors(mol))
    rotb = int(Lipinski.NumRotatableBonds(mol))
    rings = int(rdMolDescriptors.CalcNumAromaticRings(mol))
    tpsa = round(float(rdMolDescriptors.CalcTPSA(mol)), 2)
    hatom = int(mol.GetNumHeavyAtoms())       # FIX: no argument needed
    qed_val = round(float(QED.qed(mol)), 3)

    # --- Apply all filters ---
    v = ro5_violations(mw, logp, hbd, hba)
    ro5_ok = pass_ro5(v)
    ro3_ok = pass_ro3(mw, logp, hbd, hba, hatom, rotb)
    lead_ok = pass_lead_like(mw, logp, hbd, hba, rotb, tpsa)
    veber_ok = pass_veber(rotb, tpsa)
    bio_ok = pass_bioavailability(ro5_ok, veber_ok)

    return Parameters(
        mw=mw, logp=logp, hbd=hbd, hba=hba,
        tpsa=tpsa, rotb=rotb, rings=rings, hatom=hatom,
        qed=qed_val,
        ro5_violations=v,
        ro5_pass=int(ro5_ok),
        ro3_pass=int(ro3_ok),
        leadlike_pass=int(lead_ok),
        veber_pass=int(veber_ok),
        bioavail_pass=int(bio_ok),
    )

def db_row(p: Parameters) -> Dict[str, object]:
    """
    convert dataclass to a  dict compatible with compound properties columns.
    """
    return {
        "mw":   p.mw,
        "logp": p.logp,
        "hbd":  p.hbd,
        "hba":  p.hba,
        "tpsa": p.tpsa,
        "rotb": p.rotb,
        "rings":    p.rings,
        "hatom":    p.hatom,
        "qed":  p.qed,
        "ro5_violations":   p.ro5_violations,
        "ro5_pass": p.ro5_pass,
        "ro3_pass": p.ro3_pass,
        "leadlike_pass":    p.leadlike_pass,
        "veber_pass":   p.veber_pass,
        "bioavail_pass":    p.bioavail_pass,
    }

# Filter for GUI 
# Maps display name -> database column name for pre-computed filter flags
FILTER_COLUMNS = {
    "Lipinski Ro5":         "ro5_pass",
    "Rule of 3 (Fragment)": "ro3_pass",
    "Lead-likeness":        "leadlike_pass",
    "Veber Rules":          "veber_pass",
    "Bioavailability":      "bioavailability_pass",
}

FILTER_DESCRIPTIONS = {
    "Lipinski Ro5": (
        "Lipinski's Rule of 5 (Lipinski et al., 1997)\n"
        "At most 1 violation of:\n"
        "  MW  <= 500 Da\n"
        "  LogP <= 5\n"
        "  HBD <= 5 (sum of NH and OH)\n"
        "  HBA <= 10 (sum of N and O)"
    ),
    "Rule of 3 (Fragment)": (
        "Rule of 3 — Fragment-likeness (Congreve et al., 2003)\n"
        "ALL must be satisfied:\n"
        "  MW    <= 300 Da\n"
        "  LogP  <= 3\n"
        "  HBD   <= 3\n"
        "  HBA   <= 3\n"
        "  Heavy atoms < 20\n"
        "  Rotatable bonds <= 3"
    ),
    "Lead-likeness": (
        "Lead-likeness Filter\n"
        "ALL must be satisfied:\n"
        "  250 <= MW <= 350 Da\n"
        "  LogP  <= 3\n"
        "  HBD   <= 3\n"
        "  HBA   <= 6\n"
        "  Rotatable bonds <= 7\n"
        "  TPSA  <= 120 A^2"
    ),
    "Veber Rules": (
        "Veber's Rules (Veber et al., 2002)\n"
        "BOTH must be satisfied:\n"
        "  Rotatable bonds <= 10\n"
        "  TPSA <= 140 A^2"
    ),
    "Bioavailability": (
        "Bioavailability Filter (Ro5 + Veber combined)\n"
        "Must pass BOTH:\n"
        "  Lipinski Ro5 (<=1 violation)\n"
        "  Veber Rules (RotB<=10, TPSA<=140)"
    ),
}