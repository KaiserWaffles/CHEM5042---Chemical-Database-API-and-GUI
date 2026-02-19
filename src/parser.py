import io
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw, QED
from rdkit.Chem import rdDepictor

# Compounds heavier than this are excluded before DB insertion.
# Many bioactive natural products (e.g. macrolides, cyclic peptides) exceed
# Lipinski-style limits and would distort statistical summaries.
MW_MAX_THRESHOLD = 1000.0

# API
def parse_sdf(filepath: str) -> list[dict]:
    """
    Parse .SDF file and return a list of compound parameters dicts; Each valid compound is tidied by RDKit, properties calculated, and an optional molecular weight (MW) ceiling filter applied.  Molecules unable to be parsed or exceed MW_MAX_THRESHOLD skipped

    Returns: list[dict]
        One dict per compound, with KEYS matching the SQLite schema
    """
