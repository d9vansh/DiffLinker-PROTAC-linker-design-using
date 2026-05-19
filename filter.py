from rdkit import Chem
from rdkit.Chem import QED, Crippen
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import RDConfig
import os, glob, shutil, sys

# ---- SA scorer ----
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

input_dir = "/home/rajnish/DiffLinker/DiffLinker/RUN/final_R_NA_P_1"
output_dir = "/home/rajnish/DiffLinker/DiffLinker/RUN/R1_f_NA_P"
os.makedirs(output_dir, exist_ok=True)

# ---- RELAXED thresholds (start loose) ----
QED_MIN = 0.10
SA_MAX = 7.0
LOGP_MIN, LOGP_MAX = -2, 5

stats = {
    "total": 0,
    "loaded": 0,
    "connected": 0,
    "sanitized": 0,
    "single_frag": 0,
    "descriptor_ok": 0,
    "passed": 0
}

for f in glob.glob(os.path.join(input_dir, "*.xyz")):
    stats["total"] += 1

    mol = Chem.MolFromXYZFile(f)
    if mol is None:
        continue
    stats["loaded"] += 1

    try:
        rdDetermineBonds.DetermineConnectivity(mol)
    except:
        continue
    stats["connected"] += 1

    # try sanitization, but don't die if it fails
    try:
        Chem.SanitizeMol(mol)
        stats["sanitized"] += 1
    except:
        pass  # keep going anyway

    if len(Chem.GetMolFrags(mol)) != 1:
        continue
    stats["single_frag"] += 1

    try:
        mol = Chem.AddHs(mol)  # IMPORTANT
        qed = QED.qed(mol)
        sa = sascorer.calculateScore(mol)
        logp = Crippen.MolLogP(mol)
        stats["descriptor_ok"] += 1
    except:
        continue

    if (qed >= QED_MIN and
        sa <= SA_MAX and
        LOGP_MIN <= logp <= LOGP_MAX):

        base = os.path.basename(f)
        new_name = f"Q{qed:.2f}_SA{sa:.2f}_logP{logp:.2f}_" + base
        shutil.copy(f, os.path.join(output_dir, new_name))
        stats["passed"] += 1

# ---- PRINT DEBUG ----
print("\n=== PIPELINE STATS ===")
for k, v in stats.items():
    print(f"{k}: {v}")