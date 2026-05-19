from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
import glob, shutil, os

input_dir = "/home/rajnish/DiffLinker/DiffLinker/RUN/run_zinc_noanc_nopock_linksiz"
output_dir = "/home/rajnish/DiffLinker/DiffLinker/RUN/final_R_NA_NP_link"
prefix = "2"

os.makedirs(output_dir, exist_ok=True)

for f in sorted(glob.glob(os.path.join(input_dir, "*.xyz"))):
    mol = Chem.MolFromXYZFile(f)

    if mol is None:
        continue

    try:
        rdDetermineBonds.DetermineConnectivity(mol)  # <-- key step
    except:
        continue

    if len(Chem.GetMolFrags(mol)) == 1:
        base_name = os.path.basename(f)
        new_name = prefix + base_name
        shutil.copy(f, os.path.join(output_dir, new_name))

print("Done.")