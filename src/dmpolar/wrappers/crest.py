from rdkit import Chem
import subprocess
import numpy as np
import multiprocessing as mp
from ..io import read_xyz


def crest_conf_search(rdmol: Chem.rdchem.Mol, max_e: float = 3.0):
    nproc = mp.cpu_count()
    with open("struct.xyz", "w") as f:
        geom = []
        pos = rdmol.GetConformer().GetPositions()
        for atom in rdmol.GetAtoms():
            geom.append(
                f"{atom.GetSymbol()} {pos[atom.GetIdx()][0]} {pos[atom.GetIdx()][1]} {pos[atom.GetIdx()][2]}\n")
        f.write(f"{len(geom)}\n\n")
    subprocess.run(
        ["crest", "struct.xyz", "--gfn2", "-T", f"{nproc}", "--mquick"]
    )
    coords = read_xyz("crest_conformers.xyz")
    with open("crest.energies", "r") as f:
        ene = np.array([float(i.strip().split()[1]) for i in f])
    coords = coords[ene < max_e]
    return coords