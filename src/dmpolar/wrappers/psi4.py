import psi4
import sys
import numpy as np
from tempfile import TemporaryDirectory
from pathlib import Path
from rdkit import Chem


def convert_harmonic_to_cartesian(harmonic_info):
    ret = {}
    Q00, Q10, Q11c, Q11s, Q20, Q21c, Q21s, Q22c, Q22s = harmonic_info[:9]
    ret["c"] = Q00
    ret["mu_z"] = Q10
    ret["mu_x"] = Q11c
    ret["mu_y"] = Q11s
    ret["theta_xx"] = -0.5 * Q20 + 0.5 * np.sqrt(3.0) * Q22c
    ret["theta_yy"] = -0.5 * Q20 - 0.5 * np.sqrt(3.0) * Q22c
    ret["theta_zz"] = Q20
    ret["theta_xy"] = 0.5 * np.sqrt(3.0) * Q22s
    ret["theta_xz"] = 0.5 * np.sqrt(3.0) * Q21c
    ret["theta_yz"] = 0.5 * np.sqrt(3.0) * Q21s
    return ret


def read_gdma_result(gdma_matrix):
    ret = []
    for multipole in gdma_matrix:
        ret.append(convert_harmonic_to_cartesian(multipole))
    return ret

def compute_gdma(wfn):
    out = psi4.gdma(wfn)
    gdma_mat = out.variable("DMA DISTRIBUTED MULTIPOLES").np
    return read_gdma_result(gdma_mat)

def check_psi4_installation():
    with TemporaryDirectory() as dirname:
        psi4.set_output_file(f'{dirname}/psi4_check.log', False)

        h2o = psi4.geometry("""
        O
        H 1 0.96
        H 1 0.96 2 104.5
        """)

        energy = psi4.energy('scf/cc-pvdz', molecule=h2o)
        return energy

    
def compute_wfn(rdmol: Chem.Mol, pos=None, folder=None, method="hf/cc-pvdz"):
    charge = Chem.GetFormalCharge(rdmol)
    density = "SCF"
    if "mp2" in method.split("/")[0].lower() or "cc" in method.split("/")[0].lower():
        density = "MP2"
    if folder is not None:
        dirname = folder
        Path(dirname).mkdir(parents=True, exist_ok=True)
    else:
        dirname = TemporaryDirectory()
    psi4.set_output_file(f'{dirname}/gdma.log', False)
    geom = [f"{int(charge)} 1"]
    for atom in rdmol.GetAtoms():
        if pos is None:
            pos = rdmol.GetConformer().GetPositions() * 0.1
        geom.append(f"{atom.GetSymbol()} {pos[atom.GetIdx()][0]} {pos[atom.GetIdx()][1]} {pos[atom.GetIdx()][2]}")
    geom.append("units angstrom")
    geom.append("noreorient")
    geom.append("nocom")
    geom = "\n".join(geom)
    mol = psi4.geometry(geom)
    energy, wfn = psi4.energy(method, molecule=mol, return_wfn=True)
    return wfn

def compute_esp(wfn, points):
    with open("grid.dat", "w") as f:
        for xx, yy, zz in points:
            f.write(f"{xx:.8f} {yy:.8f} {zz:.8f}\n")
    ret = psi4.oeprop(wfn, "GRID_ESP")
    with open("grid_esp.dat", "r") as f:
        esp = [float(l.strip()) for l in f]
    return np.array(esp)