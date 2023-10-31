from dmpolar.wrappers.psi4 import check_psi4_installation, check_gdma_installation, compute_gdma
import numpy as np
import os


def test_psi4_check():
    import psi4
    psi4.set_memory('2000 MB')
    psi4.set_num_threads(os.cpu_count())
    energy = check_psi4_installation()
    np.testing.assert_allclose(energy, -76.0266327350922779, atol=1e-4)

def test_gdma_check():
    import psi4
    psi4.set_memory('2000 MB')
    psi4.set_num_threads(os.cpu_count())
    gdma = check_gdma_installation()
    print(gdma)
    assert len(gdma) == 3

def test_compute_gdma():
    from rdkit import Chem
    rdmol = Chem.MolFromMolFile("tests/data/mol1.mol", removeHs=False, sanitize=True)
    import psi4
    psi4.set_memory('2000 MB')
    psi4.set_num_threads(os.cpu_count())

    gdma = compute_gdma(rdmol, folder="tests/data/out/mol1")
    print(gdma)