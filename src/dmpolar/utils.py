import subprocess
from rdkit import Chem


def run_command(command: str, **kwargs) -> subprocess.CompletedProcess:
    command_list = command.split()
    return subprocess.run(command_list, **kwargs)


def generate_points(rdmol: Chem.rdmol.Mol, method: str = "CHLEGP"):
    anums = [atom.GetAtomicNum() for atom in rdmol.GetAtoms()]
    coords = rdmol.GetConformer().GetPositions() # coordinates in angstrom