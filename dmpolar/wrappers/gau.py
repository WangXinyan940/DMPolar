import subprocess
import tempfile
from pathlib import Path


def check_gaussian_installation() -> str:
    # try to run simple qm job to check installation
    with tempfile.TemporaryDirectory() as tempdir:
        temp_path = Path(tempdir) / "test.com"
        temp_out = Path(tempdir) / "test.log"
        with open(temp_path, 'w') as f:
            f.write("%nproc=1\n")
            f.write("%mem=1000MW\n")
            f.write("#p HF 3-21G sp\n")
            f.write("\n")
            f.write("Title\n")
            f.write("\n")
            f.write("0 1\n")
            f.write("O    0.0    0.0    0.0\nH    1.0    0.0    0.0\nH    0.0    1.0    0.0\n\n\n\n")
        subprocess.run(
            ['g16', 'test.com'],
            cwd=tempdir,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        if temp_out.exists():
            return "g16"
        
        subprocess.run(
            ['g09', 'test.com'],
            cwd=tempdir,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        if temp_out.exists():
            return "g09"
        
        return "none"