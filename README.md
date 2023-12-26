# DMPolar
Differentiable molecular polarizable force field builder based on DMFF. Can automatically generate MPID type force field files for OpenMM with MPID plugin.

## Installation

### Install DMFF

DMFF is a python package. Please follow the instructions in [DMFF](https://github.com/deepmodeling/DMFF).

### Install requirements
```bash
mamba create -n dmpolar python=3.10 psi4 gdma pygdma openmm rdkit scipy numba xtb-python crest -c conda-forge
```

### Install DMPolar
```bash
# Move to the project root directory
pip install .
```

## Usage

### Example input file
```yaml
```

### Run DMPolar force field builder

```bash
dmpolar -i input.sdf -n MOL --config config.yaml
```