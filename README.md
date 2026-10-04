# drg: DRG reduction of combustion reaction models with detonation and PSR targets

A modified version of [pyMARS](https://github.com/Niemeyer-Research-Group/pyMARS) (Python Mechanism Automatic Reduction Software, Niemeyer Research Group) for reducing detailed combustion reaction models with the directed relation graph (DRG) method. Species are removed as long as the reduced model reproduces the detailed model within an error limit on a set of sampled targets.

## Changes from pyMARS

- **ZND detonation targets** (`InputZND`, `Simulation_znd`): Chapman-Jouguet speed from SDToolbox and the ZND induction length.
- **Perfectly stirred reactor (PSR) targets** (`InputPSR`, `Simulation_psr`): extinction residence time, found by following the S-curve.
- `get_plot_data()` in `sampling.py` exports the target values of the detailed and reduced models for comparison plots.
- Chemkin writer (`soln2ck.py`): third-body efficiencies are written only for reactions with a generic third body (`+ M`), not for reactions with an explicit collision partner.

## Files

| File | Description |
| --- | --- |
| `drg.py` | DRG interaction matrix, graph search from the target species, and the reduction loop (`run_drg`) |
| `sampling.py` | Target definitions (autoignition, laminar flame speed, ZND, PSR), parallel sampling, and error calculation |
| `simulation.py` | Simulation classes for each target type |
| `reduce_model.py` | Removes species and their reactions to write the reduced model |
| `soln2ck.py`, `soln2cti.py` | Write a Cantera solution in Chemkin or Cantera CTI format |
| `parse_yaml.py`, `drgep.py`, `pfa.py`, `sensitivity_analysis.py`, `tools.py` | Kept from pyMARS (input parsing, DRGEP and PFA methods, sensitivity analysis, format conversion); `pfa.py` still uses package-relative imports |

## Usage

Call `run_drg()` in `drg.py` from Python with the model file, the target conditions (lists of `InputIgnition`, `InputPSR`, `InputLaminarFlame`, `InputZND`), the error limit, and the target and retained species.

Dependencies: Cantera, NumPy, NetworkX, PyYAML, PyTables, [SDToolbox](https://shepherd.caltech.edu/EDL/PublicResources/sdt/) for the detonation targets, and the MATLAB Engine API for Python. The ZND induction length is computed by a MATLAB function `znd_solve_matlab`, which is not included in this repository.

## License

pyMARS is distributed under the MIT License, reproduced in `LICENSE-pyMARS`.
