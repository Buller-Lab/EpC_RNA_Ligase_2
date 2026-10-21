# MD trajectories (input data — placeholder)

This folder is a **placeholder** for the Molecular Dynamics trajectory inputs
required to regenerate the MD post-analysis (README section 4). The trajectory
files are not included in the repository (too large); add them here before
running the MD post-analysis scripts.

Expected inputs:

- `production_0.h5` … `production_4.h5` — per-replicate MD trajectories
  (MDTraj/HDF5), consumed by `scripts/identify_top100_frames.py` to select the
  top frames by ligation-site↔catalytic-site distance.
- `productive_frames_replicate_0.pdb` … `productive_frames_replicate_5.pdb` —
  the productive-frame ensembles, consumed by
  `scripts/calculate_binding_energy.py` (PyRosetta) to compute binding energies.

Once the input files are in place, run (in the `rnalig` environment, which must
also have PyRosetta installed — see `environment.yml`):

```bash
cd molecular_analysis/MD_trajectories
python ../../scripts/identify_top100_frames.py
python ../../scripts/calculate_binding_energy.py
```
