# Distal Radius FEA (Python)

This branch now supports **both 3D and 2D** homogeneous FEA, plus optional visualizations.

## Highlights

- 3D workflow: segmentation -> surface mesh -> tetrahedral mesh -> linear elasticity solve
- 2D workflow: select slice -> 2D mesh -> linear elasticity solve
- Off-axis loading in both modes with reported target location
- Optional PNG visualizations (`--save-plots`)

## Install

```bash
pip install -r python_fea/requirements.txt
```

## Quick run (3D)

```bash
python python_fea/radius_fea.py path/to/scan.nii.gz --mode 3d --save-plots
```

## Quick run (2D)

```bash
python python_fea/radius_fea.py path/to/scan.nii.gz --mode 2d --slice-axis 2 --save-plots
```

## Common options

- `--force FX FY FZ` total force vector in N
- `--offset OX OY OZ` off-axis load offset in mm
- `--axis {0,1,2}` long axis index for proximal/distal
- `--voxel SX SY SZ` voxel size in mm
- `--threshold T` segmentation threshold override
- `--slice-axis {0,1,2}` slice axis for 2D mode
- `--slice-index N` explicit slice index for 2D mode
- `--mesh-step-2d N` coarse/fine control for 2D mesh (smaller = finer)
- `--save-plots` save visualization images
- `--viz-prefix NAME` prefix for output plot filenames

## Output files

- JSON summary (`--out-json`, default `radius_fea_python_result.json`)
- NPZ arrays (`--out-npz`, default `radius_fea_python_result.npz`)
- Optional PNG plots when `--save-plots` is set

Key output field:

- `fea.load_point_target_mm`

## Notes

- 2D mode uses in-plane components of force and offset.
- Units are SI-consistent in solver (geometry from mm is scaled to meters).
