# Distal Radius FEA (Python)

Python implementation of a homogeneous distal-radius FEA workflow from converted image volumes.

## What this does

- Reads converted scan volumes (`.mat`, `.nii`/`.nii.gz`, `.tif`/`.tiff`)
- Segments bone from intensities (Otsu by default)
- Builds a surface mesh from the segmentation (`marching_cubes`)
- Tetrahedralizes the surface mesh
- Solves linear-elastic homogeneous FEA
- Applies off-axis loading on distal nodes near an offset load target
- Reports load application location and displacement metrics

Primary requested output:

- `fea.load_point_target_mm` (force application target in mm)

## Requirements

- Python 3.10+
- Recommended: virtual environment
- Native build tools may be required for meshing dependencies on some systems

Install dependencies:

```bash
pip install -r python_fea/requirements.txt
```

## Files

- `python_fea/radius_fea.py` - main pipeline and CLI entry
- `python_fea/requirements.txt` - dependencies
- `python_fea/README.md` - this guide

## Quick start

From repository root:

```bash
python python_fea/radius_fea.py path/to/scan.nii.gz
```

Default outputs:

- `radius_fea_python_result.json`
- `radius_fea_python_result.npz`

## Typical run with parameters

```bash
python python_fea/radius_fea.py path/to/radius.nii.gz \
  --voxel 0.082 0.082 0.082 \
  --force 0 0 -800 \
  --offset 3 0 0 \
  --axis 2 \
  --out-json results/radius_case1.json \
  --out-npz results/radius_case1.npz
```

For `.mat` files, optionally specify the 3D variable name:

```bash
python python_fea/radius_fea.py path/to/radius.mat --mat-variable volume3d
```

## CLI arguments

- positional `image_path` - input volume path (`.mat`, `.nii/.nii.gz`, `.tif/.tiff`)
- `--out-json` - summary output JSON (default `radius_fea_python_result.json`)
- `--out-npz` - arrays output NPZ (default `radius_fea_python_result.npz`)
- `--mat-variable` - variable name inside `.mat` file (optional)
- `--voxel SX SY SZ` - voxel size in mm
- `--force FX FY FZ` - total force vector in N
- `--offset OX OY OZ` - off-axis load offset in mm
- `--axis {0,1,2}` - long-axis index (0=x, 1=y, 2=z)
- `--threshold` - segmentation threshold override (float)

## Coordinate and unit conventions

- Geometry input/outputs are in **mm**
- FEA solve uses **meters** internally for SI consistency
- Material model defaults:
  - `E = 15e9 Pa`
  - `nu = 0.30`

## Output structure

### JSON (`--out-json`)

Includes:

- input metadata (`input_path`, `volume_shape`)
- segmentation threshold used
- mesh size stats (`num_tet_nodes`, `num_tets`, etc.)
- parameters used
- FEA summary:
  - `load_point_target_mm`
  - `distal_centroid_mm`
  - `force_total_n`
  - `num_load_nodes`, `num_fixed_nodes`
  - `max_displacement_m`, `mean_displacement_m`

### NPZ (`--out-npz`)

Stores arrays:

- `nodes_mm`
- `tets`
- `displacements_m`
- `load_nodes`
- `fixed_nodes`

## Programmatic usage

```python
from python_fea.radius_fea import RadiusFEAParams, run_radius_fea_python

params = RadiusFEAParams(
    voxel_size_mm=(0.082, 0.082, 0.082),
    force_total_n=(0.0, 0.0, -500.0),
    load_offset_mm=(2.0, 0.0, 0.0),
    axis_for_length=2,
)

result = run_radius_fea_python(
    image_path="path/to/scan.nii.gz",
    params=params,
    output_json="result.json",
    output_npz="result.npz",
)
print(result["fea"]["load_point_target_mm"])
```

## Troubleshooting

- `Python was not found`:
  - Install Python and ensure `python` is on PATH.
- Meshing/tetrahedralization errors:
  - Check segmentation quality and scan orientation.
  - Try different threshold or preprocessing.
- No distal/proximal node sets found:
  - Verify `--axis` matches your long bone orientation.
- Unrealistic displacement magnitudes:
  - Verify voxel spacing and force units.
  - Tune material properties for your study.

## Notes

- This is a homogeneous isotropic model by design.
- Validate against your protocol before research or clinical interpretation.
