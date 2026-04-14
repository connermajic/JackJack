# Python Distal Radius FEA

This branch adds a Python equivalent of the MATLAB homogeneous distal-radius workflow.

## What it does

- Loads converted volume data (`.mat`, `.nii`/`.nii.gz`, `.tif`/`.tiff`)
- Segments bone with Otsu threshold (or user threshold)
- Extracts a surface mesh from the segmented volume
- Tetrahedralizes the surface
- Solves homogeneous linear elasticity (small strain)
- Applies off-axis load by selecting nodes in a distal patch around an offset target
- Outputs the load target location and displacement metrics

## Install

```bash
pip install -r python_fea/requirements.txt
```

## Run

```bash
python python_fea/radius_fea.py path/to/volume.nii.gz \
  --voxel 0.082 0.082 0.082 \
  --force 0 0 -500 \
  --offset 2 0 0 \
  --axis 2
```

Outputs:

- `radius_fea_python_result.json`
- `radius_fea_python_result.npz`

Main reported off-axis load location:

- `fea.load_point_target_mm`
