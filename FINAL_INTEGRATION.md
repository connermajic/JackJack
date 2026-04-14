# Distal Radius FEA Project - Final Integration Branch

This branch combines practical upgrades so you can review one final request.

## Included upgrades

1. **Python 2D support**
   - New `--mode 2d` path in `python_fea/radius_fea.py`
   - Slice-based finite element solve from converted volumes

2. **Python visualizations**
   - `--save-plots` to export summary PNGs
   - 3D displacement scatter with load target marker
   - 2D segmentation + displacement summary figure

3. **MATLAB visualization utility**
   - New `visualizeRadiusFEAResult.m` helper for mesh, displacement, stress, and load target display

## Recommended run sequence

### MATLAB

```matlab
addpath(pwd);
result = radiusFEA_run('path/to/scan.nii.gz');
visualizeRadiusFEAResult(result);
```

### Python 3D

```bash
python python_fea/radius_fea.py path/to/scan.nii.gz --mode 3d --save-plots --viz-prefix case3d
```

### Python 2D

```bash
python python_fea/radius_fea.py path/to/scan.nii.gz --mode 2d --slice-axis 2 --slice-index 120 --save-plots --viz-prefix case2d
```

## Acceptance strategy

You can treat this branch as the single final review/acceptance candidate.
