# Distal Radius FEA (MATLAB)

Homogeneous finite element analysis pipeline for distal-radius scans from converted volume images.

## What this does

- Loads converted scan volumes (`.mat`, `.nii`/`.nii.gz`, `.tif`/`.tiff`)
- Segments bone from intensity data
- Builds a surface mesh (`.stl`) from the segmentation
- Runs 3D homogeneous linear-elastic FEA using MATLAB PDE Toolbox
- Applies off-axis loading on the distal region
- Reports force application location (`loadPointTarget_mm`) and solution fields

## Requirements

- MATLAB R2023a+ recommended
- Partial Differential Equation Toolbox (**required**)
- Image Processing Toolbox (recommended for best segmentation)

## Files

- `radiusFEA_run.m` - end-to-end entry point
- `defaultRadiusFEAParams.m` - default settings
- `loadConvertedVolume.m` - input loading
- `segmentBoneVolume.m` - segmentation
- `binaryVolumeToSTL.m` - mask to STL
- `structuralRadiusFEA.m` - model setup + solve
- helper files: `findFacesInAxisBand.m`, `findFacesNearPoint.m`, `faceAreasFromGeometry.m`, `projectOffsetOrthogonalToForce.m`

## Quick start

1. Open MATLAB in this repository folder.
2. Add folder to path:

```matlab
addpath(pwd);
```

3. Run with defaults:

```matlab
result = radiusFEA_run('path/to/your/scan.nii.gz');
```

4. Read load location:

```matlab
result.loadPointTarget_mm
```

## Typical custom run

```matlab
params = struct( ...
    'voxelSize_mm', [0.082 0.082 0.082], ...
    'forceTotal_N', [0 0 -800], ...
    'loadOffset_mm', [3 0 0], ...
    'axisForLength', 3, ...
    'threshold', [] ...
);

result = radiusFEA_run('path/to/radius.mat', params);
```

## Key parameters

Edit in `defaultRadiusFEAParams.m` or pass as overrides.

- `voxelSize_mm` - physical voxel size
- `volumeOrigin_mm` - world origin of voxel grid
- `threshold` - segmentation threshold (`[]` uses Otsu)
- `invertMask` - invert contrast if needed
- `youngModulus_Pa`, `poissonRatio` - homogeneous material model
- `axisForLength` - 1/2/3 for X/Y/Z long axis
- `forceTotal_N` - total applied force vector
- `loadOffset_mm` - off-axis offset from distal centroid
- `loadPatchRadius_mm` - radius of load patch
- `proximalBand_mm`, `distalBand_mm` - BC/load band thickness
- `resultMatFile` - output MAT file name

## Outputs

`result` includes:

- `loadPointTarget_mm` - target off-axis load location (primary requested output)
- `forceTotal_N` - applied total force
- `surfaceTraction_Pa` - equivalent traction vector used on load faces
- `faceIdsLoad` - loaded face IDs
- `faceIdsProximal` - fixed face IDs
- `solution` - PDE solution object (displacement/stress)
- `meshStats` - node/element counts

If `writeResultMat = true`, output MAT file is saved (default `radiusFEA_result.mat`).

## Coordinate and unit conventions

- Inputs from image geometry are in **mm**
- Internal solve is scaled to **meters** for SI consistency with `E` in Pa and force in N
- `loadPointTarget_mm` is reported in mm in the same axis convention as your volume

## Troubleshooting

- **No faces found for proximal/distal bands**
  - Check `axisForLength` and scan orientation
  - Increase `proximalBand_mm` / `distalBand_mm`
- **Empty surface or zero-area load patch**
  - Adjust `threshold`, `invertMask`, and segmentation settings
- **Unrealistic stiffness/displacement**
  - Verify voxel size and material properties (`youngModulus_Pa`)

## Notes

- This is a homogeneous model by design (single isotropic material).
- Validate against your lab loading protocol before using for decision-making.
