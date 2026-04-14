function p = defaultRadiusFEAParams()
%DEFAULTRADIUSFEAPARAMS Default parameters for distal-radius homogeneous FEA.
%
% Recommended defaults (you can override any field when calling radiusFEA_run):
% - Converted CT/MicroCT volumes are assumed isotropic voxels in mm unless
%   you set voxelSize_mm otherwise.
% - Bone is modeled as homogeneous isotropic linear elastic solid.
% - Long axis for "proximal/distal" is taken as the third image axis (Z).
%   If your scan is oriented differently, rotate the volume before export or
%   extend the code to pass a rotation matrix.

    p.voxelSize_mm = [1 1 1];          % [sx sy sz] physical size per voxel
    p.volumeOrigin_mm = [0 0 0];       % world coords of voxel (1,1,1) corner

    % Segmentation: Hounsfield-like CT assumes bone brighter than soft tissue.
    % For microCT or inverted contrast, set invertMask = true.
    p.threshold = [];                % empty => Otsu on masked region
    p.invertMask = false;
    p.morphOpenRadius = 1;           % voxels, 0 to disable
    p.fillHoles = true;

    % Surface / STL
    p.isovalue = 0.5;                % for binary mask in [0,1]
    p.stlFeatureAngle = 30;          % degrees; merges coplanar STL facets
    p.tempStlName = 'radius_segmented.stl';

    % Mesh
    p.meshHmax_mm = [];              % empty => 0.15 * min(bounding box extent)
    p.meshHmin_mm = [];
    p.meshHgrad = 1.5;

    % Material (homogeneous cortical-bone order of magnitude; tune to your study)
    p.youngModulus_Pa = 15e9;
    p.poissonRatio = 0.3;

    % Boundary / loading (global Cartesian coords in mm, aligned with volume axes)
    % Fixed: faces sampled near the "proximal" end (minimum Z by default).
    % Load:  traction patch on the "distal" end (maximum Z by default).
    p.axisForLength = 3;             % 1=X, 2=Y, 3=Z (default: Z = column 3)
    p.proximalBand_mm = [];          % empty => 3% of model extent along axis
    p.distalBand_mm = [];

    % Total force vector (N). Example: axial compression along -Z
    p.forceTotal_N = [0 0 -500];

    % Off-axis: offset of load patch center from distal-surface centroid, in the
    % plane orthogonal to the force direction (when possible), units mm.
    p.loadOffset_mm = [2 0 0];

    % Radius of the loaded patch on the distal surface (mm), measured geodesically
    % via Euclidean proxy on surface samples (adequate for small offsets).
    p.loadPatchRadius_mm = 3;

    % Outputs
    p.writeResultMat = true;
    p.resultMatFile = 'radiusFEA_result.mat';
    p.showPlots = true;
end
