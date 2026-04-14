function visualizeRadiusFEAResult(result, varargin)
%VISUALIZERADIUSFEARESULT Quick visualization helper for MATLAB results.
%
% visualizeRadiusFEAResult(result)
% visualizeRadiusFEAResult(result, 'showMesh', true, 'showLoadPoint', true)

    p = inputParser;
    addParameter(p, 'showMesh', true, @islogical);
    addParameter(p, 'showLoadPoint', true, @islogical);
    parse(p, varargin{:});

    model = result.model;
    R = result.solution;

    if p.Results.showMesh
        figure('Name', 'FE Mesh');
        pdemesh(model);
        title('Finite Element Mesh');
    end

    figure('Name', 'Displacement Magnitude (m)');
    pdeplot3D(model.Mesh, 'ColorMapData', R.Displacement.Magnitude);
    title('Displacement Magnitude');
    colorbar;

    if isfield(R, 'VonMisesStress')
        vm = R.VonMisesStress;
    elseif isfield(R, 'vonMisesStress')
        vm = R.vonMisesStress;
    else
        vm = [];
    end

    if ~isempty(vm)
        figure('Name', 'von Mises Stress (Pa)');
        pdeplot3D(model.Mesh, 'ColorMapData', vm);
        title('von Mises Stress');
        colorbar;
    end

    if p.Results.showLoadPoint && isfield(result, 'loadPointTarget_mm')
        lp = result.loadPointTarget_mm;
        figure('Name', 'Geometry + Load Target');
        pdegplot(model, 'FaceLabels', 'on', 'FaceAlpha', 0.35);
        hold on;
        plot3(lp(1) * 1e-3, lp(2) * 1e-3, lp(3) * 1e-3, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
        title('Faces and Off-axis Load Target (red x)');
        axis equal;
        hold off;
    end
end
