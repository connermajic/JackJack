function [mask, info] = segmentBoneVolume(V, p)
%SEGMENTBONEVOLUME Simple bone segmentation for homogeneous FEA mask.

    info = struct('threshold', [], 'usedOtsu', false);

    Vn = V;
    if p.invertMask
        Vn = -Vn;
    end

    if isempty(p.threshold)
        if exist('graythresh', 'file') == 2
            t = graythresh(rescaleTo01(Vn));
            thr = min(Vn(:)) + t * (max(Vn(:)) - min(Vn(:)));
            info.usedOtsu = true;
        else
            thr = median(Vn(:));
            info.usedOtsu = false;
        end
    else
        thr = p.threshold;
    end

    info.threshold = thr;
    mask = Vn >= thr;

    if p.fillHoles
        mask = imfill(mask, 'holes');
    end

    if p.morphOpenRadius > 0
        se = strel('sphere', p.morphOpenRadius);
        mask = imopen(mask, se);
    end
end

function u = rescaleTo01(V)
    a = min(V(:));
    b = max(V(:));
    if b <= a
        u = zeros(size(V));
    else
        u = (V - a) ./ (b - a);
    end
end
