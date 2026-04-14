function off = projectOffsetOrthogonalToForce(loadOffset_m, force_N)
%PROJECTOFFSETORTHOGONALTOFORCE Remove component of offset along force direction (same units as geometry).

    f = force_N(:);
    nf = norm(f);
    if nf < eps
        error('projectOffsetOrthogonalToForce:ZeroForce', 'forceTotal_N must be non-zero.');
    end
    u = f / nf;
    o = loadOffset_m(:);
    o = o - dot(o, u) * u;
    off = o.';
end
