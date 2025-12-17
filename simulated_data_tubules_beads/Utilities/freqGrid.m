function [fx, fy] = freqGrid(N, px_um)
    df = 1/(N*px_um);
    f = (-floor(N/2):ceil(N/2)-1)*df;
    [fx, fy] = meshgrid(f, f);
end