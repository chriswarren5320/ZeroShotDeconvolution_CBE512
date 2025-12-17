function I = drawCurvedTubule(N, widthPx, varargin)
% DRAWCURVEDTUBULE Draw a single curved tubule on an N×N canvas.
% The centerline is a cubic spline through random control points,
% then rasterized as many short line segments via insertShape.
%
% Usage:
%   I = drawCurvedTubule(N, widthPx);
%   I = drawCurvedTubule(N, widthPx, 'K',5,'Margin',10);
%
% Inputs:
%   N        - image size (NxN)
%   widthPx  - line width (tubule thickness in pixels)
%
% Optional name-value parameters:
%   'K'       - number of spline control points (>=3), default 4
%   'Margin'  - margin from borders for control points, default 8
%   'nSamples'- number of samples along spline, default 400
%   'SmoothSigma' - Gaussian blur sigma to soften edges, default 0.7
%
% Output:
%   I - single, N×N, intensity in [0,1]

    % ----- Parse optional parameters -----
    p = inputParser;
    addParameter(p, 'K', 8);
    addParameter(p, 'Margin', 50);
    addParameter(p, 'nSamples', 400);
    addParameter(p, 'SmoothSigma', 0.7);
    parse(p, varargin{:});

    K           = p.Results.K;
    margin      = p.Results.Margin;
    nSamples    = p.Results.nSamples;
    smoothSigma = p.Results.SmoothSigma;

    % ----- Random control points (keep away from borders) -----
    % K control points with random positions
    xr = randi([1+margin, N-margin], [1 K]);
    yr = randi([1+margin, N-margin], [1 K]);

    % Use a monotonic "parameter" along the chain to avoid self-crossing
    L = [0, cumsum(sqrt(diff(xr).^2 + diff(yr).^2))];
    if L(end) == 0
        L(end) = 1; % avoid division by zero if all control points happen to coincide
    end
    L = L / L(end);  % normalize to [0,1]
    tSample = linspace(0, 1, nSamples);

    % Cubic spline interpolation for x(t), y(t)
    xs = spline(L, xr, tSample);
    ys = spline(L, yr, tSample);

    % ----- Build polyline segments for insertShape -----
    % Each row: [x1 y1 x2 y2]
    xs_round = xs;
    ys_round = ys;

    % Optionally clip to image bounds for stability
    % (we still build segments only if both endpoints are inside)
    segments = [];
    for i = 1:(numel(xs_round)-1)
        x1 = xs_round(i);
        y1 = ys_round(i);
        x2 = xs_round(i+1);
        y2 = ys_round(i+1);

        % keep only segments whose endpoints are inside the image
        if  x1>=1 && x1<=N && y1>=1 && y1<=N && ...
            x2>=1 && x2<=N && y2>=1 && y2<=N

            segments(end+1,:) = [x1, y1, x2, y2]; %#ok<AGROW>
        end
    end

    % If no valid segments (rare edge case), return empty image
    if isempty(segments)
        I = zeros(N,N,'single');
        return;
    end

    % ----- Rasterize with insertShape -----
    % insertShape requires uint8/double; we use uint8 RGB canvas
    rgb = zeros(N, N, 3, 'uint8');
    rgb = insertShape(rgb, 'Line', segments, ...
                      'Color', 'white', ...
                      'LineWidth', widthPx, ...
                      'SmoothEdges', true);

    % Convert to grayscale single in [0,1]
    Ig = rgb2gray(rgb);          % uint8
    I  = im2single(Ig);          % convert to single 0..1
    figure;imagesc(I)
    % Optionally soften edges a bit
    if smoothSigma > 0
        I = imgaussfilt(I, smoothSigma);
    end

    % Normalize to [0,1] just in case
    mx = max(I(:));
    if mx > 0
        I = I ./ mx;
    end
end
