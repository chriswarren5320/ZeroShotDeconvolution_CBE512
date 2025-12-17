function I = drawCurvedTubuleControlled(N, widthPx, bendAmp)
% Draw a *controlled-curvature* curved tubule on an N×N canvas using insertShape.
%
% N        : image size (NxN)
% widthPx  : tubule thickness (line width in pixels)
% bendAmp  : bending amplitude in pixels (0 = perfectly straight)
%
% Idea:
%   - Start from a straight line (centerline) of length L along a random angle.
%   - Add a smooth lateral offset: offset(s) = bendAmp * sin(2*pi*s + phi)
%     where s in [0,1] is normalized arclength along the centerline.
%   - Rasterize as many short segments via insertShape.

    if nargin < 3
        bendAmp = 20;   % default moderate bending
    end

    % --- basic geometry parameters ---
    margin   = 20;          % keep tubule away from borders
    Lmin     = 0.8*N;       % min length of tubule
    Lmax     = 0.9*N;       % max length of tubule
    nSamples = 400;         % number of points along the curve

    % --- choose random center, length, and main orientation ---
    cx = randi([1+margin, N-margin]);
    cy = randi([1+margin, N-margin]);
    L  = Lmin + (Lmax - Lmin)*rand();      % tubule length in pixels
    theta = 2*pi*rand();                   % main axis angle

    ax = cos(theta);                       % main axis direction
    ay = sin(theta);
    px = -sin(theta);                      % perpendicular direction
    py =  cos(theta);

    % --- parameter along the centerline ---
    % t in [-0.5, 0.5] * L (centered at (cx,cy))
    s = linspace(0, 1, nSamples);          % normalized arclength
    t = (s - 0.5) * L;                     % signed distance along main axis (pixels)

    % --- controlled lateral offset (this controls curvature) ---
    % bendAmp is in pixels; smaller -> almost straight, larger -> more curved
    phi = 2*pi*rand();                     % random phase to vary the shape
    offset = bendAmp * sin(2*pi*s + phi);  % lateral displacement in pixels

    % --- final coordinates of the centerline ---
    xs = cx + t*ax + offset.*px;
    ys = cy + t*ay + offset.*py;

    % --- build segments for insertShape [x1 y1 x2 y2] ---
    segments = [];
    for i = 1:(numel(xs)-1)
        x1 = xs(i);   y1 = ys(i);
        x2 = xs(i+1); y2 = ys(i+1);
        % keep segments inside image (simple clipping)
        if  x1>=1 && x1<=N && y1>=1 && y1<=N && ...
            x2>=1 && x2<=N && y2>=1 && y2<=N
            segments(end+1,:) = [x1, y1, x2, y2]; %#ok<AGROW>
        end
    end

    % if no valid segments (very rare), return empty
    if isempty(segments)
        I = zeros(N,N,'single');
        return;
    end

    % --- rasterize with insertShape ---
    rgb = zeros(N, N, 3, 'uint8');
    rgb = insertShape(rgb, 'Line', segments, ...
                      'Color', 'white', ...
                      'LineWidth', widthPx, ...
                      'SmoothEdges', true);

    Ig = rgb2gray(rgb);     % uint8
    I  = im2single(Ig);     % [0,1] single
    % figure; imagesc(I);
    % small blur to soften jaggies (optional)
    I = imgaussfilt(I, 0.7);
    mx = max(I(:));
    if mx > 0
        I = I ./ mx;
    end
end
