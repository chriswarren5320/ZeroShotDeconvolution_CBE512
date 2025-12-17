function [fig, ax] = plot_10(ori, RL, Sp, ZS, ZS_R, x1, y1, x2, y2, W, varargin)
% PLOT_10  Five full images + five W×W zooms, with one consistent display range.
%   [fig, ax] = plot_10(ori, RL, Sp, ZS, ZS_R, x1, y1, x2, y2, W, ...)
%
% Inputs (each of ori, RL, Sp, ZS, ZS_R can be an array or a file path):
%   x1,y1,x2,y2 : endpoints of the line (pixels; x=col, y=row; 1-based)
%   W           : zoom box width in pixels
%
% Name-Value options:
%   'Mode'        : 'auto_all' (default) | 'auto_ori' | 'manual'
%   'Percentiles' : [pLow pHigh] for auto modes (default [0.5 99.5])
%   'Range'       : [vmin vmax] used only if Mode='manual'
%
% Example (fixed range 0–100):
%   plot_10(A, B, C, D, E, x1, y1, x2, y2, 128, 'Mode','manual','Range',[0 100]);

    p = inputParser;
    addParameter(p, 'Mode', 'auto_all', @(s)ischar(s)||isstring(s));
    addParameter(p, 'Percentiles', [0.5 99.5], @(v)isnumeric(v)&&numel(v)==2);
    addParameter(p, 'Range', [], @(v)isnumeric(v)&&numel(v)==2);
    parse(p, varargin{:});
    Mode = lower(string(p.Results.Mode));
    pr   = p.Results.Percentiles;
    Rng  = p.Results.Range;

    % --- load images as 2D (no normalization) ---
    imgs = {ori, RL, Sp, ZS, ZS_R};
    for k = 1:5
        imgs{k} = i_load_2d(imgs{k});
    end
    titles = {'Original','RL','Sparse','Net','Net with Hessian'};

    % --- compute one global display range ---
    switch Mode
        case "manual"
            if isempty(Rng), error('Mode=manual requires Range=[vmin vmax].'); end
            vmin = Rng(1); vmax = Rng(2);
        case "auto_ori"
            data = double(imgs{1}(:));
            vmin = prctile(data, pr(1)); vmax = prctile(data, pr(2));
        otherwise % "auto_all"
            data = [];
            for k = 1:5, data = [data; double(imgs{k}(:))]; end %#ok<AGROW>
            vmin = prctile(data, pr(1)); vmax = prctile(data, pr(2));
    end
    if ~(vmax > vmin)
        vmin = min(data); vmax = max(data);
        if ~(vmax > vmin), vmax = vmin + 1; end
    end

    % --- layout ---
    cx = 0.5*(x1 + x2); cy = 0.5*(y1 + y2); W = round(W);
    fig = figure('Color','w');
    tl = tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
    ax = gobjects(2,5);

    % --- top row: full images ---
    for i = 1:5
        ax(1,i) = nexttile(tl, i);
        imagesc(imgs{i}); axis image off; colormap(ax(1,i), gray);
        set(ax(1,i),'YDir','reverse'); title(titles{i});
        caxis(ax(1,i), [vmin vmax]);
        hold(ax(1,i), 'on');
        plot(ax(1,i), [x1 x2], [y1 y2], 'LineWidth', 1,'Color', 'r');
        % zoom box
        x0 = floor(cx - W/2); y0 = floor(cy - W/2);
        rectangle(ax(1,i), 'Position', [x0 y0 W W], 'EdgeColor','w', 'LineStyle','-');
        hold(ax(1,i), 'off');
    end

    % --- bottom row: zooms ---
    for i = 1:5
        [crop, x0, y0] = i_crop_center_with_pad(imgs{i}, cx, cy, W);
        ax(2,i) = nexttile(tl, 5+i);
        imagesc(crop); axis image off; colormap(ax(2,i), gray);
        set(ax(2,i),'YDir','reverse'); % title(['Zoom: ' titles{i}]);
        caxis(ax(2,i), [vmin vmax]);
        hold(ax(2,i), 'on');
        plot(ax(2,i), [x1 - x0, x2 - x0], [y1 - y0, y2 - y0], 'LineWidth', 1.2,'Color', 'r');
        hold(ax(2,i), 'off');
    end

    % annotate the used display range
    annotation(fig,'textbox',[0.75 0.01 0.24 0.05], 'String', ...
        sprintf('range=[%.3g, %.3g]  (%s, p=[%.2g, %.2g])', vmin, vmax, Mode, pr(1), pr(2)), ...
        'EdgeColor','none','HorizontalAlignment','right','FontSize',8);
end

% ----- helpers -----
function I = i_load_2d(x)
    % Load from path or accept array; convert to grayscale if RGB; pick mid-slice if stack.
    if ischar(x) || (isstring(x) && isscalar(x))
        path = char(x);
        info = imfinfo(path);
        if numel(info) > 1
            I = imread(path, ceil(numel(info)/2)); % middle slice
        else
            I = imread(path);
        end
    else
        I = x;
    end
    if ndims(I) == 3
        if size(I,3) == 3 || size(I,3) == 4
            if size(I,3) == 4, I = I(:,:,1:3); end
            I = rgb2gray(I);
        else
            I = I(:,:, ceil(size(I,3)/2)); % middle slice of volume
        end
    end
    I = double(I); % keep raw scale; display uses caxis
end

function [crop, x0, y0] = i_crop_center_with_pad(I, cx, cy, W)
    % Crop W×W around (cx,cy) with replicate padding. Returns crop and origin (x0,y0) in original coords.
    pad = ceil(W/2) + 2;
    Ip  = padarray(I, [pad pad], 'replicate', 'both');
    cxp = cx + pad; cyp = cy + pad;
    x0p = floor(cxp - W/2); y0p = floor(cyp - W/2);
    crop = Ip(y0p:(y0p+W-1), x0p:(x0p+W-1));
    x0   = x0p - pad; y0 = y0p - pad; % origin in original coordinates
end

function [d, P] = plot_profiles(ori, RL, Sp, ZS, ZS_R, x1, y1, x2, y2, N, varargin)
% PLOT_PROFILES  Plot intensity profiles for five images along the same line.
%   [d, P] = plot_profiles(ori, RL, Sp, ZS, ZS_R, x1, y1, x2, y2, N)
%       ori, RL, Sp, ZS, ZS_R : image arrays or file paths
%       (x1,y1)->(x2,y2) : endpoints (pixels; x=col, y=row; 1-based)
%       N : number of samples along the line (e.g., 512; if empty -> auto)
%
%   Returns:
%       d : 1×N distances (pixels) from (x1,y1)
%       P : N×5 matrix of profiles [ori RL Sp ZS ZS_R]
%
%   Options (Name-Value):
%       'Normalize' (false) : if true, each profile is scaled to [0,1] for display

    p = inputParser;
    addParameter(p, 'Normalize', false, @(b)islogical(b)||ismember(b,[0,1]));
    parse(p, varargin{:});
    doNorm = p.Results.Normalize;

    imgs = {ori, RL, Sp, ZS, ZS_R};
    names = {'Original','RL','Sparse','Net','Net with Hessian'};
    for k = 1:5
        imgs{k} = i_load_2d(imgs{k});        % grayscale 2D, raw scale
    end

    % clamp coordinates to be inside image 1 (you can change target if needed)
    [H,W] = size(imgs{1});
    x1 = min(max(x1,1),W); x2 = min(max(x2,1),W);
    y1 = min(max(y1,1),H); y2 = min(max(y2,1),H);

    % number of samples along the line
    if nargin < 11 || isempty(N)
        N = max(2, round(hypot(x2-x1, y2-y1))+1);   % ~1 sample per pixel
    end

    % sample profiles
    P = zeros(N,5);
    for k = 1:5
        % If images differ in size, clamp coords to each image
        [Hk, Wk] = size(imgs{k});
        x1k = min(max(x1,1),Wk); x2k = min(max(x2,1),Wk);
        y1k = min(max(y1,1),Hk); y2k = min(max(y2,1),Hk);
        pk = improfile(imgs{k}, [x1k x2k], [y1k y2k], N, 'bilinear');
        pk = pk(:);
        if doNorm
            mn = min(pk); mx = max(pk); 
            if mx>mn, pk = (pk-mn)/(mx-mn); end
        end
        P(:,k) = pk;
    end

    % distance axis in pixels (based on the requested line length)
    L = hypot(x2-x1, y2-y1);
    d = linspace(0, L, N);

    % plot
    figure('Color','w');
    plot(d, P, 'LineWidth', 1.4);
    grid on; box on;
    xlabel('Distance (pixels)');
    ylabel('Intensity (a.u.)');
    legend(names, 'Location','best');
    title('Line intensity profiles');
end

function [fig, ax] = plot_6(A, B, C, x1, y1, x2, y2, W)
% 3 images on top + 3 W×W zooms below.
% Red line between (x1,y1)-(x2,y2). White dashed zoom box.
    I = {i_load2d(A), i_load2d(B), i_load2d(C)};
    cx = 0.5*(x1+x2);  cy = 0.5*(y1+y2);  W = round(W);

    fig = figure('Color','w');  ax = gobjects(2,3);

    % --- FULL VIEWS ---
    for k = 1:3
        [H,Wimg] = size(I{k});
        ax(1,k) = subplot(2,3,k);
        imagesc(I{k}); axis image; axis ij; axis off; colormap(gray); hold on;
        % white dashed zoom box (clamped)
        bx = max(1, min(round(cx - W/2), Wimg - W + 1));
        by = max(1, min(round(cy - W/2), H    - W + 1));
        rectangle('Position',[bx by W W], 'EdgeColor','w', 'LineStyle','-', 'LineWidth',1);
        % thick red line (clamped)
        x1c = min(max(x1,1),Wimg);  x2c = min(max(x2,1),Wimg);
        y1c = min(max(y1,1),H   );  y2c = min(max(y2,1),H   );
        line([x1c x2c],[y1c y2c],'Color',[1 0 0],'LineWidth',1.2);
        hold off; 
        if k == 1 
            title(sprintf('Original'));
        end
        if k == 2
            title(sprintf('RL'));
        end
        if k == 3
            title(sprintf('Deconvolved'));
        end
    end

    % --- ZOOMS ---
    for k = 1:3
        [H,Wimg] = size(I{k});
        x0 = max(1, min(round(cx - W/2), Wimg - W + 1));
        y0 = max(1, min(round(cy - W/2), H    - W + 1));
        crop = I{k}(y0:y0+W-1, x0:x0+W-1);

        ax(2,k) = subplot(2,3,3+k);
        imagesc(crop); axis image; axis ij; axis off; colormap(gray); hold on;
        % same red line in crop coordinates
        line([x1 - x0 + 1, x2 - x0 + 1], [y1 - y0 + 1, y2 - y0 + 1], ...
             'Color',[1 0 0],'LineWidth',1.2);
        hold off; 
    end
end

% ---- helper: load as 2D double, middle slice if stack, grayscale if RGB ----
function I = i_load2d(x)
    if ischar(x) || (isstring(x) && isscalar(x))
        info = imfinfo(char(x));
        if numel(info) > 1, I = imread(char(x), ceil(numel(info)/2));
        else,               I = imread(char(x));
        end
    else
        I = x;
    end
    if ndims(I)==3, I = mean(double(I),3); else, I = double(I); end
end

function [d, P] = plot_profiles3(A, B, C, x1, y1, x2, y2, N)
% PLOT_PROFILES3  Plot intensity profiles for three images along one line.
%   [d, P] = plot_profiles3(A,B,C, x1,y1,x2,y2, N)
% Inputs:
%   A,B,C : images (array or file path). RGB/stack OK.
%   (x1,y1)->(x2,y2) : endpoints (pixels; x=col, y=row; 1-based)
%   N : number of samples along the line (e.g., 512). If empty/omitted,
%       it defaults to ~1 sample/pixel along the segment.
%
% Outputs:
%   d : 1×N distance (pixels) from (x1,y1)
%   P : N×3 profiles [A B C] (raw intensities)

    if nargin < 9 || isempty(N)
        N = max(2, round(hypot(x2-x1, y2-y1)) + 1);
    end

    % load/grayscale/mid-slice → double
    imgs = {i_load2d(A), i_load2d(B), i_load2d(C)};

    % sample profiles with simple clamping to each image
    P = zeros(N,3);
    for k = 1:3
        [H,W] = size(imgs{k});
        x1k = min(max(x1,1),W);  x2k = min(max(x2,1),W);
        y1k = min(max(y1,1),H);  y2k = min(max(y2,1),H);
        pk  = improfile(imgs{k}, [x1k x2k], [y1k y2k], N, 'bilinear');
        P(:,k) = pk(:);
    end

    % distance axis in pixels
    d = linspace(0, hypot(x2-x1, y2-y1), N);

    % plot
    figure('Color','w');
    plot(d, P, 'LineWidth', 1.8); grid on; box on;
    xlabel('Distance (pixels)'); ylabel('Intensity (a.u.)');
    legend({'Original','RL','Deconvolved'}, 'Location','best');
    title('Line intensity profiles');
end


ori = double(imread('ori.tif')); ori = ori / max(max(ori));
RL = double(imread('cell01_RL10.tif')); RL = RL /max(max(RL));
Sp = double(imresize(imread('cell01_Sparse.tif'),0.5));  Sp = Sp / max(max(Sp));
ZS = double(imread('Dec.tif'));  ZS = ZS / max(max(ZS));
ZS_R = double(imread('Hess_Dec.tif')); ZS_R = ZS_R / max(max(ZS_R));
plot_10(ori, RL, Sp, ZS, ZS_R, 270, 380, 320, 430, 100, 'Mode','manual','Range',[0 1]);
plot_profiles(ori, RL, Sp, ZS, ZS_R, 270, 380, 320, 430, 512);
vs = double(imread('Ves.tif')); vs = vs / max(max(vs));
dc = double(imread('HNZSDec.tif')); dc = dc / max(max(dc));
dn = double(imread('HNZSDen.tif')); dn = dn / max(max(dn));
vrl = double(imread('VRL.tif')); vrl = vrl / max(max(vrl));
plot_6(vs, vrl, dc, 165, 226, 190, 251, 50);
plot_profiles3(vs, vrl, dc, 165, 226, 190, 251, 1)
%imshow(ori);

