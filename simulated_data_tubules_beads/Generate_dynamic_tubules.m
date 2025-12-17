clc; clear;close all
filePath = mfilename('fullpath');
[currentDir,~,~] = fileparts(filePath);
addpath(genpath(currentDir));
%% Global parameters
imSize      = 512;         % image size (NxN)
pixelsize   = 62 ;       % pixelsize of PSF (lateral sampling)
nFrames     = 5;
rng(1);                    % reproducibility

% Tubule geometry
nTubules    = 40;          % number of tubules
tubuleWidth = 3;           % line thickness (pixels)
addNotch    = true;        % whether to simulate incomplete labeling
notchSigma  = 2;           
notchProb   = 0.12;

% Translation velocity (pixels/frame)
% If you know physical speed v_um_per_s and frame rate fps,
% use: v_pix = v_um_per_s / (px_um * fps)
v_pix       = 1;          % ~1 µm/s @ 1 Hz and 31.3 nm/px

% Bending (local elastic deformation)
gridN_bend  = 4;           % control-grid size for random flow
bendAmp     = 4;           % amplitude of bending (pixels)
temporalScl = 0.15;        % temporal smoothness of the flow
maskDilateR = 4;           % dilation radius for mask expansion
maskSigma   = 2.0;         % Gaussian blur sigma for smooth mask

% Wide-field imaging
RI = 1.3;          % refractive index
NA          = 1.3;       % numerical aperture
lambda   = 525;       % emission wavelength, in (nm)
alphaPhot   = 800;         % photon scaling
bgCamera    = 700;         % camera background, Gaussian noise's mean value
sigmaRead   = 500;         % Gaussian standard deviation (std)

%% --------- Static tubule layers ----------
layers = zeros(imSize, imSize, nTubules, 'single');
pts    = randi([1 imSize], nTubules, 4);  % [x1 y1 x2 y2] endpoints

for k = 1:nTubules
    % layers(:,:,k) = drawSingleTubule(imSize, pts(k,:), tubuleWidth);
    layers(:,:,k) = drawCurvedTubuleControlled(imSize,tubuleWidth);
    if addNotch
        layers(:,:,k) = addNotches(layers(:,:,k), notchSigma, notchProb);
    end
end
base0 = compositeLayers(layers);  % first static composite

% figure; imagesc(base0);axis square

%% --------- PSF for WF simulation ----------
% info = imfinfo('PSF_simulated.tif');
% nSlices = numel(info);
% psf_stack = zeros(info(1).Height, info(1).Width, nSlices, 'single');
% for k = 1:nSlices
%     psf_stack(:,:,k) = im2single(imread('PSF_simulated.tif', k));
% end
% psf1 = psf_stack(:,:,ceil(nSlices/2));

% psf = PSFGenerator_2D(pixelsize,imSize,lambda,NA);
% figure; imagesc(psf)

if pixelsize==31
    PSF_filename = 'PSF_imageJ\PSF-dxy_31.tif';
elseif pixelsize==62
    PSF_filename = 'PSF_imageJ\PSF-dxy_62.tif';
end
    
info = imfinfo(PSF_filename);
nSlices = numel(info);
psf_stack = zeros(info(1).Height, info(1).Width, nSlices, 'single');
for k = 1:nSlices
    psf_stack(:,:,k) = im2single(imread(PSF_filename, k));
end
psf = psf_stack(:,:,ceil(nSlices/2));
figure; imagesc(psf);axis square

%% --------- Initial bending fields for each tubule ----------
U0 = cell(nTubules,1); V0 = cell(nTubules,1);
U1 = cell(nTubules,1); V1 = cell(nTubules,1);
for i = 1:nTubules
    [U0{i}, V0{i}] = randomGridFlow(gridN_bend, bendAmp);
    [U1{i}, V1{i}] = randomGridFlow(gridN_bend, bendAmp);
end

%% --------- Output containers ----------
GT = zeros(imSize, imSize, nFrames, 'single');
WF = zeros(imSize, imSize, nFrames, 'single');
cleanWF_stack = zeros(imSize, imSize, nFrames, 'single');

%% --------- Coordinate grid ----------
[X, Y] = meshgrid(1:imSize, 1:imSize);
se = strel('disk', maskDilateR);  % for dilating masks

%% --------- Time-lapse generation ----------
for t = 1:nFrames
    frameGT = zeros(imSize, imSize, 'single');
    tau = 0.5 * (1 + sin(2*pi*temporalScl*(t-1)));  % smooth interpolation factor

    for i = 1:nTubules
        % --- (1) Translation per paper Eq. 45 ---
        ui_t  = rand();                       % random in [0,1]
        Di_t  = ((1 + 2*ui_t)/3) * v_pix;     % displacement magnitude
        theta = 2*pi*rand();                  % random direction
        dx_t  = Di_t * cos(theta);
        dy_t  = Di_t * sin(theta);

        % --- (2) Bending deformation field ---
        [Ub, Vb] = blendAndResizeFlow(U0{i}, V0{i}, U1{i}, V1{i}, tau, imSize);
        mask = layers(:,:,i) > 0.2;           % tubule mask
        mask = imdilate(mask, se);            % expand influence region
        mask = single(imgaussfilt(single(mask), maskSigma)); % smooth mask
        Ub = Ub .* mask;  Vb = Vb .* mask;    % apply locally

        % --- (3) Combine translation + bending ---
        xw = X + dx_t + Ub; 
        yw = Y + dy_t + Vb;

        % --- (4) Non-rigid warp using interp2 ---
        moved = interp2(X, Y, layers(:,:,i), xw, yw, 'linear', 0);

        % --- (5) Composite current tubule layer ---
        frameGT = max(frameGT, moved);
    end

    % Normalize GT frame
    mx = max(frameGT(:)); if mx>0, frameGT = frameGT/mx; end
    GT(:,:,t) = frameGT;

    % --- WF simulation: PSF convolution + noise ---
    cleanWF = real(ifftshift(ifft2(fft2(fftshift(frameGT)).*fft2(fftshift(psf)))));
    cleanWF_stack(:,:,t) = cleanWF * alphaPhot;
    % % 2× down-sample via bicubic (match experimental pixel size)
    % cleanWF = imresize(cleanWF, 0.5, 'bicubic'); % pixelsize_image = 31.3*2=62.6 nm
    p = prctile(cleanWF(:), 99.9);
    I_WF_clear = cleanWF / max(p, eps);
    I_WF_clear = I_WF_clear * alphaPhot;
    WF_poisson =  poissrnd(I_WF_clear);
    gaussNoise = bgCamera + sigmaRead* randn(imSize, 'single');
    WF(:,:,t) = WF_poisson+gaussNoise;
    % Periodically refresh bending fields for smooth temporal changes
    if mod(t, round(1/temporalScl))==0
        U0 = U1; V0 = V1;
        for i = 1:nTubules
            [U1{i}, V1{i}] = randomGridFlow(gridN_bend, bendAmp);
        end
    end
end

if nFrames==1
    WF = squeeze(WF);
    GT = squeeze(GT);
    cleanWF_stack = squeeze(cleanWF_stack);
    figure;imagesc(cleanWF_stack);title('cleanWF');axis square
    figure;imagesc(WF);title('WF');axis square
    figure;imagesc(GT);title('GT');axis square
else
    imageslicer(WF);title('WF');axis square
    imageslicer(GT);title('GT');axis square
    imageslicer(cleanWF_stack);title('cleanWF_stack');axis square
end

WF = uint16(WF);
GT = uint16(GT* 65535);
cleanWF_stack = uint16(cleanWF_stack);


%% Save data
% folderName = sprintf('Simulated_tubules_62_nosample/bg_%d_std_%d', bgCamera,sigmaRead);  
% GTfolderName = fullfile(folderName,'GT');
% WFfolderName = fullfile(folderName,'WF');
% CleanWFfolderName = fullfile(folderName,'CleanWF');
% 
% if ~exist(folderName, 'dir')
%     mkdir(folderName);
%     mkdir(fullfile(folderName,'GT'));
%     mkdir(fullfile(folderName,'WF'));
%     mkdir(fullfile(folderName,'CleanWF'));
%     fprintf('Created folder: %s\n', folderName);
% else
%     fprintf('Folder already exists: %s\n', folderName);
% end
% 
% if size(WF,3)==1
%     WFfilename = [WFfolderName,'\',sprintf('WF.tif')];
%     GTfilename = [GTfolderName,'\',sprintf('GT.tif')];
%     CleanWFfilename = [CleanWFfolderName,'\',sprintf('CleanWF.tif')];
%     imwrite(WF, WFfilename, 'tif', 'Compression', 'none');
%     imwrite(GT, GTfilename, 'tif', 'Compression', 'none');
%     imwrite(cleanWF_stack, CleanWFfilename, 'tif', 'Compression', 'none');
% else
%     for k = 1:size(WF,3)
%         % Convert to uint16 or uint8 if needed for display
%         WFfilename = [WFfolderName,'\',sprintf('WF_%d.tif',k)];
%         GTfilename = [GTfolderName,'\',sprintf('GT_%d.tif',k)];
%         CleanWFfilename = [CleanWFfolderName,'\',sprintf('CleanWF_%d.tif',k)];
%         WFslice = WF(:,:,k);
%         GTslice = GT(:,:,k);
%         CleanWFslice = cleanWF_stack(:,:,k);
%         % Write first frame with 'overwrite' mode
%         imwrite(WFslice, WFfilename, 'tif', 'Compression', 'none');
%         imwrite(GTslice, GTfilename, 'tif', 'Compression', 'none');
%         imwrite(CleanWFslice, CleanWFfilename, 'tif', 'Compression', 'none');
%     end
% end