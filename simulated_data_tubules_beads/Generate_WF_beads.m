close all
clear;
clc;
addpath('Utilities\')

%% Global parameters
imSize      = 512;         % image size (NxN)
xpixelsize   = 62 ;       % pixelsize of psf (lateral sampling)
zpixelsize  =  50;
R = 100;   % radius of the bead (unit: nm)
rng(1);
% Wide-field imaging
RI = 1.3;          % refractive index
NA       = 1.3;       % numerical aperture
lambda   = 525;       % emission wavelength, in (nm)
alphaPhot   = 800;         % photon scaling
bgCamera    = 700;         % camera background, Gaussian noise's mean value
sigmaRead   = 200;         % Gaussian standard deviation (std)

Rx_pixel = R/xpixelsize;
Rz_pixel = R/zpixelsize;
sz_x = 2*(floor(Rx_pixel)+1);  % size of the template single bead stack
sz_z = 2*(floor(Rz_pixel)+1);
if mod(sz_x,2)==0
    sz_x =sz_x + 1;
end
if mod(sz_z,2)==0
    sz_z =sz_z + 1;
end

Up_rate = 9;  % upsamping to get a more accurate bead template (should be an odd number)
% 1.1 generate a bead template
Upsz_x = sz_x * Up_rate;
Upsz_z = sz_z * Up_rate;
Up_midx = (Upsz_x-1)/2;
Up_midz = (Upsz_z-1)/2;
[xx, yy, zz] = meshgrid(-Up_midx:Up_midx, -Up_midx:Up_midx, -Up_midz:Up_midz);
Up_template = zeros(size(xx));
xx = xx*xpixelsize/Up_rate;
yy = yy*xpixelsize/Up_rate;
zz = zz*zpixelsize/Up_rate;
Up_template(xx.^2+yy.^2+zz.^2<R^2)=1;
% isolim = 0.9;
% figure;
% isosurface(xx, yy, zz, Up_template, isolim);
% axis equal
% title('single bead')
template = zeros(sz_x, sz_x, sz_z);
for i = 1:sz_x
    for j = 1:sz_x
        for k = 1:sz_z
            tmp = Up_template( (i-1)*Up_rate+1: i*Up_rate, (j-1)*Up_rate+1: j*Up_rate, (k-1)*Up_rate+1: k*Up_rate);
            template(i,j,k) = sum(tmp(:));
        end
    end
end
template = template./max(template(:));

%% generate 2D beads
num_bead = 150;
dis = 3;
range = 1+floor(dis*(Rx_pixel)):imSize-floor(dis*(Rx_pixel));
xlist = randi([range(1) range(end)], [num_bead,1]);
ylist = randi([range(1) range(end)], [num_bead,1]);
R1 = (sz_x -1)/2;
R2 = (sz_z -1)/2;

template_z0 = template(:,:,(size(template,3)-1)/2);
data_xy = zeros(imSize,imSize);
for i = 1:num_bead
    xloc = xlist(i);
    yloc = ylist(i);
    data_xy(xloc-R1:xloc+R1, yloc-R1:yloc+R1) = data_xy(xloc-R1:xloc+R1, yloc-R1:yloc+R1) + template_z0;
end
figure; imagesc(data_xy); axis square
GT = data_xy;

%% --------- PSF for WF simulation ----------
% psf = PSFGenerator_2D(xpixelsize,imSize,lambda,NA);
if xpixelsize==31
    PSF_filename = 'PSF_imageJ\PSF-dxy_31.tif';
elseif xpixelsize==62
    PSF_filename = 'PSF_imageJ\PSF-dxy_62.tif';
end
    
info = imfinfo(PSF_filename);
nSlices = numel(info);
psf_stack = zeros(info(1).Height, info(1).Width, nSlices, 'single');
for k = 1:nSlices
    psf_stack(:,:,k) = im2single(imread(PSF_filename, k));
end
psf = psf_stack(:,:,ceil(nSlices/2));
figure; imagesc(psf)

%% --- WF simulation: PSF convolution + noise ---
cleanWF = real(ifftshift(ifft2(fft2(fftshift(data_xy)).*fft2(fftshift(psf)))));
figure; imagesc(cleanWF); axis square
% 2× down-sample via bicubic (match experimental pixel size)
% cleanWF = imresize(cleanWF, 0.5, 'bicubic'); % pixelsize_image = 31.3*2=62.6 nm

p = prctile(cleanWF(:), 99.9);
I_WF_clear = cleanWF / max(p, eps);
I_WF_clear = I_WF_clear * alphaPhot;
WF_poisson =  poissrnd(I_WF_clear);
gaussNoise = bgCamera + sigmaRead* randn(imSize, 'single');
WF = WF_poisson+gaussNoise;
figure; imagesc(WF); axis square

WF = uint16(WF);
GT = uint16(GT* 65535);
PSF = uint16(psf* 65535);
cleanWF = uint16(cleanWF);

folderName = sprintf('Simulated_beads_62_nosample/bg_%d_std_%d', bgCamera,sigmaRead);  
GTfolderName = fullfile(folderName,'GT');
WFfolderName = fullfile(folderName,'WF');
cleanWFfolderName = fullfile(folderName,'cleanWF');
if ~exist(folderName, 'dir')
    mkdir(folderName);
    mkdir(fullfile(folderName,'GT'));
    mkdir(fullfile(folderName,'WF'));
    mkdir(fullfile(folderName,'cleanWF'));
    fprintf('Created folder: %s\n', folderName);
else
    fprintf('Folder already exists: %s\n', folderName);
end
imwrite(PSF, fullfile(folderName, 'PSF.tif'), 'tif', 'Compression', 'none');

WFfilename = [WFfolderName,'\',sprintf('WF.tif')];
GTfilename = [GTfolderName,'\',sprintf('GT.tif')];
cleanWFfilename = [cleanWFfolderName,'\',sprintf('cleanWF.tif')];
% Write first frame with 'overwrite' mode
imwrite(WF, WFfilename, 'tif', 'Compression', 'none');
imwrite(GT, GTfilename, 'tif', 'Compression', 'none');
imwrite(cleanWF, cleanWFfilename, 'tif', 'Compression', 'none');


