clear;clc;
close all
addpath('Utilities\')
base = 'Experiments_beads';
base_folder = fullfile(base,'SR_1024');% folder that contains 7 subfolders
% reg_folders = { ...
%     'Adapt', ...
%     'Energy', ...
%     'Hess', ...
%     'l1', ...
%     'sparse', ...
%     'TV', ...
%     'Wavelet'};     % <-- change to your real folder names

d = dir(base_folder);
subfolders = d([d.isdir] & ~ismember({d.name}, {'.','..'}));
reg_folders = {subfolders.name};

gt_path  = fullfile(base,'CleanWF.tif');   % ground truth image
psf_path = fullfile(base,'PSF-dxy_62.tif');      % file that contains variable "psf"

% Fluorescence signal labels (for legend)
signal_labels = {'50','100','150','200'};
gt_img = im2double(imread(gt_path));
psf = im2double(imread(psf_path ));
figure;imagesc(gt_img);axis square

nReg = numel(reg_folders);
nSig = numel(signal_labels);  % 4
psnr_vals = nan(nReg, nSig); % 7 x 4 matrix

%% Loop over regularization terms and signal levels
for i = 1:nReg
    folder = fullfile(base_folder, reg_folders{i});
    files = dir(fullfile(folder, '*.tif'));
    % Sort files by name so that signal levels have a consistent order
    nFiles = numel(files);
    signal_values = zeros(nFiles, 1);
    for k = 1:nFiles
        fname = files(k).name;
        
        % Extract all numbers from filename (as strings)
        num_str = regexp(fname, '\d+', 'match');
        
        % Convert to number (take the first number found)
        signal_values(k) = str2double(num_str{1});
    end
    [~, idx] = sort(signal_values, 'ascend');
    files = files(idx);

    if numel(files) < nSig
        warning('Folder %s has only %d tif files (expected %d).', ...
                folder, numel(files), nSig);
    end
    
    for s = 1:min(nSig, numel(files))
        fname = fullfile(folder, files(s).name);
        sr_img = im2double(imread(fname));
        
        % Compute PSNR for this reg term & signal level
        [psnr_vals(i, s), ~, ~, ~] = compute_psnr_sr(sr_img, gt_img, psf);
    end
end
%% Optional: create a table with Reg, Signal, PSNR
Reg  = strings(nReg*nSig,1);
Sig  = strings(nReg*nSig,1);
PSNR = nan(nReg*nSig,1);
idx = 1;

for i = 1:nReg
    for s = 1:nSig
        Reg(idx)  = reg_folders{i};
        Sig(idx)  = signal_labels{s};
        PSNR(idx) = psnr_vals(i,s);
        idx = idx + 1;
    end
end

T = table(Reg, Sig, PSNR);
disp(T);


figure; 
b = bar(psnr_vals);   % 7x4 matrix -> grouped bar

% Set colors for the 4 signal levels (each group has 4 bars)
colors = lines(nSig);  % or turbo(nSig), parula(nSig), etc.
for s = 1:nSig
    b(s).FaceColor = colors(s,:);
end

% Axis labels and legend
set(gca, 'XTick', 1:nReg, 'XTickLabel', reg_folders, 'FontSize', 11);
xtickangle(30);  % rotate x labels a bit if they are long

ylabel('PSNR (dB)');
xlabel('Regularization term');
title('PSNR for 7 regularization terms and 4 fluorescence signal levels');

% lgd=legend(signal_labels, 'Location', 'northwest');
% title(lgd, 'Fluorescence signal');
grid on; box on;