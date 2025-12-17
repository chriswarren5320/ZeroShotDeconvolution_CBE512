function [psnr_val, a, b, I_trans] = compute_psnr_sr(sr_img, gt_img, psf)
% sr_img, gt_img, psf 都为 2D 矩阵

    %----------- 1) SR 与 PSF 卷积并下采样到 GT 尺寸 ------------
    psf = double(psf);
    psf = psf / (sum(psf(:)) + eps);      % 归一化 PSF

    I = conv2(double(sr_img), psf, 'same');  % 或者 imfilter
    figure; imagesc(I);axis square
    if ~isequal(size(I), size(gt_img))
        % 根据需要选择插值/步长
        I = imresize(I, size(gt_img), 'bicubic');
    end

    %----------- 2) 归一化 GT 到 [0,1]，线性拟合 I -> x -----------
    x = double(gt_img);
    x = (x - min(x(:))) / (max(x(:)) - min(x(:)) + eps);

    I = double(I);
    N = numel(I);

    % 构造最小二乘 [I(:)  ones] * [a; b] ≈ x(:)
    A = [I(:), ones(N,1)];
    theta = A \ x(:);   % 最小二乘解
    a = theta(1);
    b = theta(2);

    I_trans = a * I + b;

    %----------- 3) 计算 PSNR（MAX_I = 1） -----------------------
    mse = mean( (x(:) - I_trans(:)).^2 );
    if mse == 0
        psnr_val = Inf;
    else
        psnr_val = 10 * log10(1 / mse);
    end
end
