function out = addNotches(in, sigma, prob)
    N = size(in,1);
    noise = randn(N,N,'single');
    g = fspecial('gaussian', 9, sigma);
    blobs = imfilter(noise, g, 'replicate');
    blobs = (blobs - min(blobs(:))) / max(1e-6, (max(blobs(:))-min(blobs(:))));
    notch = blobs > (1-prob);
    out = in; out(notch & in>0.2) = 0;
end