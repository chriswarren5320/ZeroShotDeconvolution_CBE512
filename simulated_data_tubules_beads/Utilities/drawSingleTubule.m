function I = drawSingleTubule(N, p, w)
    I = zeros(N,N,'single');
    try
        rgb = repmat(I,[1 1 3]);
        rgb = insertShape(rgb,'Line',[p(1) p(2) p(3) p(4)], ...
                          'Color','white','LineWidth',w,'SmoothEdges',true);
        I = rgb2gray(rgb); I = single(I)/255;
    catch
        I = drawThickLine(I, p(1), p(2), p(3), p(4), w);
    end
    I = imfilter(I, fspecial('gaussian',3,0.7), 'replicate');
    I = min(max(I,0),1);
    I = I./max(I,[],'all');
end