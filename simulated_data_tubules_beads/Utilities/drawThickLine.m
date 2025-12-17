function I = drawThickLine(I, x1,y1,x2,y2, w)
    n = max(abs([x2-x1, y2-y1])) + 1;
    xs = round(linspace(x1,x2,n));
    ys = round(linspace(y1,y2,n));
    [X,Y] = meshgrid(1:size(I,2), 1:size(I,1));
    rad = max(1, w/2);
    for i=1:n
        mask = ((X-xs(i)).^2 + (Y-ys(i)).^2) <= rad^2;
        I(mask) = 1;
    end
end