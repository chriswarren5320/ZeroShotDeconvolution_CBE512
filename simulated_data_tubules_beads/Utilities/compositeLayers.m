function I = compositeLayers(layers)
    I = max(layers,[],3);
    mx = max(I(:)); if mx>0, I = I/mx; end
end