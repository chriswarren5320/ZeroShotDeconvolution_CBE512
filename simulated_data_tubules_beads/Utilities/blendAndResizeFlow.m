function [Ut, Vt] = blendAndResizeFlow(U0,V0,U1,V1,tau,N)
    Uc = (1-tau)*U0 + tau*U1;
    Vc = (1-tau)*V0 + tau*V1;
    Ut = imresize(Uc, [N N], 'bicubic');
    Vt = imresize(Vc, [N N], 'bicubic');
    h = fspecial('gaussian', 11, 3);
    Ut = imfilter(Ut, h, 'replicate');
    Vt = imfilter(Vt, h, 'replicate');
end