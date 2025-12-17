function PSF = PSFGenerator_2D(dxy,SizeXY,lamd,NA)

dk = 2*pi/dxy/SizeXY;
kx = (-(SizeXY-1)/2:1:(SizeXY-1)/2)*dk;
[kx,ky] = meshgrid(kx,kx);
kr_sq = kx.^2+ky.^2;

PupilMask = (kr_sq<=(2*pi/lamd*NA)^2);
tmp = fftshift(fft2(ifftshift(PupilMask)));
PSF = (abs(tmp)).^2;
PSF = PSF/max(PSF(:));