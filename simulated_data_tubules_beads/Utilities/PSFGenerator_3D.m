function PSF = PSFGenerator_3D(dxy,dz,SizeXY,SizeZ,lamd,NA,RI)

dk = 2*pi/dxy/SizeXY;
kx = (-(SizeXY-1)/2:1:(SizeXY-1)/2)*dk;
[kx,ky] = meshgrid(kx,kx);
kr_sq = kx.^2+ky.^2;
z = (-(SizeZ-1)/2:1:(SizeZ-1)/2)*dz;

PupilMask = (kr_sq<=(2*pi/lamd*NA)^2);
kz = sqrt((2*pi/lamd*RI)^2-kr_sq).*PupilMask;
for ii = 1:SizeZ
    tmp = PupilMask.*exp(1i*kz*z(ii));
    tmp = fftshift(fft2(ifftshift(tmp)));
    PSF(:,:,ii) = (abs(tmp)).^2;
end

PSF = PSF/max(PSF(:));
% PSF = uint16(PSF);