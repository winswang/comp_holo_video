function image = Forward4D(obj, ref, nx, ny, nz, nt, phase, mask, flag_intensity)
obj4D = V2C4D(obj,ny,nx,nz,nt);
Volume = reshape(obj4D, ny, nx, nz, nt);
ef = zeros(ny,nx,nz);
image = zeros(ny,nx);
for t = 1:nt
    obj = Volume(:,:,:,t);
    e = obj.*ref;
   
    for i=1:nz
        ef(:,:,i)=fftshift(fft2(ifftshift(e(:,:,i))));
    end

    cEsp=sum(ef.*phase,3);

    E=fftshift((ifft2(ifftshift(cEsp)))).*mask(:,:,t);
    image = E + image;
end




if nargin < 9
    flag_intensity = 0;
end
if flag_intensity == 0
    image = real(image);
elseif flag_intensity == 1
    image = 1 + 0.1*image;
    image = image.*conj(image);
end
image = MyC2V(image(:));
end