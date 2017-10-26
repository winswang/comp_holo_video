function E = Forward3D(obj, ref, nx, ny, nz, phase, mask)
obj = reshape(MyV2C(obj), ny, nx, nz);
e = obj.*ref;
ef = zeros(ny,nx,nz);

for i=1:nz
    ef(:,:,i)=fftshift(fft2(ifftshift(e(:,:,i))));
end

cEsp=sum(ef.*phase,3);

E=real(fftshift((ifft2(ifftshift(cEsp)))));

if nargin <7
    mask = ones(ny,nx);
end
E = E.*mask;

E = MyC2V(E(:));
end