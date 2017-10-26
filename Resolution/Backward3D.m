function e = Backward3D(sense,ref,nx,ny,nz,phase)

Sense = reshape(MyV2C(sense), ny,nx);

f=fftshift(fft2(ifftshift(Sense)));

cEs=conj(phase).*repmat(f,[1 1 nz]);

e=zeros(ny,nx,nz);
for i=1:nz
    e(:,:,i)=fftshift(ifft2(ifftshift(cEs(:,:,i))));
end

e=conj(ref).*e;

e = MyC2V(e(:));
end