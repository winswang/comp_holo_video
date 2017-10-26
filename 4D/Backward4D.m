function e = Backward4D(sense,ref,nx,ny,nz,nt,phase,mask)

Sense = reshape(MyV2C(sense), ny,nx);
e=zeros(ny,nx,nz,nt);
for t = 1:nt
    image = Sense.*mask(:,:,t);
    f=fftshift(fft2(ifftshift(image)));

    cEs=conj(phase).*repmat(f,[1 1 nz]);

    for i=1:nz
        e(:,:,i,t)=fftshift(ifft2(ifftshift(cEs(:,:,i))));
    end

    e(:,:,:,t)=conj(ref).*e(:,:,:,t);

end
e = C2V4D(e);
end