function Ref = getRef3D(X0, Y0, col, row, nz, zmin, zmax, lambda)
E = ones(row, col);
Ref = zeros(row, col, nz);
phase3D = getPhase3D(X0, Y0, col, row, nz, zmin, zmax, lambda);
for i=1:nz
    cE0=fftshift(fft2(E));

    cE=cE0.*conj(phase3D(:,:,i));

    Ref(:,:,i)=ifft2(ifftshift(cE));
end
end