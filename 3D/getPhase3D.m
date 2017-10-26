function phase = getPhase3D(X0, Y0, col, row, nz, zmin, zmax, lambda)
% 2016.5.9 this function gets a pre-computed phase factor

k = 2*pi/lambda;

dx = X0/col; dy = Y0/row; % spatial frequency
Kx = 2*pi/dx; Ky = 2*pi/dy;
kx = linspace((-Kx/2), (Kx/2), col);
ky = linspace((-Ky/2), (Ky/2), row);

[kxgrid, kygrid] = meshgrid(kx, ky);
z = linspace(zmin,zmax,nz);

% construct the circle function
circ = sqrt(kxgrid.^2 + kygrid.^2)/k;
circ(circ>1) = 0;
circ(circ<=1) = 1;

for i = 1:nz
    phase(:,:,i) = exp(1i*z(i)*sqrt(k^2 - kxgrid.^2 - kygrid.^2)).*circ;
end