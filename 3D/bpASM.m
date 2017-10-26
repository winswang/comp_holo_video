function E0 = bpASM(E, X, Y, z, lambda)

k = 2*pi/lambda;
[row, col] = size(E);
Er = ones(row, col);
dx = X/col; dy = Y/row; % spatial frequency
Kx = 2*pi/dx; Ky = 2*pi/dy;
kx = linspace((-Kx/2), (Kx/2), col);
ky = linspace((-Ky/2), (Ky/2), row);

[kxgrid, kygrid] = meshgrid(kx, ky);

% construct the circle function
circ = sqrt(kxgrid.^2 + kygrid.^2)/k;
circ(circ>1) = 0;
circ(circ<=1) = 1;

F = fftshift(fft2(ifftshift(E)));
factor = exp(1i*z*sqrt(k^2 - kxgrid.^2 - kygrid.^2));
E0 = fftshift(ifft2(ifftshift(F.*conj(factor.*circ))));

end