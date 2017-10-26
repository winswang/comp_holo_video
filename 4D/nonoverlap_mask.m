function [mask, perc] = nonoverlap_mask(nx,ny,nt)
% this function generates a temporal mask
% nx: horizontal pixel num
% ny: vertical pixel num
% nt: time step

a = rand(ny,nx);
bin = linspace(0,1,nt+1);
mask = zeros(ny,nx,nt);
perc = zeros(1,nt);
for i = 1:nt
    b = zeros(ny,nx);
    b(a>bin(i)&a<bin(i+1)) = 1;
    perc(i) = nnz(b)/nx/ny;
    mask(:,:,i) = b;
end
perc = mean(perc);