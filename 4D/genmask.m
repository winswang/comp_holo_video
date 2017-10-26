function mask = genmask(nx, ny, nt, percentage)
% this function generates a temporal mask
% nx: horizontal pixel num
% ny: vertical pixel num
% nt: time step
% percentage: percentage of 1 exposure unit
randd = rand(ny,nx,nt);
mask = zeros(ny,nx,nt);
idx = randd <= percentage;
mask(idx) = 1;
% realperc = sum(idx(:))/(nx*ny*nt);
end