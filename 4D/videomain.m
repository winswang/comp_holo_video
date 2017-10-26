% this script simulates the process of a 4-D video
%% draw the figure
close all;clc;
clear all;
% addpath('./Functions');
%% generate a 4D scene
nx = 128;
ny = 128;
nz = 5;
nt = 5;
f = zeros(ny, nx, nz, nt);
obj1 = zeros(ny,nx);
obj1(64-3:64+3,:) = 1;
obj2 = zeros(ny,nx);
obj2(:,64-3:64+3) = 1;
obj3 = zeros(ny, nx);
obj3(20:23,20:23) = 1;
obj4 = zeros(ny, nx);
obj4(90:93,90:93) = 1;
f(:,:,1,1) = circshift(obj1, [3, 0]);
f(:,:,1,2) = circshift(obj1, [6, 0]);
f(:,:,1,3) = circshift(obj1, [9, 0]);
f(:,:,1,4) = circshift(obj1, [12, 0]);
f(:,:,1,5) = circshift(obj1, [15, 0]);
f(:,:,4,1) = circshift(obj2, [0, 1]);
f(:,:,4,2) = circshift(obj2, [0, 2]);
f(:,:,4,3) = circshift(obj2, [0, 3]);
f(:,:,4,4) = circshift(obj2, [0, 4]);
f(:,:,4,5) = circshift(obj2, [0, 5]);
f(:,:,2,1) = circshift(obj3, [0, 1]);
f(:,:,2,2) = circshift(obj3, [1, 1]);
f(:,:,2,3) = circshift(obj3, [1, 2]);
f(:,:,2,4) = circshift(obj3, [2, 2]);
f(:,:,2,5) = circshift(obj3, [2, 3]);
f(:,:,3,1) = circshift(obj4, [0, 3]);
f(:,:,3,2) = circshift(obj4, [3, 3]);
f(:,:,3,3) = circshift(obj4, [3, 6]);
f(:,:,3,4) = circshift(obj4, [6, 6]);
f(:,:,3,5) = circshift(obj4, [6, 9]);

%% generate mask
percentage = 0.5;
mask = genmask(nx, ny, nt, percentage);

%% Paprameters
dx=5.86e1; dy = 5.86e1; % pixel patch (um)
X = dx*nx; Y = dy*ny;  % detector size (um)

lambda=0.532;  % wavelength (um)
deltaZ=20e3;  % axial spacing (um) 20k
offsetZ=300e3;  % distance from detector to first reconstructed plane (um)
zmin = offsetZ;
zmax = offsetZ + (nz-1)*deltaZ;

N1=ny;
N2=nx*nz*2*nt;
N3=1;

%% Kernel propagation
phase3D = getPhase3D(X, Y, nx, ny, nz, zmin, zmax, lambda);
ref = getRef3D(X, Y, nx, ny, nz, zmin, zmax, lambda);

%% Field propagation
image = Forward4D(C2V4D(f), ref, nx, ny, nz, nt, phase3D, mask);
back = Backward4D(image,ref,nx,ny,nz,nt,phase3D,mask);

%% Add noise
sigma = 1e-3;
image_noise = image + sigma*randn(size(image));

%% TwIST
%% Propagation operator (5)
A = @(f_twist) Forward4D(f_twist,ref,nx,ny,nz,nt,phase3D,mask);  % forward propagation operator
AT = @(g) Backward4D(g,ref,nx,ny,nz,nt,phase3D,mask);  % backward propagation operator


%% TwIST algorithm (6)
% twist parameters
tau = 0.01; 
piter = 4;
tolA = 1e-6;
iterations = 3000;

% g = MyC2V(g(:));

Psi = @(f,th) MyTVpsi(f,th,0.05,piter,N1,N2,N3);
Phi = @(f) MyTVphi(f,N1,N2,N3);

[f_reconstruct,dummy,obj_twist,...
    times_twist,dummy,mse_twist]= ...
    TwIST(image_noise,A,tau,...
    'AT', AT, ...
    'Psi', Psi, ...
    'Phi',Phi, ...
    'Initialization',0,...
    'Monotone',1,...
    'StopCriterion',1,...
    'MaxIterA',iterations,...
    'MinIterA',iterations,...
    'ToleranceA',tolA,...
    'Verbose', 1);


% evaluation

f_reconstruct=reshape(V2C4D(f_reconstruct,ny,nx,nz,nt),ny,nx,nz,nt);
frec = abs(f_reconstruct);
frec = frec/max(frec(:));
for i = 1:5
    fig(:,:,i) = f(:,:,i,1);
    fig(:,:,i+5) = frec(:,:,i,1);
    fig(:,:,i+10) = f(:,:,i,5);
    fig(:,:,i+15) = frec(:,:,i,5);
end
figure(1);imagesc(plotdatacube(fig),5);axis off;colormap gray;
implay(frec(:,:,1,:));
v = VideoWriter('6.1');
open(v);
rec(:,:,1,:) = frec(:,:,1,:);
writeVideo(v,rec);
close(v);
