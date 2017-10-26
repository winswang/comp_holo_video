clear all; clear;

%% Construct the object
nx = 128; ny = 128; nz = 5; nt = 5;
spd = 6;
obj1 = zeros(ny,nx);
obj1(21:40,21:40) = 0.5;
obj2 = zeros(ny,nx);
obj2(21:30,46:55) = 0.7;
obj3 = zeros(ny,nx);
obj3(21:30,46:55) = 0.7;
obj4 = zeros(ny,nx);
obj4(21:25,61:65) = 0.9;
obj5 = zeros(ny,nx);
obj5(21:23,81:83) = 1;
obj = zeros(ny,nx,nz,nt);
obj(:,:,1,1) = obj1;
obj(:,:,1,2) = circshift(obj1,[1*spd,0]);
obj(:,:,1,3) = circshift(obj1,[2*spd,0]);
obj(:,:,1,4) = circshift(obj1,[3*spd,0]);
obj(:,:,1,5) = circshift(obj1,[4*spd,0]);
obj(:,:,2,1) = obj2;
obj(:,:,2,2) = circshift(obj2,[1*spd,0]);
obj(:,:,2,3) = circshift(obj2,[2*spd,0]);
obj(:,:,2,4) = circshift(obj2,[3*spd,0]);
obj(:,:,2,5) = circshift(obj2,[4*spd,0]);
obj(:,:,3,1) = obj3;
obj(:,:,3,2) = circshift(obj3,[2*spd,0]);
obj(:,:,3,3) = circshift(obj3,[4*spd,0]);
obj(:,:,3,4) = circshift(obj3,[6*spd,0]);
obj(:,:,3,5) = circshift(obj3,[8*spd,0]);
obj(:,:,4,1) = obj4;
obj(:,:,4,2) = circshift(obj4,[2*spd,0]);
obj(:,:,4,3) = circshift(obj4,[4*spd,0]);
obj(:,:,4,4) = circshift(obj4,[6*spd,0]);
obj(:,:,4,5) = circshift(obj4,[8*spd,0]);
obj(:,:,5,1) = obj5;
obj(:,:,5,2) = circshift(obj5,[3*spd,0]);
obj(:,:,5,3) = circshift(obj5,[6*spd,0]);
obj(:,:,5,4) = circshift(obj5,[9*spd,0]);
obj(:,:,5,5) = circshift(obj5,[12*spd,0]);


complex = 0;    % 1: add random phase

if complex == 1
    rand_phase = rand(size(obj));
    obj = obj.*exp(1i*rand_phase*2);
end

%% display the object (1)
figure(1);
imagesc(plotdatacube(abs(reshape(obj,nx,ny,nz*nt)),nt)); axis equal; axis off; colormap hot;
title('Constructed object field');
%% parameters
lambda = 532e-3;        % in um
dx = 1; dy = 1;
X = nx*dx; Y = ny*dy;
z_offset = 500;
z_interval = 5;

%% Forward model
phase3D = getPhase3D(X, Y, nx, ny, nz, z_offset, z_offset + (nz-1)*z_interval, lambda);
ref3D = getRef3D(X, Y, nx, ny, nz, z_offset, z_offset + (nz-1)*z_interval, lambda);
flag_mask = 1;  % 0--uniform random masks; 1--non-overlap random masks
percentage = 0.5;
if flag_mask == 0
    mask = genmask(nx,ny,nt,0.5);
elseif flag_mask == 1
    [mask,percentage] = nonoverlap_mask(nx,ny,nt);
end
intensity = 0; % intensity formation: 0--real; 1--square term
out_forward = Forward4D(MyC2V(obj(:)), ref3D, nx, ny, nz, nt, phase3D, mask, intensity);
%% display the sensed image
figure(2);
% reshape the output
sensor_image = reshape(MyV2C(out_forward), nx, ny);
imagesc(abs(sensor_image));colormap hot; axis off; axis equal;
title('Sensed image');

%% Back propagation (3)
out_backward = Backward4D(out_forward, ref3D, nx, ny, nz, nt, phase3D, mask);
obj_backward = reshape(MyV2C(out_backward), nx, ny, nz, nt);
figure(3);
imagesc(plotdatacube(abs(reshape(obj_backward,nx,ny,nz*nt)),nz)); axis equal; axis off; colormap hot;
title('Back propagation');

%% TwIST Reconstruction
% function handlers
A = @(f_twist) Forward4D(f_twist,ref3D,nx,ny,nz,nt,phase3D,mask,intensity);  % forward propagation operator
AT = @(g) Backward4D(g,ref3D,nx,ny,nz,nt,phase3D,mask);  % backward propagation operator
% twist parameters
tau_string = 10.^(-linspace(2,4,8)); 
para_string = 10.^(-linspace(1,3,8));
para_space = zeros(length(tau_string),length(para_string));
best_tau = tau_string(1);
best_para = para_string(1);
best_psnr = 0;

piter = 4;
tolA = 1e-6;
iterations = 100;
N1 = nx; N2 = ny*nz*2*nt; N3 = 1;

% init: 0--all zeros; 1--random; 2--A'y
g = out_forward;

for i = 1:length(tau_string)
    tau = tau_string(i);
    for j = 1:length(para_string)
        para = para_string(j);
        Psi = @(f,th) MyTVpsi(f,th,para,piter,N1,N2,N3);
        Phi = @(f) MyTVphi(f,N1,N2,N3);

        obj_rec= ...
            TwIST(g,A,tau,...
            'AT', AT, ...
            'Psi', Psi, ...
            'Phi',Phi, ...
            'Initialization',2,...
            'Monotone',1,...
            'StopCriterion',1,...
            'MaxIterA',iterations,...
            'MinIterA',iterations,...
            'ToleranceA',tolA,...
            'Verbose', 1);

        obj_rec = reshape(V2C4D(obj_rec,nx,ny,nz,nt), nx, ny, nz, nt);
        obj_rec_real = real(obj_rec);
        obj_rec_imag = imag(obj_rec);
        psnr_obj = psnr(abs(obj_rec),abs(obj));
        para_space(i,j) = psnr_obj;
        if psnr_obj > best_psnr
            best_psnr = psnr_obj;
            best_tau = tau;
            best_para = para;
            obj_rec_best = obj_rec;
            psnr_obj_real = psnr(obj_rec_real,real(obj));
        psnr_obj_imag = psnr(obj_rec_imag,imag(obj));
        end
    end
end