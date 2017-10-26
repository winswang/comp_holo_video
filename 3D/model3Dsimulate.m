clear all; clear;

%% Construct the object
nx = 128; ny = 128; nz = 5;
obj = zeros(ny,nx,nz);
obj(21:40,21:40,1) = 1;
obj(41:50,41:50,2) = 0.8;
obj(:,100:105,3) = 0.6;
obj(60:65,120:125,4) = 1;
obj(80:82,20:60,5) = 0.3;

complex = 0;    % 1: add random phase

if complex == 1
    rand_phase = rand(size(obj));
    obj = obj.*exp(1i*rand_phase*2);
end

%% display the object (1)
figure(1);
imagesc(plotdatacube(abs(obj),nz)); axis equal; axis off; colormap hot;
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
mask = ones(ny,nx);intensity = 0; % intensity formation: 0--real; 1--square term
out_forward = Forward3D(MyC2V(obj(:)), ref3D, nx, ny, nz, phase3D, mask, intensity);
%% display the sensed image
figure(2);
% reshape the output
sensor_image = reshape(MyV2C(out_forward), nx, ny);
imagesc(abs(sensor_image));colormap hot; axis off; axis equal;
title('Sensed image');

%% Back propagation (3)
out_backward = Backward3D(out_forward, ref3D, nx, ny, nz, phase3D);
obj_backward = reshape(MyV2C(out_backward), nx, ny, nz);
figure(3);
imagesc(plotdatacube(abs(obj_backward),nz)); axis equal; axis off; colormap hot;
title('Back propagation');

%% TwIST Reconstruction
% function handlers
A = @(f_twist) Forward3D(f_twist,ref3D,nx,ny,nz,phase3D,mask,intensity);  % forward propagation operator
AT = @(g) Backward3D(g,ref3D,nx,ny,nz,phase3D);  % backward propagation operator
% twist parameters
tau_string = 10.^(-linspace(0,4,20)); 
para_string = 10.^(-linspace(0,3,10));
para_space = zeros(length(tau_string),length(para_string));
best_tau = tau_string(1);
best_para = para_string(1);
best_psnr = 0;

piter = 4;
tolA = 1e-6;
iterations = 300;
N1 = nx; N2 = ny*nz*2; N3 = 1;

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

        obj_rec = reshape(MyV2C(obj_rec), nx, ny, nz);
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
%% display reconstructed objects (4)
figure(4);
% subplot(3,1,1);title('Reconstructed objects');colormap hot;
imagesc(plotdatacube(abs(obj_rec_best),nz));axis equal; axis off;colormap hot;
title('Reconstructed objects');
% subplot(3,1,2);title('Reconstructed (real part)');colormap hot;
% imagesc(imagesc(plotdatacube(abs(obj_rec_real),nz)));axis equal; axis off;
% subplot(3,1,3);title('Reconstructed (imaginary part)');colormap hot;
% imagesc(imagesc(plotdatacube(abs(obj_rec_imag),nz)));axis equal; axis off;



