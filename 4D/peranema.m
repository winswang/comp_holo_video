% clear all; clear;
pera = im2double(rgb2gray(imread('peranema.png')));
pera = 1-pera;
pera(pera<0.5) = 0;
pera(pera>0.5) = 1;
%% parameters
lambda = 532e-3;        % in um
X = 912*10.8; Y = 912*10.8;
nx = 256; ny = 256;
z_offset = 70e3;
NA = (X/2)/sqrt((X/2)^2 + z_offset^2);
limit = lambda/(2*NA^2);
spacing = 50*limit;
ratio = 0.25;
nz = 2;

t_sq = 5;
t_num = length(t_sq);
d_num = 3;

for num = 1:t_num
    t = t_sq(num);
    for d = 1:d_num
%         if t < 3 || d < 3
%             break;
%         end
        %% construct object
        nt = t; spacing = 600*limit/d; dz(d) = spacing;nz = 2;
        obj = place_obj(pera,nx,ny,nz,nt,ratio);
        %% Forward model
        nz = size(obj,3);
        phase3D = getPhase3D(X, Y, nx, ny, nz, z_offset, z_offset + spacing, lambda);
        ref3D = getRef3D(X, Y, nx, ny, nz, z_offset, z_offset + spacing, lambda);
        [mask,percentage] = nonoverlap_mask(nx,ny,nt);
        out_forward = Forward4D(C2V4D(obj), ref3D, nx, ny, nz, nt, phase3D, mask, 0);
        backward = Backward4D(out_forward, ref3D, nx,ny,nz,nt,phase3D,mask);
        backward = reshape(V2C4D(backward,nx,ny,nz,nt), ny, nx, nz, nt);
        backward = backward/max(abs(backward(:)));
        out_reshape = reshape(MyV2C(out_forward),ny,nx);
%         out_noise = max(out_forward)*imnoise(out_reshape/max(out_forward),...
%             'gaussian',0,1e-20);
%         out_noise_res = reshape(out_noise, ny,nx);
%         imagesc(out_noise_res);
%         figure;
%         imagesc(out_noise_res - out_reshape);
%         out_snr(t,d) = snr(out_reshape,out_noise-out_reshape);
        
        %% TwIST
        g = out_forward;
        %% function handlers
        A = @(f_twist) Forward4D(f_twist,ref3D,nx,ny,nz,nt,phase3D,mask,0);  % forward propagation operator
        AT = @(g) Backward4D(g,ref3D,nx,ny,nz,nt,phase3D,mask);  % backward propagation operator
        
        N1 = ny; N2 = nx*2*nz*nt; N3 = 1;
        piter = 4;
        % tau = 0.001;
        % para = 0.1;
        %% try different parameters
        tau_string = 10.^(-linspace(1.3,3.3,6));
        para_string = 10.^(-linspace(0.4,1.4,6));
        best_tau = tau_string(1);best_para = para_string(1);best_mse = 1e10;
        for i = 1:length(tau_string)
            for j = 1:length(para_string)
                tau = tau_string(i);
                para = para_string(j);
                Psi = @(f,th) MyTVpsi(f,th,para,piter,N1,N2,N3);
                Phi = @(f) MyTVphi(f,N1,N2,N3);

                tolA = 1e-6;
                iterations = 20;

                [try_rec, debias, object, times, debiass, meanse]= ...
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
                    'Verbose', 1,...
                    'True_x', C2V4D(obj));
                try_rec = reshape(V2C4D(try_rec,nx,ny,nz,nt), ny, nx, nz, nt);
                try_rec = try_rec/max(try_rec(:));
%                 try_data(:,:,:,:,i,j) = try_rec;
                mserror(num,j) = meanse(end);
%                 try_psnr(i,j) = psnr(try_rec,obj);
                if mserror(num,j) <best_mse
                    best_tau = tau;
                    best_para = para;
                    best_mse = mserror(num,j);
                end
            end
        end
%         best_tau = 0.02;best_para = 0.1;
        tau = best_tau;
        para = best_para;
        
        Psi = @(f,th) MyTVpsi(f,th,para,piter,N1,N2,N3);
        Phi = @(f) MyTVphi(f,N1,N2,N3);

        tolA = 1e-6;
        iterations = 2000;

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
        obj_rec = reshape(V2C4D(obj_rec,nx,ny,nz,nt), ny, nx, nz, nt);
        obj_rec = obj_rec/max(obj_rec(:));
        psnr_bp(t,d) = psnr(backward,obj);
        psnr_cs(t,d) = psnr(obj_rec,obj);

        switch t
            case 1
                switch d
                    case 1
                        rec_t1_d1 = obj_rec;
                        obj_t1_d1 = obj;
                        img_t1_d1 = reshape(MyV2C(out_forward),ny,nx);
                    case 5
                        rec_t1_d5 = obj_rec;
                        obj_t1_d5 = obj;
                        img_t1_d5 = reshape(MyV2C(out_forward),ny,nx);
                end
            case 2
                switch d
                    case 1
                        rec_t2_d1 = obj_rec;
                        obj_t2_d1 = obj;
                        img_t2_d1 = reshape(MyV2C(out_forward),ny,nx);
                    case 3
                        rec_t2_d3 = obj_rec;
                        obj_t2_d3 = obj;
                        img_t2_d3 = reshape(MyV2C(out_forward),ny,nx);
                end
            case 5
                switch d
                    case 1
                        rec_t5_d1 = obj_rec;
                        obj_t5_d1 = obj;
                        img_t5_d1 = reshape(MyV2C(out_forward),ny,nx);
                    case 3
                        rec_t5_d3 = obj_rec;
                        obj_t5_d3 = obj;
                        img_t5_d3 = reshape(MyV2C(out_forward),ny,nx);
                end
            case 10
            switch d
                case 1
                    rec_t10_d1 = obj_rec;
                    obj_t10_d1 = obj;
                    img_t10_d1 = reshape(MyV2C(out_forward),ny,nx);
                case 3
                    rec_t10_d3 = obj_rec;
                    obj_t10_d3 = obj;
                    img_t10_d3 = reshape(MyV2C(out_forward),ny,nx);
            end
        end
    name = strcat('joint_reso_t_',num2str(t),'_d_',num2str(d),'.mat');
    save(name,'-v7.3');
    clear obj_rec try_rec
    end
end

% figure(1);
% plot(dz, psnr_cs(1,:),'o-');hold on
% plot(dz, psnr_cs(2,:),'o-');hold on
% plot(dz, psnr_cs(3,:),'o-');hold on
% plot(dz, psnr_cs(4,:),'o-');hold on
% plot(dz, psnr_cs(5,:),'o-');hold on
% plot(dz, psnr_bp(1,:),'o-');hold on
% 
% % xlabel('Object spacing (\mum)');ylabel('Temporal factor');zlabel('PSNR');legend('CS','BP');

