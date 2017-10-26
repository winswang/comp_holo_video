% clear all; clear;
pera = im2double(rgb2gray(imread('peranema.png')));
pera = 1-pera;
%% parameters
lambda = 532e-3;        % in um
X = 912*10.8; Y = 912*10.8;
nx = 256; ny = 256;
z_offset = 70e3;
NA = (X/2)/sqrt((X/2)^2 + z_offset^2);
limit = lambda/(2*NA^2);
spacing = 5*limit;
ratio = 0.25;
nz = 51;
t_num = 1;
d_max = 1;

for t = 1:t_num
    for d = 1:d_max
%         if t < 3 || d < 3
%             break;
%         end
        %% construct object
        nt = t;nz = 2;
        spacing = 50*limit/d;
        obj = place_obj(pera,nx,ny,nz,nt,ratio);
        %% Forward model
        nz = size(obj,3);
        phase3D = getPhase3D(X, Y, nx, ny, nz, z_offset, z_offset + spacing, lambda);
        ref3D = getRef3D(X, Y, nx, ny, nz, z_offset, z_offset + spacing, lambda);
        [mask,percentage] = nonoverlap_mask(nx,ny,nt);
        for tt = 1:t
            maskt = mask(:,:,tt);
            out_forward = Forward3D(MyC2V(obj), ref3D, nx, ny, nz, nt, phase3D, maskt);
            backward = Backward4D(out_forward, ref3D, nx,ny,nz,nt,phase3D);
            backward = reshape(V2C4D(backward,nx,ny,nz,nt), ny, nx, nz, nt);
    %         out_reshape = reshape(MyV2C(out_forward),ny,nx);
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
            A = @(f_twist) Forward3D(f_twist,ref3D,nx,ny,nz,nt,phase3D,maskt);  % forward propagation operator
            AT = @(g) Backward3D(g,ref3D,nx,ny,nz,nt,phase3D);  % backward propagation operator

            N1 = ny; N2 = nx*2*nz; N3 = 1;

            tau = 0.001;
            para = 0.1;
            piter = 4;
            Psi = @(f,th) MyTVpsi(f,th,para,piter,N1,N2,N3);
            Phi = @(f) MyTVphi(f,N1,N2,N3);

            tolA = 1e-6;
            iterations = 800+(t+d)*50;

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
            rec_psnr(t,d) = psnr(obj_rec,obj);
        end

        switch t
            case 10
                switch d
                    case 5
                        rec_t10_d5 = obj_rec;
                        obj_t10_d5 = obj;
                        img_t10_d5 = reshape(MyV2C(out_forward),ny,nx);
                    case 10
                        rec_t10_d5 = obj_rec;
                        obj_t10_d5 = obj;
                        img_t10_d5 = reshape(MyV2C(out_forward),ny,nx);
                end
            case 15
            switch d
                case 5
                    rec_t15_d5 = obj_rec;
                    obj_t15_d5 = obj;
                    img_t15_d5 = reshape(MyV2C(out_forward),ny,nx);
                case 10
                    rec_t15_d5 = obj_rec;
                    obj_t15_d5 = obj;
                    img_t15_d5 = reshape(MyV2C(out_forward),ny,nx);
            end
        end
    end
end
save('8.29_joint_reso2.mat','-v7.3');
