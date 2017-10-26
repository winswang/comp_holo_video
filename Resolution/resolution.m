clear all;close all;
%% Lateral resolution
% nx = 1000;ny = 1000;
% dx = 7.6; dy = 7.6;  %um
% Img = zeros(nx,ny);
% Img(498,500) = 1;
% Img(500,500) = 1;
% figure(1);imagesc(abs(Img));
% lambda = 532e-3;    %um green laser
% X0 = nx*dx; Y0 = ny*dy;
% z = 1000e3;
% I = ASM(Img, X0, Y0, z, lambda);
% figure(2);imagesc(abs(I));
% Ib = bpASM(I, X0, Y0, z, lambda);
% figure(3);imagesc(abs(Ib));
% Icrop = I(101:900,101:900);
% Ibb = bpASM(Icrop, 0.8*X0, 0.8*Y0, z, lambda);
% figure(4);imagesc(abs(Ibb));
%% Depth resolution
lambda = 532e-3;    %um
X0 = 128; Y0 = 128; % 
nx = 128; ny = 128;
z0 = 500;   % propagation distance 
dx = X0/nx; dy = Y0/ny; 
NA = (X0/2)/sqrt((X0/2)^2 + z0^2);
limit = lambda/(2*NA^2);    %
for k = 1:3
    switch k
        case 1  
            dz = limit;% separation of the points
        case 2  
            dz = limit/2;
        case 3  
            dz = limit/3;
    end
    canvas = zeros(nx, ny);
    canvas(nx/2,ny/2) = 1;
    Obj_volume(:,:,1) = canvas;Obj_volume(:,:,2) = canvas;
    % Kernel propagation
    phase = getPhase3D(X0, Y0, nx, ny, 2, z0, z0+dz, lambda);
    ref = getRef3D(X0, Y0, nx, ny, 2, z0, z0+dz, lambda);
    forward_out = Forward3D(MyC2V(Obj_volume(:)), ref, nx, ny, 2, phase);
    % reshape the output
    Image = reshape(MyV2C(forward_out), nx, ny);
    % simply backpropagate
    num = 100; period = 1.5;
    offset = 477;   di = 65/num;
    for i = 1:num
        f_back = bpASM(Image,X0,Y0,offset + di*(i-1),lambda);
        x(i) = offset + di*(i-1);
        I(i) = abs(f_back(nx/2,ny/2));
    end

    % with mask
    mask1 = genmask(nx,ny,1,0.5);
    Image_mask1 = Image.*mask1;
    for i = 1:num
        f_back = bpASM(Image_mask1,X0,Y0,offset + di*(i-1),lambda);
        x_mask1(i) = offset + di*(i-1);
        I_mask1(i) = abs(f_back(nx/2,ny/2));
    end
    mask2 = genmask(nx,ny,1,0.2);
    Image_mask2 = Image.*mask2;
    for i = 1:num
        f_back = bpASM(Image_mask2,X0,Y0,offset + di*(i-1),lambda);
        x_mask2(i) = offset + di*(i-1);
        I_mask2(i) = abs(f_back(nx/2,ny/2));
    end
    % draw backpropagation
    figure;
    plot(x,I,'linewidth',2);hold on
    plot(x_mask1,I_mask1,'linewidth',2);hold on
    plot(x_mask2,I_mask2,'linewidth',2);hold on
    x_obj = z0*ones(1,20);y_obj = linspace(0,1,20);line(x_obj,y_obj,'color','k');
    x_obj = (z0 + dz)*ones(1,20);y_obj = linspace(0,1,20);line(x_obj,y_obj,'color','k');
    axis([min(x) max(x) 0 1]); title(strcat('back-propagation (dz =',num2str(dz),'\mum)'));
    xlabel('Distance'); legend('No mask','50% random mask', '20% random mask', 'Object locations');
    drawnow;

    %% compressed sensing reconstruction
    nz = 100; % construct 30 layers
    zmin = 477; zmax = 477+65;
    phase = getPhase3D(X0, Y0, nx, ny, nz, zmin, zmax, lambda);
    ref = getRef3D(X0, Y0, nx, ny, nz, zmin, zmax, lambda);
    % Propagation operator
    mask = ones(nx,ny);
    A = @(f_twist) Forward3D(f_twist,ref,nx,ny,nz,phase,mask);  % forward propagation operator
    AT = @(g) Backward3D(g,ref,nx,ny,nz,phase);  % backward propagation operator

    %% TwIST algorithm (no mask)
    % twist parameters
    tau = 0.001; 
    piter = 4;
    tolA = 1e-6;
    iterations = 800;

    g = MyC2V(Image(:));

    Psi = @(f,th) MyTVpsi(f,th,0.05,piter,nx,ny*nz*2,1);
    Phi = @(f) MyTVphi(f,nx,ny*nz*2,1);

    [f_rec,dummy,obj_twist,...
        times_twist,dummy,mse_twist]= ...
        TwIST(g,A,tau,...
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

    f_rec = reshape(MyV2C(f_rec), nx, ny, nz);
    x_rec = linspace(zmin, zmax, nz);
    I_rec = squeeze(abs(f_rec(nx/2,ny/2,:)));

    % Propagation operator
    A = @(f_twist) Forward3D(f_twist,ref,nx,ny,nz,phase,mask1);  % forward propagation operator
    AT = @(g) Backward3D(g,ref,nx,ny,nz,phase);  % backward propagation operator

    %% TwIST algorithm (50% mask)

    g = MyC2V(Image_mask1(:));

    Psi = @(f,th) MyTVpsi(f,th,0.05,piter,nx,ny*nz*2,1);
    Phi = @(f) MyTVphi(f,nx,ny*nz*2,1);

    [f_rec,dummy,obj_twist,...
        times_twist,dummy,mse_twist]= ...
        TwIST(g,A,tau,...
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

    f_rec = reshape(MyV2C(f_rec), nx, ny, nz);
    I_rec1 = squeeze(abs(f_rec(nx/2,ny/2,:)));

    % Propagation operator
    A = @(f_twist) Forward3D(f_twist,ref,nx,ny,nz,phase,mask2);  % forward propagation operator
    AT = @(g) Backward3D(g,ref,nx,ny,nz,phase);  % backward propagation operator

    %% TwIST algorithm (20% mask)

    g = MyC2V(Image_mask2(:));

    Psi = @(f,th) MyTVpsi(f,th,0.05,piter,nx,ny*nz*2,1);
    Phi = @(f) MyTVphi(f,nx,ny*nz*2,1);

    [f_rec,dummy,obj_twist,...
        times_twist,dummy,mse_twist]= ...
        TwIST(g,A,tau,...
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

    f_rec = reshape(MyV2C(f_rec), nx, ny, nz);
    I_rec2 = squeeze(abs(f_rec(nx/2,ny/2,:)));

    % draw reconstruction
    figure;
    plot(x_rec,I_rec,'linewidth',2);hold on
    plot(x_rec,I_rec1,'linewidth',2);hold on
    plot(x_rec,I_rec2,'linewidth',2);hold on
    x_obj = z0*ones(1,20);y_obj = linspace(0,1,20);line(x_obj,y_obj,'color','k');
    x_obj = (z0 + dz)*ones(1,20);y_obj = linspace(0,1,20);line(x_obj,y_obj,'color','k');
    axis([min(x_rec) max(x_rec) 0 1]); title(strcat('CS (dz =',num2str(dz),'\mum)'));
    xlabel('Distance'); legend('No mask','50% random mask', '20% random mask');
end




