psnr_cs = [42.1	33.3689	31.8771	31.3584	30.4557	30.1042	29.7486
32.3443	29.8878	28.3839	28.1587	27.8517	27.4479	27.0594
30.7733	28.3389	27.4439	26.783	26.1941	25.9897	25.7327
28.774	26.9572	25.8392	25.8269	25.7834	25.779	25.7231];
psnr_bp = [21.94	21.52	21.5198	21.5181	21.5148	21.507	21.496
22.568	21.887	21.87	21.59	21.568	21.495	21.24
21.8	21.626	21.243	21.148	21.111	21.0585	20.687
21.75	21.75	21.67	21.37	21.33	21.25	21.04];
dz = [108000	54000	27000	13500	6750	3375	1687.5]*1e-3;
close all;
figure(1);
plot(dz,psnr_cs(1,:),'-o','linewidth',1.5,'color',[0.7 0 0],'markerfacecolor',[0.9 0 0]);hold on
plot(dz,psnr_cs(2,:),'-p','linewidth',1.5,'color',[0 0.7 0],'markerfacecolor',[0 0.8 0]);hold on
plot(dz,psnr_cs(3,:),'-d','linewidth',1.5,'color',[0 0 0.7],'markerfacecolor',[0 0 0.9]);hold on
plot(dz,psnr_cs(4,:),'-s','linewidth',1.5,'color',[1 0.5 0],'markerfacecolor',[1 0.7 0]);hold on

plot(dz,psnr_bp(1,:),'--o','linewidth',1.5,'color',[0.7 0 0]);hold on
plot(dz,psnr_bp(2,:),'--p','linewidth',1.5,'color',[0 0.7 0]);hold on
plot(dz,psnr_bp(3,:),'--d','linewidth',1.5,'color',[0 0 0.7]);hold on
plot(dz,psnr_bp(4,:),'--s','linewidth',1.5,'color',[1 0.5 0]);hold on
legend('100% CS','50%   CS','20%   CS','10%   CS','100% BP','50%   BP','20%   BP','10%   BP');
xlabel('Object spacing (mm)');
ylabel('PSNR');
axis([dz(7) dz(1) 20 45]);
axis square
for i = 1:5
    fig_obj1(:,:,i) = obj(:,:,1,i);
    fig_obj2(:,:,i) = obj(:,:,6,i);
    fig_obj1(:,:,i+5) = abs(backward(:,:,1,i));
    fig_obj2(:,:,i+5) = abs(backward(:,:,6,i));
    fig_obj1(:,:,i+10) = abs(rec_t5_d1(:,:,1,i));
    
end
figure(2);
colormap hot;imagesc(plotdatacube(fig_obj1,5));axis equal;axis off;
figure(3);
colormap gray;imagesc(plotdatacube(fig_obj2,5));axis equal;axis off;
figure(4);
colormap summer;imagesc(plotdatacube(fig_obj1(:,:,11:15),5));axis equal;axis off;
