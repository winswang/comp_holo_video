function volume = place_obj(object,nx,ny,nz,nt,ratio)
% this function place the object to each layer and 
nz = 5*nz-4;
volume = zeros(ny,nx,nz,nt);
offset = 40;
for i = 1:5:nz
    obj = imresize(object,ratio*0.9^(i-1));
    [oy,ox] = size(obj);
    offset = offset + ox*(i-2);
    canvas = zeros(ny,nx);
    canvas(1:oy,1:ox) = obj;
    canvas = circshift(canvas,[20 round(nx*0.3) + offset]);
    speed = floor((ny-oy)/nt*0.7^(i-1));
    for j = 1:nt
        move = circshift(canvas,[speed*(j-1) 0]);
        volume(:,:,i,j) = move;
    end
end
end