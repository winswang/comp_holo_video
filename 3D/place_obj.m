function volume = place_obj(object,nx,ny,nz,nt,ratio)
% this function place the object to each layer and 
volume = zeros(ny,nx,nz,nt);
offset = 0;
for i = 1:nz
    obj = imresize(object,ratio*0.7^(i-1));
    [oy,ox] = size(obj);
    offset = offset + ox*(i-1);
    canvas = zeros(ny,nx);
    canvas(1:oy,1:ox) = obj;
    canvas = circshift(canvas,[20 ox*0.3 + offset]);
    speed = floor((ny-oy)/nt*0.7^(i-1));
    for j = 1:nt
        move = circshift(canvas,[speed*(j-1) 0]);
        volume(:,:,i,j) = move;
    end
end
end