function y=MyTVphi(x,Nvx,Nvy,Nvz)
% x = x(1:length(x)/2) + 1i*x(length(x)/2+1:end);
X=reshape(x,Nvx,Nvy,Nvz);

[y,dif]=MyTVnorm(X);
% re = real(y); im = imag(y);
% y = [re;im];
function [y,dif]=MyTVnorm(x)

TV=MyTV3D_conv(x);

dif=sqrt(sum(TV.*conj(TV),4));

y=sum(dif(:));
end
end