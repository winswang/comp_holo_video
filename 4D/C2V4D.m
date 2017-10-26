function y = C2V4D(x)
% x is a complex 4D data
Real = real(x);
Imag = imag(x);
[n1,n2,n3,n4] = size(x);
rr = reshape(Real,n1*n2*n3,n4);
ii = reshape(Imag,n1*n2*n3,n4);
for i = 1:n4
    vector(:,:,i) = [rr(:,i);ii(:,i)];
end
y = vector(:);
end