function y = V2C4D(x,n1,n2,n3,n4)
leng = n1*n2*n3;
vector = reshape(x,leng*2,n4);
re = vector(1:leng,:);
im = vector(leng+1:end,:);
for i = 1:n4
    com(:,i) = re(:,i) + 1i*im(:,i);
end
y = com(:);
end