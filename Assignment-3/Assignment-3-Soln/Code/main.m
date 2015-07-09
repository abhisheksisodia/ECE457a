n=6;
x = zeros(1,n);
for i=1:n
    [bs,~,~] = gasimple(i);
    x(i) = bs(end);
    
end