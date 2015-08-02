Lx = 500;
Ly = 500;
D = 70;
Pwt = 2.3;
xkrow = zeros(1,1);
xkcol = zeros(1,1);
N = zeros(1,1);

for i=1:1
    [xkrow(i,1),xkcol(i,1),N(i,1)] = aco(Lx,Ly,Pwt,D);
end

xkx = mean(xkrow(:))
xky = mean(xkcol(:))
num = mean(N(:))