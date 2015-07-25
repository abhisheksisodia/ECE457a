Lx = 2000;
Ly = 2000;
D = 70;
Pwt = 2.3;
xkrow = zeros(10,1);
xkcol = zeros(10,1);
N = zeros(10,1);

for i=1:100
    [xkrow(i,1),xkcol(i,1),N(i,1)] = pso(Lx,Ly,Pwt,D);
end

xkx = mean(xkrow(:))
xky = mean(xkcol(:))
num = mean(N(:))
