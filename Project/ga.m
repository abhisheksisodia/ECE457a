function[xkx,xky,numofTurb] = ga(Lx,Ly,Pwt,D)
% Initalize GA parameters
popsize = 10; % Population size
MaxGen = 3; % Max number of generations/iterations
count = 0; % counter
nsite = 2; % number of mutation sites
pc = 0.95; % Crossover probability
pm = 0.05; % Mutation probability

% start off initail population with random values in the given range
xmin = 4.5;
xmax = 5.5;
xkrow = xmin+(xmax-xmin)*rand(1,popsize);
xkcol = xmin+(xmax-xmin)*rand(1,popsize);

p1 = xkrow;
p2 = xkcol; 

for j=1:MaxGen
p1;
p2;

count = count + 1;
% crossover parents
[child1, child2] = crossover(p1, p2);


% create generation
gen = [[child1, child2], p1, p2];
gen;
[r c] = size(gen);

% calculate objective function value of each
for k=1:c
  fx(1, k) = calcObjfun(gen(1,k),gen(1, k), 500 , 500, 70, 2.3);
end

[val1 idx1] = max(fx(:));
fx(1, idx1) = 0;
[val2 idx2] = max(fx(:));
fx(1, idx2) = 0;


newp1 = gen(1, idx1);
newp2 = gen(1, idx2);

p1 = newp1;
p2 = newp2;

end 

%final values for krow, kcol and numofTurb
xkx = p1;
xky = p2;
numofTurb = ((Lx/(xkx*D))+1)*((Ly/(xky*D))+1);


end


%function [xkrow, xkcol] = crossover(a)
function [c1, c2] = crossover(xkrow, xkcol)

    c1 = xkrow * 0.4 + xkcol * 0.6;
    c2 = xkrow * 0.6 + xkcol * 0.4;

end


% calculate objective function value
function [y] = calcObjfun(xkrow, xkcol, Lx, Ly, D, Pwt)
    N = ((Lx/(xkrow*D))+1)*((Ly/(xkcol*D))+1); %function to calculate number of turbines
    p = 8760*0.3*N*Pwt; %Power function = hy*nominal_power_util_factor*N*Pwt
    c = N*((2/3)+(1/3)*exp(-0.00174*(N^2))); %cost function
    y = p/c; %objective function
end