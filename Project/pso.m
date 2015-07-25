function[xkx,xky,numofTurb] = pso(Lx,Ly,Pwt,D)
%Initialize the c1, c2, w values
c1 = 1.4944;
c2 = 0.25;
w = 0;
numOfIt = 10000;

%Initially generate random values for krow and kcol positions
xmin = 4.5;
xmax = 5.5;
xkrow = xmin+(xmax-xmin)*rand(4,1);
xkcol = xmin+(xmax-xmin)*rand(4,1);

%Initialize velocity of particles to zero
vkrow = zeros(4,1);
vkcol = zeros(4,1);
fx = zeros(4,1);

for j=1:numOfIt
    %generate r1 and r2 between 0 and 1 for each iteration
    r1 = rand();
    r2 = rand();
    %Calculate the objective function value for each of particle's position values
    for k=1:4
        fx(k,1) = calcObjfun(xkrow(k,1),xkcol(k,1),Lx,Ly,D,Pwt);
    end

    %find highest objective function row
    [num idx] = max(fx(:));
    [a b] = ind2sub(size(fx),idx);
    
    %update pbest and gbest
    pbestkrow = xkrow;
    pbestkcol = xkcol;
    
    %find corresponding particle positions which give the highest obj
    %function value
    gbestkrow = xkrow(a,1);
    gbestkcol = xkcol(a,1);
    
    % update velocity and position
    for i=1:4
       vkrow(i,1) = (w*vkrow(i,1))+(c1*r1*(pbestkrow(i,1)-xkrow(i,1)))+(c2*r2*(gbestkrow-xkrow(i,1)));
       vkcol(i,1) = (w*vkcol(i,1))+(c1*r1*(pbestkcol(i,1)-xkcol(i,1)))+(c2*r2*(gbestkcol-xkcol(i,1)));
       xkrow(i,1) = xkrow(i,1) + vkrow(i,1);
       xkcol(i,1) = xkcol(i,1) + vkcol(i,1);
    end
end

    xkx = mean(xkrow(:));
    xky = mean(xkcol(:));
    numofTurb = ((Lx/(xkx*D))+1)*((Ly/(xky*D))+1);
end

function[y] = calcObjfun(xkrow,xkcol,Lx,Ly,D,Pwt)
    N = ((Lx/(xkrow*D))+1)*((Ly/(xkcol*D))+1); %function to calculate number of turbines
    p = 8760*0.3*N*Pwt; %Power function = hy*nominal_power_util_factor*N*Pwt
    c = N*((2/3)+(1/3)*exp(-0.00174*(N^2))); %cost function
    y = p/c; %objective function
end

