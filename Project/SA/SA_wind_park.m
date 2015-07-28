Lx = 2000;
Ly = 2000;
D = 82;
Pwt = 2;
T=1000;
alpha = 0.3; 
T_init =1000; % Initial temperature
T_min = 1e-10; % Final stopping temperature
F_min = -1e+100; % Min value of the function
max_rej=2500; % Maximum number of rejections
max_run=500; % Maximum number of runs
max_accept = 15; % Maximum number of accept
k = 1; % Boltzmann constant
i= 0; j = 0; accept = 0; totaleval = 0;

%Initial 
k_row=4.8+rand(1)*(5.2-4.8)
k_col=4.8+rand(1)*(5.2-4.8)

N=(Lx/(k_row*D)+1)*(Ly/(k_col*D)+1)
Power=2628*N*Pwt
Cost=N*(2/3+1/3*exp(-0.00174*(N^2)))
f_x=Power/Cost

while ((T > T_min) & (j <= max_rej))
    cur_tmp
    i = i+1;
    % Check if max numbers of run/accept are met
    if (i >= max_run) | (accept >= max_accept)
        % Cooling according to a cooling schedule
        T = T_init/(1 + (alpha * i));
        totaleval = totaleval + i;

        % reset the counters
        i = 1; accept = 1;
    end
    
    new_k_row=4.5+rand(1)*(5.5-4.5)
    new_k_col=4.5+rand(1)*(5.5-4.5)

    new_N=(Lx/(k_row*D)+1)*(Ly/(k_col*D)+1)
    new_Power=2628*N*Pwt
    new_Cost=N*(2/3+1/3*exp(-0.00174*(N^2)))
    new_f_x=Power/Cost 

    if new_f_x > f_x
        k_row = new_k_row
        k_col = new_k_col
        N = new_N 
        Power =new_Power
        Cost = new_Cost
        f_x = new_f_x
        accept=accept+1;
    else
        r = rand(1) 
        cur_tmp 
        Probability = exp(-1*(f_x-new_f_x)/T)
        if Probability > r
            k_row = new_k_row
            k_col = new_k_col
            N = new_N 
            Power =new_Power
            Cost = new_Cost
            f_x = new_f_x
            accept=accept+1;
        else
            j=j+1 
        end
    end
end
f_x
N
k_row
k_col