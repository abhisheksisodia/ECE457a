function[xkx,xky,numofTurb] = sa(Lx,Ly,Pwt,D)
    alpha = 0.3; %Geometric cooling factor
    T_init =100; % Initial temperature
    T_min = 1e-10; % Final stopping temperature
    T=T_init; %Iterator variable
    max_rej=2500; % Maximum number of rejections
    max_run=500; % Maximum number of runs
    max_accept = 15; % Maximum number of accept
    k = 1; % Boltzmann constant
    i= 0; j = 0; accept = 0; totaleval = 0;

    %Initial 
    xmin = 4.5;
    xmax = 5.5;
    k_row=xmin+rand(1)*(xmax-xmin);
    k_col=xmin+rand(1)*(xmax-xmin);

    N=(Lx/(k_row*D)+1)*(Ly/(k_col*D)+1);
    Power=2628*N*Pwt;
    Cost=N*(2/3+1/3*exp(-0.00174*(N^2)));
    f_x=Power/Cost;

    while ((T > T_min) & (j <= max_rej))
        i = i+1;
        % Check if max numbers of run/accept are met
       
        %the maximum accept value increases as the temperature decreases so
        %as to increase intensification. This makes the algorithm adaptive
        if (i >= max_run) | (accept >= max_accept*(T_init+1-T)) 
            % Cooling according to a cooling schedule
            T = T_init/(1 + (alpha * i));
            totaleval = totaleval + i;

            % reset the counters
            i = 1; accept = 1;
        end

        new_k_row=xmin+rand(1)*(xmax-xmin);
        new_k_col=xmin+rand(1)*(xmax-xmin);

        new_N=(Lx/(k_row*D)+1)*(Ly/(k_col*D)+1);
        new_Power=2628*N*Pwt;
        new_Cost=N*(2/3+1/3*exp(-0.00174*(N^2)));
        new_f_x=Power/Cost;

        if new_f_x > f_x
            k_row = new_k_row;
            k_col = new_k_col;
            N = new_N ;
            Power =new_Power;
            Cost = new_Cost;
            f_x = new_f_x;
            accept=accept+1;
        else
            r = rand(1);
            Probability = exp(-1*(f_x-new_f_x)/T);
            if Probability > r
                k_row = new_k_row;
                k_col = new_k_col;
                N = new_N;
                Power =new_Power;
                Cost = new_Cost;
                f_x = new_f_x;
                accept=accept+1;
            else
                j=j+1;
            end
        end
    end
    f_x
    numofTurb = N
    xkx = k_row
    xky = k_col
end 