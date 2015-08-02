function [xkx,xky,numofTurb] = aco(Lx,Ly,Pwt,D)
%aco using ACS algorithm to optimize the placements of wind turbines

    % getting information
    [x, y, xgrid, ygrid, Ndiv, n, dim, iter, m, t, h, alpha, beta, e, el] = ants_information;

    % create colonies
    Ncolony = 3;
    for i = 1:Ncolony
      colony(i).m = m;
      colony(i).t = t;
      colony(i).h = h;
      colony(i).alpha = alpha;
      colony(i).beta = beta;
      colony(i).e = e;
      colony(i).el = el;
    end
    
    for i = 1:iter
      for j = 1:Ncolony
        colony(j).app = ants_primaryplacing(colony(j).m, dim, n);  
        colony(j).at = ants_cycle(colony(j).app, colony(j).m, n, dim, ...
            colony(j).h, colony(j).t, colony(j).alpha, colony(j).beta);

        % calculate the cost
        [colony(j).cost, colony(j).f] = ants_cost(colony(j).m, xgrid, ygrid, Ndiv, ...
            colony(j).at, el, Lx, Ly, D);

        colony(j).t = ants_traceupdating(n,colony(j).m, dim, colony(j).t, ...
            colony(j).at, colony(j).f, colony(j).e);

        colony(j).costoa(i) = mean(colony(j).cost);
        [colony(j).maxcost(i),  colony(j).number] = max(colony(j).cost);    
        colony(j).besttour(i, :) =  colony(j).at(colony(j).number, :);
        colony(j).iteration(i) = i;
      end
    end
    
%     subplot(2, 1, 1); 
%     plot(colony(1).iteration,  mean(vertcat(colony.maxcost),1));
%     title('Best Power/Cost vs. Iteration');
%     xlabel('Iteration');
%     ylabel('Power/Cost');

%     [k, l] = max([colony.maxcost]);
%     for i = 1:(dim+2)
%       besttour = vertcat(colony.besttour);
%       indx = besttour(l, i);
%       X(i) = x(indx);
%       Y(i) = y(indx);
%     end

%     subplot(2, 1, 2); 
%     plot(X, Y, '--rs', 'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10)
%     xlabel('X'); ylabel('Y'); axis('equal');
    
%     for i = 1:(dim+2)
%       besttour = vertcat(colony.besttour);
%       indx = besttour(l, i);
%       if i == 2 
%         value = xgrid(indx - 1);        
%         text(X(i) + .5, Y(i), ['\leftarrow node ', num2str(indx), ', value ', num2str(value)]);
%       elseif i == 3
%         value = ygrid(indx - 1 - Ndiv);
%         text(X(i) + .5, Y(i), ['\leftarrow node ', num2str(indx), ', value ', num2str(value)]);
%       else
%         text(X(i) + .5, Y(i), ['\leftarrow node ', num2str(indx)]);
%       end
%     end
%     title(['optimum course -- MAX(Power/Cost) = ', num2str(k)]);
    
    
    [k, l] = max([colony.maxcost]);
    besttour = vertcat(colony.besttour);
    xkx = xgrid(besttour(l, 2) - 1);
    xky = ygrid(besttour(l, 3) - 1 - Ndiv);
    numofTurb = ((Ly/(xky*D))+1)*((Lx/(xkx*D))+1);
end

% Ant information
function [x, y, xgrid, ygrid, Ndiv, n, dim, iter, m, t, h, alpha, beta, e, el] = ants_information;

  iter = 500; % number of cycles.
  m = 20; % number of ants.

  % points of x and y
  %x = [8 0 -1 2 4 6 3 10 2.5 -5 7 9 11 13];
  %y = [2 4 6 -1 -2 0.5 0 3.7 1.8 1 0 4 3 2];
  %n = length(x); % number of nodes.
  
  % range= [xmin xmax ymin ymax] ;
  range=[4.5 5.5 4.5 5.5];
  dim = 2;
  % Grid values
  Ndiv=100;
  dx=(range(2)-range(1))/(Ndiv - 1);
  dy=(range(4)-range(3))/(Ndiv - 1);
  xgrid=range(1):dx:range(2);
  ygrid=range(3):dy:range(4);
  n = 2 * Ndiv + 2;
  
  % adjacency matrix
  % root x_1 x_2 .. x_n y_1 y_2 .. y_n dest
  % n x n  
  
  % for graphical reasons
  x = horzcat(0, 0:Ndiv-1, 0:Ndiv-1, 0);
  y = horzcat(0, ones(1,Ndiv), ones(1,Ndiv)*2, 3);

  e = 0.000001; % evaporation coefficient.
  alpha = 1; % order of trace?s effect.
  beta = 0.1; % order of effect of ants? sight.

  % h - desirability (start node, end node)
  for i = 1:n % generating sight matrix.
    for j = 1:n
      if (i == 1 && j == 1) || ...
         (i == n && j == n) || ...
         (i == j)
            h(i, j) = 0;         
      elseif (i == 1 && j >= 2 && j < (2+Ndiv)) || ...
             (i >= 2 && i < (2+Ndiv) && j >= (2+Ndiv) && j < (2+2*Ndiv)) || ...
             (i >= (2+Ndiv) && i < (2+2*Ndiv) && j == n)
                h(i, j) = 1;
      else
          h(i, j) = 0;
      end
    end
  end

  % t - pheronomes
  t = 0.0001 * ones(n); % primary tracing.
  el = .96; % coefficient of common cost elimination.
end
  
  %------------------End of ant_information------------
  
function [app] = ants_primaryplacing(m, dim, n);
  % seed the rand
  rand('state', sum(100 * clock));
  %rand('state', 0);
  
  % m - number of ants
  % app - 2D Matrix (ant number, node number)
  % choose a random first node
  for i = 1:m
    % random node between 1 to n    
    % n - number of node
    % app(i, 1) = fix(1 + rand * (n - 1));
    
    % all ants starts at node 1
    app(i, 1) = 1;
    app(i, 2 + dim) = n;
  end
end
  %-------------- End of ants_primaryplacing ------------
  
function [at]=ants_cycle(app,m,n,dim,h,t,alpha,beta);
  % m - number of ants
  % app - 2D Matrix (ant number, node number)
  for i = 1:m
    % h - desirability
    % save the desirability to mh
    mh = h;
    
    % for the rest of the nodes 1 to n - 1
    % n - number of node
    %for j = 1:n - 1
    
    % select number of dimensions + end node    
    for j = 1:(dim + 1)
      % get the current node
      c = app(i, j);
      % set the desirability to 0 for that node so that it will not go
      % there again
      mh(:, c) = 0;
      % calculate the numerator of the transition prob
      temp = (t(c, :) .^ alpha) .* (mh(c, :) .^ beta);
      % calculate the denominator of the transition prob
      s = (sum(temp));
      % calculate the transition probability
      p = (1 / s) .* temp;
      
      % select a rand number for at roulette wheel
      r = rand;
      % select a rand number for ACS
      r_o = 0.5;
      
      r_acs = rand;
      if (r_acs <= r_o)
        max_value = max(temp);
        max_indices = find(temp == max_value);
        k = max_indices(round(rand*(length(max_indices)-1)+1));
        app(i, j + 1) = k;
      end
      
      % if the max cant find the best indices or the rand is greater
      if (r_acs > r_o) 
        cum = 0;
        % roulette probability
        for k = 1:n
          cum = cum + p(k);
          % if the cumulative sum is greater than rand, then that's the node
          % if p(k) == 0, then that s will not be selected
          if cum > r          
            app(i, j + 1) = k;
            break
          end
        end
      end
    end
  end
  at = app; % generation of ants tour matrix during a cycle.
end
  %---------------- End of ants_cycle -----------------
  
function [cost, f] = ants_cost(m, xgrid, ygrid, Ndiv, at, el, L_x, L_y, D);
  % calculating the cost of all ants
  % m - number of ants
  for i = 1:m    
    % sum up the total distance of the paths that each and every ant
    %s = 0;
    %for j = 1:n
    %  s = s + d(at(i, j), at(i, j + 1));
    %end
    
    % calculate the values
    % first convert the paths to values
    % root x_1 x_2 .. x_n y_1 y_2 .. y_n dest to values = [x, y]
    values(1) = xgrid(at(i,2) - 1);    
    values(2) = ygrid(at(i,3) - 1 - Ndiv);
    
    % calculate the value
    % func(h_y, eta, L_x, L_y, D, P_wt, k_row, k_col)
    s = func(8760,0.3,L_x,L_y,D,2.31, values(1), values(2));
    % fitness of the ant is equal to s
    f(i) = s;
  end

  cost = f;
  
  % each fitness is subtract by the minimum of all fitness * a constant
  f = f - el * min(f); % elimination of common cost.
end
  
function [dec] = func(h_y, eta, L_x, L_y, D, P_wt, k_row, k_col);
  N = ((L_x/(k_row*D))+1)*((L_y/(k_col*D))+1);
  dec = h_y*eta*P_wt;
  dec = dec/(2/3 + 1/3 * exp(-0.00174*N^2));
end
%------------------ End of ants_cost----------------

function [t] = ants_traceupdating(n, m, dim, t, at, f, e);
  % find best ant
  [max_value, max_ant] = max(f);
  
  % evaporate all the paths
  for i = 1:n % generating sight matrix.
    for j = 1:n
      t(i, j) = (1 - e) * t(i, j);    
    end
  end
  
  % m - number of ants
  for i = 1:m    
    % n - number of nodes
    for j = 1:(dim+1)
      if (i == max_ant)
        % Q*L_k -, L_k = f(i)
        Q = 1/2500000;
        dt = Q * f(i);
      else
        dt = 0;
      end
      % evaporate and deposit dt
      t(at(i, j), at(i, j + 1)) = t(at(i, j), at(i, j + 1)) + dt; % updating traces.      
    end
  end
end
  %-------------- End of ants traceupdating --------------