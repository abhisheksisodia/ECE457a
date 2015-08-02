function[xkx,xky,numofTurb] = ts(Lx,Ly,Pwt1,D1)

	global L_x L_y D Pwt rows columns;
	D=D1;
    Pwt=Pwt1;
	L_x = Lx;
	L_y = Ly;
	
    tabu_tenure = 3;
    num_iterations = 20;
	
	rows = linspace(4.5,5.5,1000);
	columns = linspace(4.5,5.5,1000);
	
	[initial_solution, initial_cost] = getInitialSolution();
	
	[~, num_columns] = size(rows);
	num_nodes = num_columns;
    
    % Tabu list intialization to zero.
    tabu_list = initializeTabuList(num_nodes);
    
    current_solution = initial_solution; % The current solution is the intial solution.
    best_solution = initial_solution;    % The best solution is the intial solution.
    best_cost = initial_cost;            % The min. cost is the cost of the intial solution. 

    for i=1:num_iterations
        % Obtaining best neighboring solution.
        [best_neighborhood_solution, best_neighborhood_cost, tabu_list] = getNeighborhood (current_solution, best_cost, tabu_list, tabu_tenure);
        
        % The obtained neighboring solution is the current solution.
        current_solution = best_neighborhood_solution;
        
    	% Comparing if it's the best solution so far.
        if (best_neighborhood_cost < best_cost)
            best_solution = best_neighborhood_solution;
            best_cost = best_neighborhood_cost;
        end 
    end
    
    xkx = best_solution(1);
    xky = best_solution(2);
    k_row=xkx;
    k_col=xky;
    N_row = L_x/(k_row*D) + 1;
	N_col = L_y/(k_col*D) + 1;
	
	numofTurb = N_row * N_col;
end

function [best_neighborhood_solution, best_neighborhood_cost, tabu_list] = getNeighborhood (current_solution, best_cost, tabu_list, tabu_tenure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that obtains the best neighboring solution by realzing a set of
% modifications in the current solution. 
% Generate different neighboring solutions by swapping and keep the solution with min. cost.
% Input parameters:
%   - current_solution: current solution.
%   - best_cost: cost of the best solution obatined so far.
%   - tabu_tenure: number of iterations a swapping will be in the tabu list.
%   - tabu_list: Tabu list.
%   - cost_matrix: Mcost matrix.
% Output parameters:
%   - best_neighborhood_solution: best neighboring solution obtained.
%   - best_neighborhood_cost: cost of the best neighboring solution obatined.
%   - tabu_list: Tabu list.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	best_neighborhood_solution = [];
    best_neighborhood_cost = 0;
    best_node1 = 0;
    best_node2 = 0;
	min_cost = best_cost;
    
    for i=1:5;
        %obtain 5 neighbouring solutions.
        [neighborhood_solution, neighborhood_cost] = getInitialSolution();
        
        %swapping_frecuency = getSwappingFrecuency(neighborhood_solution(1),neighborhood_solution(2),tabu_list);
        %diversification_cost = neighborhood_cost + swapping_frecuency;
        
        % Checking that the solution is not in the tabu.
        if (isTabu(neighborhood_solution(1),neighborhood_solution(2),tabu_list) == false)
            % Is not a tabued solution.
            if ((neighborhood_cost < min_cost))
                % best solution obatined so far
                min_cost = neighborhood_cost;
                best_neighborhood_solution = neighborhood_solution;
                best_neighborhood_cost = neighborhood_cost;
                best_node1 = neighborhood_solution(1);
                best_node2 = neighborhood_solution(2);
            else
                best_neighborhood_solution=current_solution;
                best_neighborhood_cost=best_cost;
            end
        else
            % Is a tabued solution:
            % To avoid stagnation, apply asipration criteria
            % If the cost of the solution is less that the cost of best solution obatined so far, accept this solution
            if (diversification_cost < best_cost)
                % best solution obatined so far.
                min_cost = diversification_cost;
                best_neighborhood_solution = neighborhood_solution;
                best_neighborhood_cost = neighborhood_cost;
                best_node1 = neighborhood_solution(1);
                best_node2 = neighborhood_solution(2);
            else
                % Tabued solution is not permited.
            end
        end
    end
    
	% after obtaining the best neighborhood solution, update the tabu list
    tabu_list = updateTabuList(tabu_list);
    tabu_list = addSwappingTabuList(best_node1,best_node2,tabu_list,tabu_tenure);
	
end

function [initial_solution, initial_cost] = getInitialSolution()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to obtain the intial solution selecting the most nearest nodes.
% Input parameters:
%   - D: diameter.
%   - Pwt: power coefficient.
% Output parameters:
%   - initial_path: initial path.
%   - total_cost: total cost.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	global D Pwt;
	
    k_row = 4.5+rand();
	k_col = 4.5+rand();
    
	initial_solution = [k_row k_col D Pwt];
    initial_cost = CalculateCostFunc(initial_solution);

end

% calculate objective function
function obj_func = CalculateCostFunc(sol)
	global L_x L_y D Pwt;
	
	k_row = sol(1);
	k_col = sol(2);
	
	N_row = L_x/(k_row*D) + 1;
	N_col = L_y/(k_col*D) + 1;
	
	N = N_row * N_col;
	
	C = N*((2/3)+(1/3)*exp(-0.00174*N^2));
	P = 8760*0.3*N*Pwt;
	
	obj_func = C/(1+P) * 10000;
end


function [tabu_list] = initializeTabuList (num_nodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funciton to intialize the memory structure of the tabu list to zero 
% Input parameters:
%   - num_nodes: Number of nodes.
% Output parametersa:
%   - tabu_list: tabu list intialized to zero.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    for i=1:num_nodes
        for j=1:num_nodes
            tabu_list(i,j) = 0;
        end
    end
end

function [tabu_list] = updateTabuList (tabu_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to update tabu list in every iteration
% Input parameters:
%   - tabu_list: tabu list.
% Output parameters:
%   - tabu_list: Updated tabu list.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [num_rows,num_columns] = size(tabu_list);
    num_nodes = num_rows;
    
    for i=1:num_nodes
        for j=i+1:num_nodes
            if (tabu_list(i,j) > 0)
                tabu_list(i,j) = tabu_list(i,j) - 1;
            end
        end
    end
end

function [tabu_list] = addSwappingTabuList (node1, node2, tabu_list, tabu_tenure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to add a swapping to the tabu list controlling the frequency of this swapp.
% Input parameters:
%   - node1, node2: swapping nodes.
%   - tabu_list: tabu list.
%   - tabu_tenure: tabu tenure.
% Output parameters:
%   - tabu_list: Updated tabu list.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global rows columns;
    
	[row1 column1] = find(rows==node1);
	[row2 column2] = find(columns==node2);
	
    if (column1 < column2)
        % Inicializamos al valor de tabu_tenure el tiempo que va a estar el
        % swapping en la lista tabú.
        tabu_list(column1,column2) = tabu_tenure;
        % Incrementamos el número de veces que hemos hecho este swapping.
        tabu_list(column1,column2) = tabu_list(column2,column1) + 1;
    else
        % Inicializamos al valor de tabu_tenure el tiempo que va a estar el
        % swapping en la lista tabú.
        tabu_list(column2,column1) = tabu_tenure;
        % Incrementamos el número de veces que hemos hecho este swapping.
        tabu_list(column1,column2) = tabu_list(column1,column2) + 1;
    end
end

function [is_tabu_solution] = isTabu (node1, node2, tabu_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to check the solution is listed as tabued. 
% Input parameters:
%   - node1, node2: swapping nodes.
%   - tabu_list: tabu list.
% Output parameters:
%   - is_tabu_solution: boolean variable to indicate the solution is taubed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    is_tabu_solution = false;
    
	global rows columns;
	
	[row1 column1] = find(rows==node1);
	[row2 column2] = find(columns==node2);
    
    if (column1 < column2)
        % Checking the swap (node1,node2) in tabu list.
        if (tabu_list(column1,column2) > 0)
            % Is a tabued solution.
            is_tabu_solution = true;
        end
    else
        % Checking the swap (node2,node1) in tabu list.
        if (tabu_list(column2,column1) > 0)
            % Is a tabued solution.
            is_tabu_solution = true;
        end
    end
end

function [swapping_frecuency] = getSwappingFrecuency (node1, node2, tabu_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to obtain the sawpping frequency
% Input parameters:
%   - node1, node2: swapping nodes.
%   - tabu_list: tabu list.
% Output parameters:
%   - is_tabu_solution: boolean variable to indicate the solution is taubed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global rows columns;
	[row1 column1] = find(rows==node1);
	[row2 column2] = find(columns==node2);
	node1 = column1;
	node2 = column2;
    if (node1 < node2)
        % Obtaining the value of the swap (node2,node1) from the tabu list.
        swapping_frecuency = tabu_list(node2,node1);
    else
        % Obtaining the value of the swap (node1,node2) from the tabu list.
        swapping_frecuency = tabu_list(node1,node2);
    end
end
