% Genetic Algorithm (Simple Demo) Matlab/Octave Program
% Written by X S Yang (Cambridge University)
% Usage: gasimple or gasimple('x*exp(-x)');

function [bestsol, bestfun, count] = gasimple()
global solnew sol pop popnew fitness fitold f range dimension;

dimension = 10; % Dimensions
range = [-1 1]; % Range/Domain
f = @func; % Function

% Initializing the parameters
rand('state', 0'); % Reset the random generator
popsize = 20; % Population size
MaxGen = 100; % Max number of generations
count = 0; % counter
nsite = 2; % number of mutation sites
pc = 0.95; % Crossover probability
pm = 0.05; % Mutation probability
nsbit = 16; % String length (bits)

% Generating the initial population
% rows: dimension
% columns: string length
% height: instance
% ie. 
%  a(:,:,1) =
%
%     0     0     0
%     1     0     0
%
%
% a(:,:,2) =
%
%     1     0     0
%     1     1     1
%
popnew = init_gen(dimension, nsbit, popsize);
fitness = zeros(1, popsize); % fitness array

% Display the shape of the function
x = range(1):0.1:range(2); % row vector between ranges with 0.1 steps
plot(x, f(x));

% Initialize solution <- initial population
% Convert the 16 bit string to decimal
for i = 1:size(popnew, 3),
    for j = 1:size(popnew, 1),
        solnew(j,i) = bintodec(popnew(j,:,i));
    end
end

% Start the evolution loop
for i = 1:MaxGen,
    % Record as the history
    fitold = fitness; 
    pop = popnew; 
    sol = solnew;
    
    for j = 1:popsize,
        % Crossover pair
        ii = floor(popsize * rand) + 1; 
        jj = floor(popsize * rand) + 1;
        
        % Cross over
        if pc > rand,
            [popnew(:, :, ii), popnew(:, :, jj)] = crossover(pop(:, :, ii), pop(:, :, jj));
            % Evaluate the new pairs
            count = count + 2;
            evolve(ii); 
            evolve(jj);
        end
     
        % Mutation at n sites
        if pm > rand,
            kk = floor(popsize * rand) + 1; count = count + 1;
            popnew(:, :, kk) = mutate(pop(:, :, kk), nsite);
            evolve(kk);
        end
    end % end for j
    
    % Record the current best
    bestfun(i) = min(fitness);
    bestsol(i) = mean(sol(bestfun(i) == fitness));
end

% Display results
set(gcf, 'color', 'w');
subplot (2, 1, 1); plot(bestsol); title('Best estimates');
subplot(2, 1, 2); plot(bestfun); title('Fitness');

% All the sub functions

% function func
function [dec] = func(x)
    dec = 0;
    for i = 1:size(x),
        dec = dec + abs(x(i))^(i+1);
    end;

% generation of the initial population
% create 
function pop = init_gen(ndimension, nsbit, npop)
    % String length=nsbit+l with pop(:,l) for the Sign
    pop = rand(ndimension, nsbit + 1, npop) > 0.5;

% Evolving the new generation
function evolve(j)
    global solnew popnew fitness fitold pop sol f;
    
    % loop through all the dimensions
    for k = 1:size(popnew, 1),
        solnew(k, j) = bintodec(popnew(k, :, j));
    end
    
    fitness(j) = 1/(1 +f(solnew(:, j)));
    
    if fitness(j) > fitold(j),
        pop(:, :, j) = popnew(:, :, j);
        sol(:, j) = solnew(:, j);
    end
    
% Convert a binary string into a decimal number
function [dec] = bintodec(bin)
    global range;
    % Length of the string without sign
    nn = length(bin) - 1;
    num = bin(2:end); % get the binary
    % Sign=+1 if bin(l)=0; Sign=-l if bin(l)=l.
    Sign = 1 - 2 * bin(1);
    dec = 0;
    % floating point/decimal place in a binary string
    dp = floor(log2(max(abs(range))));
    for i = 1:nn,
        dec = dec + num(i) * 2 ^ (dp - i);
    end
    dec = dec * Sign;
    
% Crossover operator
function [c, d] = crossover(a, b)
    for i = 1:size(a, 1),
        row_a = a(i,:);
        row_b = b(i,:);
        nn = length(row_a) - 1;
        % generating a random crossover point
        cpoint = floor(nn * rand) + 1;
        c(i,:) = [row_a(1:cpoint) row_b(cpoint + 1:end)];
        d(i,:) = [row_b(1:cpoint) row_a(cpoint + 1:end)];       
    end

% Mutatation operator
function anew = mutate(a, nsite)
    for i = 1:size(a,1),
        row_a = a(i,:);
        nn = length(row_a);
        anew = row_a;
        for j = 1:nsite,
            k = floor(rand * nn) + 1;
            anew(i,k) = mod(row_a(k) + 1, 2);
        end
    end
