function score = evaluateLine(vec)
% Evaluates the score of each horizontal, vertical, or diagonal line
x = 1;
o = -1;
b = 0;


score = 0;
% vec = vector of three elements
% Score three Xs as 100
if isequal(vec, [x x x])
    score = score + 100;
% Score two Xs and one blank as 10
elseif isequal(vec, [x x b]) ...
        || isequal(vec, [x b x]) ...
        || isequal(vec, [b x x])
    score = score + 10;
% Score one X and two blanks as 1
elseif isequal(vec, [x b b]) ...
        || isequal(vec, [b x b]) ...
        || isequal(vec, [b b x])
    score = score + 1;
    
% Score three Os as 100
elseif isequal(vec, [o o o])
    score = score - 100;
% Score two Os and one blank as 10
elseif isequal(vec, [o o b]) ...
        || isequal(vec, [o b o]) ...
        || isequal(vec, [b o o])
    score = score - 10;
% Score one O and two blanks as 1
elseif isequal(vec, [o b b]) ...
        || isequal(vec, [b o b]) ...
        || isequal(vec, [b b o])
    score = score - 1;
end