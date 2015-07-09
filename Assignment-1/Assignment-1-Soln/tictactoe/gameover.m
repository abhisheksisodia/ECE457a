function over = gameover(board)
% Returns 1 or 0

over = false; % Default position.


x = 1;
o = -1;
b = 0;


bvec = board(:);% Bvec is the node written as a vector
numBlanks = sum(bvec == b);% Show how many bs there are in the node


if numBlanks == 0
    over = true;
%     'The game is over!'
    return;
end


diag1 = [board(1,1) board(2,2) board(3,3)];
diag2 = [board(1,3) board(2,2) board(3,1)];
for i = 1:3
    if allSame(board(i,:)) || allSame(board(:,i)')... % if the rows or columns are all Xs or all Os
        || allSame(diag1) || allSame(diag2) % if either diagonal is all Xs or all Os
        over = true; % The game is over.
%         'The game is over!'
        return;
    end
end
end