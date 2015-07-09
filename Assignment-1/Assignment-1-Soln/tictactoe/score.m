function score = valuation(board,player)
% Determines score of board for the current player
x = 1; %#ok
o = -1;
b = 0; %#ok


% Heuristic Evaluation Function (heuristic = rule of thumb)
% There are 8 lines that can win: 3 rows, 3 columns, and 2 diagonals.
score = 0;


for i = 1:3
    % Rows
    score = score + evaluateLine(board(i,:));
    % Columns
    score = score + evaluateLine(board(:,i)');
end


% Diagonals
diag1 = [board(1,1) board(2,2) board(3,3)];
score = score + evaluateLine(diag1);
diag2 = [board(1,3) board(2,2) board(3,1)];
score = score + evaluateLine(diag2);
end