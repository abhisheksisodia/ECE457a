x = 1;
o = -1;
b = 0;

currBoard = [
    b b b
    b b b
    b b b];


depth = 4;

turncount = 0;


% player = o;
player = x;


% While the game is still going, get minimax score and make next move
while gameover(currBoard) == 0


%     turncount
%     player
    if player == x
        [score, currBoard] = minimaxAlphaBeta(currBoard, depth, player, -inf, inf);
    else
        % Show board
        %  Display 1 as X, -1 as O, 0 as -
        str = board2str(currBoard)
        i = input('Enter the row: ');
%         while ~ismember(i,1:3) || ~ismember(j,1:3); only allow legal
%         moves
        j = input('Enter the column: ');
        currBoard(i,j) = o; % How to put a string in a matrix?
    end
   turncount = turncount + 1;
   player = -player;
end