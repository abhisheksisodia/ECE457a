function [score, bestChild] = minimaxAlphaBeta(node,depth,player,alpha,beta)
x = 1; 
node = reshape(node, 3, 3);
if depth == 0 || gameover(node) == 1;
    score = valuation(node, player);
    bestChild = node;
%     str = board2str(bestChild)
    return;
end


if depth == 3 && gameover(node) == 1;
% if gameover(node) == 1
    score = valuation(node, player);
    bestChild = node;
    str = board2str(bestChild)
    'The game is over!'
    return;
end


children = generateChildren(node,player);
width = size(children,2);% Each column is a node
   
if player == x
    % Find maximum score of all children.
    alpha = -inf;
    for i = 1:width% Try each node
        childMatrix = reshape(children(:,i),3,3);% reshapes children back into 3x3 matreces
        score = minimaxAlphaBeta(childMatrix,depth-1,-player, alpha, beta);% Returns the score of the best child.
        if (score == alpha)% If tie, choose the child with three in a row if there is one
            diag1 = [childMatrix(1,1) childMatrix(2,2) childMatrix(3,3)];
            diag2 = [childMatrix(1,3) childMatrix(2,2) childMatrix(3,1)];
            for j = 1:3
                if allSame(childMatrix(j,:)) || allSame(childMatrix(:,j)')... % if the rows or columns are all Xs or all Os
                        || allSame(diag1) || allSame(diag2) % if either diagonal is all Xs or all Os
                    alpha = score;
                    bestChild = reshape(children(:,i),3,3); % Returns the best child.
                end
            end
        elseif (score > alpha)
            alpha = score;
            bestChild = reshape(children(:,i),3,3); % Returns the best child.
        end
        if alpha >= beta
            break % beta cut-off
        end
    end
    score = alpha;
    return;
else
    % Player == o
    % Find minimum score of all children.
    beta = inf;
    for i = 1:width
        childMatrix = reshape(children(:,i),3,3);% reshapes children back into 3x3 matrix
        score = minimaxAlphaBeta(childMatrix,depth-1,-player, alpha, beta);% Returns the score of that node
        if (score == beta)% If tie, choose the child with three in a row if there is one
            diag1 = [childMatrix(1,1) childMatrix(2,2) childMatrix(3,3)];
            diag2 = [childMatrix(1,3) childMatrix(2,2) childMatrix(3,1)];
            for j = 1:3
                if allSame(childMatrix(j,:)) || allSame(childMatrix(:,j)')... % if the rows or columns are all Xs or all Os
                        || allSame(diag1) || allSame(diag2) % if either diagonal is all Xs or all Os
                    beta = score;
                    bestChild = reshape(children(:,i),3,3); % Returns the best child.
                end
            end
        elseif (score < beta)
                beta = score;
                bestChild = reshape(children(:,i),3,3); % Returns the best child.
        end
    end
    score = beta;
    return;
end
end