function str = board2str(board)
%UNTITLED3 Produces simple picture of board
% E.g. 
% board = 
% 1 0 1
% 1 1 0  
% -1 0 0
%
% then make str
%
% X-O
% XX
% O--
%
[n,m] = size(board);
str = reshape(arrayfun(@player2char, board(:)), n, m);
end


function c = player2char(player)
symbols = 'O-X';
c = symbols(player+2);
end