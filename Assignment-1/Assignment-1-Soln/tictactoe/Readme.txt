ECE457A: Cooperative and Adaptive Algorithms
Spring 2014
Code developed by: Keyvan Golestan

In this package, the main file is 'tictactoe.m' in which collaboration with the Graphical User Interface, and switching turns between player and computer is performed. 
The good news is that you do not need to be worried about most of the functions defined in this file. 
The only one that is so important is called 'decision' that contains the AI behind computer's brain.
'decision' calculates the next move to be played by computer and returns the new board created after performing computer's action. 
This action is then reflected on the game board by calling the 'picksquare' function. As you can see, 'minimaxAlphaBeta' is called in 'decision' function.

The other main functions in this package are:
1. 'minimaxAlphaBeta.m' that runs the core of alpha beta pruning algorithm
2. 'valuation' that determines the values of leaf nodes
3. 'generateChildren.m' that gets the current board and returns a list of children
4. 'allSame.m' that checks wether a line contains all Xs or all Os
