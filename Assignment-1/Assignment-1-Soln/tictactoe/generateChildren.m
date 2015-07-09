function children = generateChildren(node,player)
% Find blanks, put x in first blank, store that, then put x in next blank and store that
x = 1;
o = -1;
b = 0;


bvec = node(:);% Bvec is the node written as a vector
locations = find(bvec == b);% shows the locations of b
numBlanks = sum(bvec == b);% Show how many bs there are in the node
% children = zeros (1,numBlanks);% Shows how many children there are.
children = zeros(9,numBlanks);% Columns are the number of new nodes; rows are the number of blanks


for i = 1:length(locations)
    % for j = 1:length(locations)
    newvec = bvec;% Defines newvec as bvec
    newvec(locations(i)) = player;% puts an x or o in the place of the blank
    % newboard = reshape(newvec,3,3);% reshapes newvec back into 3x3 matrix
    children(:,i) = newvec;% Each column is a child
end
end