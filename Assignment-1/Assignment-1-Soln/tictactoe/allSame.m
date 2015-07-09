function same = allSame(line)
% A line is either all Xs or all Os.


x = 1;
o = -1;
b = 0;
same = all(line == x) || all(line == o);
end