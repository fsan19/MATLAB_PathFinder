% This script uses the dfs function to solve a maze, provided a maze (as a
% text file), a start location and a target location.
% The results m and s of the dfs function are then displayed on a figure
% using the plotmap function

% set the co-ordinates of the start and target locations
startloc = [14,1];
targetloc = [3,8];

[m,v,s]=dfs('map_8.txt', startloc, targetloc);
plotmap(m,s);