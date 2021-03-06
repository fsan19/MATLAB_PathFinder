% About 

% This is an implementation of Djakstras Algorithm to find the best path between 2 points on the map 
% Open this project using MATLAB and run the commands below to calculate a path.
% An image is generated using the plot function
% The map can be altered by changing the map_1.txt file




%% Function descriptions
%% To view a map, given a map file

viewmap('map_1.txt',0); % shows walls
%viewmap('map_1.txt',1); % shows track

%or else

[m]=map_convert('map_1.txt');
plotmap(m);
%Here   'm' is the map file returned as a matrix


% Your simulations will use this structure
[m,v,s]=dfs_1('map_1.txt',[2,2],[14,14]);
%Here   'm' is the map file returned as a matrix
%       'v' is a matrix that shows which cells have been visited, '0' means
%       visited, '1' means not visited
%       's' is the vector of steps to reach the target,
%       [startloc] is a vector i.e. [2,2]
%       [targetloc] is also a vector ie [4,18]

%To view the path determined above use 
plotmap(m,s);

