%% This is a shell that you will have to follow strictly.
% You will use the plotmap() and viewmap() to display the outcome of your algorithm.

% Load sample_data_map_8, three variables will be created in your workspace. These were created as a
% result of [m,v,s]=dfs('map_8.txt',[14,1],[1,18]);
% The solution can be viewed using
% plotmap(m,s)

% write your own function for the DFS algorithm.
function [retmap,retvisited,retsteps] = dfs_1( mapfile,startlocation,targetlocation)
%#codegen
retmap = map_convert(mapfile);
nodeID = retmap;

currentloc = startlocation;

% relations = zeros(200,200);
relations = [];


retvisited = [];
retsteps = [];


startnode = 0;
targetnode = 0;

% Initialisation for Dijkstra's algorithm.
% Go through the maze and give each valid co-ordinate location a unique ID
% value, starting from 2 (as all walls have a value of 1). The ID values of
% the startlocation and targetlocation are given to startnode and
% targetnode respectively
count = 2;
for i=1:15
    for j=1:19
        if retmap(i,j) == 0
            nodeID(i,j) = count;
            if startlocation(1) == i && startlocation(2) == j
                startnode = nodeID(i,j);
            elseif targetlocation(1) == i && targetlocation(2) == j
                targetnode = nodeID(i,j);
            end
            count = count + 1;
        end
    end
end

% Exit the function if the startnode or targetnode have an ID of 0. This
% means that either the startlocation and targetlocation are either the
% same(only 1 unique ID value) or that either the startlocation or
% targetlocation are invalid(not in the maze, not in a free space)
if startnode == 0 || targetnode == 0
    disp("invalid start/target location or start and target nodes are the same");
    return;
end

% Initialise the set, unvisited, retvisited and prevnode matrices. The set
% matrix is all the nodes in the maze by their unique ID values. unvisited
% is initially set to be the same as set since initially, no nodes have
% been visited yet. retvisited is also set to set since by the end of this
% function, everycell should have been visited and evaluated. The prevnode
% matrix is the prevnode taken to reach the current node being evaluated.
set = 1:1:count-1;
unvisited = set;
retvisited = set;
prevnode = ones(size(unvisited));

% Initialise the distance values associated with each node to Inf, then set
% the distance value of the startnode to 0
distance = Inf(size(unvisited));
distance(startnode) = 0;

% This block is responsible for creating the relations between each valid
% cell on the map. This is evaluated by checking all the valid directions
% that a cell can move (up, right, down, left) and inserting that relation
% into the relaticons matrix
for i=1:15
    for j=1:19
        if retmap(i,j) == 0
            currentloc(1) = i;
            currentloc(2) = j;
            if currentloc(1) - 1 ~= 0 && retmap(currentloc(1)-1,currentloc(2)) == 0
                relations = [relations;[(nodeID(currentloc(1),currentloc(2))),nodeID(currentloc(1)-1,currentloc(2))]];
            end
            if currentloc(2) + 1 ~= 20 && retmap(currentloc(1),currentloc(2)+1) == 0
                relations = [relations;[(nodeID(currentloc(1),currentloc(2))),nodeID(currentloc(1),currentloc(2)+1)]];
            end
            if currentloc(1) + 1 ~= 16 && retmap(currentloc(1)+1,currentloc(2)) == 0
                relations = [relations;[(nodeID(currentloc(1),currentloc(2))),nodeID(currentloc(1)+1,currentloc(2))]];
            end
            if currentloc(2) - 1  ~= 0 && retmap(currentloc(1),currentloc(2)-1) == 0
                relations = [relations;[(nodeID(currentloc(1),currentloc(2))),nodeID(currentloc(1),currentloc(2)-1)]];
            end
        end
    end
end

% This block is responsible for the implementation of Dijkstra's algorithm.
% This is done by evaluating the cost to reach each node from the starting
% node, and doing this until every node has been visited and every cost has
% been evaluated. If the cost found to reach a node is less than the
% current cost to reach that node, then that cost is updated to the
% lesser value.
while ~isempty(unvisited)
    currentnode = nextnode(set, unvisited, distance);
    unvisited = unvisited(unvisited~=currentnode);
    [relationSize, ~] = size(relations);
    for i=1:relationSize
        relanode = relations(i,:);
        if currentnode == relanode(1)
            if distance(currentnode) ~= Inf
                newdist = distance(currentnode) + 1;
            else
                newdist = Inf;
            end
            j = relanode(2);
            if newdist < distance(j)
                distance(j) = newdist;
                prevnode(j) = currentnode;
            end
        end
    end
end



% This block is used to create the lowest cost path to reach the targetnode
% from the startnode, P. This is done by adding the prevnode associated
% with the targetnode until the startnode is reached, and reversing this
% order.
shortestdist = distance(targetnode);
makepath = targetnode;
P = zeros(1,shortestdist + 1);
P(shortestdist + 1) = targetnode;
for j=1:shortestdist
    for i=1:length(set)
        if set(i) == makepath
            P(shortestdist - j + 1) = prevnode(i);
            makepath = prevnode(i);
            break
        end
    end
end

% This block creates the co-ordinate steps taken from the startnode to the
% targetnode by converting the nodeID values stored in P and converting
% them back to their cell values.
retsteps = zeros(1,2,length(P));
for k=1:length(P)
    for i=1:15
        for j=1:19
            if nodeID(i,j) == P(k)
                retsteps(k,1) = i;
                retsteps(k,2) = j;
            end
        end
    end
end

end

% This helper function is used to evaluate the nextnode to be checked in
% Dijkstra's algorithm. This is done by checking which node has not been
% visited that has the lowest cost value associated with it and returning
% that index value.
function [index] = nextnode(set, unvisited, distance)
index = 1;
for i=2:length(set)
    if(~(ismember(set(index), unvisited)))
        index = i;
    end
    if distance(i) <= distance(index)
        if (ismember(set(i), unvisited))
            index = i;
        end
    end
end

end

% Function provided with code, don't thin it actually gets used
function placestep(position,i)
% This function will plot a insert yellow rectangle and also print a number in this rectangle. Use with plotmap/viewmap.
position = [16-position(1) position(2)];
position=[position(2)+0.1 position(1)+0.1];
rectangle('Position',[position,0.8,0.8],'FaceColor','y');
c=sprintf('%d',i);
text(position(1)+0.2,position(2)+0.2,c,'FontSize',10);
end