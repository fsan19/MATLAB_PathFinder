%% This is a shell that you will have to follow strictly. 
% You will use the plotmap() and viewmap() to display the outcome of your algorithm.

% Load sample_data_map_8, three variables will be created in your workspace. These were created as a 
% result of [m,v,s]=dfs('map_8.txt',[14,1],[1,18]);
% The solution can be viewed using 
% plotmap(m,s) 

% write your own function for the DFS algorithm.
function [retmap,retvisited,retsteps] = dfs( mapfile,startlocation,targetlocation)
 %import map
 map =map_convert(mapfile);
 retmap=map;
 %take in locations
 start=startlocation;
 finish=targetlocation;
 %save a copy of the unchanged map
 original_map= map;
 %get the size of the map
 s= size(map);
 rows=s(1);
 cols=s(2);
 %iterate through map and find the nodes on the path
 %each node has an ID start at 5 ( value which is not 0 (path) or 1(wall)
 offset =1;
 ID=offset+1;
 %create a set of nodes to be populated with an ID
 nodes=[];
 locations=[];
 for i=1:1:rows,
     for j=1:1:cols,
         %if location is a node
         if(map(i,j)==0)
             %printing all the locations on the path
             %{
             fprintf("(");
             fprintf('%d',i);
             fprintf(",");
             fprintf('%d',j);
             fprintf(")"+'\n');
             %}
             %put the ID 
             map(i,j)=ID;
             % add the node to a Set
             nodes=[nodes;ID];
             %save location of the set (row number corresponds to node)
             locations=[locations; i j];
             %increment ID
             ID=ID+1;
         end;
     end;
 end;
 
%find node id of start point and end point
 startNode=map(start(1),start(2));
 endNode=map(finish(1),finish(2));
 
%create an array of unvisited nodes
unVisited=nodes;
%the distance of all nodes is initialised as Infinity 
distN= Inf(size(nodes));
%the parent node is null for all nodes at beginning
prevN=zeros(size(nodes));
% matrix for the edges or relations between nodes
relations=[];
%init start node as 0 distance offset is 1
distN(startNode-1)=0; 
%Find all relations
relations=find_relations(nodes,original_map,rows,cols,locations);
%DJAKSTRAS ALGORITHM see pseudo code for explanations
while(~(isempty(unVisited)))
        %the node with the minimum distance
        u=find_minDist_unvisited(unVisited,distN,nodes);
        %remove node from unVisited nodes
        unVisited=setdiff(unVisited,u);
        r=size(relations);
        %iterate through all the relations
        for relation=1:1:r(1),
            %for all edges u,v
            if ( relations(relation,1)==u ),
                %calculate the new distN(v) and compare to old
                new = distN(u-1)+1;

                v=relations(relation,2);
                %relations(2)=v but have to subtract offset 1
                if( new < distN(v-1) ),
                    distN(v-1)=new;
                    prevN(v-1)=u;
                end;   
            end;    
        end;
 end
 disp("djakstras done");
path=[];
currentNode=endNode;
%perform relaxation 
while(~ismember(startNode,path))

%add the node to the path
path=[path;currentNode];
%find the previous node
currentNode=prevN(currentNode-1);    
end    
%flip the path array so that is from start to finish
retpath=flip(path);
psize =size(retpath);
pathsize=psize(1);
%store x y coords of the nodes in the path
retsteps=zeros(pathsize,2);
for i=1:pathsize,
    
    %find the path by nodeID
    for j=1:rows,
        for k=1:cols,
            if ( map(j,k)==retpath(i) ),
                %format of retstep is row,col 
                %matrix of 2 cols for co-ords -rows correspond to each node  
                retsteps(i,2)=k;
                retsteps(i,1)=j;
            end;    
        end;
    end;    

end;
disp("done");

%all nodes have been visited 
retvisited=nodes;


end

function index = findNode_index(Node,nodes)
    
    index=NaN;
    for i=1:1:size(nodes),

        if (nodes(i)==Node),
            index=i;
            return 
        end;

    end;

    return

end

 %find an node and show its location on the map
 function [x,y] = findNode_xy(Node,nodes,locations)
     %Node=10;
     for i=1:1:size(nodes),
         
         if (nodes(i)==Node),
             disp(i);
             %index value i corresponds to row of locations matrix
             x=locations(i,1);
             y=locations(i,2);
         end;
         
     end;
 
 end
 
 %finds all the set relations for the nodes
 function relations=find_relations(nodes,original_map,rows,cols,locations)
    %for all the nodes:
    %1.Check for all the squares around 
    rel=0;
    relations=[];
    a=size(nodes);
    
    for i= 1:1:rows,
        for j=1:1:cols,
            if original_map(i,j) == 0
                currentloc(1) = i;
                currentloc(2) = j;
                if currentloc(1) - 1 ~= 0 && original_map(currentloc(1)-1,currentloc(2)) == 0
                    relations = [relations;[(Node(currentloc(1),currentloc(2),locations)),Node(currentloc(1)-1,currentloc(2),locations)]];
                end
                if currentloc(2) + 1 ~= 20 && original_map(currentloc(1),currentloc(2)+1) == 0
                    relations = [relations;[(Node(currentloc(1),currentloc(2),locations)),Node(currentloc(1),currentloc(2)+1,locations)]];
                end
                if currentloc(1) + 1 ~= 16 && original_map(currentloc(1)+1,currentloc(2)) == 0
                    relations = [relations;[(Node(currentloc(1),currentloc(2),locations)),Node(currentloc(1)+1,currentloc(2),locations)]];
                end
                if currentloc(2) - 1  ~= 0 && original_map(currentloc(1),currentloc(2)-1) == 0
                    relations = [relations;[(Node(currentloc(1),currentloc(2),locations)),Node(currentloc(1),currentloc(2)-1,locations)]];
                end
            end;
            
        end;
    end;    
    
 end
 
 %uses locations index to find node ID +1 as offset by 1
 function node = Node(x,y,locations)
    
    node = 0;
    s =size(locations);
    s=s(1);
    %loop through all locations and find if this point is in the path
    for i=1:1:s,
        if( (locations(i,1)==x) && (locations(i,2)==y)),
            node=i+1;%+1 is offset as we start at id 2
            %disp("found node index");
            return 
        end;    
    end;    
 
 
 end
 
 
 %returns the node with the min dist
 function node_min=find_minDist_unvisited(unVisited,distN,nodes)
    s=size(nodes);
    %store the min distance
    min_distance=Inf;
    %return 0 by default we know this is wrong
    node_min=0;
    %loop through the nodes 
    for node=1:s(1),
        %check if node is unvisited
        if ismember(nodes(node),unVisited),
            %find the distance of this node and compare
            if(distN(node) <= min_distance),
                %update the new distance and min_node
                min_distance=distN(node);
                node_min=nodes(node);
            end;
        end;   
    end;  
 end

function placestep(position,i)
% This function will plot a insert yellow rectangle and also print a number in this rectangle. Use with plotmap/viewmap. 
position = [16-position(1) position(2)];
position=[position(2)+0.1 position(1)+0.1];
rectangle('Position',[position,0.8,0.8],'FaceColor','y');
c=sprintf('%d',i);
text(position(1)+0.2,position(2)+0.2,c,'FontSize',10);
end
