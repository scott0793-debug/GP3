function path_dijkstraPP=func_dijkstraPP(A,b,mu,edge_index)

[~,num_edges] = size(A);
%% 将地图改为每条路径的起点和终点的集合
for i=1:num_edges
    start_node(i)=find(A(:,i)==1);
    end_node(i)=find(A(:,i)==-1);
end
%% 找到path的起点和终点
s=find(b==1);
r=find(b==-1);

%% 找到X值最小的那条路径的起点和终点
edge_start=find(A(:,edge_index)==1);
edge_end=find(A(:,edge_index)==-1);
g = digraph(start_node,end_node,mu');
path_dijkstraPP=zeros(num_edges,1);
%% 如果路径的起点与路的起点相同
if edge_start==s && edge_end==r
    path_dijkstraPP=edge_index;
elseif edge_start==s && edge_end~=r
    [path_point,distance,edgepath] = shortestpath(g,edge_end,r);
    path_dijkstraPP=[edge_index,edgepath];
%% 如果路径的终点与路的终点相同
elseif edge_end==r && edge_start~=s
    [path_point,distance,edgepath] = shortestpath(g,s,edge_start);
    path_dijkstraPP=[edgepath,edge_index];
else
    [path_point1,distance1,edgepath1] = shortestpath(g,s,edge_start);
    [path_point2,distance2,edgepath2] = shortestpath(g,edge_end,r);
    path_dijkstraPP=[edgepath1,edge_index,edgepath2];
end

end