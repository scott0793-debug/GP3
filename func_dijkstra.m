function x_let=func_dijkstra(A,b,mu)
%% calculate #nodes and #edges out of the map
[~,num_edges]=size(A); % num_nodes is the number of nodes.num_edges is the number of edges

%% construct the Dijkstra understandable graph
start_node = zeros(num_edges,1);
end_node = zeros(num_edges,1);

for i=1:num_edges
    start_node(i)=find(A(:,i)==1,1); % find the first element whose value is 1
    end_node(i)=find(A(:,i)==-1,1); % find the first element whose value is -1
end
r=find(b==1); % r is origin
s=find(b==-1);% s is the destinatioin
start_node_real=start_node';
end_node_real=end_node';
mu_real=mu';
graph_dijkstra = digraph(start_node_real,end_node_real,mu_real);

%% use dijkstra's algorithm to find the shortest path
 [~,~,x_let] =shortestpath(graph_dijkstra,r,s);




end