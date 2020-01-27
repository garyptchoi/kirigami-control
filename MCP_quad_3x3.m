% MCP_quad_3x3: Plot an exmaple of minimum connecting link patterns for 3x3
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences USA, 2020.

%% Parameters
M = 3; %width
N = 3; %height
nquad = M*N; %Number of quads
nlink = M*N-1;% theoretical lower bound for number of links

% inner links
linkpairs=[
    4*2, 4*5-3;
    4*5-2, 4*6-3;
    4*5, 4*8-3;
    4*6, 4*9-3;
    4*1-2, 4*2-3;
    4*2-2, 4*3-3;
    4*1, 4*4-3;
    4*4, 4*7-3;
    ];

linkpairsc = ceil(linkpairs/4);
adjacencyMatrix = sparse([linkpairsc(:,1); linkpairsc(:,2)], [linkpairsc(:,2); linkpairsc(:,1)], ones(length(linkpairsc)*2,1), nquad, nquad);
G = graph(adjacencyMatrix);
bins = conncomp(G);
num_component = length(unique(bins));
disp(['# connected components = ',num2str(num_component)])

%% generate plot
v = zeros(4*M*N,2);
f = [];
for i = 0:N-1 
    for j = 0:M-1
        n = M*i + j + 1;
        v(4*n-3,:) = [2*j,2*i];
        v(4*n-2,:) = [2*j+1.3,2*i];
        v(4*n-1,:) = [2*j+1.3,2*i+1.3];
        v(4*n,:) = [2*j,2*i+1.3];
        f = [f; 4*n-3 4*n-2 4*n-1 4*n];
    end
end

% plot the quads
figure; hold on;
for i = 1:length(linkpairs)
    plot(v(linkpairs(i,:),1), v(linkpairs(i,:),2),'Color',[255 51 51]/255,'LineWidth',3);
end
patch('Faces',f,'Vertices',v,'FaceColor',[89 197 255]/255,'EdgeColor','k','LineWidth',3);
axis equal tight off
   

    