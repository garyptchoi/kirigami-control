% MRP_quad_7x7: Construct a minimum rigidifying link pattern for 30x30 quad
% kirigami by hierarchical construction
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences, 117(9), 4511-4517, 2020.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Prepare the basic link patterns for 2x2, 3x3, 5x5 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MRP for 2x2
L = 2;
linkpairs2 = [];
% Boundary links
for i=1:L-1
    linkpairs2=[linkpairs2;4*i-2,4*(i+1)-3];
end
for i=L^2-L+1:L^2-1
    linkpairs2=[linkpairs2;4*i-1,4*(i+1)];
end
for i=1:L:L^2-L
    linkpairs2=[linkpairs2;4*i,4*(i+L)-3];
end
for i=L:L:L^2-L
    linkpairs2=[linkpairs2;4*i-1,4*(i+L)-2];
end
% inner links
linkpairs2=[linkpairs2;
    4*1-1, 4*2;
    ];

if length(linkpairs2) ~= ceil((3*L^2-3)/2)
    error('Incorrect number of link pairs!');
end

% vertex coordinates for the 3x3 quads
v2 = zeros(4*L^2,2);
f2 = [];
for i = 0:L-1 
    for j = 0:L-1
        n = L*i + j + 1;
        v2(4*n-3,:) = [2*j,2*i];
        v2(4*n-2,:) = [2*j+1,2*i];
        v2(4*n-1,:) = [2*j+1,2*i+1];
        v2(4*n,:) = [2*j,2*i+1];
        f2 = [f2; 4*n-3 4*n-2 4*n-1 4*n];
    end
end

%% MRP for 3x3
L = 3;
linkpairs3 = [];
% Boundary links
for i=1:L-1
    linkpairs3=[linkpairs3;4*i-2,4*(i+1)-3];
end
for i=L^2-L+1:L^2-1
    linkpairs3=[linkpairs3;4*i-1,4*(i+1)];
end
for i=1:L:L^2-L
    linkpairs3=[linkpairs3;4*i,4*(i+L)-3];
end
for i=L:L:L^2-L
    linkpairs3=[linkpairs3;4*i-1,4*(i+L)-2];
end
% inner links
linkpairs3=[linkpairs3;
    4*2, 4*5-3;
    4*5-2, 4*6-3;
    4*5-1, 4*8-2;
    4*4-1, 4*5;
    ];

if length(linkpairs3) ~= ceil((3*L^2-3)/2)
    error('Incorrect number of link pairs!');
end

% vertex coordinates for the 3x3 quads
v3 = zeros(4*L^2,2);
f3 = [];
for i = 0:L-1 
    for j = 0:L-1
        n = L*i + j + 1;
        v3(4*n-3,:) = [2*j,2*i];
        v3(4*n-2,:) = [2*j+1,2*i];
        v3(4*n-1,:) = [2*j+1,2*i+1];
        v3(4*n,:) = [2*j,2*i+1];
        f3 = [f3; 4*n-3 4*n-2 4*n-1 4*n];
    end
end

%% MRP for 5x5
L = 5;
linkpairs5 = [];
% Boundary links
for i=1:L-1
    linkpairs5=[linkpairs5;4*i-2,4*(i+1)-3];
end
for i=L^2-L+1:L^2-1
    linkpairs5=[linkpairs5;4*i-1,4*(i+1)];
end
for i=1:L:L^2-L
    linkpairs5=[linkpairs5;4*i,4*(i+L)-3];
end
for i=L:L:L^2-L
    linkpairs5=[linkpairs5;4*i-1,4*(i+L)-2];
end
% inner links
linkpairs5=[linkpairs5;
    4*11-2,4*12-3;
    4*14-1,4*15;
    4*3,4*8-3;
    4*8-1,4*13-2;
    4*13,4*18-3;
    4*18-1,4*23-2;
    
    4*2, 4*7-3;
    4*7-2,4*8-3;
    4*7-1,4*12-2;
    
    4*12-1,4*17-2;
    4*16-2,4*17-3;
    4*17-1,4*18;
    4*17,4*22-3;
    
    4*4-1,4*9-2;
    4*8-2,4*9-3;
    4*9-1,4*10;
    4*9,4*14-3;
    
    4*14, 4*19-3;
    4*18-1,4*19;
    4*19-1,4*24-2;
    ];

if length(linkpairs5) ~= ceil((3*L^2-3)/2)
    error('Incorrect number of link pairs!');
end

% vertex coordinates for the 5x5 quads
v5 = zeros(4*L^2,2);
f5 = [];
for i = 0:L-1 
    for j = 0:L-1
        n = L*i + j + 1;
        v5(4*n-3,:) = [2*j,2*i];
        v5(4*n-2,:) = [2*j+1,2*i];
        v5(4*n-1,:) = [2*j+1,2*i+1];
        v5(4*n,:) = [2*j,2*i+1];
        f5 = [f5; 4*n-3 4*n-2 4*n-1 4*n];
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Construct the link patterns for 30x30 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 30; 
nquad = L^2; % Number of quads
nlink = ceil((3*L^2-3)/2);% theoretical lower bound for number of links

v = zeros(4*L^2,2);
f = [];
for i = 0:L-1 
    for j = 0:L-1
        n = L*i + j + 1;
        v(4*n-3,:) = [2*j,2*i];
        v(4*n-2,:) = [2*j+1,2*i];
        v(4*n-1,:) = [2*j+1,2*i+1];
        v(4*n,:) = [2*j,2*i+1];
        f = [f; 4*n-3 4*n-2 4*n-1 4*n];
    end
end

% Construct link pattern from the small patterns
linkpairs=[];
% fill in the 3x3 links
width3 = range(v3(:,1)) + 1;
for i = 0:(L/3-1)
    for j = 0:(L/3-1)
        % shift v3 
        v3_shifted = [v3(:,1) + i*width3, v3(:,2) + j*width3];
        id3_list = zeros(length(v3),1);
        % find corresponding vertices in the 60x60 pattern
        for k = 1:length(v3)
            id = find(v(:,1) == v3_shifted(k,1) & v(:,2) == v3_shifted(k,2));
            id3_list(k) = id;
        end
        linkpairs=[linkpairs;id3_list(linkpairs3)];
    end
end
       
% fill in the 5x5 links
v5_magnified = zeros(size(v5));
for i = 1:length(v5)
    if mod(v5(i,1),2) == 0
        v5_magnified(i,1) = width3*(v5(i,1)/2);
    else
        v5_magnified(i,1) = width3*((v5(i,1)+1)/2)-1;
    end
    if mod(v5(i,2),2) == 0
        v5_magnified(i,2) = width3*(v5(i,2)/2);
    else
        v5_magnified(i,2) = width3*((v5(i,2)+1)/2)-1;
    end
end
width5 = range(v5_magnified(:,1)) + 1;
for i = 0:(L/3/5-1)
    for j = 0:(L/3/5-1)
        % shift v5 with magnification 
        v5_shifted = [v5_magnified(:,1) + i*width5, v5_magnified(:,2) + j*width5];
        id5_list = zeros(length(v5),1);
        % find corresponding vertices in the 60x60 pattern
        for k = 1:length(v5)
            id = find(v(:,1) == v5_shifted(k,1) & v(:,2) == v5_shifted(k,2));
            id5_list(k) = id;
        end
        linkpairs=[linkpairs;id5_list(linkpairs5)];
    end
end

% fill in the 2x2 links
v2_magnified = zeros(size(v2));
for i = 1:length(v2)
    if mod(v2(i,1),2) == 0
        v2_magnified(i,1) = width5*(v2(i,1)/2);
    else
        v2_magnified(i,1) = width5*((v2(i,1)+1)/2)-1;
    end
    if mod(v2(i,2),2) == 0
        v2_magnified(i,2) = width5*(v2(i,2)/2);
    else
        v2_magnified(i,2) = width5*((v2(i,2)+1)/2)-1;
    end
end
id2_list = zeros(length(v2),1);
% find corresponding vertices in the 60x60 pattern
for k = 1:length(v2)
    id = find(v(:,1) == v2_magnified(k,1) & v(:,2) == v2_magnified(k,2));
    id2_list(k) = id;
end
linkpairs=[linkpairs;id2_list(linkpairs2)];

% Check the total number of links
disp(['Theoretical minimum number of links = ', num2str(nlink)])
disp(['Number of links chosen = ', num2str(length(linkpairs))])

%% generate plot
figure;
patch('Faces',f,'Vertices',v,'FaceColor',[89 197 255]/255,'EdgeColor','k');
axis equal tight off

hold on;
for t = 1:length(linkpairs)
    plot([v(linkpairs(t,1),1), v(linkpairs(t,2),1)], [v(linkpairs(t,1),2), v(linkpairs(t,2),2)], 'r-');
end
