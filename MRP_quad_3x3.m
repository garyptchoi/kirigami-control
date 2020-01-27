% MRP_quad_3x3: Construct a minimum rigidifying link pattern for 3x3
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences USA, 2020.

%% Parameters
L = 3;
nquad = L^2; %Number of quads
nlink = ceil((3*L^2-3)/2);% theoretical lower bound for number of links
mat = zeros(nquad*5+nlink*2,nquad*4); 

%% Edge length constraints
% 4 quad boundary constraints, and 1 no shear constraints (Direction fixed
% for now, from bottom left to top right)
for i=1:nquad
    mat(i*5-4,i*8-7)=-1;
    mat(i*5-4,i*8-5)=1;
    mat(i*5-3,i*8-7)=1;
    mat(i*5-3,i*8-3)=-1;
    mat(i*5-3,i*8-6)=1;
    mat(i*5-3,i*8-2)=-1;
    mat(i*5-2,i*8-6)=1;
    mat(i*5-2,i*8-0)=-1;
    mat(i*5-1,i*8-4)=1;
    mat(i*5-1,i*8-2)=-1;
    mat(i*5-0,i*8-3)=1;
    mat(i*5-0,i*8-1)=-1;
end

%% Boundary links
rown=nquad*5+1;
for i=1:L-1
    [mat,rown]=constrain(mat,rown,4*i-2,4*(i+1)-3);
end
for i=L*L-L+1:L*L-1
    [mat,rown]=constrain(mat,rown,4*i-1,4*(i+1));
end
for i=1:L:L*L-L
    [mat,rown]=constrain(mat,rown,4*i,4*(i+L)-3);
end
for i=L:L:L*L-L
    [mat,rown]=constrain(mat,rown,4*i-1,4*(i+L)-2);
end

%% inner links
linkpairs=[
    4*2, 4*5-3;
    4*5-2, 4*6-3;
    4*5-1, 4*8-2;
    4*4-1, 4*5;
    ];

for t=1:size(linkpairs,1)
    [mat,rown]=constrain(mat,rown,linkpairs(t,1),linkpairs(t,2));
end

disp(['rank = ',num2str(rank(mat))])
disp(['DoF = ',num2str(nquad*8-rank(mat))])

%% generate plot
v = zeros(4*L^2,2);
f = [];
for i = 0:L-1 
    for j = 0:L-1
        n = L*i + j + 1;
        v(4*n-3,:) = [2*j,2*i];
        v(4*n-2,:) = [2*j+1.3,2*i];
        v(4*n-1,:) = [2*j+1.3,2*i+1.3];
        v(4*n,:) = [2*j,2*i+1.3];
        f = [f; 4*n-3 4*n-2 4*n-1 4*n];
    end
end

% plot the quads
figure; hold on;
% plot the boundary links
for i=1:L-1
    plot(v([4*i-2,4*(i+1)-3] ,1), v([4*i-2,4*(i+1)-3],2),'Color',[255 51 51]/255,'LineWidth',3);
end
for i=L^2-L+1:L^2-1
    plot(v([4*i-1,4*(i+1)],1), v([4*i-1,4*(i+1)],2),'Color',[255 51 51]/255,'LineWidth',3);
end
for i=1:L:L^2-L
    plot(v([4*i,4*(i+L)-3],1), v([4*i,4*(i+L)-3],2),'Color',[255 51 51]/255,'LineWidth',3);
end
for i=L:L:L^2-L
    plot(v([4*i-1,4*(i+L)-2],1), v([4*i-1,4*(i+L)-2],2),'Color',[255 51 51]/255,'LineWidth',3);
end
% plot the internal links
for i = 1:length(linkpairs)
    plot(v(linkpairs(i,:),1), v(linkpairs(i,:),2),'Color',[255 51 51]/255,'LineWidth',3);
end
patch('Faces',f,'Vertices',v,'FaceColor',[89 197 255]/255,'EdgeColor','k','Linewidth',3);
axis equal tight off

%%
function [mat, rown] = constrain(mat,rown,i,j)
    mat(rown,i*2-1)=1;
    mat(rown,j*2-1)=-1;
    mat(rown+1,i*2)=1;
    mat(rown+1,j*2)=-1;
    rown = rown+2;
end
