% MRP_quad_3x5: Construct a minimum rigidifying link pattern for 3x5 quad
% kirigami
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences, 117(9), 4511-4517, 2020.

%% Parameters
M = 3; %width
N = 5; %height
nquad = M*N; %Number of quads
nlink = ceil((3*M*N-3)/2);% theoretical lower bound for number of links
mat=zeros(nquad*5+nlink*2,nquad*4); 

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
for i=1:M-1
    [mat,rown]=constrain(mat,rown,4*i-2,4*(i+1)-3);
end
for i=M*N-M+1:M*N-1
    [mat,rown]=constrain(mat,rown,4*i-1,4*(i+1));
end
for i=1:M:M*N-M
    [mat,rown]=constrain(mat,rown,4*i,4*(i+M)-3);
end
for i=M:M:M*N-M
    [mat,rown]=constrain(mat,rown,4*i-1,4*(i+M)-2);
end

%% inner links
linkpairs=[
    4*2-1, 4*5-2;
    4*4-2, 4*5-3;
    4*5-1, 4*6;
    4*5, 4*8-3;
    4*8-1, 4*9;
    4*8-1, 4*11-2;
    4*10-2, 4*11-3;
    4*11-1, 4*12;
    4*11, 4*14-3;
    ];

for t=1:size(linkpairs,1)
    [mat,rown]=constrain(mat,rown,linkpairs(t,1),linkpairs(t,2));
end

disp(['rank = ',num2str(rank(mat))])
disp(['DoF = ',num2str(nquad*8-rank(mat))])

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
% plot the boundary links
for i=1:M-1
    plot(v([4*i-2,4*(i+1)-3] ,1), v([4*i-2,4*(i+1)-3],2),'Color',[255 51 51]/255,'LineWidth',3);
end
for i=M*N-M+1:M*N-1
    plot(v([4*i-1,4*(i+1)],1), v([4*i-1,4*(i+1)],2),'Color',[255 51 51]/255,'LineWidth',3);
end
for i=1:M:M*N-M
    plot(v([4*i,4*(i+M)-3],1), v([4*i,4*(i+M)-3],2),'Color',[255 51 51]/255,'LineWidth',3);
end
for i=M:M:M*N-M
    plot(v([4*i-1,4*(i+M)-2],1), v([4*i-1,4*(i+M)-2],2),'Color',[255 51 51]/255,'LineWidth',3);
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
