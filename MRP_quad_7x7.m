% MRP_quad_7x7: Construct a minimum rigidifying link pattern for 7x7 quad
% kirigami
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences, 117(9), 4511-4517, 2020.

%% Parameters
L = 7; 
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
for i=L^2-L+1:L^2-1
    [mat,rown]=constrain(mat,rown,4*i-1,4*(i+1));
end
for i=1:L:L^2-L
    [mat,rown]=constrain(mat,rown,4*i,4*(i+L)-3);
end
for i=L:L:L^2-L
    [mat,rown]=constrain(mat,rown,4*i-1,4*(i+L)-2);
end

%% inner links
linkpairs=[
    4*2, 4*9-3;
    4*4, 4*11-3;
    4*5-1, 4*12-2; 
    4*13-2, 4*14-3;
    4*27-2, 4*28-3;
    4*34-1, 4*35;
    4*41-1, 4*48-2;
    4*39-1, 4*46-2;
    4*38, 4*45-3;
    4*36-1, 4*37;
    4*22-1, 4*23;
    4*15-2, 4*16-3;
    
    4*9-2, 4*10-3;
    4*10-2, 4*11-3;
    4*11-2, 4*12-3;
    4*12-2, 4*13-3;
    
    4*37-1, 4*38;
    4*38-1, 4*39;
    4*39-1, 4*40;
    4*40-1, 4*41;
    
    4*13-1, 4*20-2;
    4*20-1, 4*27-2;
    4*27-1, 4*34-2;
    4*34-1, 4*41-2;
    
    4*9, 4*16-3;
    4*16, 4*23-3;
    4*23, 4*30-3;
    4*30, 4*37-3;
    
    4*10, 4*17-3;
    4*11, 4*18-3;
    4*12-1, 4*19-2;
    4*17-2, 4*18-3;
    4*18-2, 4*19-3;
    4*19-1, 4*20;
    4*17-1, 4*24-2;
    4*18-1, 4*25-2;
    4*19, 4*26-3;
    4*23-2, 4*24-3;
    4*26-1, 4*27;
    4*24-1, 4*31-2;
    4*25, 4*32-3;
    4*26, 4*33-3;
    4*30-2, 4*31-3;
    4*31-1, 4*32;
    4*32-1, 4*33;
    4*31, 4*38-3;
    4*32-1, 4*39-2;
    4*33-1, 4*40-2;
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
