% MRP_kagome_3x5: Construct a minimum rigidifying link pattern for 7x7
% kagome kirigami
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences USA, 2020.

M = 3; %width
N = 5; %height
ntri = M*N; %Number of triangles
nlink = ceil((3*M*N-3)/2);% theoretical lower bound for number of links


Linkpairs = [];
% bdy links
for n = 0
    for m = 0:M-2
        i = M*n + m + 1;
        if mod(n,2) ==  mod(m,2)
            Linkpairs = [Linkpairs;
            3*i-1, 3*(i+1)-2];
        
        else
            Linkpairs = [Linkpairs;
            3*i-2, 3*(i+1)-2];
        end
    end
end
for n = N-1
    for m = 0:M-2
        i = M*n + m + 1;
        if mod(n,2) ==  mod(m,2)
            Linkpairs = [Linkpairs;
            3*i, 3*(i+1)];
        
        else
            Linkpairs = [Linkpairs;
            3*i-1, 3*(i+1)];
        end
    end
end

for n = 0:N-2
    for m = 0
        i = M*n + m + 1;
        if mod(n,2) ==  mod(m,2)
            Linkpairs = [Linkpairs;
            3*i, 3*(i+M)-2];
        
        else
            Linkpairs = [Linkpairs;
            3*i, 3*(i+M)-2];
        end
    end
    
    for m = M-1
        i = M*n + m + 1;
        if mod(n,2) ==  mod(m,2)
            Linkpairs = [Linkpairs;
            3*i, 3*(i+M)-2];
        
        else
            Linkpairs = [Linkpairs;
            3*i-1, 3*(i+M)-1];
        end
    end
end
num_bdy_links = length(Linkpairs);

% horizontal links
for n = 0:N-1
    for m = 0:M-2
        i = M*n + m + 1;
        if mod(n,2) ==  mod(m,2)
            Linkpairs = [Linkpairs;
            3*i-1, 3*(i+1)-2;
            3*i, 3*(i+1)];
        
        else
            Linkpairs = [Linkpairs;
            3*i-2, 3*(i+1)-2;
            3*i-1, 3*(i+1)];
        end
    end
end

% vertical links
for n = 0:N-2
    for m = 0:M-1
        i = M*n + m + 1;
        if mod(n,2) ==  mod(m,2)
            Linkpairs = [Linkpairs;
            3*i, 3*(i+M)-2];
        
        else
            Linkpairs = [Linkpairs;
            3*i-1, 3*(i+M)-1;
            3*i, 3*(i+M)-2];
        end
    end
end

% remove the duplicated boundary links
Linkpairs = unique(Linkpairs,'rows','stable');

mat=zeros(ntri*3+nlink*2,ntri*3*2);

linkpairs = Linkpairs([1:num_bdy_links, 15 18 20 24 25 29 32 37 22],:); 

% Edge length constraints
    for n = 0:N-1 
        for m = 0:M-1
            i = M*n + m + 1;
            if mod(n,2) ==  mod(m,2)
                mat(i*3-2,i*6-5)=-1;
                mat(i*3-2,i*6-3)=1;

                mat(i*3-1,i*6-5)=-1;
                mat(i*3-1,i*6-1)=1;
                mat(i*3-1,i*6-4)=-sqrt(3);
                mat(i*3-1,i*6-0)=sqrt(3);

                mat(i*3,i*6-0)=sqrt(3);
                mat(i*3,i*6-2)=-sqrt(3);
                mat(i*3,i*6-1)=-1;
                mat(i*3,i*6-3)=1;

            else
                mat(i*3-2,i*6-1)=-1;
                mat(i*3-2,i*6-3)=1;

                mat(i*3-1,i*6-5)=-1;
                mat(i*3-1,i*6-3)=1;
                mat(i*3-1,i*6-4)=-sqrt(3);
                mat(i*3-1,i*6-2)=sqrt(3);

                mat(i*3  ,i*6-1)=-1;
                mat(i*3  ,i*6-5)=1;
                mat(i*3  ,i*6-4)=-sqrt(3);
                mat(i*3  ,i*6-0)=sqrt(3);
            end
        end
    end

rown=ntri*3+1;

% link constraints
for t=1:size(linkpairs,1)
    [mat,rown]=constrain(mat,rown,linkpairs(t,1),linkpairs(t,2));
end

disp(['DoF = ',num2str(ntri*6-rank(mat))]);

%%
% generate plot
v = zeros(3*M*N,2);
f = [];
edgelength = 2.5;
hs = 0.25;
vs = 1;
for i = 0:N-1 
    for j = 0:M-1
        n = M*i + j + 1;
        if mod(i,2) ==  mod(j,2)
        
            v(3*n-2,:) = [2*j-edgelength/2, 2*i+i*vs];
            v(3*n-1,:) = [2*j+edgelength/2, 2*i+i*vs];
            v(3*n,:)   = [2*j , 2*i+edgelength*sqrt(3)/2+i*vs];
        else
            v(3*n-2,:) = [2*j, 2*i+i*vs];
            v(3*n-1,:) = [2*j+edgelength/2, 2*i+edgelength*sqrt(3)/2+i*vs];
            v(3*n,:) = [2*j-edgelength/2, 2*i+edgelength*sqrt(3)/2+i*vs];
        end
        f = [f; 3*n-2 3*n-1 3*n];
    end
end
% plot the tri
figure;
hold on
patch('Faces',f,'Vertices',v,'FaceColor',[89 197 255]/255,'EdgeColor','k','LineWidth',3);
axis equal tight off
   
% plot the links
for i = 1:size(linkpairs,1)
    plot(v(linkpairs(i,:),1), v(linkpairs(i,:),2),'Color',[255 51 51]/255,'LineWidth',3);
end
    
%%
function [mat, rown] = constrain(mat,rown,i,j)
    mat(rown,i*2-1)=1;
    mat(rown,j*2-1)=-1;
    mat(rown+1,i*2)=1;
    mat(rown+1,j*2)=-1;
    rown = rown+2;
end
