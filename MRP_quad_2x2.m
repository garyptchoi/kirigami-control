% MRP_quad_2x2: Construct all minimum rigidifying link patterns for 2x2
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences USA, 2020.

%% 
L = 2;
nquad = L^2; %Number of quads
nlink = ceil((3*L^2-3)/2);% theoretical lower bound for number of links

% inner links
Linkpairs=[
    4*1-1, 4*2;
    4*1-1, 4*3-2;
    4*2, 4*4-3;
    4*3-2, 4*4-3;
    
    4*1-2, 4*2-3;
    4*2-1, 4*4-2;
    4*1, 4*3-3;
    4*3-1, 4*4;
    ];

combinations = combnk(1:length(Linkpairs),nlink);
DoF = zeros(length(combinations),1);
tic;
for k = 1:length(combinations)
    
    mat=zeros(nquad*5+nlink*2,nquad*4*2); 
    linkpairs = Linkpairs(combinations(k,:),:);
    
    %% Edge length constraints
    % 4 quad boundary constraints, and 1 no shear constraints
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
    %%
    rown=nquad*5+1;

    for t=1:size(linkpairs,1)
        [mat,rown]=constrain(mat,rown,linkpairs(t,1),linkpairs(t,2));
    end

    DoF(k) = nquad*8-rank(mat);
%     disp(['rank = ',num2str(rank(mat))])
%     disp(['DoF = ',num2str(nquad*8-rank(mat))])

end
toc;
solution = find(DoF == 3);
disp(['# optimal patterns = ',num2str(length(solution))])

%% generate plot
v = zeros(4*L*L,2);
f = [];
for i = 0:L-1 
    for j = 0:L-1
        n = L*i + j + 1;
        v(4*n-3,:) = [2*j,2*i];
        v(4*n-2,:) = [2*j+1.5,2*i];
        v(4*n-1,:) = [2*j+1.5,2*i+1.5];
        v(4*n,:) = [2*j,2*i+1.5];
        f = [f; 4*n-3 4*n-2 4*n-1 4*n];
    end
end

for k = 1:length(solution)
    %%
    linkpairs = Linkpairs(combinations(solution(k),:),:);
    % plot the quads
    figure;
    patch('Faces',f,'Vertices',v,'FaceColor',[89 197 255]/255,'EdgeColor','k','LineWidth',3);
    axis equal tight off
    hold on
    % plot the links
    for i = 1:size(linkpairs,1)
        plot(v(linkpairs(i,:),1), v(linkpairs(i,:),2),'r-','LineWidth',3);
    end
end

%%
function [mat, rown] = constrain(mat,rown,i,j)
    mat(rown,i*2-1)=1;
    mat(rown,j*2-1)=-1;
    mat(rown+1,i*2)=1;
    mat(rown+1,j*2)=-1;
    rown = rown+2;
end
