% MRP_quad_3x3_all_with_bdy: Construct all minimum rigidifying link
% patterns for 3x3
%
% Remark: The computation may take several minutes
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences USA, 2020.

%% 
L = 3;
nquad = L^2; %Number of quads
nlink = ceil((3*L^2-3)/2);% theoretical lower bound for number of links

% inner links
Linkpairs=[
    4*1-1, 4*2;
    4*1-1, 4*4-2;
    4*2, 4*5-3;
    4*2-1, 4*3;
    4*2-1, 4*5-2;
    4*3, 4*6-3;
    4*4-2, 4*5-3;
    4*5-2, 4*6-3;
    4*4-1, 4*7-2;
    4*4-1, 4*5;
    4*5, 4*8-3;
    4*5-1, 4*8-2;
    4*5-1, 4*6;
    4*6, 4*9-3;
    4*7-2, 4*8-3;
    4*8-2, 4*9-3;
    
    % bdy link
    4*1-2, 4*2-3;
    4*2-2, 4*3-3;
    4*1, 4*4-3;
    4*4, 4*7-3;
    4*7-1, 4*8;
    4*8-1, 4*9;
    4*3-1, 4*6-2;
    4*6-1, 4*9-2;
    ];

combinations = combnk(1:length(Linkpairs),nlink);
DoF = zeros(length(combinations),1);
tic;
parfor k = 1:length(combinations)
    
    mat=zeros(nquad*5+nlink*2,nquad*4); 
    linkpairs = Linkpairs(combinations(k,:),:);
    
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
    %%
    rown=nquad*5+1;

    for t=1:size(linkpairs,1)
        [mat,rown]=constrain(mat,rown,linkpairs(t,1),linkpairs(t,2));
    end

    DoF(k) = nquad*8-rank(mat);

end
toc;
solution = find(DoF == 3);
disp(['# MRPs = ',num2str(length(solution))])

%%
function [mat, rown] = constrain(mat,rown,i,j)
    mat(rown,i*2-1)=1;
    mat(rown,j*2-1)=-1;
    mat(rown+1,i*2)=1;
    mat(rown+1,j*2)=-1;
    rown = rown+2;
end
