% MRP_kagome_all_with_bdy: Construct all minimum rigidifying link
% patterns for quad kirigami, assuming all boundary links are included
% 
% Remark: The computation may take very long time for 5x5 or larger
% 
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences, 117(9), 4511-4517, 2020.

% M = 2; %width 
% N = 2; %height
M = 3; %width 
N = 3; %height
% M = 4; %width 
% N = 4; %height
ntri = M*N; %Number of quads
nlink = ceil((3*M*N-3)/2);% theoretical lower bound for number of links

Linkpairs = [];
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

%%

combinations = combnk(1:length(Linkpairs),nlink);
DoF = zeros(length(combinations),1);

tic;
parfor k = 1:length(combinations)
    
    mat=zeros(ntri*3+nlink*2,ntri*3*2); 
    linkpairs = Linkpairs(combinations(k,:),:);

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

    DoF(k) = ntri*6-rank(mat);
    
end
toc;
solution = find(DoF == 3);
disp(['# MRPs assuming all boundary links = ',num2str(length(solution))])

%%
function [mat, rown] = constrain(mat,rown,i,j)
    mat(rown,i*2-1)=1;
    mat(rown,j*2-1)=-1;
    mat(rown+1,i*2)=1;
    mat(rown+1,j*2)=-1;
    rown = rown+2;
end