% The codes below calculate the Degree of Freedom (DoF) and NCC (number of
% connected components in a kagome kirigami structure, with varying number 
% of links added
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences USA, 2020.

%%
M = 30; %width
N = 30; %height
ntri = M*N; %Number of quads
nlink = rigidity_bound(M,N);% theoretical lower bound for number of links


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

n_simu=100;
n_bin=20;
n_maxlink=round((3*M*floor((N-1)/2))+(1.5*M-0.5*rem(M,2))*rem(N-1,2)+2*(M-1)*N);%Number of links
% combinations = combnk(1:length(Linkpairs),nlink);


link_list = round(linspace(0,n_maxlink,n_bin+1));
dof_all = zeros(n_bin+1,n_simu);
num_conncomp = zeros(n_bin+1,n_simu);
size_conncomp = zeros(n_bin+1,n_simu);
dof_rot = zeros(n_bin+1,n_simu);

tic;
for link_i = 1:length(link_list)
    
    n_link = link_list(link_i);
    mat=zeros(ntri*10+n_link*4,3);
    
    % Edge length constraints

    for n = 0:N-1 
        for m = 0:M-1
            i = M*n + m + 1;
            if mod(n,2) ==  mod(m,2)
                mat(i*10-9,:)=[i*3-2,i*6-5,-1];
                mat(i*10-8,:)=[i*3-2,i*6-3,1];

                mat(i*10-7,:)=[i*3-1,i*6-5,-1];
                mat(i*10-6,:)=[i*3-1,i*6-1,1];
                mat(i*10-5,:)=[i*3-1,i*6-4,-sqrt(3)];
                mat(i*10-4,:)=[i*3-1,i*6-0,sqrt(3)];

                mat(i*10-3,:)=[i*3,  i*6-0,sqrt(3)];
                mat(i*10-2,:)=[i*3,  i*6-2,-sqrt(3)];
                mat(i*10-1,:)=[i*3,  i*6-1,-1];
                mat(i*10-0,:)=[i*3,  i*6-3,1];

            else
                mat(i*10-9,:)=[i*3-2,i*6-1,-1];
                mat(i*10-8,:)=[i*3-2,i*6-3,1];

                mat(i*10-7,:)=[i*3-1,i*6-5,-1];
                mat(i*10-6,:)=[i*3-1,i*6-3,1];
                mat(i*10-5,:)=[i*3-1,i*6-4,-sqrt(3)];
                mat(i*10-4,:)=[i*3-1,i*6-2,sqrt(3)];

                mat(i*10-3,:)=[i*3  ,i*6-1,-1];
                mat(i*10-2,:)=[i*3  ,i*6-5,1];
                mat(i*10-1,:)=[i*3  ,i*6-4,-sqrt(3)];
                mat(i*10-0,:)=[i*3  ,i*6-0,sqrt(3)];
            end
        end
    end
    for jjj=1:n_simu
        linkpairs = Linkpairs(randsample(size(Linkpairs,1),n_link),:);
        disp([num2str(M),' ',num2str(n_link),' ',num2str(jjj)]);
        rown=0;
        newmat=mat;
        
        % Add link constraints
        for t=1:size(linkpairs,1)
            [newmat,rown]=constrain(newmat,rown,linkpairs(t,1),linkpairs(t,2),ntri);
        end
        
        % Calculate DoF
        [r,rgd_Matrix]=calc_rank(newmat,ntri,n_link,M);%%% if M=N
        dof=ntri*6-r;
        dof_all(link_i,jjj) = dof;
        
        % Calculate the NCC and the size of the largest connected component
        linkpairs_adj=ceil(linkpairs/3);
        adjacencyMatrix = sparse([linkpairs_adj(:,1); linkpairs_adj(:,2)], [linkpairs_adj(:,2); linkpairs_adj(:,1)], ones(size(linkpairs,1)*2,1), ntri, ntri);
        G = graph(adjacencyMatrix);
        bins = conncomp(G);
        num_conncomp(link_i,jjj)=max(bins);
        [MM,FF]=mode(bins);
        size_conncomp(link_i,jjj)=FF;
        
        % Calculate the number of internal rotational DoF
        dof_rot(link_i,jjj)=dof-3*max(bins);
    end

end
toc;

save(['DoF_and_ConnComp_Kagome_rect_L',num2str(M),'_',datestr(datetime('now'),'yyyymmdd'),'.mat'],...
    'dof_all','dof_rot','num_conncomp','size_conncomp');


%%
function [mat, rown] = constrain(mat,rown,i,j,ntri)
    mat(10*ntri+rown*2+1,:)=[3*ntri+rown+1,i*2-1,1];
    mat(10*ntri+rown*2+2,:)=[3*ntri+rown+1,j*2-1,-1];
    mat(10*ntri+rown*2+3,:)=[3*ntri+rown+2,i*2,1];
    mat(10*ntri+rown*2+4,:)=[3*ntri+rown+2,j*2,-1];
    rown = rown+2;
end
%%
function n = rigidity_bound(M,N)
    if nargin == 1
        n = ceil((3*M^2-3)/2);
    else 
        n = ceil((3*M*N-3)/2);
    end
end