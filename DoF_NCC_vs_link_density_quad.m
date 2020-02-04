% The codes below calculate the Degree of Freedom (DoF) and NCC (number of
% connected components in a quad kirigami structure, with varying number of
% links added.
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences USA, 2020.

%% Parameters
for LL=1:10 %for L from 10 to 100


L=LL*10;
nquad=L^2; %Number of quads

%% Construct all the links
%% Boundary links

linkpairs=[
];
for i=1:L-1
    linkpairs(end+1,:)=[4*i-2,4*(i+1)-3];
end
for i=L^2-L+1:L^2-1
    linkpairs(end+1,:)=[4*i-1,4*(i+1)];
end
for i=1:L:L^2-L
    linkpairs(end+1,:)=[4*i,4*(i+L)-3];
end
for i=L:L:L^2-L
    linkpairs(end+1,:)=[4*i-1,4*(i+L)-2];
end


%% inner links

% horizontal
for jj=1:L-1
    for i=jj*L-(L-1):jj*L-1
        linkpairs(end+1,:)=[4*i-1,4*(i+1)];
    end
end
for jj=1:L-1
    for i=jj*L+1:jj*L+L-1
        linkpairs(end+1,:)=[4*i-2,4*(i+1)-3];
    end
end
% Vertical
for jj=1:L-1
    for i=jj:L:nquad-L
        linkpairs(end+1,:)=[4*i-1,4*(i+L)-2];
    end
end
for jj=1:L-1
    for i=jj+1:L:nquad-L
        linkpairs(end+1,:)=[4*i,4*(i+L)-3];
    end
end
%% Add Edge Length constraint

n_simu=100;
n_bin=60;
n_maxlink=4*L*(L-1);

if n_maxlink ~= size(linkpairs,1)
    error('Link Number Incorrect!')
end

% Link constraints
mat=zeros(nquad*12,3);

% Edge length constraints
% 4 quad boundary constraints, and 1 no shear constraints (Direction fixed
% for now, from bottom left to top right)
for i=1:nquad
    mat(i*12-11,:)=[i*5-4,i*8-7,-1];
    mat(i*12-10,:)=[i*5-4,i*8-5,1];
    mat(i*12-9,:)=[i*5-3,i*8-7,1];
    mat(i*12-8,:)=[i*5-3,i*8-3,-1];
    mat(i*12-7,:)=[i*5-3,i*8-6,1];
    mat(i*12-6,:)=[i*5-3,i*8-2,-1];
    mat(i*12-5,:)=[i*5-2,i*8-6,1];
    mat(i*12-4,:)=[i*5-2,i*8-0,-1];
    mat(i*12-3,:)=[i*5-1,i*8-4,1];
    mat(i*12-2,:)=[i*5-1,i*8-2,-1];
    mat(i*12-1,:)=[i*5-0,i*8-3,1];
    mat(i*12-0,:)=[i*5-0,i*8-1,-1];
end


link_list = round(linspace(0*n_maxlink,1*n_maxlink,n_bin+1));

dof_all=zeros(length(link_list),n_simu);
dof_rot=zeros(length(link_list),n_simu);
num_conncomp=zeros(length(link_list),n_simu);
size_conncomp=zeros(length(link_list),n_simu);

for link_i = 1:length(link_list)
    n_link = link_list(link_i);

    parfor jjj=1:n_simu
        disp([num2str(L),' ',num2str(n_link),' ',num2str(jjj)]);
        ids = randsample(1:n_maxlink,n_link);
        
        % Add the link constraint
        newmat=[mat;zeros(4*length(ids),size(mat,2))];
        newrown=0;%nquad*5;
        for t=1:length(ids)
            [newmat,newrown]=constrain(newmat,newrown,linkpairs(ids(t),1),linkpairs(ids(t),2),nquad);
        end
        
        % Calculate the DoF
        rgd_Matrix = sparse(newmat(:,1),newmat(:,2), ...
            newmat(:,3), nquad*5+n_link*2, 8*L^2);
        [r]=calc_rank(rgd_Matrix);
        dof=nquad*8-r;
        dof_all(link_i,jjj)=dof;
        
        % Calculate the NCC and size of the largest connected component
        linkpairs_adj=ceil(linkpairs/4);
        adjacencyMatrix = sparse([linkpairs_adj(ids,1); linkpairs_adj(ids,2)], [linkpairs_adj(ids,2); linkpairs_adj(ids,1)], ...
            ones(size(linkpairs(ids,:),1)*2,1), nquad, nquad);
        G = graph(adjacencyMatrix);
        bins = conncomp(G);
        num_conncomp(link_i,jjj)=max(bins);
        [MM,FF]=mode(bins);
        size_conncomp(link_i,jjj)=FF;
        
        % Calculate the number of internal rotatoinal modes
        dof_rot(link_i,jjj)=dof-3*max(bins);
        
        
    end
    
end

ppsave(L,link_list, dof_all, dof_rot, num_conncomp, size_conncomp, svd_vals);

end


function [mat, rown] = constrain(mat,rown,i,j,nquad)
    % Add constraint to the rigitidy matrix
    mat(12*nquad+rown*2+1,:)=[5*nquad+rown+1,i*2-1,1];
    mat(12*nquad+rown*2+2,:)=[5*nquad+rown+1,j*2-1,-1];
    mat(12*nquad+rown*2+3,:)=[5*nquad+rown+2,i*2,1];
    mat(12*nquad+rown*2+4,:)=[5*nquad+rown+2,j*2,-1];
    rown = rown+2;
end

function ppsave(L,link_list,dof_all,dof_rot,num_conncomp,size_conncomp,svd_vals) 
    % Save data in parallel for
    save(['DoF_and_numbin_rho030to060_cst_finer_',num2str(L),'_',datestr(datetime('now'),'yyyymmdd'),'.mat'],...
        'link_list','dof_all','dof_rot','num_conncomp','size_conncomp','svd_vals');

end