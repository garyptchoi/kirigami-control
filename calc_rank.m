% Calculate rank of a sparse non-square matrix
%
% Reference:
% S. Chen, G. P. T. Choi, L. Mahadevan, 
% ``Deterministic and stochastic control of kirigami topology.''
% Proceedings of the National Academy of Sciences USA, 2020.


function [r] = calc_rank(rgd_Matrix)

    [m,n] = size(rgd_Matrix);
    q = colamd(rgd_Matrix, [n m]);
    R=qr(rgd_Matrix(:, q), 0);

    tol = max(size(rgd_Matrix))*eps*abs(R(1,1)); %Adaptive tolerence for nonzero element
    % Count nonzero diagonal entries

    r = 0;
    j = 0;
    while 1
        if r>=m || j>=n
            break
        else
            if abs(R(r+1,j+1)) >= tol % If the element is larger
                r = r+1;
                j = j+1;
            else % If the element is smaller than the tolerence
                s = 0;
                while 1 % Check the next element to the right
                    j=j+1;
                    if j>=n
                        break
                    end
                    if abs(R(r+1, j+1))>=tol
                        r=r+1;
                        j=j+1;
                        s = 1;
                        break
                    end
                end
                if s == 0
                    break
                end
            end
        end
    end
end

