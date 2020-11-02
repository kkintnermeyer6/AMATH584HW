function [Q,R] = factorToQR(A)
% Test code to debug ----
% Q = A;
% R = A;
% -----------------------

A_NumRows = size(A,1);
A_NumColumns = size(A,2);
R = zeros(A_NumColumns,A_NumColumns);
Q = A;



for j = 1:A_NumColumns
    for k = 1:j-1
        if j == 1
            break;
        end
       R(k,j) = (Q(:,k)'*Q(:,j));
       Q(:,j) =  Q(:,j) - (Q(:,k)'*Q(:,j))*Q(:,k);
    end
    R(j,j) = norm(Q(:,j),2);
    Q(:,j) = Q(:,j)/norm(Q(:,j),2);
end


end