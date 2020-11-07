function [P,L,U] = LUfactor(A)




P = eye(size(A,1));
L = eye(size(A,1));
U = A;
rowsize = size(A,1);
columnsize = size(A,2);

if rowsize > columnsize
   stopValOffset = 0;
elseif rowsize == columnsize
    stopValOffset = 1;
else
    stopValOffset = columnsize - rowsize + 1;
end


P_i_storage = zeros(size(A,1),size(A,1),size(A,2)-stopValOffset);
inv_L_i_storage = P_i_storage;

for k = 1:size(A,2)-stopValOffset % looping through columns of A, except last column
    % Permute first
    [~,I] = max(abs(U(k:end,k)));
    I = I+k-1;
    P_i = eye(size(A,1));
    P_i(k,k) = 0;
    P_i(I,I) = 0;
    P_i(k,I) = 1;
    P_i(I,k) = 1;  
    U = P_i*U;
    P = P_i*P;
    P_i_storage(:,:,k) = P_i;
    
    

    
    % Then row reduce one column
    L_i = eye(size(A,1));
    for k2 = k+1:size(L_i,1)
        L_i(k2,k) = -U(k2,k)/U(k,k);
    end
    U = L_i*U;
    
%     if k < size(A,2)-stopValOffset
%         % P_i_lookahead calc
%         k_bar = k+1;
%         [~,I] = max(abs(U(k_bar:end,k_bar)));
%         I = I+k_bar-1;
%         P_i_lookahead = eye(size(A,1));
%         P_i_lookahead(k_bar,k_bar) = 0;
%         P_i_lookahead(I,I) = 0;
%         P_i_lookahead(k_bar,I) = 1;
%         P_i_lookahead(I,k_bar) = 1;
%     end
    
    
    % Current Inverse
    inv_L_i = 2*eye(size(A,1))-L_i;
    inv_L_i_storage(:,:,k) = inv_L_i;
    
    
   
    
%     if k < size(A,2)-stopValOffset
%         inv_L_i = P_i_lookahead*inv_L_i*P_i_lookahead;
%     end

    
    
    
end


for k = 1:size(A,2)-stopValOffset
    L_inv_temp = inv_L_i_storage(:,:,k);
    
    if k ~= size(A,2)-stopValOffset
        for j = k+1:size(A,2)-stopValOffset
            L_inv_temp = P_i_storage(:,:,j)*L_inv_temp*P_i_storage(:,:,j);
        end
    end
    
    
    L=L*L_inv_temp;
    
end



end