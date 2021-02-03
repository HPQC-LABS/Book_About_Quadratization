%% DC2
[LHS,RHS] = lhs2rhs('zzz',1,'P(3->2)DC2');
[V_RHS, E_RHS] = eig(RHS); E_RHS=diag(E_RHS)';
[V_LHS, E_LHS] = eig(LHS); E_LHS=diag(E_LHS)';

k=3;
E_RHS = round(E_RHS(:,1:2^k)); % rounding just to make 1 and 1.0000 take the same amount of space, but beware that 0.5 will erroneously look the same as 1.
V_RHS = V_RHS(:,1:2^k);        % looking at the eigenvectors with the lowest 8 eigenvalues.

dist=zeros(2^k,2^(k-1));
for ind_R = 1:2^(k-1)          % We go up to 4 because any single 3-body term of Pauli operators will have 4 eigenvalues of -1 and 4 eigenvalues of +1.
  for ind_L = 1:2^k
    dist(ind_L,ind_R) = norm(V_RHS(:,ind_R)-V_LHS(:,ind_L));
  end
end
[row, column]=find(dist<1e-13); % Rows and columns of dist, which have elements smaller than 1e-6 (i.e., rows tells us which columns of V_LHS are the same as the first 4 columns of V_RHS)
sorted_V_LHS=V_LHS(:,row);
sorted_E_LHS=E_LHS(row);
isequal(sorted_V_LHS,V_RHS(:,1:2^(k-1))) % works for zzz, zzx, zzy, etc. but doesn't work for xxx, zxx, etc. Also, seems NOT to work for Delta > 1.
isequal(sorted_E_LHS,E_RHS(:,1:2^(k-1))) % only the first two elements match, if Delta is only 1. 

%% DC1
%% KKR