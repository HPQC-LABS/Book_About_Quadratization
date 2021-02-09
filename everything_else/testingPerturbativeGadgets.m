[LHS,RHS] = lhs2rhs('zzz',2e13,'P(3->2)DC2');
tol = 1e-03;
[V_RHS,E_RHS] = eig(RHS);
[V_LHS,E_LHS] = eig(LHS);
E_RHS = diag(E_RHS);
E_LHS = diag(E_LHS);
[ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
R = V_RHS(:,ind_evals_R);
[ind_evecs_L, ind_evecs_R] = find( abs(L'*R) >= tol );                     % indices of L with matiching eigenvectors in R
sorted_L = L(:,ind_evecs_L);
isequal(sorted_L(:,1:32:end),round(abs(V_RHS(:,1:8))))                     % gives 1

%% Test results

% DC1: zzz still failed, max(min(dist)) always equals to 1.4142 even when all 16 energies match
% KKR: failed for vectors comparison as well
% DC2: We have tested zzz, zxx (4e12), zxz (4.5e12) etc., but need loop to find appropriate and more precise Delta value for zzz, zzx etx.
