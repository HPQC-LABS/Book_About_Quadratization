[LHS,RHS] = lhs2rhs('zzz',2e13,'P(3->2)DC2');
tol = 1e-03;
[V_RHS,E_RHS] = eig(RHS);
[V_LHS,E_LHS] = eig(LHS);
E_RHS = diag(E_RHS);
E_LHS = diag(E_LHS);
[ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
R = V_RHS(:,ind_evals_R);
[ind_evecs_L, ind_evecs_R] = find( abs(L'*R)-1 > -tol );                     % indices of L with matiching eigenvectors in R
sorted_L = L(:,ind_evecs_L);
repeat_times = sum(ind_evals_R == ind_evals_R(1,1))*sum(ind_evecs_R == ind_evals_R(1,1));    % times of eigenvectors repeat in sorted_L
sorted_L = sorted_L(:,1:repeat_times:end);
sorted_R = V_RHS(:,unique(ind_evals_R));
max(sqrt( sum( (abs(sorted_L)-abs(sorted_R)).^2 ) )) < tol;   % maximum distance between two matching eigenvectors is within the bound

% Another potential solution: compute the distance of two eigenvectors

dist = sqrt( sum( (abs(L)-abs(R)).^2 ) );
ind_evecs = find(dist < tol);              % indices of L with mathcing eigenvectors in R
sorted_L = L(:,ind_evecs);
isequal(sorted_L,round(abs(V_RHS(:,unique(ind_evals_R)))))  % gives 1
max(dist(ind_evecs))                        % gives the maximum distance of two matching eigenvectors

% Because iseuqal() function cannot return true when two eigenvectors are only 'close' to each other but not the same,
% we can use sorted_L and maximum distance to determine if the result is expected (i.e. look at the number of matching eigenvectors and their distance)


%% Test results

% DC1: zzz still failed, max(min(dist)) always equals to 1.4142 even when all 16 energies match
% KKR: failed for vectors comparison as well
% DC2: We have tested zzz, zzx (3.58e12), zxz (4.5e12) etc., but need loop to find appropriate and more precise Delta value for zzz, zzx etx.

% Successful example: it has 7 eigenvectors matching ( V_RHS(1-7) ), but still returns false due to small errors in V_RHS. We need to find a way to ignore those errors.

[LHS,RHS] = lhs2rhs('zzx',4.5e12,'P(3->2)DC2');
tol = 1e-03;
[V_RHS,E_RHS] = eig(RHS);
[V_LHS,E_LHS] = eig(LHS);
E_RHS = diag(E_RHS);
E_LHS = diag(E_LHS);
[ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
R = V_RHS(:,ind_evals_R);
[ind_evecs_L, ind_evecs_R] = find( abs(L'*R)-1 > -tol );                     % indices of L with matiching eigenvectors in R
sorted_L = L(:,ind_evecs_L);
repeat_times = sum(ind_evals_R == ind_evals_R(ind_evecs_R(1)))*sum(ind_evecs_R == ind_evecs_R(1));    % times of eigenvectors repeat in sorted_L
sorted_L = sorted_L(:,1:repeat_times:end);  sorted_R = V_RHS(:,unique(ind_evals_R));
dist = sqrt( sum( (abs(sorted_L)-abs(sorted_R)).^2 ) );
max(dist) < tol;   % gives 1

% Failed example: it has 5 eigenvalues matching, but only has 1 matching eigenvectors. If there were more than 1 eigenvectors matching, it gives an error.

[LHS,RHS] = lhs2rhs('xxz',4.3e12,'P(3->2)DC2');
tol = 1e-03;
[V_RHS,E_RHS] = eig(RHS);
[V_LHS,E_LHS] = eig(LHS);
E_RHS = diag(E_RHS);
E_LHS = diag(E_LHS);
[ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
R = V_RHS(:,ind_evals_R);
[ind_evecs_L, ind_evecs_R] = find( abs(L'*R)-1 > -tol );                     % indices of L with matiching eigenvectors in R
sorted_L = L(:,ind_evecs_L);
repeat_times = sum(ind_evals_R == ind_evals_R(ind_evecs_R(1)))*sum(ind_evecs_R == ind_evecs_R(1));    % times of eigenvectors repeat in sorted_L
sorted_L = sorted_L(:,1:repeat_times:end);  sorted_R = V_RHS(:,unique(ind_evals_R));
dist = sqrt( sum( (abs(sorted_L)-abs(sorted_R)).^2 ) );
max(dist) < tol;   % gives 1
