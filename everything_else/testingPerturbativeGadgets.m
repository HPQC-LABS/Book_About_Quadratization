[LHS,RHS] = lhs2rhs('zzz',2e13,'P(3->2)DC2'); tol = 1e-03;
[V_RHS,E_RHS] = eig(RHS); [V_LHS,E_LHS] = eig(LHS); 
E_RHS = diag(E_RHS); E_LHS = diag(E_LHS);
E_diff = abs(bsxfun(@minus,E_RHS',E_LHS));       % find matching eigenvalues by calculating differences between every elements in E_RHS and E_LHS[
[ind_LHS, ind_RHS] = find(E_diff < tol);         % that tells us which elements of E_LHS and E_RHS match, respectively

% unsolved problem: we use abs() here to limit the impact of a vector's sign, but it does not help for those need scalar multiplication (e.g. [0 0 2 0] and [0 0 1 0])
dist = 2*ones(max(ind_LHS),max(ind_RHS));        % set dist to 2 by default so that we can find matching eigenvectors easily in next step
for ind = 1:size(ind_LHS,1)                      % calculate the distance between two eigenvectors whose eigevalues match, and return the result to corresponding cell
    dist(ind_LHS(ind),ind_RHS(ind)) = norm(abs(V_RHS(:,ind_RHS(ind)))-abs(V_LHS(:,ind_LHS(ind)))); % rows of dist -> columns of V_LHS, columns of dist -> columns of V_RHS
end

[dist_ind_LHS, dist_ind_RHS] = find(dist < tol); % dist_ind_LHS and dist_ind_RHS refers to rows and columns of dist, which have elements smaller than 1e-3 (i.e., dist_ind_LHS and dist_ind_RHS tells us which columns of V_LHS match the columns of V_RHS)
max(min(dist(:,min(ind_RHS):max(ind_RHS))));     % farthest distance between two matching eigenvectors is 7.3682e-05
sorted_V_LHS = V_LHS(:,dist_ind_LHS);
sorted_E_LHS = E_LHS(dist_ind_LHS);

%% Test results

% DC1: zzz still failed, max(min(dist)) always equals to 1.4142 even when all 16 energies match
% KKR: failed for vectors comparison as well
% DC2: We have tested zzz, zxx (4e12), zxz (4.5e12) etc., but need loop to find appropriate and more precise Delta value for zzz, zzx etx.