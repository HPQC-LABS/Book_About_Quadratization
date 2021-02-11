delta_required = zeros(7,1); tol = 1e-03; k = 1;    % fundamental settings where k is the number of auxiliary qubit
for Delta = 1:1e12:1e13
[LHS,RHS] = lhs2rhs('zzz',Delta,'P(3->2)DC2');
[V_RHS,E_RHS] = eig(RHS);
[V_LHS,E_LHS] = eig(LHS);
E_RHS = diag(E_RHS);
E_LHS = diag(E_LHS);
[ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
  if isempty(ind_evals_L) == 0        % matching eigenvalues exist
    sorted_E_RHS = sort(E_RHS);
    if (delta_required(1) == 0) & (sum( abs(E_LHS - sorted_E_RHS(1)) < tol ) >= 2^k) % value of Delta that let ground energy match
        delta_required(1) = Delta;
    end
    if (delta_required(3) == 0) & (sum( abs(E_LHS - sorted_E_RHS(2^k + 1)) < tol ) >= 2^k) % value of Delta that let first excited energy match
        delta_required(3) = Delta;
    end
    if (delta_required(5) == 0) & (numel(unique(ind_evals_R)) == 2^(k+2))     % value of Delta that let all 8 energies match
        delta_required(5) = Delta;
    end

    L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
    R = V_RHS(:,ind_evals_R);
    [ind_evecs_L, ind_evecs_R] = find( abs(L'*R)-1 > -tol );                     % indices of L with matiching eigenvectors in R

    if isempty(ind_evecs_L) == 0      % matching eigenvectors exist
      sorted_L = L(:,ind_evecs_L);
      repeat_times = sum(ind_evals_R == ind_evals_R(1,1))*sum(ind_evecs_R == ind_evals_R(1,1));    % times of eigenvectors repeat in sorted_L
      sorted_L = sorted_L(:,1:repeat_times:end);
      sorted_R = V_RHS(:,unique(ind_evals_R));
      if size(sorted_L,2) == size(sorted_R,2)
          max(sqrt( sum( (abs(sorted_L)-abs(sorted_R)).^2 ) )) < tol;   % maximum distance between two matching eigenvectors is within the bound
      end
    end
  end
  if ne(delta_required(5),0) & (delta_required(7) == 0) & (numel(unique(ind_evals_R)) == 0)     % value of Delta that is too large to keep energies matching
      delta_required(7) = Delta;
  end
end
delta_required(delta_required == 0) = nan;

%% Test results

% DC1: zzz still failed, max(min(dist)) always equals to 1.4142 even when all 16 energies match
% KKR: failed for vectors comparison as well
% DC2: We have tested xxz (2.05e12), xxx (2.488e12), zyz (1.833e12), yyy (2.114e12) etc. successfully for eigenvalues matching.

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
