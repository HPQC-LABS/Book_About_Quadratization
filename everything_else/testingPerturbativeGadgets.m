delta_required = zeros(7,1); tol = 1e-03; k = 1;    % fundamental settings where k is the number of auxiliary qubit
for Delta = 1:1e9:1e13
[LHS,RHS] = lhs2rhs('zzz',Delta,'P(3->2)DC2');
[V_RHS,E_RHS] = eig(RHS);
[V_LHS,E_LHS] = eig(LHS);
E_RHS = diag(E_RHS);
E_LHS = diag(E_LHS);
[ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
  if isempty(ind_evals_L) == 0        % matching eigenvalues exist
    sorted_E_RHS = sort(E_RHS);
    if (sum( abs(E_LHS - sorted_E_RHS(1)) < tol ) >= 2^k) % value of Delta that let ground energy match
        if (delta_required(1) == 0)
            delta_required(1) = Delta;
        end
        if (sum( abs(E_LHS - sorted_E_RHS(2^k + 1)) < tol ) >= 2^k) % value of Delta that let first excited energy match
            if (delta_required(3) == 0)
                delta_required(3) = Delta;
            end
        end
    end
    if (delta_required(5) == 0) & (numel(unique(ind_evals_R)) == 2^(k+2))     % value of Delta that let all 8 energies match
        delta_required(5) = Delta;
    end

    L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
    R = V_RHS(:,ind_evals_R);
    [ind_evecs_L, ind_evecs_R] = find( abs(L'*R)-1 > -tol );                     % indices of L with matiching eigenvectors in R
    sorted_L = L(:,ind_evecs_L);
    if isempty(ind_evecs_L) == 0      % matching eigenvectors exist
      sorted_L = unique(sorted_L','rows','stable')';   % remove repeated eigenvectors while maintaining the order
      sorted_R = V_RHS(:,unique(ind_evals_R));
      if (size(sorted_L,2) >= 2^k) & (delta_required(2) == 0) & (sum( sqrt( sum((abs(sorted_L(:,1:2^k))-abs(sorted_R(:,1:2^k))).^2) ) < tol ) == 2^k)  % value of Delta that let ground state match
          delta_required(2) = Delta;
      end
      if (size(sorted_L,2) >= 2^(k+1)) & (delta_required(4) == 0) & (sum( sqrt( sum((abs(sorted_L(:,1:2^(k+1)))-abs(sorted_R(:,1:2^(k+1)))).^2) ) < tol ) == 2^(k+1))   % value of Delta that let first excited state match
          delta_required(4) = Delta;
      end
      if isequal(size(sorted_L,2),size(sorted_R,2),2^(k+2)) & (delta_required(6) == 0) & (max(sqrt( sum((abs(sorted_L)-abs(sorted_R)).^2) )) < tol)  % value of Delta that let all 8 states match
          delta_required(6) = Delta;
      end
    end
  end
  if ne(delta_required(5),0) & ne(delta_required(6),0) & (delta_required(7) == 0) & (numel(unique(ind_evals_R)) == 0)  % value of Delta that is too large to keep energies matching
      delta_required(7) = Delta;
  end
end
delta_required(delta_required == 0) = nan;

%% Test results

% DC1: zzz still failed, max(min(dist)) always equals to 1.4142 even when all 16 energies match
% KKR: failed for vectors comparison as well
% DC2: We have tested zzz (3.46e12), zzx(1.52e12), zzy (1.22e12), zxz (1.84e12), zyz (1.833e12), xzz (1.735e12), yzz (1.668e12)
%     xxx(nan), xxz (2.06375e12), xxy(nan), xzx(nan), xyx(nan), zxx (6.769e12), yxx(nan)
%     yyy(nan), zyy (6.5480e12), xyy(nan), yzy(nan), yxy(nan), yyx(nan), yyz(nan)
%     xzy(nan), xyz(nan), zxy (6.0378e12), zyx (5.46549e12), yxz(nan), yzx(nan) for all eigenvectors and eigenvalues to match.
