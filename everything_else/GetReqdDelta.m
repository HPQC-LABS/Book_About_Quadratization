function delta_required = GetReqdDelta(operators,Delta_lower_bound,Delta_increment,Delta_upper_bound,name_of_quadratization,tol,k)

delta_required = Delta_lower_bound*ones(1,6);                   % create a default row for the table which contains the result of testDelta()
test_times = (Delta_upper_bound - Delta_lower_bound)/Delta_increment;
parfor i = 0:test_times
  Delta = Delta_lower_bound + Delta_increment*i;                 % set the value of Delta indirectly because i can only increase by 1
  delta_required_for1 = testDelta(operators,Delta,name_of_quadratization,tol,k);    % test the requirements that each Delta meets
  delta_required = [delta_required;delta_required_for1];
end
delta_required(1,:) = [];           % remove the default row
delta_required = min(delta_required)';

end

function result = testDelta(operators,Delta,name_of_quadratization,tol,k)  % test a Delta value to see if it meets the requirements stated below
requriements_met = zeros(1,6);
[LHS,RHS] = lhs2rhs(operators,Delta,name_of_quadratization);
[V_RHS,E_RHS] = eig(RHS);
[V_LHS,E_LHS] = eig(LHS);
[E_RHS,index] = sort(diag(E_RHS)); V_RHS = V_RHS(:,index);
[E_LHS,index] = sort(diag(E_LHS)); V_LHS = V_LHS(:,index);
[ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
if isempty(ind_evals_L) == 0        % matching eigenvalues exist
    L = V_LHS(:,ind_evals_L);                                               % L = LHS eigenvectors for eigenvalues matching RHS
    R = V_RHS(:,ind_evals_R);
    ind_evecs = find( sqrt(sum( (abs(L)-abs(R)).^2 ) ) < tol);
    L = L(:,ind_evecs);
    [~,index] = unique(L','rows','first');
    sorted_L = L(:,sort(index));                                            % remove repeated eigenvectors while maintaining the order
    R = R(:,ind_evecs);
    sorted_R = R(:,sort(index));
    if (sum( abs(E_LHS(1:2^k) - E_RHS(1)) < tol ) == 2^k)
            requriements_met(1) = Delta;  % can let ground energy match
            if (isempty(ind_evecs) == 0)
                requriements_met(2) = Delta;  % can let ground state match
            end

            if (sum( abs(E_LHS((2^k + 1):2^(k+1)) - E_RHS(2)) < tol ) == 2^k)
                requriements_met(3) = Delta;  % can let first excited energy match
                if (size(sorted_L,2) >= 2)
                    requriements_met(4) = Delta;   % can let first excited state match
                end
            end

            if (numel(unique(ind_evals_R)) == 8)
                requriements_met(5) = Delta;   % can let all 8 energies match
                if isequal(size(sorted_L,2),size(sorted_R,2),8)
                    requriements_met(6) = Delta;  % can let all 8 states match
                end
            end
    end
end
requriements_met(requriements_met == 0) = nan;
result = requriements_met;                  % records the value of Delta if it meets the requiremnt and NaN if it doesn't
end
