delta_required = zeros(7,1); tol = 1e-03; k = 1;    % fundamental settings where k is the number of auxiliary qubit
for Delta = 1e9:1e9:1e13
[LHS,RHS] = lhs2rhs('xxx',Delta,'P(3->2)DC2');
[V_RHS,E_RHS] = eig(RHS);
[V_LHS,E_LHS] = eig(LHS);
E_RHS = sort(diag(E_RHS));
E_LHS = sort(diag(E_LHS));
[ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
if isempty(ind_evals_L) == 0        % matching eigenvalues exist
    L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
    R = V_RHS(:,ind_evals_R);
    [ind_evecs_L, ind_evecs_R] = find( abs(L'*R)-1 > -tol );                     % indices of L with matiching eigenvectors in R
    sorted_L = unique(L(:,ind_evecs_L)','rows','stable')';                       % remove repeated eigenvectors while maintaining the order
    sorted_R = V_RHS(:,unique(ind_evals_R));

    if (sum( abs(E_LHS(1:2^k) - E_RHS(1)) < tol ) == 2^k)
        if (delta_required(1) == 0)
            delta_required(1) = Delta;  % value of Delta that let ground energy match
        end
        if (isempty(ind_evecs_L) == 0) & (delta_required(2) == 0) & (norm( abs(sorted_L(:,1))-abs(sorted_R(:,1)) ) < tol)
            delta_required(2) = Delta;  % value of Delta that let ground state match
        end

        if (sum( abs(E_LHS((2^k + 1):2^(k+1)) - E_RHS(2)) < tol ) == 2^k)
            if (delta_required(3) == 0)
                delta_required(3) = Delta;  % value of Delta that let first excited energy match
            end

            if (size(sorted_L,2) >= 2) & (delta_required(4) == 0) & (norm( abs(sorted_L(:,2))-abs(sorted_R(:,2))) < tol )
                delta_required(4) = Delta;   % value of Delta that let first excited state match
            end
        end

        if (numel(unique(ind_evals_R)) == 8)
            if (delta_required(5) == 0)
            delta_required(5) = Delta;   % value of Delta that let all 8 energies match
            end

            if isequal(size(sorted_L,2),size(sorted_R,2),8) & (delta_required(6) == 0) & (max( sqrt(sum( (abs(sorted_L)-abs(sorted_R)).^2 )) ) < tol)
                delta_required(6) = Delta;  % value of Delta that let all 8 states match
            end
        end
    end
end
if ne(delta_required(5),0) & ne(delta_required(6),0) & (delta_required(7) == 0) & (numel(unique(ind_evals_R)) == 0)
    delta_required(7) = Delta;  % value of Delta that is too large to keep energies matching
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


%% Alternative: use euclidean distance to find matching eigenvectors

delta_required = zeros(7,1); tol = 1e-03; k = 1;    % fundamental settings where k is the number of auxiliary qubit
for Delta = 1e9:1e9:1e13
[LHS,RHS] = lhs2rhs('xxx',Delta,'P(3->2)CBBK');
[V_RHS,E_RHS] = eig(RHS);
[V_LHS,E_LHS] = eig(LHS);
E_RHS = sort(diag(E_RHS));
E_LHS = sort(diag(E_LHS));
[ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
if isempty(ind_evals_L) == 0        % matching eigenvalues exist
    L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
    R = V_RHS(:,ind_evals_R);
    ind_evecs = find( sqrt(sum( (abs(L)-abs(R)).^2 ) ) < tol);
    sorted_L = unique(L(:,ind_evecs)','rows','stable')';                       % remove repeated eigenvectors while maintaining the order
    sorted_R = unique(R(:,ind_evecs)','rows','stable')';
    if (sum( abs(E_LHS(1:2^k) - E_RHS(1)) < tol ) == 2^k)
            if (delta_required(1) == 0)
                delta_required(1) = Delta;  % value of Delta that let ground energy match
            end
            if (isempty(ind_evecs) == 0) & (delta_required(2) == 0)
                delta_required(2) = Delta;  % value of Delta that let ground state match
            end

            if (sum( abs(E_LHS((2^k + 1):2^(k+1)) - E_RHS(2)) < tol ) == 2^k)
                if (delta_required(3) == 0)
                    delta_required(3) = Delta;  % value of Delta that let first excited energy match
                end

                if (size(sorted_L,2) >= 2) & (delta_required(4) == 0)
                    delta_required(4) = Delta;   % value of Delta that let first excited state match
                end
            end

            if (numel(unique(ind_evals_R)) == 8)
                if (delta_required(5) == 0)
                delta_required(5) = Delta;   % value of Delta that let all 8 energies match
                end

                if isequal(size(sorted_L,2),size(sorted_R,2),8) & (delta_required(6) == 0)
                    delta_required(6) = Delta;  % value of Delta that let all 8 states match
                end
            end
    end
end
if ne(delta_required(5),0) & ne(delta_required(6),0) & (delta_required(7) == 0) & (numel(unique(ind_evals_R)) == 0)
    delta_required(7) = Delta;  % value of Delta that is too large to keep energies matching
end
end
delta_required(delta_required == 0) = nan;
