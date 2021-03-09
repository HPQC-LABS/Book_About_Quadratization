function delta_required = FindReqdDelta(tol,coefficient,input_choice,minDelta,maxDelta,name_of_quadratization)
% This function takes in 6 arguments and gives a double matrix as output.
% It finds the Delta values that meet requriements stated below, in the range provide, and prints the output to a text file.
%
% e.g.    delta_required = FindReqdDelta(1e-3,-1,'xyz',1e10,1e13,'P(3->2)DC2')
%         refers to find the Delta required within 1e10-1e13 for -xyz using DC2

if strcmp(input_choice,"all_cubics") || strcmp(input_choice,"27_comb")     % test all 27 cubic terms
  delta_required = zeros(7,27); combinations = cell(1,27);
  S{1} = 'x'; S{2} = 'y'; S{3} = 'z'; n_combination = 1;

  for s1 = 1:3
  for s2 = 1:3
  for s3 = 1:3
  operators = [S{s1} S{s2} S{s3}];
  if coefficient == -1
  combinations(n_combination) = {['-' operators]};
  else
  combinations(n_combination) = {operators};
  end

  Delta = minDelta;
  while Delta <= maxDelta
  Delta = Delta + 10^(floor(log10(Delta))-4);
  [LHS,RHS] = lhs2rhs(coefficient,operators,Delta,name_of_quadratization);
  if isnan(RHS) == 0
    m = log2(size(RHS,2)) - length(operators);             % the number of auxiliary qubits
    [V_RHS,E_RHS] = eig(RHS);
    [V_LHS,E_LHS] = eig(LHS);
    [E_RHS,index] = sort(diag(E_RHS)); V_RHS = V_RHS(:,index);
    [E_LHS,index] = sort(diag(E_LHS)); V_LHS = V_LHS(:,index);
    [ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
    if isempty(ind_evals_L) == 0        % matching eigenvalues exist
      L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
      R = V_RHS(:,ind_evals_R);
      ind_evecs = find( sqrt(sum( (abs(L)-abs(R)).^2 ) ) < tol);
      L = L(:,ind_evecs);
      [~,index] = unique(L','rows','first');
      sorted_L = L(:,sort(index));                                            % remove repeated eigenvectors while maintaining the order
      R = R(:,ind_evecs);
      sorted_R = R(:,sort(index));
      if (sum( abs(E_LHS(1:2^m) - E_RHS(1)) < tol ) == 2^m)
          if (delta_required(1,n_combination) == 0)
              delta_required(1,n_combination) = Delta;  % value of Delta that let ground energy match
          end
          if (isempty(ind_evecs) == 0) && (delta_required(2,n_combination) == 0)
              delta_required(2,n_combination) = Delta;  % value of Delta that let ground state match
          end

          if (sum( abs(E_LHS((2^m + 1):2^(m+1)) - E_RHS(2)) < tol ) == 2^m)
              if (delta_required(3,n_combination) == 0)
                  delta_required(3,n_combination) = Delta;  % value of Delta that let first excited energy match
              end

              if (size(sorted_L,2) >= 2) && (delta_required(4,n_combination) == 0)
                  delta_required(4,n_combination) = Delta;   % value of Delta that let first excited state match
              end
          end

          if (numel(unique(ind_evals_R)) == 8)
              if (delta_required(5,n_combination) == 0)
                delta_required(5,n_combination) = Delta;   % value of Delta that let all 8 energies match
              end

              if isequal(size(sorted_L,2),size(sorted_R,2),8) && (delta_required(6,n_combination) == 0)
                  delta_required(6,n_combination) = Delta;  % value of Delta that let all 8 states match
              end
          end
        end
      if ne(delta_required(5,n_combination),0) && (delta_required(7,n_combination) == 0) && (numel(unique(ind_evals_R)) == 0)
        delta_required(7,n_combination) = Delta;  % value of Delta that is too large to keep energies matching
      end
    end
  end
  end       % use endwhile in Octave
  n_combination = n_combination + 1;
  end
  end
  end
  delta_required(delta_required == 0) = nan;
  UpdateOutput(tol,combinations,name_of_quadratization);

elseif length(input_choice) == 3        % test a single term
  delta_required = zeros(7,27); combinations = cell(1,27);
  operators = input_choice; n_combination = 1;   % to match the description in 'all_cubics' case

  if coefficient == -1
  combinations = {['-' operators]};
  else
  combinations = {operators};
  end

  Delta = minDelta;
  while Delta <= maxDelta
  Delta = Delta + 10^(floor(log10(Delta))-4);
  [LHS,RHS] = lhs2rhs(coefficient,operators,Delta,name_of_quadratization);
  if isnan(RHS) == 0
    m = log2(size(RHS,2)) - length(operators);             % the number of auxiliary qubits
    [V_RHS,E_RHS] = eig(RHS);
    [V_LHS,E_LHS] = eig(LHS);
    [E_RHS,index] = sort(diag(E_RHS)); V_RHS = V_RHS(:,index);
    [E_LHS,index] = sort(diag(E_LHS)); V_LHS = V_LHS(:,index);
    [ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
    if isempty(ind_evals_L) == 0        % matching eigenvalues exist
      L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
      R = V_RHS(:,ind_evals_R);
      ind_evecs = find( sqrt(sum( (abs(L)-abs(R)).^2 ) ) < tol);
      L = L(:,ind_evecs);
      [~,index] = unique(L','rows','first');
      sorted_L = L(:,sort(index));                                            % remove repeated eigenvectors while maintaining the order
      R = R(:,ind_evecs);
      sorted_R = R(:,sort(index));
      if (sum( abs(E_LHS(1:2^m) - E_RHS(1)) < tol ) == 2^m)
          if (delta_required(1,n_combination) == 0)
              delta_required(1,n_combination) = Delta;  % value of Delta that let ground energy match
          end
          if (isempty(ind_evecs) == 0) && (delta_required(2,n_combination) == 0)
              delta_required(2,n_combination) = Delta;  % value of Delta that let ground state match
          end

          if (sum( abs(E_LHS((2^m + 1):2^(m+1)) - E_RHS(2)) < tol ) == 2^m)
              if (delta_required(3,n_combination) == 0)
                  delta_required(3,n_combination) = Delta;  % value of Delta that let first excited energy match
              end

              if (size(sorted_L,2) >= 2) && (delta_required(4,n_combination) == 0)
                  delta_required(4,n_combination) = Delta;   % value of Delta that let first excited state match
              end
          end

          if (numel(unique(ind_evals_R)) == 8)
              if (delta_required(5,n_combination) == 0)
                delta_required(5,n_combination) = Delta;   % value of Delta that let all 8 energies match
              end

              if isequal(size(sorted_L,2),size(sorted_R,2),8) && (delta_required(6,n_combination) == 0)
                  delta_required(6,n_combination) = Delta;  % value of Delta that let all 8 states match
              end
          end
        end
      if ne(delta_required(5,n_combination),0) && (delta_required(7,n_combination) == 0) && (numel(unique(ind_evals_R)) == 0)
        delta_required(7,n_combination) = Delta;  % value of Delta that is too large to keep energies matching
      end
    end
  end
  end       %   use endwhile in Octave
  delta_required(delta_required == 0) = nan;
  UpdateOutput(tol,combinations,name_of_quadratization);

else
  disp("Invalid Input");
end

end


function UpdateOutput(tol,combinations,name_of_quadratization)      % Update the output text file

if exist('Required_Delta.txt','file') == 0    % if the file doesn't exist, we will create a new one
RowNames_1 = strcat(name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(1));      % gives the name of row as 'P(3->2)DC2_1e-3_1'
RowNames_2 = strcat(name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(2));
T = table([{'combinations'};{RowNames_1};{RowNames_2}],[combinations;num2cell(delta_required(1,:));num2cell(delta_required(2,:))],cell(1,size(T,2)));
T = table2array(T);

else      % need to be further improved so that we can get it updated during executing
InputT = readtable('Required_Delta.txt','ReadVariableNames',false);
InputT = table2array(InputT)';
RowNames_1 = strcat(name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(1));      % gives the name of row as 'P(3->2)DC2_1e-3_1'
RowNames_2 = strcat(name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(2));
T = table([{'combinations'};{RowNames_1};{RowNames_2}],[combinations;num2cell(delta_required(1,:));num2cell(delta_required(2,:))]);
T = table2array(T);
T = [InputT;T;cell(1,size(T,2))];

end
T = array2table(T');
writetable(T,'Required_Delta.txt','Delimiter','\t','WriteVariableNames',false);

end
