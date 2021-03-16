function FindReqDelta(tol,coefficient,decimal_place,input_choice,minDelta,maxDelta,name_of_quadratization,test_times)
% This function takes in 6 arguments and gives a double matrix as output.
% It finds the Delta values that meet requriements stated below, in the range provide, and prints the output to a text file.
%
% e.g.    delta_required = FindReqdDelta(1e-3,-1,3,'xyz',1e10,1e13,'P(3->2)DC2',1)
%         refers to find the Delta required within 1e10-1e13 (3 decimal places) for -xyz using DC2 the 1st time

if strcmp(input_choice,"all_cubics") || strcmp(input_choice,"27_comb")     % test all 27 cubic terms
  delta_required = zeros(7,27); n_combination = 1;
  Op{1} = 'x'; Op{2} = 'y'; Op{3} = 'z';   % Op{:} is the array of operators

  if coefficient == -1
  FileName = strcat('negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
  else
  FileName = strcat('positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
  end

  for s1 = 1:3
  for s2 = 1:3
  for s3 = 1:3
  operators = strcat(Op{s1},Op{s2},Op{s3});

  % calculate S1,S2,S3 outside of the function lhs2rhs
    x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1]; S = cell(length(operators));
    m = GetAuxNum(name_of_quadratization);                % the number of auxiliary qubits
    for ind = 1:length(operators)
        if operators(ind) == 'x'
            S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(length(operators)+m-ind)));
        elseif operators(ind) == 'y'
            S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(length(operators)+m-ind)));
        elseif operators(ind) == 'z'
            S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(length(operators)+m-ind)));
        end
    end

  Delta = minDelta;
  checkpoint = minDelta;
  while Delta <= maxDelta
  Delta = Delta + 10^(floor(log10(Delta))-decimal_place);
  [LHS,RHS] = lhs2rhs(coefficient,S,Delta,name_of_quadratization);
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
              delta_required(1,n_combination) = floor(log10(Delta));  % value of Delta that let ground energy match
          end
          if (isempty(ind_evecs) == 0) && (delta_required(2,n_combination) == 0)
              delta_required(2,n_combination) = floor(log10(Delta));  % value of Delta that let ground state match
          end

          if (sum( abs(E_LHS((2^m + 1):2^(m+1)) - E_RHS(2)) < tol ) == 2^m)
              if (delta_required(3,n_combination) == 0)
                  delta_required(3,n_combination) = floor(log10(Delta));  % value of Delta that let first excited energy match
              end

              if (size(sorted_L,2) >= 2) && (delta_required(4,n_combination) == 0)
                  delta_required(4,n_combination) = floor(log10(Delta));   % value of Delta that let first excited state match
              end
          end

          if (numel(unique(ind_evals_R)) == 8)
              if (delta_required(5,n_combination) == 0)
                delta_required(5,n_combination) = floor(log10(Delta));   % value of Delta that let all 8 energies match
              end

              if isequal(size(sorted_L,2),size(sorted_R,2),8) && (delta_required(6,n_combination) == 0)
                  delta_required(6,n_combination) = floor(log10(Delta));  % value of Delta that let all 8 states match
              end
          end
        end
      if ne(delta_required(5,n_combination),0) && (delta_required(7,n_combination) == 0) && (numel(unique(ind_evals_R)) == 0)
        delta_required(7,n_combination) = floor(log10(Delta));  % value of Delta that is too large to keep energies matching
      end
    end
  end
  if (Delta - checkpoint) >= 10^(floor(log10(Delta))-(decimal_place - 2))    % update the output file and the checkpoint
  dlmwrite(FileName,delta_required','delimiter','\t','newline','unix');   % the delta required
  dlmwrite(FileName,Delta,'-append');     % largest Delta tested
  checkpoint = Delta;
  end
  endwhile      % use end in MATLAB
  n_combination = n_combination + 1;
  end
  end
  end
  delta_required(delta_required == 0) = 308;
  dlmwrite(FileName,delta_required','delimiter','\t','newline','unix');
  dlmwrite(FileName,maxDelta,'-append','delimiter','\t');     % range of Delta tested

elseif length(input_choice) == 3        % test a single term
  delta_required = zeros(7,1); operators = input_choice; n_combination = 1;   % to match the description in 'all_cubics' case

  if coefficient == -1     % set the name of the output file
    FileName = strcat('negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
    input_choice = strcat('-',input_choice);
  else
    FileName = strcat('positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
  end

% calculate S1,S2,S3 outside of the function lhs2rhs
  x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1]; S = cell(length(operators));
  m = GetAuxNum(name_of_quadratization);                % the number of auxiliary qubits
  for ind = 1:length(operators)
      if operators(ind) == 'x'
          S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(length(operators)+m-ind)));
      elseif operators(ind) == 'y'
          S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(length(operators)+m-ind)));
      elseif operators(ind) == 'z'
          S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(length(operators)+m-ind)));
      end
  end

  Delta = minDelta;
  checkpoint = minDelta;
  while Delta <= maxDelta
  Delta = Delta + 10^(floor(log10(Delta))-decimal_place);
  [LHS,RHS] = lhs2rhs(coefficient,S,Delta,name_of_quadratization);
  if isnan(RHS) == 0
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
              delta_required(1,n_combination) = floor(log10(Delta));  % value of Delta that let ground energy match
          end
          if (isempty(ind_evecs) == 0) && (delta_required(2,n_combination) == 0)
              delta_required(2,n_combination) = floor(log10(Delta));  % value of Delta that let ground state match
          end

          if (sum( abs(E_LHS((2^m + 1):2^(m+1)) - E_RHS(2)) < tol ) == 2^m)
              if (delta_required(3,n_combination) == 0)
                  delta_required(3,n_combination) = floor(log10(Delta));  % value of Delta that let first excited energy match
              end

              if (size(sorted_L,2) >= 2) && (delta_required(4,n_combination) == 0)
                  delta_required(4,n_combination) = floor(log10(Delta));   % value of Delta that let first excited state match
              end
          end

          if (numel(unique(ind_evals_R)) == 8)
              if (delta_required(5,n_combination) == 0)
                delta_required(5,n_combination) = floor(log10(Delta));   % value of Delta that let all 8 energies match
              end

              if isequal(size(sorted_L,2),size(sorted_R,2),8) && (delta_required(6,n_combination) == 0)
                  delta_required(6,n_combination) = floor(log10(Delta));  % value of Delta that let all 8 states match
              end
          end
        end
      if ne(delta_required(5,n_combination),0) && (delta_required(7,n_combination) == 0) && (numel(unique(ind_evals_R)) == 0)
        delta_required(7,n_combination) = floor(log10(Delta));  % value of Delta that is too large to keep energies matching
      end
    end
  end
  if (Delta - checkpoint) >= 10^(floor(log10(Delta))-(decimal_place - 2))        % update the output file and the checkpoint
  dlmwrite(FileName,delta_required','delimiter','\t','newline','unix');   % the delta required
  dlmwrite(FileName,Delta,'-append');     % largest Delta tested
  dlmwrite(FileName,input_choice,'-append');     % input choice
  checkpoint = Delta;
  end
  endwhile      %   use end in MATLAB
  delta_required(delta_required == 0) = 308;
  dlmwrite(FileName,delta_required','delimiter','\t','newline','unix');
  dlmwrite(FileName,maxDelta,'-append','delimiter','\t');     % range of Delta tested
  dlmwrite(FileName,input_choice,'-append');     % input choice

else
  disp("Invalid Input");
end

end

function number_of_auxiliary = GetAuxNum(name_of_quadratization)    % need to be updated whenever lhs2rhs.m gets updated
    if strcmp(name_of_quadratization, 'P(3->2)-DC2') || strcmp(name_of_quadratization, 'P(3->2)DC2')
        number_of_auxiliary = 1;                   % known from The Book About Quadratization
    elseif strcmp(name_of_quadratization, 'P(3->2)-KKR') || strcmp(name_of_quadratization, 'P(3->2)KKR')
        number_of_auxiliary = 3;
    elseif strcmp(name_of_quadratization, 'P(3->2)KKR-A') % no coefficient needed
        number_of_auxiliary = 3;
    elseif strcmp(name_of_quadratization, 'P(3->2)-DC1') || strcmp(name_of_quadratization, 'P(3->2)DC1')
        number_of_auxiliary = 3;
   elseif strcmp(name_of_quadratization, 'ZZZ-TI-CBBK')
        number_of_auxiliary = 1;
    elseif strcmp(name_of_quadratization, 'PSD-CBBK')
        number_of_auxiliary = 1;
    elseif strcmp(name_of_quadratization, 'PSD-OT')
        number_of_auxiliary = 1;
    elseif strcmp(name_of_quadratization, 'P(3->2)CBBK') || strcmp(name_of_quadratization, 'P(3->2)-CBBK')
         number_of_auxiliary = 1;
    elseif strcmp(name_of_quadratization, 'P(3->2)-OT') || strcmp(name_of_quadratization, 'P(3->2)OT')
        number_of_auxiliary = 1;
    elseif strcmp(name_of_quadratization, 'P1B1-CBBK')
        number_of_auxiliary = 1;
    elseif strcmp(name_of_quadratization, 'P1B1-OT')
        number_of_auxiliary = 1;
    elseif strcmp(name_of_quadratization, 'PSD-CN')
        number_of_auxiliary = 2;
    else
        disp('cannot find this method');
        number_of_auxiliary = nan;
    end
end
