function FindReqdDelta(tol,coefficient,decimal_place,input_choice,minDelta,maxDelta,name_of_quadratization,test_times)
% This function takes in 6 arguments and gives a double matrix as output.
% It finds the Delta values that meet requriements stated below, in the range provide, and prints the output to a text file.
%
% e.g.    delta_required = FindReqdDelta(1e-3,-1,3,'xyz',1e10,1e13,'P(3->2)DC2',1)
%         refers to find the Delta required within 1e10-1e13 (3 decimal places) for -xyz using DC2 the 1st time

if strcmp(input_choice,"all_cubics") || strcmp(input_choice,"27_comb")     % test all 27 cubic terms
  delta_required = zeros(7,27); n_combination = 1;
  Op{1} = 'x'; Op{2} = 'y'; Op{3} = 'z';   % Op{:} is the array of operators

  if coefficient == -1
  FileName = strcat(input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
  else
  FileName = strcat(input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
  end

  for s1 = 1:3
  for s2 = 1:3
  for s3 = 1:3
  operators = strcat(Op{s1},Op{s2},Op{s3});

  % calculate S1,S2,S3 outside of the function lhs2rhs
    x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
    n = length(operators);
    S = cell(n);
    [m,NeededM] = GetAuxNum(size(operators,2),name_of_quadratization);                %  m is the number of auxiliary qubits nad A is the array of auxiliary operators
    for ind = 1:n
        if operators(ind) == 'x'
            S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(n+m-ind)));
        elseif operators(ind) == 'y'
            S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(n+m-ind)));
        elseif operators(ind) == 'z'
            S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+m-ind)));
        end
    end
  LHS = coefficient*S{1}*S{2}*S{3};
  NeededM{end} = LHS;              % save LHS in the last cell of NeededM
  Delta = minDelta;
  checkpoint = minDelta;
  while Delta <= maxDelta
  Delta = Delta + 10^(floor(log10(Delta))-decimal_place);
  [LHS,RHS] = lhsrhs(coefficient,S,NeededM,Delta,name_of_quadratization);
  if isnan(RHS) == 0
    m = log2(size(RHS,2)) - n;             % the number of auxiliary qubits
    [V_RHS,E_RHS] = eig(RHS);
    [V_LHS,E_LHS] = eig(LHS);
    [E_RHS,index] = sort(diag(E_RHS));
    V_RHS = V_RHS(:,index);
    [E_LHS,index] = sort(diag(E_LHS));
    V_LHS = V_LHS(:,index);
    [ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
    if isempty(ind_evals_L) == 0        % matching eigenvalues exist
      L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
      R = V_RHS(:,ind_evals_R);
      ind_evecs = find( sqrt(sum( (abs(L)-abs(R)).^2 ) ) < tol);
      if isempty(ind_evecs) == 0        % matching eigenvectors exist
        L = L(:,ind_evecs);
        % [~,index] = unique(L','rows','first');
        [index] = UniqueRows(L');                                               % faster unique function
        sorted_L = L(:,sort(index));                                            % remove repeated eigenvectors while maintaining the order
        R = R(:,ind_evecs);
        sorted_R = R(:,sort(index));
      else
        sorted_L = [];
        sorted_R = [];
      end
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

          if (numel(unique(ind_evals_R)) == 2^n)
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
  end      % use end in MATLAB / endwhile in Octave
  ind = find(delta_required(:,n_combination) == 0);
  delta_required(ind,n_combination) = Inf;
  n_combination = n_combination + 1;
  end
  end
  end
  dlmwrite(FileName,delta_required','delimiter','\t','newline','unix');
  dlmwrite(FileName,maxDelta,'-append','delimiter','\t');     % range of Delta tested

elseif strcmp(input_choice,'all_quartics')
  delta_required = zeros(7,81); n_combination = 1;
  Op{1} = 'x'; Op{2} = 'y'; Op{3} = 'z';   % Op{:} is the array of operators

  if coefficient == -1
  FileName = strcat(input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
  else
  FileName = strcat(input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
  end

  for s1 = 1:3
  for s2 = 1:3
  for s3 = 1:3
  for s4 = 1:3
  operators = strcat(Op{s1},Op{s2},Op{s3},Op{s4});

  % calculate S1,S2,S3 outside of the function lhs2rhs
    x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
    n = length(operators);
    S = cell(n);
    [m,NeededM] = GetAuxNum(size(operators,2),name_of_quadratization);                %  m is the number of auxiliary qubits nad A is the array of auxiliary operators
    for ind = 1:n
        if operators(ind) == 'x'
            S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(n+m-ind)));
        elseif operators(ind) == 'y'
            S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(n+m-ind)));
        elseif operators(ind) == 'z'
            S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+m-ind)));
        end
    end
  A = S{1}*S{2};
  B = S{3}*S{4};
  LHS = A*B;
  NeededM{end - 2} = A;
  NeededM{end - 1} = B;
  NeededM{end} = LHS;              % save LHS in the last cell of NeededM
  Delta = minDelta;
  checkpoint = minDelta;
  while Delta <= maxDelta
  Delta = Delta + 10^(floor(log10(Delta))-decimal_place);
  [LHS,RHS] = lhsrhs(coefficient,S,NeededM,Delta,name_of_quadratization);
  if isnan(RHS) == 0
    m = log2(size(RHS,2)) - n;             % the number of auxiliary qubits
    [V_RHS,E_RHS] = eig(RHS);
    [V_LHS,E_LHS] = eig(LHS);
    [E_RHS,index] = sort(diag(E_RHS)); V_RHS = V_RHS(:,index);
    [E_LHS,index] = sort(diag(E_LHS)); V_LHS = V_LHS(:,index);
    [ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );              % indices of LHS and RHS where eigenvalues match within tol
    if isempty(ind_evals_L) == 0        % matching eigenvalues exist
      L = V_LHS(:,ind_evals_L);                                                   % L = LHS eigenvectors for eigenvalues matching RHS
      R = V_RHS(:,ind_evals_R);
      ind_evecs = find( sqrt(sum( (abs(L)-abs(R)).^2 ) ) < tol);
      if isempty(ind_evecs) == 0        % matching eigenvectors exist
        L = L(:,ind_evecs);
        % [~,index] = unique(L','rows','first');
        [index] = UniqueRows(L');                                               % faster unique function
        sorted_L = L(:,sort(index));                                            % remove repeated eigenvectors while maintaining the order
        R = R(:,ind_evecs);
        sorted_R = R(:,sort(index));
      else
        sorted_L = [];
        sorted_R = [];
      end
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

          if (numel(unique(ind_evals_R)) == 2^n)
              if (delta_required(5,n_combination) == 0)
                delta_required(5,n_combination) = floor(log10(Delta));   % value of Delta that let all 16 energies match
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
  end      % use end in MATLAB / endwhile in Octave
  ind = find(delta_required(:,n_combination) == 0);
  delta_required(ind,n_combination) = Inf;
  n_combination = n_combination + 1;
  end
  end
  end
  end
  dlmwrite(FileName,delta_required','delimiter','\t','newline','unix');
  dlmwrite(FileName,maxDelta,'-append','delimiter','\t');     % range of Delta tested

elseif length(input_choice) == 3        % test a single term
  delta_required = zeros(7,1); operators = input_choice; n_combination = 1;   % to match the description in 'all_cubics' case

  if coefficient == -1
  FileName = strcat(input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
  else
  FileName = strcat(input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
  end

  % calculate S1,S2,S3 outside of the function lhs2rhs
    x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
    n = length(operators);
    S = cell(n);
    [m,NeededM] = GetAuxNum(size(operators,2),name_of_quadratization);                %  m is the number of auxiliary qubits nad A is the array of auxiliary operators
    for ind = 1:n
        if operators(ind) == 'x'
            S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(n+m-ind)));
        elseif operators(ind) == 'y'
            S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(n+m-ind)));
        elseif operators(ind) == 'z'
            S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+m-ind)));
        end
    end
  LHS = coefficient*S{1}*S{2}*S{3};
  NeededM{end} = LHS;
  Delta = minDelta;
  checkpoint = minDelta;
  while Delta <= maxDelta
  Delta = Delta + 10^(floor(log10(Delta))-decimal_place);
  [LHS,RHS] = lhsrhs(coefficient,S,NeededM,Delta,name_of_quadratization);
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
      if isempty(ind_evecs) == 0        % matching eigenvectors exist
        L = L(:,ind_evecs);
        % [~,index] = unique(L','rows','first');
        [index] = UniqueRows(L');                                               % faster unique function
        sorted_L = L(:,sort(index));                                            % remove repeated eigenvectors while maintaining the order
        R = R(:,ind_evecs);
        sorted_R = R(:,sort(index));
      else
        sorted_L = [];
        sorted_R = [];
      end
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

          if (numel(unique(ind_evals_R)) == 2^n)
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
  end      % use end in MATLAB / endwhile in Octave
  ind = find(delta_required(:,n_combination) == 0);
  delta_required(ind,n_combination) = Inf;
  dlmwrite(FileName,delta_required','delimiter','\t','newline','unix');
  dlmwrite(FileName,maxDelta,'-append','delimiter','\t');     % range of Delta tested
  dlmwrite(FileName,input_choice,'-append');     % input choice

else
  disp("Invalid Input");
end

end

function [number_of_auxiliary,NeededM] = GetAuxNum(number_of_logical,name_of_quadratization)     % known from Book About Quadratization
% need to be updated whenever lhs2rhs.m gets updated
    x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
    n = number_of_logical;

    if strcmp(name_of_quadratization, 'P(3->2)-DC2') || strcmp(name_of_quadratization, 'P(3->2)DC2')
        number_of_auxiliary = 1;
        xa = kron(eye(8),x); za = kron(eye(8),z); I = eye(2^4);
        NeededM = cell(4,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'P(3->2)-KKR') || strcmp(name_of_quadratization, 'P(3->2)KKR')
        number_of_auxiliary = 3;
        xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
        za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);
        I = eye(2^6);
        NeededM = cell(8,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa1; NeededM{2} = xa2; NeededM{3} = xa3;
        NeededM{4} = za1; NeededM{5} = za2; NeededM{6} = za3; NeededM{7} = I;

    elseif strcmp(name_of_quadratization, 'P(3->2)KKR-A') % no coefficient needed
        number_of_auxiliary = 3;
        xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
        za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);
        I = eye(2^6);
        NeededM = cell(8,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa1; NeededM{2} = xa2; NeededM{3} = xa3;
        NeededM{4} = za1; NeededM{5} = za2; NeededM{6} = za3; NeededM{7} = I;

    elseif strcmp(name_of_quadratization, 'P(3->2)-DC1') || strcmp(name_of_quadratization, 'P(3->2)DC1')
        number_of_auxiliary = 3;
        xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
        za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);
        I = eye(2^6);
        NeededM = cell(8,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa1; NeededM{2} = xa2; NeededM{3} = xa3;
        NeededM{4} = za1; NeededM{5} = za2; NeededM{6} = za3; NeededM{7} = I;

   elseif strcmp(name_of_quadratization, 'ZZZ-TI-CBBK')
        number_of_auxiliary = 1;
        xa = kron(eye(8),x); za = kron(eye(8),z); I = eye(2^4);
        NeededM = cell(4,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'PSD-CBBK')
        number_of_auxiliary = 1;
        xa = kron(eye(2^n),x); za = kron(eye(2^n),z); I = eye(2^(n+1));
        NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'PSD-OT')
        number_of_auxiliary = 1;
        xa = kron(eye(2^n),x); za = kron(eye(2^n),z); I = eye(2^(n+1));
        NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'P(3->2)CBBK') || strcmp(name_of_quadratization, 'P(3->2)-CBBK')
         number_of_auxiliary = 1;
         za = kron(eye(8),z); xa = kron(eye(8),x); I = eye(16);
         NeededM = cell(4,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
         NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'P(3->2)-OT') || strcmp(name_of_quadratization, 'P(3->2)OT')
        number_of_auxiliary = 1;
        xa = kron(eye(8),x); za = kron(eye(8),z); I = eye(16);
        NeededM = cell(4,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'P1B1-CBBK')
        number_of_auxiliary = 1;
        xa = kron(eye(2^n),x); za = kron(eye(2^n),z); I = eye(2^(n+1));
        NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'P1B1-OT')
        number_of_auxiliary = 1;
        xa = kron(eye(2^n),x); za = kron(eye(2^n),z); I = eye(2^(n+1));
        NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'PSD-CN')
        number_of_auxiliary = 2;
        za_11 = kron(kron(eye(2^(n+0)),z),eye(2)); xa_11 = kron(kron(eye(2^(n+0)),x),eye(2));
        za_1 = kron(kron(eye(2^(n+1)),z),eye(1)); I = eye(2^(n+2));

        NeededM = cell(7,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = za_11; NeededM{2} = xa_11; NeededM{3} = za_1; NeededM{4} = I;
    else
        disp('cannot find this method');
        number_of_auxiliary = nan;
    end
end

function index = UniqueRows(A)        % This MATLAB function provides a faster version of MATLAB's unique rows method (i.e. 'unique(points,''rows'')')
% [index] = UniqueRows(A); % equivalent is '[~,ind]=unique(A,'rows');'
% A_Unique = A(sort(ind),:);
% faster unique function (https://www.mathworks.com/matlabcentral/fileexchange/77329-fast-unique-rows-method-with-unit-and-performance-tests)

[A_sorted,idx1] = sortrows(A);
k = find([true; any(diff(A_sorted,1,1),2); true]);
idx2 = k(diff(k) >= 1);
index=idx1(idx2);

end
