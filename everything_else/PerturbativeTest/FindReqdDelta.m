function FindReqdDelta(tol,coefficient,decimal_place,input_choice,minDelta,maxDelta,number_of_eigenvalues,name_of_quadratization,test_times)
% It finds the Delta values that meet requriements stated, under the conditions provided, and save the output to a text file
%
% Arguments:
%   tol - tolerance of eigs, usually 1e-3 or 1e-6
%   coefficient - coefficient of the tested terms, usually -1 (negative) or 1 (positive)
%   decimal_place - the precision of Delta tested
%                   e.g. decimal_place = 1 and Delta within 1e10-1e13, Delta being tested are 1.0e10, 1.1e10, ... , 9.9e12, 1.0e13
%   input_choice - several options: all_cubics, all_quartics, and a specific 3-local term (for testing purpose)
%   minDelta - lower bound of Delta
%   maxDelta - upper bound of Delta
%   number_of_eigenvalues - number of eigs wanted to consider
%                           1 is to consider the ground only, 2 is to consider the first excited eigs and ground, others is to consider all
%   name_of_quadratization - the gadget used
%   test_times - the time of test, which helps identity different tests
%
% Only the product of Pauli matrices x y z can be tested (no other operations allowed)
% And all operators should act on different qubits
%
% Based on info from Book Of Quardratization Perterbative Gadgets Section
%
% e.g.
%        FindReqdDelta(1e-3,-1,3,'xyz',1e10,1e13,1,'P(3->2)DC2',1)
%        refers to find the Delta required within 1e10-1e13 (3 decimal places) for -x1*y2*z3 using DC2 for the 1st time

  if strcmp(input_choice,"all_cubics")   % test all 3^3 cubic terms, i.e. xxx, xxy, xxz, etc.
    n = 3;                                   % the number of logical operators
    n_combination = 1;
    Op{1} = 'x'; Op{2} = 'y'; Op{3} = 'z';   % Op{:} is the array of operators

    % define the filename with conditions provided so that readers can understand
    if coefficient < 0                    % testing negative terms
      if number_of_eigenvalues == 1       % only consider ground eigs
        FileName = strcat('ground','_',input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(2,3^n);
      elseif number_of_eigenvalues == 2   % consider ground and first excited eigs
        FileName = strcat('1st_excited','_',input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(4,3^n);
      else                                % consider all eigs
        FileName = strcat(input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(7,3^n);
      end
    else                                  % testing positive terms
      if number_of_eigenvalues == 1
        FileName = strcat('ground','_',input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(2,3^n);
      elseif number_of_eigenvalues == 2
        FileName = strcat('1st_excited','_',input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(4,3^n);
      else
        FileName = strcat(input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(7,3^n);
      end
    end

    % 3^3 loops for different terms in total
    for s1 = 1:3
    for s2 = 1:3
    for s3 = 1:3
    operators = strcat(Op{s1},Op{s2},Op{s3});

    % calculate S1,S2,S3 outside of the function lhsrhs to save time
    x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
    S = cell(n);
    % m is the number of auxiliary qubits and NeededM contains matrices which are not affected by Delta
    [m,NeededM] = GetAuxNum(size(operators,2),name_of_quadratization);
    for ind = 1:n
        if operators(ind) == 'x'
            S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(n+m-ind)));
        elseif operators(ind) == 'y'
            S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(n+m-ind)));
        elseif operators(ind) == 'z'
            S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+m-ind)));
        end
    end
    A = S{1}*S{2};    % A should be the product of the first half of LHS
    B = S{3};         % B should be the product of the second half of LHS
    LHS = coefficient*A*B;
    NeededM{end - 2} = A;
    NeededM{end - 1} = B;
    NeededM{end} = LHS;              % save LHS in the last cell of NeededM

    % compute the analytic eigs of LHS to compare with the result of eig(LHS)
    ANAL_E_LHS = 1;
    for ind = 1:n
        if operators(ind) == 'x'
            ANAL_E_LHS = sort(kron(ANAL_E_LHS,eig(x)));
          elseif operators(ind) == 'y'
            ANAL_E_LHS = sort(kron(ANAL_E_LHS,eig(y)));
          elseif operators(ind) == 'z'
            ANAL_E_LHS = sort(kron(ANAL_E_LHS,eig(z)));
        end
    end
    for ind = 1:m
        ANAL_E_LHS = kron(ANAL_E_LHS,eig(eye(size(x,1))));
    end

    Delta = minDelta;
    checkpoint = minDelta;
    % test Delta with the range and precision provided
    while Delta <= maxDelta
     Delta = Delta + 10^(floor(log10(abs(Delta)))-decimal_place);
     [LHS,RHS] = lhsrhs(coefficient,S,NeededM,Delta,name_of_quadratization);
     if isnan(RHS) == 0
      % m = log2(size(RHS,2)) - n;             % the number of auxiliary qubits but not needed since we defined it in GetAuxNum

      [V_RHS,E_RHS] = eig(RHS);
      [V_LHS,E_LHS] = eig(LHS);
      [E_RHS,index] = sort(diag(E_RHS));
      V_RHS = V_RHS(:,index);                  % ensure eigenvalues are matching with corresponding eigenvectors
      [E_LHS,index] = sort(diag(E_LHS));
      V_LHS = V_LHS(:,index);
      [ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );        % indices of LHS and RHS where eigenvalues match within tol

      % check the difference between analytic eigenvalues and results of eig()
      diff_E_LHS  = abs(ANAL_E_LHS - E_LHS);
      if isequal(diff_E_LHS,zeros(size(E_LHS))) == 0
        diff_E_LHS          % if analytic eigenvalues don't match with the resultd of eig, display the difference
      end

      if isempty(ind_evals_L) == 0        % matching eigenvalues exist
        L = V_LHS(:,ind_evals_L);                                                   % L is LHS eigenvectors whose eigenvalues matching with RHS
        R = V_RHS(:,ind_evals_R);
        ind_evecs = find( sqrt(sum( (abs(L)-abs(R)).^2 ) ) < tol);                  % indices of L and R where eigenvectors match within tol
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
        if (number_of_eigenvalues >= 1) && (sum( abs(E_LHS(1:2^m) - E_RHS(1)) < tol ) == 2^m)
            if (delta_required(1,n_combination) == 0)
                delta_required(1,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let ground energy match
            end
            if (isempty(ind_evecs) == 0) && (delta_required(2,n_combination) == 0)
                delta_required(2,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let ground state match
            end

            % consider below when NOT focusing on the ground eigs
            if (number_of_eigenvalues >= 2) && (sum( abs(E_LHS((2^m + 1):2^(m+1)) - E_RHS(2)) < tol ) == 2^m)
                if (delta_required(3,n_combination) == 0)
                    delta_required(3,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let first excited energy match
                end

                if (size(sorted_L,2) >= 2) && (delta_required(4,n_combination) == 0)
                    delta_required(4,n_combination) = floor(log10(abs(Delta)));   % value of Delta that let first excited state match
                end
            end

            % consider below when NOT focusing on the ground and first excited eigs
            if (number_of_eigenvalues > 2) && (numel(unique(ind_evals_R)) == 2^n)
                if (delta_required(5,n_combination) == 0)
                    delta_required(5,n_combination) = floor(log10(abs(Delta)));   % value of Delta that let all energies match
                end

                if isequal(size(sorted_L,2),size(sorted_R,2),2^n) && (delta_required(6,n_combination) == 0)
                    delta_required(6,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let all states match
                end
            end
        end
        if (number_of_eigenvalues > 2) && ne(delta_required(5,n_combination),0) && (delta_required(7,n_combination) == 0) && (numel(unique(ind_evals_R)) == 0)
            delta_required(7,n_combination) = floor(log10(abs(Delta)));  % value of Delta that is too large to keep energies matching
        end
      end
     end
     if (Delta - checkpoint) >= 10^(floor(log10(abs(Delta)))-(decimal_place - 2))    % update the output file and the checkpoint
        dlmwrite(FileName,delta_required','delimiter','\t');   % the delta required
        dlmwrite(FileName,' ','-append');       % separate line
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
    dlmwrite(FileName,delta_required','delimiter','\t');
    dlmwrite(FileName,' ','-append');                             % separate line
    dlmwrite(FileName,floor(log10(abs(minDelta))),'-append');
    dlmwrite(FileName,floor(log10(abs(maxDelta))),'-append');     % range of Delta tested

  elseif strcmp(input_choice,'all_quartics')    % test all 3^4 quartic terms, i.e. xxxx, xxxy, xxxz, etc.
    n = 4;                                   % the number of logical operators
    n_combination = 1;
    Op{1} = 'x'; Op{2} = 'y'; Op{3} = 'z';   % Op{:} is the array of operators

    % define the filename with conditions provided so that readers can understand
    if coefficient < 0                    % testing negative terms
      if number_of_eigenvalues == 1       % only consider ground eigs
        FileName = strcat('ground','_',input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(2,3^n);
      elseif number_of_eigenvalues == 2   % consider ground and first excited eigs
        FileName = strcat('1st_excited','_',input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(4,3^n);
      else                                % consider all eigs
        FileName = strcat(input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(7,3^n);
      end
    else                                  % testing positive terms
      if number_of_eigenvalues == 1
        FileName = strcat('ground','_',input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(2,3^n);
      elseif number_of_eigenvalues == 2
        FileName = strcat('1st_excited','_',input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(4,3^n);
      else
        FileName = strcat(input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(7,3^n);
      end
    end

    % 3^4 loops for different terms in total
    for s1 = 1:3
    for s2 = 1:3
    for s3 = 1:3
    for s4 = 1:3
    operators = strcat(Op{s1},Op{s2},Op{s3},Op{s4});

    % calculate S1,S2,S3,S4 outside of the function lhsrhs to save time
    x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
    S = cell(n);
    % m is the number of auxiliary qubits and NeededM contains matrices which are not affected by Delta
    [m,NeededM] = GetAuxNum(size(operators,2),name_of_quadratization);
    for ind = 1:n
        if operators(ind) == 'x'
            S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(n+m-ind)));
        elseif operators(ind) == 'y'
            S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(n+m-ind)));
        elseif operators(ind) == 'z'
            S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+m-ind)));
        end
    end
    A = S{1}*S{2};         % A should be the product of the first half of LHS
    B = S{3}*S{4};         % B should be the product of the second half of LHS
    LHS = coefficient*A*B;
    NeededM{end - 2} = A;
    NeededM{end - 1} = B;
    NeededM{end} = LHS;              % save LHS in the last cell of NeededM

    % compute the analytic eigs of LHS to compare with the result of eig(LHS)
    ANAL_E_LHS = 1;
    for ind = 1:n
        if operators(ind) == 'x'
            ANAL_E_LHS = sort(kron(ANAL_E_LHS,eig(x)));
          elseif operators(ind) == 'y'
            ANAL_E_LHS = sort(kron(ANAL_E_LHS,eig(y)));
          elseif operators(ind) == 'z'
            ANAL_E_LHS = sort(kron(ANAL_E_LHS,eig(z)));
        end
    end
    for ind = 1:m
        ANAL_E_LHS = kron(ANAL_E_LHS,eig(eye(size(x,1))));
    end

    Delta = minDelta;
    checkpoint = minDelta;
    % test Delta with the range and precision provided
    while Delta <= maxDelta
     Delta = Delta + 10^(floor(log10(abs(Delta)))-decimal_place);
     [LHS,RHS] = lhsrhs(coefficient,S,NeededM,Delta,name_of_quadratization);
     if isnan(RHS) == 0
      % m = log2(size(RHS,2)) - n;             % the number of auxiliary qubits but not needed since we defined it in GetAuxNum

      [V_RHS,E_RHS] = eig(RHS);
      [V_LHS,E_LHS] = eig(LHS);
      [E_RHS,index] = sort(diag(E_RHS));
      V_RHS = V_RHS(:,index);                  % ensure eigenvalues are matching with corresponding eigenvectors
      [E_LHS,index] = sort(diag(E_LHS));
      V_LHS = V_LHS(:,index);
      [ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );        % indices of LHS and RHS where eigenvalues match within tol

      % check the difference between analytic eigenvalues and results of eig()
      diff_E_LHS  = abs(ANAL_E_LHS - E_LHS);
      if isequal(diff_E_LHS,zeros(size(E_LHS))) == 0
        diff_E_LHS          % if analytic eigenvalues don't match with the resultd of eig, display the difference
      end

      if isempty(ind_evals_L) == 0        % matching eigenvalues exist
        L = V_LHS(:,ind_evals_L);                                                   % L is LHS eigenvectors whose eigenvalues matching with RHS
        R = V_RHS(:,ind_evals_R);
        ind_evecs = find( sqrt(sum( (abs(L)-abs(R)).^2 ) ) < tol);                  % indices of L and R where eigenvectors match within tol
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
        if (number_of_eigenvalues >= 1) && (sum( abs(E_LHS(1:2^m) - E_RHS(1)) < tol ) == 2^m)
            if (delta_required(1,n_combination) == 0)
                delta_required(1,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let ground energy match
            end
            if (isempty(ind_evecs) == 0) && (delta_required(2,n_combination) == 0)
                delta_required(2,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let ground state match
            end

            % consider below when NOT focusing on the ground eigs
            if (number_of_eigenvalues >= 2) && (sum( abs(E_LHS((2^m + 1):2^(m+1)) - E_RHS(2)) < tol ) == 2^m)
                if (delta_required(3,n_combination) == 0)
                    delta_required(3,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let first excited energy match
                end

                if (size(sorted_L,2) >= 2) && (delta_required(4,n_combination) == 0)
                    delta_required(4,n_combination) = floor(log10(abs(Delta)));   % value of Delta that let first excited state match
                end
            end

            % consider below when NOT focusing on the ground and first excited eigs
            if (number_of_eigenvalues > 2) && (numel(unique(ind_evals_R)) == 2^n)
                if (delta_required(5,n_combination) == 0)
                    delta_required(5,n_combination) = floor(log10(abs(Delta)));   % value of Delta that let all energies match
                end

                if isequal(size(sorted_L,2),size(sorted_R,2),2^n) && (delta_required(6,n_combination) == 0)
                    delta_required(6,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let all states match
                end
            end
        end
        if (number_of_eigenvalues > 2) && ne(delta_required(5,n_combination),0) && (delta_required(7,n_combination) == 0) && (numel(unique(ind_evals_R)) == 0)
            delta_required(7,n_combination) = floor(log10(abs(Delta)));  % value of Delta that is too large to keep energies matching
        end
      end
     end
     if (Delta - checkpoint) >= 10^(floor(log10(abs(Delta)))-(decimal_place - 2))    % update the output file and the checkpoint
        dlmwrite(FileName,delta_required','delimiter','\t');   % the delta required
        dlmwrite(FileName,' ','-append');       % separate line
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
    dlmwrite(FileName,delta_required','delimiter','\t');
    dlmwrite(FileName,' ','-append');                             % separate line
    dlmwrite(FileName,floor(log10(abs(minDelta))),'-append');
    dlmwrite(FileName,floor(log10(abs(maxDelta))),'-append');     % range of Delta tested

  elseif length(input_choice) == 3           % test a specific 3-local term
    n = 3;                                   % the number of logical operators
    n_combination = 1;
    operators = input_choice;                % to match the description in 'all_cubics' case

    % define the filename with conditions provided so that readers can understand
    if coefficient < 0                    % testing negative terms
      if number_of_eigenvalues == 1       % only consider ground eigs
        FileName = strcat('ground','_',input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(2,3^n);
      elseif number_of_eigenvalues == 2   % consider ground and first excited eigs
        FileName = strcat('1st_excited','_',input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(4,3^n);
      else                                % consider all eigs
        FileName = strcat(input_choice,'_','negative','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(7,3^n);
      end
    else                                  % testing positive terms
      if number_of_eigenvalues == 1
        FileName = strcat('ground','_',input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(2,3^n);
      elseif number_of_eigenvalues == 2
        FileName = strcat('1st_excited','_',input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(4,3^n);
      else
        FileName = strcat(input_choice,'_','positive','_',name_of_quadratization,'_',num2str(tol,'%1.0e'),'_',num2str(test_times),'.txt');
        delta_required = zeros(7,3^n);
      end
    end

    % calculate S1,S2,S3 outside of the function lhsrhs to save time
    x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
    S = cell(n);
    % m is the number of auxiliary qubits and NeededM contains matrices which are not affected by Delta
    [m,NeededM] = GetAuxNum(size(operators,2),name_of_quadratization);
    for ind = 1:n
        if operators(ind) == 'x'
            S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(n+m-ind)));
        elseif operators(ind) == 'y'
            S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(n+m-ind)));
        elseif operators(ind) == 'z'
            S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+m-ind)));
        end
    end
    A = S{1}*S{2};    % A should be the product of the first half of LHS
    B = S{3};         % B should be the product of the second half of LHS
    LHS = coefficient*A*B;
    NeededM{end - 2} = A;
    NeededM{end - 1} = B;
    NeededM{end} = LHS;              % save LHS in the last cell of NeededM

    % compute the analytic eigs of LHS to compare with the result of eig(LHS)
    ANAL_E_LHS = 1;
    for ind = 1:n
        if operators(ind) == 'x'
            ANAL_E_LHS = sort(kron(ANAL_E_LHS,eig(x)));
          elseif operators(ind) == 'y'
            ANAL_E_LHS = sort(kron(ANAL_E_LHS,eig(y)));
          elseif operators(ind) == 'z'
            ANAL_E_LHS = sort(kron(ANAL_E_LHS,eig(z)));
        end
    end
    for ind = 1:m
        ANAL_E_LHS = kron(ANAL_E_LHS,eig(eye(size(x,1))));
    end

    Delta = minDelta;
    checkpoint = minDelta;
    % test Delta with the range and precision provided
    while Delta <= maxDelta
     Delta = Delta + 10^(floor(log10(abs(Delta)))-decimal_place);
     [LHS,RHS] = lhsrhs(coefficient,S,NeededM,Delta,name_of_quadratization);
     if isnan(RHS) == 0
      % m = log2(size(RHS,2)) - n;             % the number of auxiliary qubits but not needed since we defined it in GetAuxNum

      [V_RHS,E_RHS] = eig(RHS);
      [V_LHS,E_LHS] = eig(LHS);
      [E_RHS,index] = sort(diag(E_RHS));
      V_RHS = V_RHS(:,index);                  % ensure eigenvalues are matching with corresponding eigenvectors
      [E_LHS,index] = sort(diag(E_LHS));
      V_LHS = V_LHS(:,index);
      [ind_evals_L, ind_evals_R] = find( abs(E_RHS'-E_LHS) < tol );        % indices of LHS and RHS where eigenvalues match within tol

      % check the difference between analytic eigenvalues and results of eig()
      diff_E_LHS  = abs(ANAL_E_LHS - E_LHS);
      if isequal(diff_E_LHS,zeros(size(E_LHS))) == 0
        diff_E_LHS          % if analytic eigenvalues don't match with the resultd of eig, display the difference
      end

      if isempty(ind_evals_L) == 0        % matching eigenvalues exist
        L = V_LHS(:,ind_evals_L);                                                   % L is LHS eigenvectors whose eigenvalues matching with RHS
        R = V_RHS(:,ind_evals_R);
        ind_evecs = find( sqrt(sum( (abs(L)-abs(R)).^2 ) ) < tol);                  % indices of L and R where eigenvectors match within tol
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
        if (number_of_eigenvalues >= 1) && (sum( abs(E_LHS(1:2^m) - E_RHS(1)) < tol ) == 2^m)
            if (delta_required(1,n_combination) == 0)
                delta_required(1,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let ground energy match
            end
            if (isempty(ind_evecs) == 0) && (delta_required(2,n_combination) == 0)
                delta_required(2,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let ground state match
            end

            % consider below when NOT focusing on the ground eigs
            if (number_of_eigenvalues >= 2) && (sum( abs(E_LHS((2^m + 1):2^(m+1)) - E_RHS(2)) < tol ) == 2^m)
                if (delta_required(3,n_combination) == 0)
                    delta_required(3,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let first excited energy match
                end

                if (size(sorted_L,2) >= 2) && (delta_required(4,n_combination) == 0)
                    delta_required(4,n_combination) = floor(log10(abs(Delta)));   % value of Delta that let first excited state match
                end
            end

            % consider below when NOT focusing on the ground and first excited eigs
            if (number_of_eigenvalues > 2) && (numel(unique(ind_evals_R)) == 2^n)
                if (delta_required(5,n_combination) == 0)
                    delta_required(5,n_combination) = floor(log10(abs(Delta)));   % value of Delta that let all energies match
                end

                if isequal(size(sorted_L,2),size(sorted_R,2),2^n) && (delta_required(6,n_combination) == 0)
                    delta_required(6,n_combination) = floor(log10(abs(Delta)));  % value of Delta that let all states match
                end
            end
        end
        if (number_of_eigenvalues > 2) && ne(delta_required(5,n_combination),0) && (delta_required(7,n_combination) == 0) && (numel(unique(ind_evals_R)) == 0)
            delta_required(7,n_combination) = floor(log10(abs(Delta)));  % value of Delta that is too large to keep energies matching
        end
      end
     end
     if (Delta - checkpoint) >= 10^(floor(log10(abs(Delta)))-(decimal_place - 2))    % update the output file and the checkpoint
        dlmwrite(FileName,delta_required','delimiter','\t');   % the delta required
        dlmwrite(FileName,' ','-append');       % separate line
        dlmwrite(FileName,Delta,'-append');     % largest Delta tested
        checkpoint = Delta;
     end
    end      % use end in MATLAB / endwhile in Octave
    ind = find(delta_required(:,n_combination) == 0);
    delta_required(ind,n_combination) = Inf;
    dlmwrite(FileName,delta_required','delimiter','\t');
    dlmwrite(FileName,' ','-append');                             % separate line
    dlmwrite(FileName,floor(log10(abs(minDelta))),'-append');
    dlmwrite(FileName,floor(log10(abs(maxDelta))),'-append');     % range of Delta tested

  else
    disp("Invalid Input");
  end

end

function [number_of_auxiliary,NeededM] = GetAuxNum(number_of_logical,name_of_quadratization)
% Given the name of quadratization and the number of logical qubits
% Get the number of auxiliary qubits and necessary matrices which are not affected by Delta
%
% Known from Book About Quadratization
% Need to be updated whenever lhsrhs.m gets updated
%
% Notice that NeededM is always defined with 3 empty cells at the end, where A,B,LHS will be

    x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
    n = number_of_logical;

    if strcmp(name_of_quadratization, 'P(3->2)-DC2') || strcmp(name_of_quadratization, 'P(3->2)DC2')
        number_of_auxiliary = 1;
        xa = kron(eye(8),x); za = kron(eye(8),z); I = eye(2^4);
        NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'P(3->2)-KKR') || strcmp(name_of_quadratization, 'P(3->2)KKR')
        number_of_auxiliary = 3;
        xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
        za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);
        I = eye(2^6);
        NeededM = cell(10,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa1; NeededM{2} = xa2; NeededM{3} = xa3;
        NeededM{4} = za1; NeededM{5} = za2; NeededM{6} = za3; NeededM{7} = I;

    elseif strcmp(name_of_quadratization, 'P(3->2)KKR-A') % no coefficient needed
        number_of_auxiliary = 3;
        xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
        za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);
        I = eye(2^6);
        NeededM = cell(10,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa1; NeededM{2} = xa2; NeededM{3} = xa3;
        NeededM{4} = za1; NeededM{5} = za2; NeededM{6} = za3; NeededM{7} = I;

    elseif strcmp(name_of_quadratization, 'P(3->2)-DC1') || strcmp(name_of_quadratization, 'P(3->2)DC1')
        number_of_auxiliary = 3;
        xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
        za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);
        I = eye(2^6);
        NeededM = cell(10,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
        NeededM{1} = xa1; NeededM{2} = xa2; NeededM{3} = xa3;
        NeededM{4} = za1; NeededM{5} = za2; NeededM{6} = za3; NeededM{7} = I;

    elseif strcmp(name_of_quadratization, 'ZZZ-TI-CBBK')
        number_of_auxiliary = 1;
        xa = kron(eye(8),x); za = kron(eye(8),z); I = eye(2^4);
        NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
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
         NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
         NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'P(3->2)-OT') || strcmp(name_of_quadratization, 'P(3->2)OT')
        number_of_auxiliary = 1;
        xa = kron(eye(8),x); za = kron(eye(8),z); I = eye(16);
        NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
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

    elseif strcmp(name_of_quadratization, 'PD-JF')
      number_of_auxiliary = n;
      for ind = 1:n
          ZA{ind} = kron(kron(eye(2^(n-1+ind)),z),eye(2^(n-ind)));
          XA{ind} = kron(kron(eye(2^(n-1+ind)),x),eye(2^(n-ind)));
      end
      I = eye(2^(2*n));

      index = 0;
      for j = 2:n
          for i = 1:j - 1
              index = index + 1;
              Zcouple{index} = ZA{i}*ZA{j};
          end
      end

      Z_total = index*I;
      for k = 1:index
          Z_total = Z_total - Zcouple{k};
      end

      NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
      NeededM{1} = Z_total; NeededM{2} = XA; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'P(3->2)-CBBK2') || strcmp(name_of_quadratization, 'P(3->2)CBBK2')
      number_of_auxiliary = 1;

      za = kron(eye(8),z);
      xa = kron(eye(8),x);
      I = eye(16);
      NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
      NeededM{1} = xa; NeededM{2} = za; NeededM{3} = I;

    elseif strcmp(name_of_quadratization, 'PD-CK')
      number_of_auxiliary = n;

      for ind = 1:n
          ZA{ind} = kron(kron(eye(2^(n-1+ind)),z),eye(2^(n-ind)));
          XA{ind} = kron(kron(eye(2^(n-1+ind)),x),eye(2^(n-ind)));
      end
      I = eye(2^(2*n));
      index = 0;
      for j = 2:n
          for i = 1:j - 1
              index = index + 1;
              Zcouple{index} = ZA{i}*ZA{j};
          end
      end
      Z_total = index*I;
      for k = 1:index
          Z_total = Z_total - Zcouple{k};
      end

      NeededM = cell(6,1);                                    % contains the auxiliary matrices, identity matrix needed and LHS in the last cell
      NeededM{1} = Z_total; NeededM{2} = XA; NeededM{3} = I;

    else
        disp('cannot find this method');
        number_of_auxiliary = nan;
    end
end

function index = UniqueRows(A)
% This MATLAB function provides a faster version of MATLAB's unique rows method (i.e. 'unique(points,''rows'')')
% [index] = UniqueRows(A); % equivalent is '[~,ind]=unique(A,'rows');'
% A_Unique = A(sort(ind),:);
% faster unique function (https://www.mathworks.com/matlabcentral/fileexchange/77329-fast-unique-rows-method-with-unit-and-performance-tests)

[A_sorted,idx1] = sortrows(A);
k = find([true; any(diff(A_sorted,1,1),2); true]);
idx2 = k(diff(k) >= 1);
index=idx1(idx2);

end
