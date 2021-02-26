%% example of testing a specific combination (test zzz for Delta = 1e12:1e9:1e13 using P(3->2)CBBK, with tolerance set to 1e-03)
delta_required = GetReqdDelta('zzz',1e12,1e9,1e13,'P(3->2)CBBK',1e-03,1);    % the last argument represents the number of auxiliary qubits being used


%% Test results

% DC1: Not tested
% KKR: Not tested
% DC2/OT/CBBK: We have tested all s1s2s3 combinations for all eigenvectors and eigenvalues to match (results in google sheet).


%% test 27 negative combinations at one time

% fundamental settings for testing multiple Delta values, where k is the number of auxiliary qubits
name_of_quadratization = 'P(3->2)DC2'; Delta_lower_bound = 1e12; Delta_increment = 1e10; Delta_upper_bound = 1e13;  tol = 1e-06; k = 1;

combinations = cell(27,1); S{1} = 'x'; S{2} = 'y'; S{3} = 'z'; n_combination = 1;  % fundamental settings for testing multiple combinations
combined_delta_required = zeros(6,27);             % the desired Delta values for different s1s2s3 combinations

% can also set a specific s1s2s3 combination here
for s1 = 1:3
for s2 = 1:3
for s3 = 1:3
    combinations(n_combination) = {['-' S{s1} S{s2} S{s3}]};         % the order of terms being tested
    % get the disired delta values for a certain combination
    delta_required = GetReqdDelta(['-' S{s1} S{s2} S{s3}],Delta_lower_bound,Delta_increment,Delta_upper_bound,name_of_quadratization,tol,k);
    combined_delta_required(:,n_combination) = delta_required;
    delta_required = [];                % clear the array
    n_combination = n_combination + 1;
end
end
end

combinations
combined_delta_required       % gives combinations and corresponding delta_required as outputs
