function [LHS, RHS] = lhs2rhs(n,s1,s2,s3,str) 
% test DC1, DC2, and KKR with only three logical qubits

% n refers to the number of logical qubits, but it is set to 3 for
% convenience
% s_i refers to operator acting on the i_th qubit, it is chosen from x,y,z

x = [0 1 ; 1 0]; z = [1 0 ; 0 -1];

if strcmp(str, 'P(3->2)-DC2')
    n = 3;
    s1 = kron(kron(eye(2^0),s1),eye(2^3));
    s2 = kron(kron(eye(2^1),s2),eye(2^2));
    s3 = kron(kron(eye(2^2),s3),eye(2^1));
    xa = kron(kron(eye(2^3),x),eye(2^0));
    za = kron(kron(eye(2^3),z),eye(2^0));
  
    delta = 1e10;
    alpha = (1/2)*delta;
    alpha_s = ((1/4)*(delta^(2/3)) - 1);
    alpha_z = (1/2)*delta;
    alpha_ss = delta^(1/3);
    alpha_sz = (1/4)*(delta^(2/3));
    alpha_sx = delta^(2/3);
    
    LHS = s1*s2*s3;
    RHS = alpha*eye(2^(n+1)) + alpha_s*s3 + alpha_z*za + alpha_ss*((s1 + s2)^2) + alpha_sx*(s1*xa + s2*xa) + alpha_sz*s3*za;

elseif strcmp(str, 'P(3->2)-DC1')
    n = 3;
    s1 = kron(kron(eye(2^0),s1),eye(2^5));
    s2 = kron(kron(eye(2^1),s2),eye(2^4));
    s3 = kron(kron(eye(2^2),s3),eye(2^3));
    xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
    za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);

    delta = 1e10;
    alpha = (1/8)*delta;
    alpha_ss = (1/6)*(delta)^(1/3);
    alpha_sx = (-1/6)*(delta)^(2/3);
    alpha_zz = (-1/24)*delta;

    LHS = s1*s2*s3;
    RHS = alpha*eye(2^(n+3)) + alpha_ss*(s1^2 + s2^2 + s3^2) + alpha_sx*s1*xa1 + alpha_sx*s2*xa2 + alpha_sx*s3*xa3 + alpha_zz*(za1*za2 + za1*za3 + za2*za3);
    
elseif strcmp(str, 'P(3->2)-KKR')
    n = 3;
    s1 = kron(kron(eye(2^0),s1),eye(2^5));
    s2 = kron(kron(eye(2^1),s2),eye(2^4));
    s3 = kron(kron(eye(2^2),s3),eye(2^3));
    xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
    za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);
    
    delta = 1e10;
    alpha = -(1/8)*(delta);    % or they should be the same as DC1 (need to be tested)
    alpha_ss = -(1/6)*(delta)^(1/3);
    alpha_sx = (1/6)*(delta)^(2/3);
    alpha_zz = (1/24)*(delta);

    LHS = s1*s2*s3;
    RHS = alpha*eye(2^(n+3)) + alpha_ss*(s1^2 + s2^2 + s3^2) + alpha_sx*(s1*xa1 + s2*xa2 + s3*xa3) + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

else
    LHS = 0; RHS = 0;
    disp('wrong method chosen');
end

end