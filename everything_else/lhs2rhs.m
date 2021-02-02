function [LHS, RHS] = lhs2rhs(operators, Delta, name_of_quadratization)
% test DC1, DC2, KKR, and ZZZ-TI-CBBK
% operators shold be in the form of 'xyz'
%
% e.g.    [LHS, RHS] = lhs2rhs('xyz',1e10,'P(3->2)-DC2')
%         refers to quadratize x1*y2*z3 using P(3->2)-DC2 with Delta = 1e10

    n = strlength(operators);
    S = cell(n);
    x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

    if strcmp(name_of_quadratization, 'P(3->2)-DC2') || strcmp(name_of_quadratization, 'P(3->2)DC2')
        assert(n == 3, 'P(3->2)-DC2 requires a 3-local term, please only give 3 operators.');
        for ind = 1:n
            if operators(ind) == 'x'
                S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(n+1-ind)));
            elseif operators(ind) == 'y'
                S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(n+1-ind)));
            elseif operators(ind) == 'z'
                S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+1-ind)));
            end
        end
        xa = kron(eye(8),x); za = kron(eye(8),z);

        alpha = (1/2)*Delta;
        alpha_s = ((1/4)*(Delta^(2/3)) - 1);
        alpha_z = (1/2)*Delta;
        alpha_ss = Delta^(1/3);
        alpha_sz = (1/4)*(Delta^(2/3));
        alpha_sx = Delta^(2/3);

        LHS = S{1}*S{2}*S{3};
        RHS = alpha*eye(16) + alpha_s*S{3} + alpha_z*za + alpha_ss*((S{1} + S{2})^2) + alpha_sx*(S{1}*xa + S{2}*xa) + alpha_sz*za*S{3};

    elseif strcmp(name_of_quadratization, 'P(3->2)-KKR') || strcmp(name_of_quadratization, 'P(3->2)KKR')
        assert(n == 3, 'P(3->2)-KKR requires a 3-local term, please only give 3 operators.');
        for ind = 1:n
            if operators(ind) == 'x'
                S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(n+3-ind)));
            elseif operators(ind) == 'y'
                S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(n+3-ind)));
            elseif operators(ind) == 'z'
                S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+3-ind)));
            end
        end
        xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
        za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);

        alpha = -(1/8)*(Delta);
        alpha_ss = -(1/6)*(Delta)^(1/3);
        alpha_sx = (1/6)*(Delta)^(2/3);
        alpha_zz = (1/24)*(Delta);

        LHS = S{1}*S{2}*S{3};
        RHS = alpha*eye(2^(n+3)) + alpha_ss*(S{1}^2 + S{2}^2 + S{3}^2) + alpha_sx*(S{1}*xa1 + S{2}*xa2 + S{3}*xa3) + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

    elseif strcmp(name_of_quadratization, 'P(3->2)KKR-A')
        assert(n == 3, 'P(3->2)KKR-A requires a 3-local term, please only give 3 operators.');
        for ind = 1:n
            if operators(ind) == 'x'
                S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(n+3-ind)));
            elseif operators(ind) == 'y'
                S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(n+3-ind)));
            elseif operators(ind) == 'z'
                S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+3-ind)));
            end
        end
        xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
        za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);

        alpha = (3/4)*(Delta);
        alpha_ss = (Delta)^(1/3);
        alpha_sx = -(Delta)^(2/3);
        alpha_zz = -(1/4)*(Delta);

        LHS = -6*S{1}*S{2}*S{3};
        RHS = alpha*eye(64) + 3*alpha_ss*eye(64) + alpha_sx*(S{1}*xa1 + S{2}*xa2 + S{3}*xa3) + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

    elseif strcmp(name_of_quadratization, 'P(3->2)-DC1') || strcmp(name_of_quadratization, 'P(3->2)DC1')
        assert(n == 3, 'P(3->2)-DC1 requires a 3-local term, please only give 3 operators.');
        for ind = 1:n
            if operators(ind) == 'x'
                S{ind} = kron(kron(eye(2^(ind-1)),x),eye(2^(n+3-ind)));
            elseif operators(ind) == 'y'
                S{ind} = kron(kron(eye(2^(ind-1)),y),eye(2^(n+3-ind)));
            elseif operators(ind) == 'z'
                S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+3-ind)));
            end
        end
        xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
        za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);

        alpha = (1/8)*Delta;
        alpha_ss = (1/6)*(Delta)^(1/3);
        alpha_sx = (-1/6)*(Delta)^(2/3);
        alpha_zz = (-1/24)*Delta;

        LHS = S{1}*S{2}*S{3};
        RHS = alpha*eye(2^(n+3)) + alpha_ss*(S{1}^2 + S{2}^2 + S{3}^2) + alpha_sx*(S{1}*xa1 + S{2}*xa2 + S{3}*xa3) + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

   elseif strcmp(name_of_quadratization, 'ZZZ-TI-CBBK')
        assert(n == 3, 'ZZZ-TI-CBBK requires a 3-local term, please only give 3 operators.');
        for ind = 1:n
            assert(operators(ind) == 'z', 'ZZZ-TI-CBBK requires a ZZZ term.');
            S{ind} = kron(kron(eye(2^(ind-1)),z),eye(2^(n+1-ind)));
        end

        xa = kron(eye(8),x);
        za = kron(eye(8),z);

        alpha = 1;
        alpha_I = (1/2)*(Delta + ((alpha/6)^(2/5))*(Delta^(3/5)) + 6*((alpha/6)^(4/5))*(Delta^(1/5)) );
        alpha_zi = (-1/2)*(( ((7/6)*alpha) + ( ((alpha/6)^(3/5))*(Delta^(2/5))) ) - ( (alpha/6)*(Delta^4) )^(1/5));
        alpha_zj = alpha_zi;
        alpha_zk = alpha_zi;
        alpha_za = (-1/2)*( Delta - (((alpha/6)^(2/5))*(Delta^(3/5))) );
        alpha_xa = ( (alpha/6)*(Delta^4) )^(1/5);
        alpha_zzia = (-1/2)*( ( (7/6)*alpha + ((alpha/6)^(3/5))*(Delta^(2/5)) ) + ( (alpha/6)*(Delta^4) )^(1/5) );
        alpha_zzja = alpha_zzia;
        alpha_zzka = alpha_zzja;
        alpha_zzij = 2*((alpha/6)^(4/5))*((Delta)^(1/5));
        alpha_zzik = alpha_zzij;
        alpha_zzjk = alpha_zzij;

        LHS = S{1}*S{2}*S{3};
        RHS = alpha_I*eye(16) + alpha_zi*S{1} + alpha_zj*S{2} + alpha_zk*S{3} + alpha_za*za + alpha_xa*xa + alpha_zzia*S{1}*za + alpha_zzja*S{2}*za + alpha_zzka*S{3}*za + alpha_zzij*S{1}*S{2} + alpha_zzik*S{1}*S{3} + alpha_zzjk*S{2}*S{3};

    else
        disp('cannot find this method');
        LHS = []; RHS = [];
    end
end

%% Below are some tests to compare eigenvalues and eigenvectors.
%% If want to run it locally, please use Ctrl+T (or Cmd+T) to uncomment it, then use Ctrl+C to copy the test.
%% Dont forget to use Ctrl+/ to comment it out again.

%% test: P(3->2)DC2    also tested ('xxx',1e13,'P(3->2)DC2') ('xzx',1e10,'P(3->2)DC2') successfully
% [LHS RHS] = lhs2rhs('zxx',1e10,'P(3->2)DC2');
% [VR, ER] = eig(RHS); [VL, EL] = eig(LHS);
%
%% Because the sign of eigenvectors does not matter, we set any non-zero elements in the vectors to 1.
% ER = round(ER(:,1:8)); VR = +(round(VR(:,1:8)) ~= 0); VL = +(VL ~= 0);
%
%% It is divided into two parts as the ground state and the first excited state to ensure the eigenvalues (energy) match.
%% Sort VL one by one through finding the eigenvector that is closest to VR.
% for ind_R = 1:4
%   for ind_L = 1:8
%     dist(ind_L,ind_R) = norm(VR(:,ind_R)-VL(:,ind_L));
%   end
%   [dist(size(VL,2)+1,ind_R), ind] = min(dist(1:8,ind_R)); S_VL(:,ind_R) = VL(:,ind); S_EL(:,ind_R) = sum(EL(:,ind));
% end
%
% for ind_R = 5:8
%   for ind_L = 9:16
%      dist(ind_L,ind_R) = norm(VR(:,ind_R)-VL(:,ind_L));
%   end
%   [dist(size(VL,2)+1,ind_R), ind] = min(dist(9:16,ind_R)); S_VL(:,ind_R) = VL(:,ind+8); S_EL(:,ind_R) = sum(EL(:,ind+8));
% end
%
% dist = transpose(dist(size(dist,1),:));
% isequal(S_VL,VR);
% isequal(S_EL,sum(ER));
%% dist returns a zero vector because the eigenvectors are reproduced nicely, and S_EL is the corresponding eigenvalues of S_VL.

%% !! test: KKR-A/P(3->2)DC1    failed every time for any Delta input (dist >= 1)
%% Delta cannot get higher than 1e15, otherwise eigenvalues will not match

% i = 1;
% for Delta = 1:1e14:1e15
%   [LHS RHS] = lhs2rhs('zzz',Delta,'P(3->2)KKR-A');
%   [VR, ER] = eig(RHS); [VL, EL] = eig(LHS);
%   ER = round(ER(:,1:16)); VR = +(round(VR(:,1:16)) ~= 0); VL = +(VL ~= 0);
%   for ind_R = 1:8
%     for ind_L = 1:32
%       dist(ind_L,ind_R) = norm(VR(:,ind_R)-VL(:,ind_L));
%     end
%     [dist(size(VL,2)+1,ind_R), ind] = min(dist(1:32,ind_R)); S_VL(:,ind_R) = VL(:,ind); S_EL(:,ind_R) = sum(EL(:,ind));
%   end
%   for ind_R = 9:16
%     for ind_L = 33:64
%         dist(ind_L,ind_R) = norm(VR(:,ind_R)-VL(:,ind_L));
%     end
%     [dist(size(VL,2)+1,ind_R), ind] = min(dist(33:64,ind_R)); S_VL(:,ind_R) = VL(:,ind+32); S_EL(:,ind_R) = sum(EL(:,ind+32));
%   end
%   total_dist(:,i) = transpose(dist(size(dist,1),:));
%   i = i + 1;
%  end

%% !! failed test: KKR   dist >= 1 (separate because ground energy is in different columns)
% i = 1;
% for Delta = 1:1e14:1e15
%   [LHS RHS] = lhs2rhs('zzz',Delta,'P(3->2)KKR');
%   [VR, ER] = eig(RHS); [VL, EL] = eig(LHS);
%   ER = round(ER(:,49:64)); VR = +(round(VR(:,49:64)) ~= 0); VL = +(VL ~= 0);
%   for ind_R = 1:8
%     for ind_L = 1:32
%        dist(ind_L,ind_R) = norm(VR(:,ind_R)-VL(:,ind_L));
%      end
%     [dist(size(VL,2)+1,ind_R), ind] = min(dist(1:32,ind_R)); S_VL(:,ind_R) = VL(:,ind); S_EL(:,ind_R) = sum(EL(:,ind));
%   end
%   for ind_R = 9:16
%    for ind_L = 33:64
%        dist(ind_L,ind_R) = norm(VR(:,ind_R)-VL(:,ind_L));
%    end
%   [dist(size(VL,2)+1,ind_R), ind] = min(dist(33:64,ind_R)); S_VL(:,ind_R) = VL(:,ind+32); S_EL(:,ind_R) = sum(EL(:,ind+32));
%   end
% total_dist(:,i) = transpose((dist(size(dist,1),:)));
% i = i + 1;
% end
