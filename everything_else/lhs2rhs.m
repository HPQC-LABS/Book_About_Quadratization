function [LHS, RHS] = lhs2rhs(operators, delta, name_of_quadratization)
% test DC1, DC2, KKR, and ZZZ-TI-CBBK
% operators shold be in the form of 'xyz'
%
% e.g.    [LHS, RHS] = lhs2rhs('xyz',10,'P(3->2)-DC2')
%         refers to quadratize x1*y2*z3 using P(3->2)-DC2 with delta = 1e10

    n = strlength(operators);
    delta = 10^(delta);
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

        alpha = (1/2)*delta;
        alpha_s = ((1/4)*(delta^(2/3)) - 1);
        alpha_z = (1/2)*delta;
        alpha_ss = delta^(1/3);
        alpha_sz = (1/4)*(delta^(2/3));
        alpha_sx = delta^(2/3);

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

        alpha = -(1/8)*(delta);    % or they should be the same as DC1 (need to be tested)
        alpha_ss = -(1/6)*(delta)^(1/3);
        alpha_sx = (1/6)*(delta)^(2/3);
        alpha_zz = (1/24)*(delta);

        LHS = S{1}*S{2}*S{3};
        RHS = alpha*eye(2^(n+3)) + alpha_ss*(S{1}^2 + S{2}^2 + S{3}^2) + alpha_sx*(S{1}*xa1 + S{2}*xa2 + S{3}*xa3) + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

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

        alpha = (1/8)*delta;
        alpha_ss = (1/6)*(delta)^(1/3);
        alpha_sx = (-1/6)*(delta)^(2/3);
        alpha_zz = (-1/24)*delta;

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
        alpha_I = (1/2)*(delta + ((alpha/6)^(2/5))*(delta^(3/5)) + 6*((alpha/6)^(4/5))*(delta^(1/5)) );
        alpha_zi = (-1/2)*(( ((7/6)*alpha) + ( ((alpha/6)^(3/5))*(delta^(2/5))) ) - ( (alpha/6)*(delta^4) )^(1/5));
        alpha_zj = alpha_zi;
        alpha_zk = alpha_zi;
        alpha_za = (-1/2)*( delta - (((alpha/6)^(2/5))*(delta^(3/5))) );
        alpha_xa = ( (alpha/6)*(delta^4) )^(1/5);
        alpha_zzia = (-1/2)*( ( (7/6)*alpha + ((alpha/6)^(3/5))*(delta^(2/5)) ) + ( (alpha/6)*(delta^4) )^(1/5) );
        alpha_zzja = alpha_zzia;
        alpha_zzka = alpha_zzja;
        alpha_zzij = 2*((alpha/6)^(4/5))*((delta)^(1/5));
        alpha_zzik = alpha_zzij;
        alpha_zzjk = alpha_zzij;
        
        LHS = S{1}*S{2}*S{3};
        RHS = alpha_I*eye(16) + alpha_zi*S{1} + alpha_zj*S{2} + alpha_zk*S{3} + alpha_za*za + alpha_xa*xa + alpha_zzia*S{1}*za + alpha_zzja*S{2}*za + alpha_zzka*S{3}*za + alpha_zzij*S{1}*S{2} + alpha_zzik*S{1}*S{3} + alpha_zzjk*S{2}*S{3};

    else
        disp('cannot find this method');
        LHS = []; RHS = [];
    end
end
