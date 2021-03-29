function [LHS, RHS] = lhsrhs(coefficient,S,NeededM,Delta,name_of_quadratization)
% This function should only be called from FindReqdDelta.m
% test P(3->2)-DC1, P-(3->2)DC2, P-(3->2)KKR, P(3->2)-OT, P(3->2)-CBBK, P(3->2)CBBK2, ZZZ-TI-CBBK,
% P1B1-OT, P1B1-CBBK, PSD-OT, PSD-CBBK, PSD-CN, PD-JF, and PD-CK
% S needs to be defined before using this function, and it should be a cell array in which S{1} is a matrix of an operator
% NeededM needs to be defined before using this function, and it should be a cell array that contains matrices which are not affected by Delta value
%
% e.g.    S = {S1; S2; S3}; NeededM = {xa;za;LHS};
%         [LHS, RHS] = lhsrhs(-1,S,NeededM,1e10,'P(3->2)-DC2')
%         refers to quadratize -S1*S2*S3 using P(3->2)-DC2 with Delta = 1e10


    n = size(S,1);      % number of operators

    if strcmp(name_of_quadratization, 'P(3->2)-DC2') || strcmp(name_of_quadratization, 'P(3->2)DC2')
        assert(n == 3, 'P(3->2)-DC2 requires a 3-local term, please only give 3 operators.');

        xa = NeededM{1}; za = NeededM{2}; I = NeededM{3};
        S1 = S{1}; S2 = S{2}; S3 = S{3};

        alpha = (1/2)*Delta;
        alpha_s = coefficient*((1/4)*(Delta^(2/3)) - 1);
        alpha_z = (1/2)*Delta;
        alpha_ss = Delta^(1/3);
        alpha_sz = coefficient*(1/4)*(Delta^(2/3));
        alpha_sx = Delta^(2/3);

        LHS = NeededM{end};
        RHS1 = alpha*I;
        RHS2 = alpha_s*S3;
        RHS3 = alpha_sx*xa*S1;
        RHS4 = alpha_sx*xa*S2;
        RHS5 = alpha_sz*za*S3;
        RHS6 = alpha_z*za;
        RHS7 = 2*alpha_ss*I;
        RHS8 = 2*alpha_ss*S1*S2;        % we only use x,y,z here, for which S^2 = I
        RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS5 + RHS6 + RHS7 + RHS8;
        RHS = alpha*I + alpha_s*S3 + alpha_sx*xa*S1 + alpha_sx*xa*S2 + alpha_sz*za*S3 + alpha_z*za + 2*alpha_ss*(I + S1*S2);

    elseif strcmp(name_of_quadratization, 'P(3->2)-KKR') || strcmp(name_of_quadratization, 'P(3->2)KKR')
        assert(n == 3, 'P(3->2)-KKR requires a 3-local term, please only give 3 operators.');

        xa1 = NeededM{1}; xa2 = NeededM{2}; xa3 = NeededM{3};
        za1 = NeededM{4}; za2 = NeededM{5}; za3 = NeededM{6}; I = NeededM{7};

        alpha = -(1/8)*(Delta);
        alpha_ss = -(1/6)*(Delta)^(1/3);
        alpha_sx = (1/6)*(Delta)^(2/3);
        alpha_zz = (1/24)*(Delta);
 %      have not add coefficient
        LHS = NeededM{end};
        RHS = alpha*I + alpha_ss*(S1^2 + S2^2 + S3^2) + alpha_sx*(S1*xa1 + S2*xa2 + S3*xa3) + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

    elseif strcmp(name_of_quadratization, 'P(3->2)KKR-A') % no coefficient needed
        assert(n == 3, 'P(3->2)KKR-A requires a 3-local term, please only give 3 operators.');
        assert(isequal(coefficient,1), 'Cannot change the coefficient of P(3->2)KKR-A');

        xa1 = NeededM{1}; xa2 = NeededM{2}; xa3 = NeededM{3};
        za1 = NeededM{4}; za2 = NeededM{5}; za3 = NeededM{6}; I = NeededM{7};

        alpha = (3/4)*(Delta);
        alpha_ss = (Delta)^(1/3);
        alpha_sx = -(Delta)^(2/3);
        alpha_zz = -(1/4)*(Delta);

        LHS = -6*NeededM{end};
        RHS = alpha*I + 3*alpha_ss*I + alpha_sx*(S1*xa1 + S2*xa2 + S3*xa3) + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

    elseif strcmp(name_of_quadratization, 'P(3->2)-DC1') || strcmp(name_of_quadratization, 'P(3->2)DC1')
        assert(n == 3, 'P(3->2)-DC1 requires a 3-local term, please only give 3 operators.');

        xa1 = NeededM{1}; xa2 = NeededM{2}; xa3 = NeededM{3};
        za1 = NeededM{4}; za2 = NeededM{5}; za3 = NeededM{6}; I = NeededM{7};

        alpha = (1/8)*Delta;
        alpha_ss = (1/6)*(Delta)^(1/3);
        alpha_sx = (-1/6)*(Delta)^(2/3);
        alpha_zz = (-1/24)*Delta;

        LHS = NeededM{end};
        RHS = alpha*eye(2^(n+3)) + alpha_ss*(S1^2 + S2^2 + S3^2) + alpha_sx*(S1*xa1 + S2*xa2 + S3*xa3) + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

   elseif strcmp(name_of_quadratization, 'ZZZ-TI-CBBK')
        assert(n == 3, 'ZZZ-TI-CBBK requires a 3-local term, please only give 3 operators.');

        xa = NeededM{1}; za = NeededM{2}; I = NeededM{3};
        S1 = S{1}; S2 = S{2}; S3 = S{3};

        alpha_I = (1/2)*(Delta + ((coefficient/6)^(2/5))*(Delta^(3/5)) + 6*((coefficient/6)^(4/5))*(Delta^(1/5)) );
        alpha_zi = (-1/2)*(( ((7/6)*coefficient) + ( ((coefficient/6)^(3/5))*(Delta^(2/5))) ) - ( (coefficient/6)*(Delta^4) )^(1/5));
        alpha_zj = alpha_zi;
        alpha_zk = alpha_zi;
        alpha_za = (-1/2)*( Delta - (((coefficient/6)^(2/5))*(Delta^(3/5))) );
        alpha_xa = ( (coefficient/6)*(Delta^4) )^(1/5);
        alpha_zzia = (-1/2)*( ( (7/6)*coefficient + ((coefficient/6)^(3/5))*(Delta^(2/5)) ) + ( (coefficient/6)*(Delta^4) )^(1/5) );
        alpha_zzja = alpha_zzia;
        alpha_zzka = alpha_zzja;
        alpha_zzij = 2*((coefficient/6)^(4/5))*((Delta)^(1/5));
        alpha_zzik = alpha_zzij;
        alpha_zzjk = alpha_zzij;

        LHS = NeededM{end};
        RHS = alpha_I*I + alpha_zi*S1 + alpha_zj*S2 + alpha_zk*S3 + alpha_za*za + alpha_xa*xa + alpha_zzia*S1*za + alpha_zzja*S2*za + alpha_zzka*S3*za + alpha_zzij*S1*S2 + alpha_zzik*S1*S3 + alpha_zzjk*S2*S3;

    elseif strcmp(name_of_quadratization, 'PSD-CBBK')
        assert(n >= 5, 'PSD-CBBK requires at least a 5-local term, please give at least 5 operators.');

        n_a = ceil(n/2);

        xa = NeededM{1}; za = NeededM{2}; I = NeededM{3};
        A = NeededM{end - 2};
        B = NeededM{end - 1};

%        for ind = 1:n_a
%            A = A*S{ind};
%        end
%        for ind = n_a + 1:n
%            B = B*S{ind};
%        end

        LHS = NeededM{end};
        RHS = (Delta)*((1*I - za)/2) + abs(coefficient)*((1*I + za)/2) + sqrt( abs(coefficient)*Delta/2 )*(sign(coefficient)*A - B)*xa;

    elseif strcmp(name_of_quadratization, 'PSD-OT')
        assert(n >= 4, 'PSD-OT requires at least a 4-local term, please give at least 4 operators.');

        n_a = ceil(n/2);

        xa = NeededM{1}; za = NeededM{2}; I = NeededM{3};
        A = NeededM{end - 2};
        B = NeededM{end - 1};
%        A = I; B = I;
%        for ind = 1:n_a
%            A = A*S{ind};
%        end

%        for ind = n_a + 1:n
%            B = B*S{ind};
%        end

        LHS = NeededM{end};
        RHS1 = 1/2*Delta*(1*I);
        RHS2 = 1/2*Delta*(-za);
        RHS31 = (coefficient/2)*(I);
        RHS3 = (coefficient/2)*(A^2);
        RHS41 = (coefficient/2)*(I);
        RHS4 = (coefficient/2)*(B^2);
        RHS5 = B - A;
        RHS6 = sqrt( coefficient*Delta/2 )*(RHS5)*xa;
        RHS = RHS1 + RHS2 + RHS3 + RHS4 + RHS6;

    elseif strcmp(name_of_quadratization, 'P(3->2)CBBK') || strcmp(name_of_quadratization, 'P(3->2)-CBBK')
      assert(n == 3, 'P(3->2)-CBBK requires a 3-local term, please only give 3 operators.');

      xa = NeededM{1}; za = NeededM{2}; I = NeededM{3};
      S1 = S{1}; S2 = S{2}; S3 = S{3};

      alpha = Delta/2 + (1/2)*(coefficient/2)^(2/3)*Delta^(1/2)*( (sign(coefficient)^2) + 1 ) - (sign(coefficient)^2)*((coefficient/2)^(4/3))* ((sign(coefficient)^2) + 1);
      alpha_s3 = (1/2)*(coefficient/2)^(1/3)*(Delta^(1/2)) - (coefficient/4)*( (sign(coefficient)^2) + 1 );
      alpha_za = (-Delta/2) + (1/2)*(coefficient/2)^(2/3)*Delta^(1/2)*( (sign(coefficient)^2) + 1 ) - (sign(coefficient)^2)*((coefficient/2)^(4/3))* ((sign(coefficient)^2) + 1);
      alpha_s1_s2 = 2*sign(coefficient)*((coefficient/2)^(2/3))*(Delta^(1/2)) - 4*sign(coefficient)*((coefficient/2)^(4/3));
      alpha_s3_za = (-1/2)*(coefficient/2)^(1/3)*(Delta^(1/2)) - (coefficient/4)*( (sign(coefficient)^2) + 1 );
      alpha_s1_xa = sign(coefficient)*(coefficient/2)^(1/3)*(Delta^(3/4));
      alpha_s2_xa = (coefficient/2)^(1/3)*(Delta^(3/4));

      LHS = NeededM{end};
      RHS = alpha*I + alpha_s3*S3 + alpha_za*za + alpha_s3_za*S3*za + alpha_s1_xa*S1*xa + alpha_s2_xa*S2*xa + alpha_s1_s2*S1*S2;

    elseif strcmp(name_of_quadratization, 'P(3->2)-OT') || strcmp(name_of_quadratization, 'P(3->2)OT')
      assert(n == 3, 'P(3->2)-OT requires a 3-local term, please only give 3 operators.');

      xa = NeededM{1}; za = NeededM{2}; I = NeededM{3};
      S1 = S{1}; S2 = S{2}; S3 = S{3};

      alpha = (Delta/2);
      alpha_s1 = ((Delta^(1/3))*(coefficient^(2/3))/2);
      alpha_s2 = ((Delta^(1/3))*(coefficient^(2/3))/2);
      alpha_s3 = -((Delta^(2/3))*(coefficient^(1/3))/2);
      alpha_za = -(Delta/2);

      alpha_s1_s2 = -((Delta^(1/3))*(coefficient^(2/3)));
      alpha_s1_s3 = (coefficient/2);
      alpha_s2_s3 = (coefficient/2);

      alpha_s3_za = (Delta^(2/3))*(coefficient^(1/3)/2);

      alpha_s1_xa = -((Delta^(2/3))*(coefficient^(1/3))/sqrt(2));
      alpha_s2_xa = ((Delta^(2/3))*(coefficient^(1/3))/sqrt(2));

      LHS = NeededM{end};
      RHS = alpha*I + alpha_s1*(S1)^2 + alpha_s2*(S2)^2 + alpha_s3*S3 ...
      + alpha_za*za + alpha_s1_s2*S1*S2 + alpha_s1_s3*(S1^2)*S3 + alpha_s2_s3*(S2^2)*S3 + alpha_s3_za*S3*za ...
      + alpha_s1_xa*S1*xa + alpha_s2_xa*S2*xa;

    elseif strcmp(name_of_quadratization, 'P1B1-CBBK')
        assert(n >= 3, 'P1B1-CBBK requires at least a 3-local term, please give at least 3 terms.');

        xa = NeededM{1}; za = NeededM{2}; I = NeededM{3};
        prod_A = I; prod_B = I;

        for k = 1:(n-2)
            prod_A = prod_A*S{k};
        end

        for k = 1:(n-1)
            prod_B = prod_B*S{k};
        end

        LHS = NeededM{end};      % NeededM{end} is the product of S{1}*S{2}*S{3}* ...
        RHS = (Delta*I + ((coefficient/2)^(1/3))*(Delta^(1/2))*S{n})*((1*I - za)/2) ...
            + ((coefficient/2)^(1/3))*(Delta^(3/4))*(sign(coefficient)*prod_A + S{n-1})*xa ...
            + ((coefficient^(2/3))/2)*( (sign(coefficient)^2) + 1 )*((2^(1/3))*(Delta^(1/2))*I - (coefficient^(1/3))*S{n} - (sign(coefficient)^2)*((2*coefficient)^(2/3))*I)*((1*I + za)/2) ...
            + sign(coefficient)*((coefficient^(2/3))*(2^(1/3))*(Delta^(1/2)) - 4*((coefficient/2)^(4/3)))*prod_B;

    elseif strcmp(name_of_quadratization, 'P1B1-OT')
        assert(n >= 3, 'P1B1-CBBK requires at least a 3-local term, please give at least 3 terms.');

        xa = NeededM{1}; za = NeededM{2}; I = NeededM{3};
        prod_A = I;

        for k = 1:(n-2)
            prod_A = prod_A*S{k};
        end

        LHS = NeededM{end};
        RHS = (Delta*I - (Delta^(2/3))*(coefficient^(1/3))*S{n})*((1*I - za)/2) ...
        + ((Delta^(2/3))*(coefficient^(1/3))/sqrt(2))*(-prod_A + S{n-1})*xa ...
        + ((Delta^(1/3))*(coefficient^(2/3))/2)*(-prod_A + S{n-1})^2 + (coefficient/2)*((prod_A)^2 + (S{n-1})^2)*S{n};

    elseif strcmp(name_of_quadratization, 'PSD-CN')
        assert(n >= 3, 'PSD-CN requires at least a 3-local term, please give at least 3 terms.')

        za_11 = NeededM{1}; xa_11 = NeededM{2}; za_1 = NeededM{3}; I = NeededM{4};
        H_1j = NeededM{end - 2}; H_2j = NeededM{end - 1};

%        for k = 1:ceil(n/2)
%           H_1j = H_1j*S{k};
%        end
%        for k = ceil(n/2)+1:n
%            H_2j = H_2j*S{k};
%        end

        R = 1; C = 1;
        beta = sqrt((coefficient*Delta)/(2*R)); alpha = Delta/(2*C);

        LHS = NeededM{end};
        RHS = alpha*( (I - za_11*za_1) ) + alpha*( (I - za_1) ) + alpha*( (I - za_1*za_1) ) ...
          + beta*( (H_1j - H_2j)*xa_11 ) + coefficient*I;

    elseif strcmp(name_of_quadratization, 'PD-JF')
        assert(n >= 3, 'PD-JF requires at least a 3-local term, please give at least 3 terms.');

        Z_total = NeededM{1}; XA = NeededM{2}; I = NeededM{3};

        for ind = 1:n
              Xcouple{ind} = S{ind}*XA{ind};
        end

        Xcouple{1} = Xcouple{1}*coefficient;
        X_total = 0*I;
        for k = 1:n
            X_total = X_total + Xcouple{k};
        end

        LHS = NeededM{end};
        RHS = (1/2)*Z_total + (1/Delta)*(X_total) - coefficient*I;

    elseif strcmp(name_of_quadratization, 'P(3->2)-CBBK2') || strcmp(name_of_quadratization, 'P(3->2)CBBK2')
        assert(n == 3, 'P(3->2)-CBBK2 requires a 3-local term, please only give 3 operators.');

        xa = NeededM{1}; za = NeededM{2}; I = NeededM{3};
        S1 = S{1}; S2 = S{2}; S3 = S{3};

        LHS = NeededM{end};
        RHS = (Delta*I + (( abs(coefficient)/2 )^(1/3))*(Delta^(1/2))*S3)*((I - za)/2) ...
          + ((abs(coefficient)/2)^(1/3))*(Delta^(3/4))*(sign(coefficient)*S1 + S2)*xa ...
          + ((abs(coefficient)/2)^(2/3))*(Delta^(1/2))*(sign(coefficient)*S1 + S2)^2 ...
          - (abs(coefficient)/2)*(sign(coefficient)^2 + 1)*S3 ...
          - ((abs(coefficient)/2)^(4/3))*(sign(coefficient)*S1 + S2)^4;

    elseif strcmp(name_of_quadratization, 'PD-CK')
        assert(n >= 3, 'PD-CK requires at least a 3-local term, please give at least 3 terms.');

        Z_total = NeededM{1}; XA = NeededM{2}; I = NeededM{3};


        for ind = 1:n
            Xcouple{ind} = S{ind}*XA{ind};
        end
        X_total = 0*I;
        for k = 1:n
            X_total = X_total + Xcouple{k};
        end

        LHS = NeededM{end};
        RHS = (Delta/(2*(n - 1)))*(Z_total) + (coefficient^(1/n))*(X_total) - coefficient*I;

    else
        disp('cannot find this method');
        LHS = []; RHS = [];
    end
  end
