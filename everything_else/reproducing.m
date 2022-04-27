%% Pg. 1, Eqs 1-2

b= dec2bin(2^4-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);
LHS=b1.*b2 + b2.*b3 + b3.*b4 - 4*b1.*b2.*b3;
RHS=b1.*b2 + b2.*b3 + b3.*b4 + 4*b1 - 4*b1.*b2 - 4*b1.*b3;

%% Pg. 1, Eqs 4-5 (Also used for abstract for 2019 AQC paper).

x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
x1 = kron(x, eye(8)); x2 = kron(kron(eye(2), x), eye(4)); x3 = kron(kron(eye(4), x), eye(2)); x4 = kron(eye(8), x);
y1 = kron(y, eye(8)); y2 = kron(kron(eye(2), y), eye(4)); y3 = kron(kron(eye(4), y), eye(2)); y4 = kron(eye(8), y);
z1 = kron(z, eye(8)); z2 = kron(kron(eye(2), z), eye(4)); z3 = kron(kron(eye(4), z), eye(2)); z4 = kron(eye(8), z);
LHS = x1*y2*z3*y4 + y1*x2*z3*y4 + x1*x2*y3;
RHS = x1*y4 + x2*y4 + x3;
max(eig(LHS)-eig(RHS))<1e-13; % gives 1

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR); % sorted eigenvalues in ascending order and corresponding eigenvectors

for col = 1:1:size(VL,2) % compare eigenvectors and eigenvalues
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

max(E_diff); % 3.1e-15, full energy spectrum reproduced
min(V_diff); % 1.28,      no states reproduced

%% Pg. 5, Deduc-Reduc: Eqs. 6-7.
%% Pg. 6, Eqs 8 and 10

LHS=b1.*b2 + b2.*b3 + b3.*b4 - 4*b1.*b2.*b3;
RHS=b1.*b2 + b2.*b3 + b3.*b4 + 4*b1 - 4*b1.*b2 - 4*b1.*b3;

%% Pg. 7, Groebner bases, Eqs. 11-12 (not easy for a first-time contributor)
%% Pg. 8, Dattani 2018, Eqs. 13-18 (not easy for a first-time contributor)
%% Pg. 9, Split-Reduc, Eqs. 19-23.
%% Pg. 10, Eqs 25-26

b= dec2bin(2^7-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);b5=b(:,5);b6=b(:,6);ba=b(:,7);

LHS=(-2)*b1.*b2.*b3.*b4.*b5.*b6 + b5.*b6;
RHS=2*(5*ba - b1.*ba - b2.*ba - b3.*ba - b4.*ba - b5.*ba - b6.*ba) + b5.*b6;

%% Pg. 10, Eq 27

b= dec2bin(2^7-1:-1:0)-'0'; % repeated from previous cell, so not technically necessary, but makes this cell more clear.
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);b5=b(:,5);b6=b(:,6);ba=b(:,7);

LHS=min(reshape(-b1.*b2.*b3.*b4.*b5.*b6,2,[]));
RHS=min(reshape((6 - 1 - b1 - b2 - b3 - b4 - b5 - b6).*ba,2,[]));
isequal(LHS,RHS);

%% Pg. 11, NTR-ABCG Eq. 30-32
%% Pg. 12, NTR-ABCG-2 Eq. 34 needs a quadratized version of the equation to be added
b = dec2bin(2^7-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);b5=b(:,5);b6=b(:,6);ba=b(:,7);

LHS = min(-2 * b1.*b2.*b3.*b4.*b5.*b6 + b5.*b6);
% so, as per before, we can rewrite b1b2...b6 in terms of equation.33, we just need to quadratize the b1...b6 portion
RHS = min(2 * ((2 * 6 - 1) * ba - 2 * (b1.*ba + b2.*ba + b3.*ba + b4.*ba + b5.*ba + b6.*ba)) + b5.*b6);
isequal(LHS, RHS);

% retesting, for equation 36, using expanded version
RHS = min(22 * ba - 4 * b1.*ba - 4 * b2.*ba - 4*b3.*ba - 4*b4.*ba - 4*b5.*ba - 4*b6.*ba + b5.*b6)
isequal(LHS, RHS);

%% Pg. 13, NTR-GBP: -b1b2b3 = min_ba(ba - b1 + b2 + b3 - b1b2 - b1b3 + b1)

b= dec2bin(2^4-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);ba=b(:,4);
LHS = min(reshape(-1*b1.*b2.*b3,2,[]));
RHSa = min(reshape(ba.*(-b1 + b2 + b3) - b1.*b2 - b1.*b3 + b1,2,[]));
RHSb = min(reshape(ba.*(-b2 + b1 + b3) - b1.*b2 - b2.*b3 + b2,2,[]));
RHSc = min(reshape(ba.*(-b3 + b1 + b2) - b2.*b3 - b1.*b3 + b3,2,[]));
isequal(LHS,RHSa,RHSb,RHSc);

%% Pg. 14, NTR-YXKK Eq. 42-43

% Page 14 NTR-YXKK equation 42
% Verification for when k=6 and alpha=0.001, 0.5, 1, 2, 15
b= dec2bin(2^7 -1:-1:0) - '0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);b5=b(:,5);b6=b(:,6);ba=b(:,7);
LHS=min(reshape(-b1.*b2.*b3.*b4.*b5.*b6,2,[]));

% alpha=0.001
RHS=min(reshape(ba - 1 + 0.001*((1-b1).*(1-ba) + (1-b2).*(1-ba) + ...
    (1-b3).*(1-ba) + (1-b4).*(1-ba) + (1-b5).*(1-ba) + ...
    (1-b6).*(1-ba)), 2, []));

isequal(LHS,RHS);

% alpha=0.5
RHS=min(reshape(ba - 1 + 0.5*((1-b1).*(1-ba) + (1-b2).*(1-ba) + ...
    (1-b3).*(1-ba) + (1-b4).*(1-ba) + (1-b5).*(1-ba) + ...
    (1-b6).*(1-ba)), 2, []));

isequal(LHS,RHS);

% alpha=1
RHS=min(reshape(ba - 1 + 1*((1-b1).*(1-ba) + (1-b2).*(1-ba) + ...
    (1-b3).*(1-ba) + (1-b4).*(1-ba) + (1-b5).*(1-ba) + ...
    (1-b6).*(1-ba)), 2, []));

isequal(LHS,RHS);

% alpha=2
RHS=min(reshape(ba - 1 + 2*((1-b1).*(1-ba) + (1-b2).*(1-ba) + ...
    (1-b3).*(1-ba) + (1-b4).*(1-ba) + (1-b5).*(1-ba) + ...
    (1-b6).*(1-ba)), 2, []));

isequal(LHS,RHS);

% alpha=15
RHS=min(reshape(ba - 1 + 15*((1-b1).*(1-ba) + (1-b2).*(1-ba) + ...
    (1-b3).*(1-ba) + (1-b4).*(1-ba) + (1-b5).*(1-ba) + ...
    (1-b6).*(1-ba)), 2, []));

isequal(LHS,RHS);

% Equation 44

b= dec2bin(2^6 -1:-1:0) - '0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);b5=b(:,5);ba=b(:,6);
LHS=min(reshape(-b1.*b2.*b3.*b4.*b5 + 5*b1.*b4 - b3, 2, []));
RHS=min(reshape(ba + 9 -10*ba  - 2*b1 - 2*b2 - 2*b3 - 2*b4 - 2*b5 + ...
    2*b1.*ba +2*b2.*ba +2*b3.*ba + 2*b4.*ba + 2*b5.*ba +5*b1.*b4 - b3,  2, []));

isequal(LHS,RHS);
%% Pg. 15, NTR-RBL
%% Pg. 16, NTR-LHZ
%% Pg. 17, PTR-BG

b = dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);ba2=b(:,6);
LHS=min(reshape(b1.*b2.*b3.*b4,4,[]));
RHS=min(reshape(ba1.*(2+ b1 -b2 - b3 - b4) + ba2.*(1 + b2 - b3 - b4) + b3.*b4,4,[]));
isequal(LHS,RHS); % Gives 1, verified by Nike on 28 Feb 2021.

%% Pg. 18, PTR-Ishikawa

b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = ba.*(3 - 2*b1 - 2*b2 - 2*b3 - 2*b4) + b1.*b2 + b1.*b3 + b1.*b4 + b2.*b3 + b2.*b4 +  b3.*b4;

%% Pg. 19, PTR-ABCG (No example given yet!)
%% Pg. 20, PTR-BCR-1

b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = 1/2*(b1 + b2 + b3 + b4 - 2*ba1).*(b1 + b2 + b3 + b4 - 2*ba1 - 1);

%% Pg. 21, PTR-BCR-2
%% Pg. 22, PTR-BCR-3 (example appears to be the same as PTR-BCR-1, and may have to be redone)
%% Pg. 23, PTR-BCR-4 
%% Pg. 24, PTR-KZ (needs an example!)
%% PTR-KZ: b1b2b3 = min_ba(1 − (ba + b1 + b2 + b3) + ba (b1 + b2 + b3) + b1b2 + b1b3 + b2b3)

LHS=min(reshape(b1.*b2.*b3,2,[]));
RHS=min(reshape(1 - (b4+b1+b2+b3) + b4.*(b1+b2+b3) +b1.*b2+b1.*b3+b2.*b3,2,[]));
isequal(LHS,RHS);

%% Pg. 25, PTR-KZ (needs in example!)
%% PTR-GBP: b1b2b3 = min_ba(ba - b2ba - b3ba +b1ba +b2b3)

LHS = min(reshape(b1.*b2.*b3,2,[]));
RHSa = min(reshape(b4 - b2.*b4 - b3.*b4 + b1.*b4 +b2.*b3,2,[]));
RHSb = min(reshape(b4 - b1.*b4 - b3.*b4 + b2.*b4 + b1.*b3,2,[]));
RHSc = min(reshape(b4 - b1.*b4 - b2.*b4 + b3.*b4 + b1.*b2,2,[]));
isequal(LHS,RHSa,RHSb,RHSc);

%% example eq. (?)
LHS = min(reshape(b1.*b2.*b3 + b1.*b3 - b2,2,[]));
RHS = min(reshape((b4 - b1.*b4 - b3.*b4 + b2.*b4 + 2*b1.*b3)-b2,2,[]));
isequal(LHS,RHS);

%% Pg. 26, PTR-GBP (mentioned here: https://proofassistants.stackexchange.com/q/1187/8)

b=dec2bin(2^4-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);ba=b(:,4);
LHS = min(reshape(b1.*b2.*b3,2,[]));
RHS = min(reshape(ba - b2.*ba - b3.*ba + b1.*ba +b2.*b3,2,[]));
isequal(LHS,RHS) % answer is 1

%% Pg. 27, PTR-RBL
%% Pg. 28, PTR-RBL-(4->2)

%% PTR-RBL, k=4, middle of LHZ
z=[1 0; 0 -1];
z1 = kron(z,eye(16));
z2 = kron(kron(eye(2),z),eye(8));
z3 = kron(kron(eye(4),z),eye(4));
z4 = kron(kron(eye(8),z),eye(2));
za = kron(eye(16),z);

LHS=z1*z2*z3*z4;
RHS=4*za*(z1+z2+z3+z4)+2*(z1*z2+z1*z3+z1*z4+z2*z3+z2*z4+z3*z4)+8;

[diag(LHS) diag(RHS)];

%% PTR-RBL, k=5, middle of LHZ (not stated to be true in paper, so testing it).

z1 = kron(z,eye(32));
z2 = kron(kron(eye(2),z),eye(16));
z3 = kron(kron(eye(4),z),eye(8));
z4 = kron(kron(eye(8),z),eye(4));
z5 = kron(kron(eye(16),z),eye(2));
za = kron(eye(32),z);

LHS=z1*z2*z3*z4*z5;
RHS=4*za*(z1+z2+z3+z4+z5)+2*(z1*z2+z1*z3+z1*z4+z1*z5+z2*z3+z2*z4+z2*z5+z3*z4+z3*z5+z4*z5)+9;

[diag(LHS) diag(RHS)];

%% PTR-RBL, k=3, end of LHZ
z1 = kron(z,eye(8));
z2 = kron(kron(eye(2),z),eye(4));
z3 = kron(kron(eye(4),z),eye(2));
za = kron(eye(8),z);

LHS=z1*z2*z3;
RHS_theirs=4*za+2*(z1+z2+z3)+4*za*(z1+z2+z3)+2*(z1*z2+z1*z3+z2*z3)+8;
RHS=2*za+(z1+z2+z3)+2*za*(z1+z2+z3)+(z1*z2+z1*z3+z2*z3)+3;

[diag(LHS) diag(RHS)];

%% PTR-RBL, k=4, end of LHZ
z1 = kron(z,eye(16));
z2 = kron(kron(eye(2),z),eye(8));
z3 = kron(kron(eye(4),z),eye(4));
z4 = kron(kron(eye(8),z),eye(2));
za = kron(eye(16),z);

LHS=z1*z2*z3*z4;
RHS=4*za+2*(z1+z2+z3+z4)+4*za*(z1+z2+z3+z4)+2*(z1*z2+z1*z3+z1*z4+z2*z3+z2*z4+z3*z4)+9;

[diag(LHS) diag(RHS)];

%% pG. 29 PTR-YXKK: b1b2b3 + b1b3 − b2 -> b1b2 + ba1 + 2ba2 + 2(1 − b1)(1 − ba1) + 2(1 − b2)(1 − ba1) + 2(1 − b3)(1 − ba2) + 2(1 − ba2)(1 − ba1) − 2(1 − b3) − 1 + b1b3 − b2
b=dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1); b2=b(:,2); b3=b(:,3); ba1=b(:,4); ba2=b(:,5);
LHS=min(reshape(b1.*b2.*b3 + b1.*b3 - b2,4,[]));
RHS=min(reshape(b1.*b2 + ba1 + 2*ba2 + 2*(1-b1).*(1-ba1) + 2*(1-b2).*(1-ba1) + 2*(1-b3).*(1-ba2) + 2*(1-ba2).*(1-ba1) - 2*(1-b3) - 1 + b1.*b3 - b2,4,[]));
isequal(LHS,RHS);

%% Pg. 30, PTR-CZW
%% Pg. 31, Bit-flipping
%% Pg. 32, SFR-ABCG-1

%% Pg. 33, SFR-BCR-1: (b1b2 + b1b3 + b1b4 + b2b3 + b2b4 + b3b4) − 3(b1b2b3 + b1b2b4 + b1b3b4 + b2b3b4) + 6b1b2b3b4 -> (−3 + b1 + b2 + b3 + b4 − ba1 + 3ba2)^2
b= dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1); b2=b(:,2); b3=b(:,3); b4=b(:,4); ba1=b(:,5); ba2=b(:,6);
LHS= min(reshape( (b1.*b2 + b1.*b3 + b1.*b4 + b2.*b3 + b2.*b4 + b3.*b4) - 3*(b1.*b2.*b3 + b1.*b2.*b4 + b1.*b3.*b4 + b2.*b3.*b4) + 6*b1.*b2.*b3.*b4 ,4,[]));
RHS= min(reshape( (-3 + b1 + b2 + b3 + b4 - ba1 + 3*ba2).^2 ,4,[]));
isequal(LHS, RHS);

%% Pg. 34, SFR-BCR-2: (b1b2 + b1b3 + b1b4 + b2b3 + b2b4 + b3b4)−3(b1b2b3 + b1b2b4 + b1b3b4 + b2b3b4) + 6b1b2b3b4 -> (1 − b1 − b2 − b3 − b4 − ba1 + 3ba2)^2
b= dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1); b2=b(:,2); b3=b(:,3); b4=b(:,4); ba1=b(:,5); ba2=b(:,6);
LHS= min(reshape( (b1.*b2 + b1.*b3 + b1.*b4 + b2.*b3 + b2.*b4 + b3.*b4) - 3*(b1.*b2.*b3 + b1.*b2.*b4 + b1.*b3.*b4 + b2.*b3.*b4) + 6*b1.*b2.*b3.*b4 ,4,[]));
RHS= min(reshape( (1 - b1 - b2 - b3 - b4 - ba1 + 3*ba2).^2 ,4,[]));
isequal(LHS, RHS);

%% Pg. 35, SFR-BCR-3
%% Pg. 36, SFR-BCR-4
%% Pg. 37, SFR-BCR-5
%% Pg. 38, SFR-BCR-6
%% Pg. 39, SFR-ABCG-2
%% Pg. 40, SFR-ABCG-3
%% Pg. 41, SFR-BCR-7
%% Pg. 42, SFR-BCR-8
%% Pg. 43, SFR-BCR-9
%% Pg. 44, SFR-ABCG-4

b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = b4.*b2 + ba.*(b3-1);

%% Pg. 45, PFR-BCR-1
%% Pg. 46, PFR-BCR-2


%% Pg. 49, RBS-Rosenberg: b1b2b3 = min_ba(b1ba + b2b3 - 2*b2ba -2*b3ba + 3*ba) (Eq 150)

b=[0 1];
b1=kron(b,ones(1,8));
b2=kron(kron(ones(1,2),b),ones(1,4));
b3=kron(kron(ones(1,4),b),ones(1,2));
ba=kron(ones(1,8),b);
LHS=min(reshape(b1.*b2.*b3,2,[]));
RHS=min(reshape(b1.*ba + b2.*b3 - 2*b2.*ba - 2*b3.*ba + 3*ba,2,[]));
isequal(LHS,RHS);

%% RBS-Rosenberg: b1b2b3 + b1b2b4 = min_ba(b3ba + b4ba + 2*b1b2 - 4*b1ba - 4*b2ba + 6ba) (Eq 151)

b=dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba=b(:,5);
LHS=min(reshape(b1.*b2.*b3 + b1.*b2.*b4,2,[]));
RHS=min(reshape(b3.*ba + b4.*ba + 2*b1.*b2 - 4*b1.*ba - 4*b2.*ba + 6*ba,2,[]));
isequal(LHS,RHS);

%% Pg. 50, FGBZ for negative: -b1b2b3 - b1b2b4 = min_ba((1 - b1b2 - b3)ba1 + (1 - b1b2 - b4)ba1)  (Eqs 153-154)

b=dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);ba2=b(:,6);
LHS=min(reshape(-b1.*b2.*b3 - b1.*b2.*b4,4,[]));
RHSa=min(reshape((1 - b1.*b2 - b3).*ba1 + (1 - b1.*b2 - b4).*ba1,4,[]));
RHSb=min(reshape(2*ba1 - b3.*ba1 - b4.*ba1 - 2*b1.*b2.*ba1,4,[]));
isequal(LHS,RHSa,RHSb);

%% FGBZ for negative: -b1b2b3 - b1b2b4 = min_ba(2*ba1 - b3ba1 - b4ba1 + 2*(2 - b1 - b2 - ba1)ba2)  (Eqs 155-156)

LHS=min(reshape(-b1.*b2.*b3 - b1.*b2.*b4,4,[]));
RHSa=min(reshape(2*ba1 - b3.*ba1 - b4.*ba1 + 2*(2 - b1 - b2 - ba1).*ba2,4,[]));
RHSb=min(reshape(2*ba1 - b3.*ba1 - b4.*ba1 + 4*ba2 - 2*b1.*ba2 - 2*b2.*ba2 - 2*ba1.*ba2,4,[]));
isequal(LHS,RHSa,RHSb);

%% Pg. 51, FGBZ for positive: b1b2b3 + b1b2b4 = min_ba(2*ba1b1 + (1 - ba1)b2b3 + (1 - ba1)b2b4)  (Eqs 158-159)

b=dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);ba2=b(:,6);
LHS=min(reshape(b1.*b2.*b3 + b1.*b2.*b4,4,[]));
RHSa=min(reshape(2*ba1.*b1 + (1 - ba1).*b2.*b3 + (1 - ba1).*b2.*b4,4,[]));
RHSb=min(reshape(2*ba1.*b1 + b2.*b3 + b2.*b4 - ba1.*b2.*b3 - ba1.*b2.*b4,4,[]));
isequal(LHS,RHSa,RHSb);

%% FGBZ for positive: b1b2b3 + b1b2b4 = min_ba(2*ba1b1 + b2b3 + b2b4 + 4*ba2 - 2*ba2ba1 - 2*ba2b2 - ba2b3 - ba2b4)  (Eqs 160-161)

LHS=min(reshape(b1.*b2.*b3 + b1.*b2.*b4,4,[]));
RHSa=min(reshape(2*ba1.*b1 + b2.*b3 + b2.*b4 + 2*ba2 - ba2.*ba1 - ba2.*b2 - ba2.*b3 + 2*ba2 - ba2.*ba1 - ba2.*b2 - ba2.*b4,4,[]));
RHSb=min(reshape(2*ba1.*b1 + b2.*b3 + b2.*b4 + 4*ba2 - 2*ba2.*ba1 - 2*ba2.*b2 - ba2.*b3 - ba2.*b4,4,[]));
isequal(LHS,RHSa,RHSb);

%% Pg. 52, Pairwise Covers: 5*b1b2b3b4 + 4*b1b2b4 - 3*b1b3 - 2*b2b3b4
%%                  = min_ba(5*ba1ba2 + 4*b1ba2 - 3*ba1 - 2*b3ba2 + 8*(ba1(3 - 2*b1 - 2*b3) + b1b3) + 11*(ba2(3 - 2*b2 - 2*b4) + b2b4))

b=dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);ba2=b(:,6);
LHS=min(reshape(5*b1.*b2.*b3.*b4 + 4*b1.*b2.*b4 - 3*b1.*b3 - 2*b2.*b3.*b4,4,[]));
RHS=min(reshape(5*ba1.*ba2 + 4*b1.*ba2 - 3*ba1 - 2*b3.*ba2 + 8*(ba1.*(3 - 2*b1 - 2*b3) + b1.*b3) + 11*(ba2.*(3 - 2*b2 - 2*b4) + b2.*b4),4,[]));
isequal(LHS,RHS);

%% Pg. 53, Flag-based SAT mapping
%% Pg. 55, SCM-BCR
%% Pg. 56, Decomposition into symmetric and anti-symmetric parts

b = dec2bin(2^4-1:-1:0)-'0';
b1 = b(:,1); b2 = b(:,2); b3 = b(:,3); b4 = b(:,4);

f = b1 + 2*b1.*b2 - 4*b1.*b2.*b3 + 2*b1.*b2.*b3.*b4;

f_sym = (1/2)*(f + (1 - f));
f_anti = (1/2)*(f - (1 - f));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QUANTUM GADGETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pg. 57, ZZZ-TI-CBBK

z = [1 0; 0 -1];
x = [0 1;1 0];
zi = kron(z,eye(8));
zj = kron(kron(eye(2),z),eye(4));
zk = kron(kron(eye(4),z),eye(2));
za = kron(eye(8),z);
xa = kron(eye(8),x);

alpha = 1;

for delta = 1:1e10:1e11   
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

LHS = alpha*zi*zj*zk;
RHS = alpha_I*eye(16) + alpha_zi*zi + alpha_zj*zj + alpha_zk*zk + alpha_za*za + alpha_xa*xa + alpha_zzia*zi*za + alpha_zzja*zj*za + alpha_zzka*zk*za + alpha_zzij*zi*zj + alpha_zzik*zi*zk + alpha_zzjk*zj*zk;
RHS_alt = ( delta*eye(16) + ((alpha*(delta^4)/6)^(1/5))*(zi + zj + zk))*((1*eye(16) - za)/2) + (alpha*delta^4/6)^(1/5)*xa + ( ((alpha/6)^(2/5)*(delta^(3/5)))*eye(16) - (((alpha/6)^(3/5)*(delta^(2/5))) + (7*alpha/6))*(zi + zj + zk) )*((1*eye(16) + za)/2)  + ((alpha/6)^(4/5)*(delta^(1/5)))*(3*eye(16) + 2*zi*zj + 2*zi*zk + 2*zj*zk);

min(eig(LHS)) - min(eig(RHS))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:8:end); EL = EL(:,1:8:end); VR = VR(:,1:8:end); ER = ER(:,1:8:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3345

%% NP-OY: Example (8->4-body, 4->2-body is a rewrite)
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

x_113  = sparse(kron(kron(speye(2^0),x),speye(2^17)));
z_113  = sparse(kron(kron(speye(2^0),z),speye(2^17)));

x_122  = sparse(kron(kron(speye(2^1),x),speye(2^16)));
z_122  = sparse(kron(kron(speye(2^1),z),speye(2^16)));

x_114  = sparse(kron(kron(speye(2^2),x),speye(2^15)));
z_114  = sparse(kron(kron(speye(2^2),z),speye(2^15)));

x_124  = sparse(kron(kron(speye(2^3),x),speye(2^14)));

x_211  = sparse(kron(kron(speye(2^4),x),speye(2^13)));
z_211  = sparse(kron(kron(speye(2^4),z),speye(2^13)));

x_221  = sparse(kron(kron(speye(2^5),x),speye(2^12)));

x_213  = sparse(kron(kron(speye(2^6),x),speye(2^11)));

x_222  = sparse(kron(kron(speye(2^7),x),speye(2^10)));

z_111  = sparse(kron(kron(speye(2^8),z),speye(2^9)));

z_112  = sparse(kron(kron(speye(2^9),z),speye(2^8)));

z_014  = sparse(kron(kron(speye(2^10),z),speye(2^7)));

z_103  = sparse(kron(kron(speye(2^11),z),speye(2^6)));

x_a111 = sparse(kron(kron(speye(2^12),x),speye(2^5)));
y_a111 = sparse(kron(kron(speye(2^12),y),speye(2^5)));
z_a111 = sparse(kron(kron(speye(2^12),z),speye(2^5)));

x_a112 = sparse(kron(kron(speye(2^13),x),speye(2^4)));
y_a112 = sparse(kron(kron(speye(2^13),y),speye(2^4)));
z_a112 = sparse(kron(kron(speye(2^13),z),speye(2^4)));

x_a211 = sparse(kron(kron(speye(2^14),x),speye(2^3)));
y_a211 = sparse(kron(kron(speye(2^14),y),speye(2^3)));
z_a211 = sparse(kron(kron(speye(2^14),z),speye(2^3)));

x_a212 = sparse(kron(kron(speye(2^15),x),speye(2^2)));
y_a212 = sparse(kron(kron(speye(2^15),y),speye(2^2)));
z_a212 = sparse(kron(kron(speye(2^15),z),speye(2^2)));

x_a121 = sparse(kron(kron(speye(2^16),x),speye(2^1)));
y_a121 = sparse(kron(kron(speye(2^16),y),speye(2^1)));
z_a121 = sparse(kron(kron(speye(2^16),z),speye(2^1)));

x_a122 = sparse(kron(kron(speye(2^17),x),speye(2^0)));
y_a122 = sparse(kron(kron(speye(2^17),y),speye(2^0)));
z_a122 = sparse(kron(kron(speye(2^17),z),speye(2^0)));

I = speye(2^18);

U = 3; t = 0.375*U; J = 0.09*U;

LHS = sparse(-J*(x_113*x_122*x_114*x_124*x_211*x_221*x_213*x_222 + z_111*z_112*z_113*z_114 + z_014*z_111 + z_112*z_103 + z_114*z_211 + z_113*z_122));

RHS = sparse(J*(-z_111*z_112*z_113*z_114 - z_014*z_111 - z_112*z_103 - z_114*z_211 - z_113*z_122 ...
    + (I - z_a111 + z_a112 + z_a111*z_a112) * (z_a121 + z_a122 + z_a121*z_a122 - I) ...
    + (I + z_a111 - z_a112 + z_a111*z_a112) * (I - z_a211 - z_a212 - z_a211*z_a212)) ...
    - (U/2)*(z_a111 + z_a112 + z_a111*z_a112 - I) ...
    - (t/2)*((x_a112 + z_a111*x_a112)*x_113*x_114 + (x_a111*x_a112 + y_a111*y_a112)*x_122*x_124 ...
    + (x_a112 - z_a111*x_a112)*x_221*x_222 + (x_a111*x_a112 - y_a111*y_a112)*x_211*x_213 )) + 5.218055837657829*I;

abs(eigs(LHS, 1, 'smallestreal')-eigs(RHS, 1, 'smallestreal')) % 2.0428e-14

[VL, EL] = eigs(LHS, 128, 'smallestreal'); [VR, ER] = eigs(RHS, 128, 'smallestreal');
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.4013

%% NP - SJ (4z5 - 3x1 + 2*z1*y2*x5 + 9*x1*x2*x3*x4 - x1*y2*z3*x5 -> 9*xa1 + 4*za2*z5 - 3*za3*x1 - za3*xa2 + 2*xa3*x5)
x = [0 1;1 0];
y = [0 -1i;1i 0];
z = [1 0;0 -1];

x1 = kron(x,eye(128));x2 = kron(kron(eye(2),x),eye(64));x3 = kron(kron(eye(4),x),eye(32));x4 = kron(kron(eye(8),x),eye(16));x5 = kron(kron(eye(16),x),eye(8));
xa1 = kron(kron(eye(32),x),eye(4));xa2 = kron(kron(eye(64),x),eye(2));xa3 = kron(eye(128),x);
y2 = kron(kron(eye(2),y),eye(64));
z1 = kron(z,eye(128));z3 = kron(kron(eye(4),z),eye(32));z5 = kron(kron(eye(16),z),eye(8));
za1 = kron(kron(eye(32),z),eye(4));
za2 = kron(kron(eye(64),z),eye(2));
za3 = kron(eye(128),z);

LHS = 4*z5 - 3*x1 + 2*z1*y2*x5 + 9*x1*x2*x3*x4 - x1*y2*z3*x5;
RHS = 9*xa1 + 4*za2*z5 - 3*za3*x1 - za3*xa2 + 2*xa3*x5;

max(eig(LHS)-eig(RHS))<1e-13 % gives 1.

U_a1 = ((eye(256) + za1)/2) + ((eye(256) - za1)/2)*x1*x2*x3*x4;
U_a2 = ((eye(256) + za2)/2) + ((eye(256) - za2)/2)*z5;
U_a3 = ((eye(256) + za3)/2) + ((eye(256) - za3)/2)*x1;

RHS_a = U_a1*U_a2*U_a3*RHS;

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS_a);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:8:end); EL = EL(:,1:8:end); VR = VR(:,1:8:end); ER = ER(:,1:8:end); % reduce repeats due to auxiliary qubits

for col = 1:1:size(VL,2) % compare eigenvectors and eigenvalues
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.2770

%% NP-Nagaj-1: Feynman Hamiltonian and H_4-local
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

x1 = kron(kron(eye(2^0),x),eye(2^3));
y1 = kron(kron(eye(2^0),y),eye(2^3));
z1 = kron(kron(eye(2^0),z),eye(2^3));

x2 = kron(kron(eye(2^1),x),eye(2^2));
y2 = kron(kron(eye(2^1),y),eye(2^2));
z2 = kron(kron(eye(2^1),z),eye(2^2));

x3 = kron(kron(eye(2^2),x),eye(2^1));
y3 = kron(kron(eye(2^2),y),eye(2^1));
z3 = kron(kron(eye(2^2),z),eye(2^1));

x4 = kron(kron(eye(2^3),x),eye(2^0));
y4 = kron(kron(eye(2^3),y),eye(2^0));
z4 = kron(kron(eye(2^3),z),eye(2^0));

U_2_local = (1/2)*(eye(2^4) + z3 + x4 - z3*x4);

E_194 = (1/4)*(x1*x2 - 1i*y1*x2 + 1i*x1*y2 + y1*y2)*U_2_local + (1/4)*(x1*x2 + 1i*y1*x2 - 1i*x1*y2 + y1*y2)*conj(transpose(U_2_local));
E_195 = (1/4)*(x1*x2 + y1*y2 + x1*x2*z3 + x1*x2*x4 + y1*y2*z3 + y1*y2*x4 - x1*x2*z3*x4 - y1*y2*z3*x4);

isequal(E_194, E_195)

%% NP-Nagaj-1: 4-local to 2-local
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
gm1 = [0 1 0 ; 1 0 0 ; 0 0 0]; gm2 = [0 -1i 0 ; 1i 0 0 ; 0 0 0]; gm3 = [1 0 0 ; 0 -1 0 ; 0 0 0]; 
gm4 = [0 0 1 ; 0 0 0 ; 1 0 0]; gm5 = [0 0 -1i ; 0 0 0 ; 1i 0 0]; gm6 = [0 0 0 ; 0 0 1 ; 0 1 0];
gm7 = [0 0 0 ; 0 0 -1i ; 0 1i 0]; gm8 = (1/sqrt(3))*[1 0 0 ; 0 1 0 ; 0 0 -2];

x1 = kron(kron(speye(2^0),x),speye(2^7));
y1 = kron(kron(speye(2^0),y),speye(2^7));
z1 = kron(kron(speye(2^0),z),speye(2^7));

x2 = kron(kron(speye(2^1),x),speye(2^6));
y2 = kron(kron(speye(2^1),y),speye(2^6));
z2 = kron(kron(speye(2^1),z),speye(2^6));

x3 = kron(kron(speye(2^2),x),speye(2^5));
y3 = kron(kron(speye(2^2),y),speye(2^5));
z3 = kron(kron(speye(2^2),z),speye(2^5));

x4 = kron(kron(speye(2^3),x),speye(2^4));
y4 = kron(kron(speye(2^3),y),speye(2^4));
z4 = kron(kron(speye(2^3),z),speye(2^4));

xa1 = kron(kron(speye(2^4),x),speye(2^3));
ya1 = kron(kron(speye(2^4),y),speye(2^3));
za1 = kron(kron(speye(2^4),z),speye(2^3));

xa2 = kron(kron(speye(2^5),x),speye(2^2));
ya2 = kron(kron(speye(2^5),y),speye(2^2));
za2 = kron(kron(speye(2^5),z),speye(2^2));

xa3 = kron(kron(speye(2^6),x),speye(2^1));
ya3 = kron(kron(speye(2^6),y),speye(2^1));
za3 = kron(kron(speye(2^6),z),speye(2^1));

xa4 = kron(kron(speye(2^7),x),speye(2^0));
ya4 = kron(kron(speye(2^7),y),speye(2^0));
za4 = kron(kron(speye(2^7),z),speye(2^0));

lam_1a5 = kron(kron(speye(3^0),gm1),speye(3^5));
lam_2a5 = kron(kron(speye(3^0),gm2),speye(3^5));
lam_4a5 = kron(kron(speye(3^0),gm4),speye(3^5));
lam_5a5 = kron(kron(speye(3^0),gm5),speye(3^5));
lam_6a5 = kron(kron(speye(3^0),gm6),speye(3^5));

lam_1a6 = kron(kron(speye(3^1),gm1),speye(3^4));
lam_2a6 = kron(kron(speye(3^1),gm2),speye(3^4));
lam_4a6 = kron(kron(speye(3^1),gm4),speye(3^4));
lam_5a6 = kron(kron(speye(3^1),gm5),speye(3^4));
lam_6a6 = kron(kron(speye(3^1),gm6),speye(3^4));

lam_1a7 = kron(kron(speye(3^2),gm1),speye(3^3));
lam_2a7 = kron(kron(speye(3^2),gm2),speye(3^3));
lam_4a7 = kron(kron(speye(3^2),gm4),speye(3^3));
lam_5a7 = kron(kron(speye(3^2),gm5),speye(3^3));
lam_6a7 = kron(kron(speye(3^2),gm6),speye(3^3));

lam_1a8 = kron(kron(speye(3^3),gm1),speye(3^2));
lam_2a8 = kron(kron(speye(3^3),gm2),speye(3^2));
lam_4a8 = kron(kron(speye(3^3),gm4),speye(3^2));
lam_5a8 = kron(kron(speye(3^3),gm5),speye(3^2));
lam_6a8 = kron(kron(speye(3^3),gm6),speye(3^2));

lam_1a9 = kron(kron(speye(3^4),gm1),speye(3^1));
lam_2a9 = kron(kron(speye(3^4),gm2),speye(3^1));
lam_4a9 = kron(kron(speye(3^4),gm4),speye(3^1));
lam_5a9 = kron(kron(speye(3^4),gm5),speye(3^1));
lam_6a9 = kron(kron(speye(3^4),gm6),speye(3^1));

lam_1a10 = kron(kron(speye(3^5),gm1),speye(3^0));
lam_2a10 = kron(kron(speye(3^5),gm2),speye(3^0));
lam_4a10 = kron(kron(speye(3^5),gm4),speye(3^0));
lam_5a10 = kron(kron(speye(3^5),gm5),speye(3^0));
lam_6a10 = kron(kron(speye(3^5),gm6),speye(3^0));

LHS = (1/4)*(kron(x1*x2,eye(3^6)) + kron(y1*y2,eye(3^6)) + kron(x1*x2*z3,eye(3^6)) + kron(x1*x2*x4,eye(3^6)) + kron(y1*y2*z3,eye(3^6)) + kron(y1*y2*x4,eye(3^6)) - kron(x1*x2*z3*x4,eye(3^6)) - kron(y1*y2*z3*x4,eye(3^6)));
RHS = (1/2)*(kron(x1,lam_1a5) + kron(y1,lam_2a5) + kron(lam_6a5,eye(2^8)) - kron(lam_6a5,z3) + kron(lam_5a5,ya1) + kron(lam_4a5,xa1) + kron(xa1,lam_1a6) + kron(ya1,lam_2a6) ...
    + 2*kron(lam_6a6,x4) + kron(lam_5a6,ya2) + kron(lam_4a6,xa2) + kron(xa2,lam_1a7) + kron(ya2,lam_2a7) + kron(lam_6a7,eye(2^8)) - kron(lam_6a7,z3) ...
    + kron(lam_5a7,y2) + kron(lam_4a7,x2) + kron(x1,lam_1a8) + kron(y1,lam_2a8) + kron(lam_6a8,eye(2^8)) + kron(lam_6a8,z3) + kron(lam_5a8,ya3) + kron(lam_4a8,xa3) ...
    + kron(xa3,lam_1a9) + kron(ya3,lam_2a9) + 2*kron(lam_6a9,eye(2^8)) + kron(lam_5a9,ya4) + kron(lam_4a9,xa4) + kron(xa4,lam_1a10) + kron(ya4,lam_2a10) ...
    + kron(lam_6a10,eye(2^8)) + kron(lam_6a10,z3) + kron(lam_5a10,y2) + kron(lam_4a10,x2) ) + 6.764818490514657*speye(186624);

eigs(LHS, 1, 'smallestreal')-eigs(RHS, 1, 'smallestreal') % -1.9984e-15

[VL, EL] = eigs(LHS, 256, 'smallestreal'); [VR, ER] = eigs(RHS, 256, 'smallestreal');
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.4114

%% NP-Nagaj-2
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

x1 = kron(kron(speye(2^0),x),speye(2^19));
y1 = kron(kron(speye(2^0),y),speye(2^19));
z1 = kron(kron(speye(2^0),z),speye(2^19));

x2 = kron(kron(speye(2^1),x),speye(2^18));
y2 = kron(kron(speye(2^1),y),speye(2^18));
z2 = kron(kron(speye(2^1),z),speye(2^18));

x3 = kron(kron(speye(2^2),x),speye(2^17));
y3 = kron(kron(speye(2^2),y),speye(2^17));
z3 = kron(kron(speye(2^2),z),speye(2^17));

x4 = kron(kron(speye(2^3),x),speye(2^16));
y4 = kron(kron(speye(2^3),y),speye(2^16));
z4 = kron(kron(speye(2^3),z),speye(2^16));

xa1 = kron(kron(speye(2^4),x),speye(2^15));
ya1 = kron(kron(speye(2^4),y),speye(2^15));
za1 = kron(kron(speye(2^4),z),speye(2^15));

xa2 = kron(kron(speye(2^5),x),speye(2^14));
ya2 = kron(kron(speye(2^5),y),speye(2^14));
za2 = kron(kron(speye(2^5),z),speye(2^14));

xa3 = kron(kron(speye(2^6),x),speye(2^13));
ya3 = kron(kron(speye(2^6),y),speye(2^13));
za3 = kron(kron(speye(2^6),z),speye(2^13));

xa4 = kron(kron(speye(2^7),x),speye(2^12));
ya4 = kron(kron(speye(2^7),y),speye(2^12));
za4 = kron(kron(speye(2^7),z),speye(2^12));

xa5 = kron(kron(speye(2^8),x),speye(2^11));
ya5 = kron(kron(speye(2^8),y),speye(2^11));
za5 = kron(kron(speye(2^8),z),speye(2^11));

xa6 = kron(kron(speye(2^9),x),speye(2^10));
ya6 = kron(kron(speye(2^9),y),speye(2^10));
za6 = kron(kron(speye(2^9),z),speye(2^10));

xa7 = kron(kron(speye(2^10),x),speye(2^9));
ya7 = kron(kron(speye(2^10),y),speye(2^9));
za7 = kron(kron(speye(2^10),z),speye(2^9));

xa8 = kron(kron(speye(2^11),x),speye(2^8));
ya8 = kron(kron(speye(2^11),y),speye(2^8));
za8 = kron(kron(speye(2^11),z),speye(2^8));

xa9 = kron(kron(speye(2^12),x),speye(2^7));
ya9 = kron(kron(speye(2^12),y),speye(2^7));
za9 = kron(kron(speye(2^12),z),speye(2^7));

xa10 = kron(kron(speye(2^13),x),speye(2^6));
ya10 = kron(kron(speye(2^13),y),speye(2^6));
za10 = kron(kron(speye(2^13),z),speye(2^6));

xa11 = kron(kron(speye(2^14),x),speye(2^5));
ya11 = kron(kron(speye(2^14),y),speye(2^5));
za11 = kron(kron(speye(2^14),z),speye(2^5));

xa12 = kron(kron(speye(2^15),x),speye(2^4));
ya12 = kron(kron(speye(2^15),y),speye(2^4));
za12 = kron(kron(speye(2^15),z),speye(2^4));

xa13 = kron(kron(speye(2^16),x),speye(2^3));
ya13 = kron(kron(speye(2^16),y),speye(2^3));
za13 = kron(kron(speye(2^16),z),speye(2^3));

xa14 = kron(kron(speye(2^17),x),speye(2^2));
ya14 = kron(kron(speye(2^17),y),speye(2^2));
za14 = kron(kron(speye(2^17),z),speye(2^2));

xa15 = kron(kron(speye(2^18),x),speye(2^1));
ya15 = kron(kron(speye(2^18),y),speye(2^1));
za15 = kron(kron(speye(2^18),z),speye(2^1));

xa16 = kron(kron(speye(2^19),x),speye(2^0));
ya16 = kron(kron(speye(2^19),y),speye(2^0));
za16 = kron(kron(speye(2^19),z),speye(2^0));

LHS = (1/4)*(x1*x2 + y1*y2 + x1*x2*z3 + x1*x2*x4 + y1*y2*z3 + y1*y2*x4 - x1*x2*z3*x4 - y1*y2*z3*x4);
RHS = (1/4)*(za5 - za6 - za5*z3 + za6*z3 + za9 - za10 - za9*z3 + za10*z3 + za11 - za12 + za11*z3 - za12*z3 + za15 - za16 + za15*z3 - za16*z3) ...
    + (1/2)*(za7*x4 - za8*x4 + za13 - za14) ...
    + (1/(2*sqrt(2)))*(x1*xa5 + y1*ya5 + x1*xa6 + y1*ya6 - xa5*xa1 - ya5*ya1 + xa6*xa1 + ya6*ya1 + xa1*xa7 + ya1*ya7 + xa1*xa8 + ya1*ya8 ...
    - xa7*xa2 - ya7*ya2 + xa8*xa2 + ya8*ya2 + xa2*xa9 + ya2*ya9 + xa2*xa10 + ya2*ya10 - xa9*x2 - ya9*y2 + xa10*x2 + ya10*y2 ...
    + x1*xa11 + y1*ya11 + x1*xa12 + y1*ya12 - xa11*xa3 - ya11*ya3 + xa12*xa3 + ya12*ya3 + xa3*xa13 + ya3*ya13 + xa3*xa14 + ya3*ya14 ...
    - xa13*xa4 - ya13*ya4 + xa14*xa4 + ya14*ya4 + xa4*xa15 + ya4*ya15 + xa4*xa16 + ya4*ya16 - xa15*x2 - ya15*y2 + xa16*x2 + ya16*y2) + 8.211261928531997*speye(2^20);

eigs(LHS, 1, 'smallestreal')-eigs(RHS, 1, 'smallestreal') % 2.9088e-14

[VL, EL] = eigs(LHS, 32, 'smallestreal'); [VR, ER] = eigs(RHS, 32, 'smallestreal');
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.4129

%% P(3->2)-DC1

x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
x1 = kron(x,eye(64)); x2 = kron(kron(eye(2),x),eye(32));
xa1 = kron(kron(eye(16),x),eye(4)); xa2 = kron(kron(eye(32),x),eye(2)); xa3 = kron(eye(64),x);
z1 = kron(z,eye(64)); z3 = kron(kron(eye(4),z),eye(16));
za1 = kron(kron(eye(16),z),eye(4)); za2 = kron(kron(eye(32),z),eye(2)); za3 = kron(eye(64),z);
y4 = kron(kron(eye(8),y),eye(8));

for delta = 1:1e9:1e10
        alpha = (1/8)*delta;
        alpha_ss = (1/6)*(delta)^(1/3);
        alpha_sx = (-1/6)*(delta)^(2/3);
        alpha_zz = (-1/24)*delta;

LHS = (x1 + 3*x2)*z3*y4 + z1*x2;
RHS = z1*x2 + alpha*eye(128) + alpha_ss*((x1 + 3*x2)^2 + z3^2 + y4^2)+ alpha_sx*(x1 + 3*x2)*xa1 + alpha_sx*z3*xa2 + alpha_sx*y4*xa3 + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

min(eig(LHS))-min(eig(RHS));
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:8:end); EL = EL(:,1:8:end); VR = VR(:,1:8:end); ER = ER(:,1:8:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3345

%% P(3->2)-DC2

x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
x1 = kron(x,eye(32)); x2 = kron(kron(eye(2),x),eye(16)); x3 = kron(kron(eye(4),x),eye(8)); x4 = kron(kron(eye(8),x),eye(4));
xa1 = kron(kron(eye(16),x),eye(2)); xa2 = kron(eye(32),x);
y3 = kron(kron(eye(4),y),eye(8)); y4 = kron(kron(eye(8),y),eye(4));
z1 = kron(z,eye(32)); z2 = kron(kron(eye(2),z),eye(16));
za1 = kron(kron(eye(16),z),eye(2)); za2 = kron(eye(32),z);

a1 = 1;
a2 = -3;

for delta = 1:1e8:1e9
alpha = (1/2)*delta;

alpha_s1 = a1*((1/4)*(delta^(2/3)) - 1);
alpha_s2 = a2*((1/4)*(delta^(2/3)) - 1);

alpha_z = (1/2)*delta;
alpha_ss = delta^(1/3);

alpha_sz1 = (a1/4)*(delta^(2/3));
alpha_sz2 = (a2/4)*(delta^(2/3));

alpha_sx = delta^(2/3);

LHS = x1*z2*y3 - 3*x1*x2*y4 + z1*x2;
RHS = z1*x2 + (2*alpha)*eye(64) + alpha_s1*y3 + alpha_s2*y4 + alpha_z*(za1 + za2) + alpha_ss*((x1 + z2)^2 + (x1 + x2)^2) + alpha_sx*(x1*xa1 + z2*xa1 + x1*xa2 + x2*xa2) + alpha_sz1*y3*za1 + alpha_sz2*y4*za2;

abs(min(eig(LHS))-min(eig(RHS)));
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.2708

%% P(3->2)-KKR, example.
x = [0 1;1 0]; y = [0 -1i;1i 0]; z = [1 0;0 -1];
z1 = kron(z,eye(512));z2 = kron(kron(eye(2),z),eye(256));
x1 = kron(x,eye(512));x2 = kron(kron(eye(2),x),eye(256));
y3 = kron(kron(eye(4),y),eye(128));y4 = kron(kron(eye(8),y),eye(64));
xa_11 = kron(kron(eye(16),x),eye(32)); xa_12 = kron(kron(eye(32),x),eye(16)); xa_13 = kron(kron(eye(64),x),eye(8));
xa_21 = kron(kron(eye(128),x),eye(4)); xa_22 = kron(kron(eye(256),x),eye(2)); xa_23 = kron(eye(512),x);
za_11 = kron(kron(eye(16),z),eye(32)); za_12 = kron(kron(eye(32),z),eye(16)); za_13 = kron(kron(eye(64),z),eye(8));
za_21 = kron(kron(eye(128),z),eye(4)); za_22 = kron(kron(eye(256),z),eye(2)); za_23 = kron(eye(512),z);

I_size = 1024;
a = -1;

for delta = 1:1e10:1e11
    
alpha = (3/4)*(delta);
alpha_ss = (delta)^(1/3);
alpha_sx = -(delta)^(2/3);
alpha_zz = -(1/4)*(delta);

LHS = z1*x2 - 5*x1*z2*y3 + 8*x1*x2*y4;
RHS = z1*x2 + 2*alpha*eye(I_size) + alpha_ss*(2*x1^2 + z2^2 + y3^2 + x2^2 + y4^2) ...
    + alpha_sx*(x1*xa_11 + z2*xa_12 + y3*xa_13 + x1*xa_21 + x2*xa_22 + y4*xa_23) ...
    + alpha_zz*(za_11*za_12 + za_11*za_13 + za_12*za_13 + za_21*za_22 + za_21*za_23 + za_22*za_23) - 0.779475460407596*eye(I_size);

abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:64:end); EL = EL(:,1:64:end); VR = VR(:,1:64:end); ER = ER(:,1:64:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3932

%% P(3->2)-KKR, Alternative Form (i.e. original from KKR paper, near Eq. 13 on the arXiv version).

x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
x1 = kron(x,eye(32));
xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
y3 = kron(kron(eye(4),y),eye(8));
z2 = kron(kron(eye(2),z),eye(16));
za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);

for delta = 1:1e10:1e11

alpha = (3/4)*(delta);
alpha_ss = (delta)^(1/3);
alpha_sx = -(delta)^(2/3);
alpha_zz = -(1/4)*(delta);

LHS = -6*x1*z2*y3;
RHS = alpha*eye(64) + 3*alpha_ss*eye(64) + alpha_sx*(x1*xa1 + z2*xa2 + y3*xa3) + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:8:end); EL = EL(:,1:8:end); VR = VR(:,1:8:end); ER = ER(:,1:8:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3576

%% P(3->2)-OT: Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

a = 3;
b = 4;

y1 = kron(y,eye(16));
z1 = kron(z,eye(16));

x2 = kron(kron(eye(2),x),eye(8));
y2 = kron(kron(eye(2),y),eye(8));
z2 = kron(kron(eye(2),z),eye(8));

x3 = kron(kron(eye(4),x),eye(4));
y3 = kron(kron(eye(4),y),eye(4));

za1 = kron(kron(eye(8),z),eye(2));
xa1 = kron(kron(eye(8),x),eye(2));

za2 = kron(eye(16),z);
xa2 = kron(eye(16),x);

for delta = 1:1e9:1e10
    alpha = (delta/2);
    alpha_y1 = ((delta^(1/3))*(a^(2/3))/2);
    alpha_z2 = ((delta^(1/3))*(a^(2/3))/2);
    alpha_x3 = -((delta^(2/3))*(a^(1/3))/2);
    alpha_za = -(delta/2);
    alpha_y1_z2 = -((delta^(1/3))*(a^(2/3)));
    alpha_y1_x3 = (a/2);
    alpha_z2_x3 = (a/2);
    alpha_x3_za = (delta^(2/3))*(a^(1/3)/2);
    alpha_y1_xa = -((delta^(2/3))*(a^(1/3))/sqrt(2));
    alpha_z2_xa = ((delta^(2/3))*(a^(1/3))/sqrt(2));
    
    alpha_z1 = ((delta^(1/3))*(b^(2/3))/2);
    alpha_x2 = ((delta^(1/3))*(b^(2/3))/2);
    alpha_y3 = -((delta^(2/3))*(b^(1/3))/2);
    alpha_z1_x2 = -((delta^(1/3))*(b^(2/3)));
    alpha_z1_y3 = (b/2);
    alpha_x2_y3 = (b/2);
    alpha_y3_za = (delta^(2/3))*(b^(1/3)/2);
    alpha_z1_xa = -((delta^(2/3))*(b^(1/3))/sqrt(2));
    alpha_x2_xa = ((delta^(2/3))*(b^(1/3))/sqrt(2));
    
    LHS = a*y1*z2*x3 + b*z1*x2*y3 - z1*y2;
    RHS = 2*alpha*eye(32) + alpha_y1*(y1)^2 + alpha_z1*(z1)^2 + alpha_z2*(z2)^2 + alpha_x2*(x2)^2 + alpha_x3*x3 + alpha_y3*y3 ...
        + alpha_za*(za1 + za2) + alpha_y1_z2*y1*z2 + alpha_z1_x2*z1*x2 + alpha_y1_x3*(y1^2)*x3 + alpha_z1_y3*(z1^2)*y3 + alpha_z2_x3*(z2^2)*x3 + alpha_x2_y3*(x2^2)*y3 ...
        + alpha_x3_za*x3*za1 + alpha_y3_za*y3*za2 + alpha_y1_xa*y1*xa1 + alpha_z1_xa*z1*xa2 + alpha_z2_xa*z2*xa1 + alpha_x2_xa*x2*xa2 - z1*y2;
    
    abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3672

%% P(3->2)-CBBK: Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

a = 3;
b = 2;
x1 = kron(x,eye(32));
y1 = kron(y,eye(32));
z1 = kron(z,eye(32));

x2 = kron(kron(eye(2),x),eye(16));
y2 = kron(kron(eye(2),y),eye(16));
z2 = kron(kron(eye(2),z),eye(16));

y3 = kron(kron(eye(4),y),eye(8));

z4 = kron(kron(eye(8),z),eye(4));

za1 = kron(kron(eye(16),z),eye(2));
xa1 = kron(kron(eye(16),x),eye(2));

za2 = kron(eye(32),z);
xa2 = kron(eye(32),x);

for Delta = 1:1e11:1e12
    alpha_I1 = Delta/2 + (1/2)*(a/2)^(2/3)*Delta^(1/2)*( (sign(a)^2) + 1 ) - (sign(a)^2)*((a/2)^(4/3))* ((sign(a)^2) + 1);
    alpha_y3 = (1/2)*(a/2)^(1/3)*(Delta^(1/2)) - (a/4)*( (sign(a)^2) + 1 );
    alpha_za1 = (-Delta/2) + (1/2)*(a/2)^(2/3)*Delta^(1/2)*( (sign(a)^2) + 1 ) - (sign(a)^2)*((a/2)^(4/3))* ((sign(a)^2) + 1);
    alpha_x1_z2 = 2*sign(a)*((a/2)^(2/3))*(Delta^(1/2)) - 4*sign(a)*((a/2)^(4/3));
    alpha_y3_za = (-1/2)*(a/2)^(1/3)*(Delta^(1/2)) - (a/4)*( (sign(a)^2) + 1 );
    alpha_x1_xa = sign(a)*(a/2)^(1/3)*(Delta^(3/4));
    alpha_z2_xa = (a/2)^(1/3)*(Delta^(3/4));
    
    alpha_I2 = Delta/2 + (1/2)*(b/2)^(2/3)*Delta^(1/2)*( (sign(b)^2) + 1 ) - (sign(b)^2)*((b/2)^(4/3))* ((sign(b)^2) + 1);
    alpha_z4 = (1/2)*(b/2)^(1/3)*(Delta^(1/2)) - (b/4)*( (sign(b)^2) + 1 );
    alpha_za2 = (-Delta/2) + (1/2)*(b/2)^(2/3)*Delta^(1/2)*( (sign(b)^2) + 1 ) - (sign(b)^2)*((b/2)^(4/3))* ((sign(b)^2) + 1);
    alpha_y1_x2 = 2*sign(b)*((b/2)^(2/3))*(Delta^(1/2)) - 4*sign(b)*((b/2)^(4/3));
    alpha_z4_za = (-1/2)*(b/2)^(1/3)*(Delta^(1/2)) - (b/4)*( (sign(b)^2) + 1 );
    alpha_y1_xa = sign(b)*(b/2)^(1/3)*(Delta^(3/4));
    alpha_x2_xa = (b/2)^(1/3)*(Delta^(3/4));
    
LHS = a*x1*z2*y3 + b*y1*x2*z4 - z1*x2;

RHS = alpha_I1*eye(64) + alpha_I2*eye(64) + alpha_y3*y3 + alpha_z4*z4 + alpha_za1*za1 + alpha_za2*za2 + alpha_x1_z2*x1*z2 + alpha_y1_x2*y1*x2 ...
    + alpha_y3_za*y3*za1 + alpha_z4_za*z4*za2 + alpha_x1_xa*x1*xa1 + alpha_y1_xa*y1*xa2 + alpha_z2_xa*z2*xa1 + alpha_x2_xa*x2*xa2 - z1*x2;

abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3694

%% PSD-OT: Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

alpha = 3;
x1 = kron(x,eye(64));
z2 = kron(kron(eye(2),z),eye(32));
y3 = kron(kron(eye(4),y),eye(16));
z4 = kron(kron(eye(8),z),eye(8));
x5 = kron(kron(eye(16),x),eye(4));
y6 = kron(kron(eye(32),y),eye(2));

A = x1*z2*y3; B = z4*x5*y6;

za = kron(eye(64),z);
xa = kron(eye(64),x);

delta_array = [];
ans_array = [];

for delta = 1:1e6:1e7
    LHS = alpha*x1*z2*y3*z4*x5*y6;
    RHS = delta*((1*eye(128) - za)/2) + (alpha/2)*(A^2 + B^2) + sqrt( alpha*delta/2 )*(-A + B)*xa;   
end


[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:16:end); EL = EL(:,1:16:end); VR = VR(:,1:16:end); ER = ER(:,1:16:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 0.7113

%% PSD-CBBK: Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

alpha = 5;
x1 = kron(x,eye(64));
z2 = kron(kron(eye(2),z),eye(32));
y3 = kron(kron(eye(4),y),eye(16));
z4 = kron(kron(eye(8),z),eye(8));
x5 = kron(kron(eye(16),x),eye(4));
y6 = kron(kron(eye(32),y),eye(2));

A = x1*z2*y3; B = z4*x5*y6;

za = kron(eye(64),z);
xa = kron(eye(64),x);

for delta = 1:1e5:1e6
    LHS = alpha*A*B;
    RHS = (delta)*((1*eye(128) - za)/2) + abs(alpha)*((1*eye(128) + za)/2) + sqrt( abs(alpha)*delta/2 )*( sign(alpha)*A - B)*xa;
    
    abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.2948

%% P1B1-OT: Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

a = 3;

y1 = kron(y,eye(32));
z1 = kron(z,eye(32));

y2 = kron(kron(eye(2),y),eye(16));
z2 = kron(kron(eye(2),z),eye(16));

x3 = kron(kron(eye(4),x),eye(8));
y3 = kron(kron(eye(4),y),eye(8));

x4 = kron(kron(eye(8),x),eye(4));

y5 = kron(kron(eye(16),y),eye(2));

za = kron(eye(32),z);
xa = kron(eye(32),x);

I_size = 64;

for delta = 1:1e3:1e4
    LHS = a*y1*z2*x3*x4*y5;
    RHS = (delta*eye(I_size) - (delta^(2/3))*(a^(1/3))*y5)*((1*eye(I_size) - za)/2) ...
        + ((delta^(2/3))*(a^(1/3))/sqrt(2))*(-y1*z2*x3 + x4)*xa ...
        + ((delta^(1/3))*(a^(2/3))/2)*(-y1*z2*x3 + x4)^2 + (a/2)*((y1*z2*x3)^2 + (x4)^2)*y5;
    
    abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.4142

%% P1B1-CBBK Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

a = 3;
x1 = kron(x,eye(32));
z2 = kron(kron(eye(2),z),eye(16));
y3 = kron(kron(eye(4),y),eye(8));
y4 = kron(kron(eye(8),y),eye(4));
x5 = kron(kron(eye(16),x),eye(2));

z1 = kron(z,eye(32));
x3 = kron(kron(eye(4),x),eye(8));

za = kron(eye(32),z);
xa = kron(eye(32),x);

s1 = x1; s2 = z2; s3 = y3; s4 = y4; s5 = x5;

I_size = 64;

for Delta = 1:1e9:1e10
LHS = a*s1*s2*s3*s4*s5 - z1*x3;
RHS = (Delta*eye(I_size) + ((a/2)^(1/3))*(Delta^(1/2))*s5)*((1*eye(I_size) - za)/2) ...
    + ((a/2)^(1/3))*(Delta^(3/4))*(sign(a)*s1*s2*s3 + s4)*xa ...
    + ((a^(2/3))/2)*( (sign(a)^2) + 1 )*((2^(1/3))*(Delta^(1/2))*eye(I_size) - (a^(1/3))*s5 - (sign(a)^2)*((2*a)^(2/3))*eye(I_size))*((1*eye(I_size) + za)/2) ...
    + sign(a)*((a^(2/3))*(2^(1/3))*(Delta^(1/2)) - 4*((a/2)^(4/3)))*s1*s2*s3*s4 - z1*x3;

abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3401

%% PSD-CN: Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

a_1 = 3;

x1 = kron(kron(eye(1),x),eye(32));
y2 = kron(kron(eye(2),y),eye(16));
z3 = kron(kron(eye(4),z),eye(8));
y4 = kron(kron(eye(8),y),eye(4));

H_1j = x1*y2;
H_2j = z3*y4;

za_11 = kron(kron(eye(16),z),eye(2));
xa_11 = kron(kron(eye(16),x),eye(2));

za_1 = kron(kron(eye(32),z),eye(1));

I_size = 64;

R = 1;
C = 1;

for delta = 1:1e5:1e6
    beta = sqrt((a_1*delta)/(2*R));
    alpha = delta/(2*C);
    
    LHS = a_1*H_1j*H_2j;
    RHS = alpha*( (eye(I_size) - za_11*za_1) ) ...
        + alpha*( (eye(I_size) - za_1) ) + alpha*( (eye(I_size) - za_1*za_1) ) ...
        + beta*( (H_1j - H_2j)*xa_11 ) + a_1*eye(I_size);
    abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.4142

%% P(3->2)-CBBK2: Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

a_1 = 3;
a_2 = 5;
x1 = kron(x,eye(32));
y1 = kron(y,eye(32));
z1 = kron(z,eye(32));

x2 = kron(kron(eye(2),x),eye(16));
y2 = kron(kron(eye(2),y),eye(16));
z2 = kron(kron(eye(2),z),eye(16));

y3 = kron(kron(eye(4),y),eye(8));

z4 = kron(kron(eye(8),z),eye(4));

za1 = kron(kron(eye(16),z),eye(2));
xa1 = kron(kron(eye(16),x),eye(2));

za2 = kron(eye(32),z);
xa2 = kron(eye(32),x);
I_size = 64;

s_11_12 = 0; s_12_12 = 0;
if ( ~isequal(x1*y1 - y1*x1,0*eye(I_size)) && isequal(z2*x2 - x2*z2, 0*eye(I_size)) ) || ( isequal(x1*y1 - y1*x1,0*eye(I_size)) && ~isequal(z2*x2 - x2*z2, 0*eye(I_size)) )
    s_11_12 = 1;
else
    s_11_12 = 0;
end

if ~isequal(x1*x2 - x2*x1, 0*eye(I_size)) || ~isequal(y1*z2 - z2*y1,0*eye(I_size))
    s_12_12 = 1;
else
    s_12_12 = 0;
end

s_1_12 = s_11_12 + s_12_12;

if ( ~isequal(x1*y1 - y1*x1,0*eye(I_size)) && ~isequal(z2*x2 - x2*z2, 0*eye(I_size)) )
    s_2_12 = 1;
else
    s_2_12 = 0;
end


for delta = 1:1e6:1e7
    LHS = a_1*x1*z2*y3 + a_2*y1*x2*z4 - z1*x2;
    RHS = (delta*eye(I_size) + (( abs(a_1)/2 )^(1/3))*(delta^(1/2))*y3)*((eye(I_size) - za1)/2) ...
        + ((abs(a_1)/2)^(1/3))*(delta^(3/4))*(sign(a_1)*x1 + z2)*xa1 ...
        + ((abs(a_1)/2)^(2/3))*(delta^(1/2))*(sign(a_1)*x1 + z2)^2 ...
        - (abs(a_1)/2)*(sign(a_1)^2 + 1)*y3 ...
        - ((abs(a_1)/2)^(4/3))*(sign(a_1)*x1 + z2)^4 ...
        + (delta*eye(I_size) + (( abs(a_2)/2 )^(1/3))*(delta^(1/2))*z4)*((eye(I_size) - za2)/2) ...
        + ((abs(a_2)/2)^(1/3))*(delta^(3/4))*(sign(a_2)*y1 + x2)*xa2 ...
        + ((abs(a_2)/2)^(2/3))*(delta^(1/2))*(sign(a_2)*y1 + x2)^2 ...
        - (abs(a_2)/2)*(sign(a_2)^2 + 1)*z4 ...
        - ((abs(a_2)/2)^(4/3))*(sign(a_2)*y1 + x2)^4 ...
        - s_1_12*(sign(a_1)^2)*(sign(a_2)^2)*((1/2)^(4/3))*(abs(a_1)^(2/3))*(abs(a_2)^(2/3))*eye(I_size) ...
        - s_2_12*( 2*(sign(a_1)^2)*(sign(a_2)^2)*((1/2)^(4/3))*(abs(a_1)^(2/3))*(abs(a_2)^(2/3))*eye(I_size) ...
        - 2*sign(a_1)*sign(a_2)*((abs(a_1)/2)^(2/3))*((abs(a_2)/2)^(2/3))*x1*y1*z2*x2 ) - z1*x2;
    
    abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3748

%% PD-JF: Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

a1 = 3;
x1 = kron(kron(eye(1),x),eye(512));
y2 = kron(kron(eye(2),y),eye(256));
z3 = kron(kron(eye(4),z),eye(128));

a2 = 5;
y1 = kron(kron(eye(1),y),eye(512));
z2 = kron(kron(eye(2),z),eye(256));
y4 = kron(kron(eye(8),y),eye(64));

za11 = kron(kron(eye(16),z),eye(32));
za12 = kron(kron(eye(32),z),eye(16));
za13 = kron(kron(eye(64),z),eye(8));

xa11 = kron(kron(eye(16),x),eye(32));
xa12 = kron(kron(eye(32),x),eye(16));
xa13 = kron(kron(eye(64),x),eye(8));

za21 = kron(kron(eye(128),z),eye(4));
za22 = kron(kron(eye(256),z),eye(2));
za23 = kron(kron(eye(512),z),eye(1));

xa21 = kron(kron(eye(128),x),eye(4));
xa22 = kron(kron(eye(256),x),eye(2));
xa23 = kron(kron(eye(512),z),eye(1));

I_size = 1024;

for delta = 1:1e7:1e8
    LHS = a1*x1*y2*z3 + a2*y1*z2*y4;
    RHS = (1/2)*(6*eye(I_size) - za11*za12 - za11*za13 - za12*za13 - za21*za22 - za21*za23 - za22*za23) ...
        + (1/delta)*(a1*x1*xa11 + y2*xa12 + z3*xa13 + a2*y1*xa21 + z2*xa22 + y4*xa23) - (a1 + a2)*eye(I_size);
    abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3646

%% PD-BFBD: Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

z_4i_1_j  = sparse(kron(kron(speye(2^0),z),speye(2^13)));
z_4i_2_j  = sparse(kron(kron(speye(2^1),z),speye(2^12)));
z_4i_3_j  = sparse(kron(kron(speye(2^2),z),speye(2^11)));
z_4i_4_j  = sparse(kron(kron(speye(2^3),z),speye(2^10)));

x_4i_3_j  = sparse(kron(kron(speye(2^2),x),speye(2^11)));
x_4i_4_j  = sparse(kron(kron(speye(2^3),x),speye(2^10)));
x_4i_6_j  = sparse(kron(kron(speye(2^4),x),speye(2^9)));
x_4i_4_j1 = sparse(kron(kron(speye(2^5),x),speye(2^8)));

x_8i_4_j  = sparse(kron(kron(speye(2^6),x),speye(2^7)));
z_8i_4_j  = sparse(kron(kron(speye(2^6),z),speye(2^7)));

x_8i_6_j  = sparse(kron(kron(speye(2^7),x),speye(2^6)));
z_8i_6_j  = sparse(kron(kron(speye(2^7),z),speye(2^6)));

x_8i_3_j1 = sparse(kron(kron(speye(2^8),x),speye(2^5)));
z_8i_3_j1 = sparse(kron(kron(speye(2^8),z),speye(2^5)));

x_8i_5_j1 = sparse(kron(kron(speye(2^9),x),speye(2^4)));
z_8i_5_j1 = sparse(kron(kron(speye(2^9),z),speye(2^4)));

x_8i_1_j  = sparse(kron(kron(speye(2^10),x),speye(2^3)));
z_8i_1_j  = sparse(kron(kron(speye(2^10),z),speye(2^3)));

x_8i_2_j  = sparse(kron(kron(speye(2^11),x),speye(2^2)));
z_8i_2_j  = sparse(kron(kron(speye(2^11),z),speye(2^2)));

x_8i_7_j  = sparse(kron(kron(speye(2^12),x),speye(2^1)));
z_8i_7_j  = sparse(kron(kron(speye(2^12),z),speye(2^1)));

x_8i_8_j  = sparse(kron(kron(speye(2^13),x),speye(2^0)));
z_8i_8_j  = sparse(kron(kron(speye(2^13),z),speye(2^0)));


LHS = sparse(-(z_4i_1_j*z_4i_2_j*z_4i_3_j*z_4i_4_j + x_4i_3_j*x_4i_4_j*x_4i_6_j*x_4i_4_j1));

for delta = 1:1e7:1e8
RHS = sparse(-(x_8i_4_j*x_8i_6_j + x_8i_3_j1*x_8i_5_j1 + z_8i_4_j*z_8i_3_j1 + z_8i_6_j*z_8i_5_j1 ...
    + (1/delta)*( x_8i_1_j*x_8i_3_j1 + x_8i_2_j*x_8i_4_j + x_8i_5_j1*x_8i_7_j + x_8i_6_j*x_8i_8_j ...
    + z_8i_1_j*z_8i_3_j1 + z_8i_2_j*z_8i_4_j + z_8i_5_j1*z_8i_7_j + z_8i_6_j*z_8i_8_j ))) + 0.828427124746191*speye(2^14);

abs(eigs(LHS, 1, 'smallestreal')-eigs(RHS, 1, 'smallestreal'))
end

%% PD-CK: Example
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

a_1 = 0.1;
a_2 = 0.2;

x1 = kron(kron(eye(1),x),eye(1024));
x2 = kron(kron(eye(2),x),eye(512));
x3 = kron(kron(eye(4),x),eye(256));
y4 = kron(kron(eye(8),y),eye(128));
z5 = kron(kron(eye(16),z),eye(64));

za_11 = kron(kron(eye(32),z),eye(32));
xa_11 = kron(kron(eye(32),x),eye(32));

za_12 = kron(kron(eye(64),z),eye(16));
xa_12 = kron(kron(eye(64),x),eye(16));

za_13 = kron(kron(eye(128),z),eye(8));
xa_13 = kron(kron(eye(128),x),eye(8));


za_21 = kron(kron(eye(256),z),eye(4));
xa_21 = kron(kron(eye(256),x),eye(4));

za_22 = kron(kron(eye(512),z),eye(2));
xa_22 = kron(kron(eye(512),x),eye(2));

za_23 = kron(kron(eye(1024),z),eye(1));
xa_23 = kron(kron(eye(1024),x),eye(1));

I_size = 2048;

for delta = 1:1e5:1e6
    mu_1 = (a_1)^(1/3);
    mu_2 = (a_2)^(1/3);
    
    LHS = a_1*x1*x2*x3 + a_2*x2*y4*z5;
    RHS = (delta/4)*(3*eye(I_size) - za_11*za_12 - za_12*za_13 - za_11*za_13) ...
        + (delta/4)*(3*eye(I_size) - za_21*za_22 - za_22*za_23 - za_11*za_13) ...
        + mu_1*(x1*xa_11 + x2*xa_12 + x3*xa_13) ...
        + mu_2*(x2*xa_21 + y4*xa_22 + z5*xa_23) - (a_1 + a_2)*eye(I_size);
    abs(min(eig(LHS))-min(eig(RHS)))
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:4:end); EL = EL(:,1:4:end); VR = VR(:,1:4:end); ER = ER(:,1:4:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3669

%% III.D 15-term, 5-variable, degree-4 function
%% blue
b=dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);
LHS=min(reshape(5*b1.*b2.*b3.*b4 - 3*b1.*b2.*b3 - b1.*b2.*b4 - b1.*b3.*b4 - 2*b2.*b3.*b4,2,[]));
RHS=min(reshape(b1.*b2 + b1.*b3 + 3*b1.*b4 + 2*b2.*b4 + 2*b3.*b4 - 5*b1.*ba1 - 4*b2.*ba1 - 4*b3.*ba1 - 6*b4.*ba1 + 8*ba1,2,[]))
isequal(LHS,RHS)
%% red
b=dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b5=b(:,4);ba2=b(:,5);
LHS=min(reshape(4*b1.*b2.*b3.*b5 - 5*b1.*b2.*b5 - b1.*b3.*b5 - b2.*b3.*b5,2,[]));
RHS=min(reshape(-3*b1 + 6*b2 - 3*b3 + 5*b5 - 5*b1.*b2 + 3*b1.*b3 - 5*b1.*b5 - b2.*b3 - b3.*b5 + 8*b1.*ba2 - 6*b2.*ba2 + 4*b3.*ba2 - 5*b5.*ba2 - 3*ba2 + 3,2,[]));
isequal(LHS,RHS)
%% purple
b=dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b4=b(:,3);b5=b(:,4);ba3=b(:,5);
LHS=min(reshape(3*b1.*b2.*b4.*b5 - b1.*b4.*b5 - 4*b2.*b4.*b5,2,[]));
RHS=min(reshape(b1 + 4*b2 + 3*b1.*b2 - b1.*b4 - b1.*b5 - 4*b2.*b4 - 4*b2.*b5 - 4*b1.*ba3 - 7*b2.*ba3 + 5*b4.*ba3 + 5*b5.*ba3 + 3*ba3,2,[]));
isequal(LHS,RHS)
%% green
b=dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b3=b(:,2);b4=b(:,3);b5=b(:,4);ba4=b(:,5);
LHS=min(reshape(2*b1.*b3.*b4.*b5 - 3*b3.*b4.*b5,2,[]));
RHS=min(reshape(2*b1.*ba4 - 3*b3.*ba4 - 3*b4.*ba4 - 3*b5.*ba4 + 6*ba4,2,[]));
isequal(LHS,RHS)
%% orange
b=dec2bin(2^5-1:-1:0)-'0';
b2=b(:,1);b2=b(:,2);b4=b(:,3);b5=b(:,4);ba5=b(:,5);
LHS=min(reshape(b2.*b3.*b4.*b5,2,[]));
RHS=min(reshape(b2.*b3 + b2.*b4 + b2.*b5 + b3.*b4 + b3.*b5 + b4.*b5 + 3*ba5 - 2*b2.*ba5 - 2*b3.*ba5 - 2*b4.*ba5 - 2*b5.*ba5,2,[]));
isequal(LHS,RHS)
