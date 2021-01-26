%% Pg. 1, Eqs 1-2

b= dec2bin(2^4-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);
LHS=b1.*b2 + b2.*b3 + b3.*b4 - 4*b1.*b2.*b3;
RHS=b1.*b2 + b2.*b3 + b3.*b4 + 4*b1 - 4*b1.*b2 - 4*b1.*b3;

%% Pg. 1, Eqs 4-5

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

for col = 1:1:length(VL) % compare eigenvectors and eigenvalues
V_diff(1,col) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(1,col) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

max(E_diff); % 3.1e-15, full energy spectrum reproduced
min(V_diff); % 1.28,      no states reproduced

%% Pg. 6, Eqs 8 and 10

LHS=b1.*b2 + b2.*b3 + b3.*b4 - 4*b1.*b2.*b3;
RHS=b1.*b2 + b2.*b3 + b3.*b4 + 4*b1 - 4*b1.*b2 - 4*b1.*b3;

%% Pg. 10, Eqs 25-26

b= dec2bin(2^7-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);b5=b(:,5);b6=b(:,6);ba=b(:,7);

LHS=(-2)*b1.*b2.*b3.*b4.*b5.*b6 + b5.*b6;
RHS=2*(5*ba - b1.*ba - b2.*ba - b3.*ba - b4.*ba - b5.*ba - b6.*ba) + b5.*b6;

%% Pg. 10, Eq 27

LHS=min(reshape(-b1.*b2.*b3.*b4.*b5.*b6,2,[]));
RHS=min(reshape((6 - 1 - b1 - b2 - b3 - b4 - b5 - b6).*ba,2,[]));
isequal(LHS,RHS);

%% PTR_BG Pg: 17, Eqn 50

b = dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);ba2=b(:,6);
LHS = b1.*b2.*b3.*b4;
RHS = ba1.*(2+ b1 -b2 - b3 - b4) + ba2.*(1 + b2 - b3 - b4) + b3.*b4;

%% PTR_Ishikawa Pg. 18, Eqn 52

b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = ba.*(3 - 2*b1 - 2*b2 - 2*b3 - 2*b4) + b1.*b2 + b1.*b3 + b1.*b4 + b2.*b3 + b2.*b4 +  b3.*b4;

%% PTR-BCR-1: Pg. 20, Eqns 59

b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = 1/2*(b1 + b2 + b3 + b4 - 2*ba1).*(b1 + b2 + b3 + b4 - 2*ba1 - 1);

%% PTR-BCR-3 Pg. 22, Eqns 65

b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = 1/2*(b1 + b2 + b3 + b4 - 2*ba).*(b1 + b2 + b3 + b4 - 2*ba - 1);

%% SFR-ABCG-4 Pg. __, Eqns __

b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = b4.*b2 + ba.*(b3-1);

%% NTR-GBP: -b1b2b3 = min_ba(ba - b1 + b2 + b3 - b1b2 - b1b3 + b1)

b= dec2bin(2^4-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);ba=b(:,4);
LHS = min(reshape(-1*b1.*b2.*b3,2,[]));
RHSa = min(reshape(ba.*(-b1 + b2 + b3) - b1.*b2 - b1.*b3 + b1,2,[]));
RHSb = min(reshape(ba.*(-b2 + b1 + b3) - b1.*b2 - b2.*b3 + b2,2,[]));
RHSc = min(reshape(ba.*(-b3 + b1 + b2) - b2.*b3 - b1.*b3 + b3,2,[]));
isequal(LHS,RHSa,RHSb,RHSc);


%%

z=[1 0; 0 -1];

%% PTR-RBL, k=4, middle of LHZ

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

%% RBS-Rosenberg: b1b2b3 = min_ba(b1ba + b2b3 - 2*b2ba -2*b3ba + 3*ba) (Eq 150)

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

%% FGBZ for negative: -b1b2b3 - b1b2b4 = min_ba((1 - b1b2 - b3)ba1 + (1 - b1b2 - b4)ba1)  (Eqs 153-154)

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

%% FGBZ for positive: b1b2b3 + b1b2b4 = min_ba(2*ba1b1 + (1 - ba1)b2b3 + (1 - ba1)b2b4)  (Eqs 158-159)

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

%% Pairwise Covers: 5*b1b2b3b4 + 4*b1b2b4 - 3*b1b3 - 2*b2b3b4
%%                  = min_ba(5*ba1ba2 + 4*b1ba2 - 3*ba1 - 2*b3ba2 + 8*(ba1(3 - 2*b1 - 2*b3) + b1b3) + 11*(ba2(3 - 2*b2 - 2*b4) + b2b4))

b=dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);ba2=b(:,6);
LHS=min(reshape(5*b1.*b2.*b3.*b4 + 4*b1.*b2.*b4 - 3*b1.*b3 - 2*b2.*b3.*b4,4,[]));
RHS=min(reshape(5*ba1.*ba2 + 4*b1.*ba2 - 3*ba1 - 2*b3.*ba2 + 8*(ba1.*(3 - 2*b1 - 2*b3) + b1.*b3) + 11*(ba2.*(3 - 2*b2 - 2*b4) + b2.*b4),4,[]));
isequal(LHS,RHS);

%% PTR-KZ: b1b2b3 = min_ba(1 − (ba + b1 + b2 + b3) + ba (b1 + b2 + b3) + b1b2 + b1b3 + b2b3)

LHS=min(reshape(b1.*b2.*b3,2,[]));
RHS=min(reshape(1 - (b4+b1+b2+b3) + b4.*(b1+b2+b3) +b1.*b2+b1.*b3+b2.*b3,2,[]));
isequal(LHS,RHS);

%% PTR-GBP: b1b2b3 = min_ba(ba - b2ba - b3ba +b1ba +b2b3)

LHS = min(reshape(b1.*b2.*b3,2,[]));
RHSa = min(reshape(b4 - b2.*b4 - b3.*b4 + b1.*b4 +b2.*b3,2,[]));
RHSb = min(reshape(b4 - b1.*b4 - b3.*b4 + b2.*b4 + b1.*b3,2,[]));
RHSc = min(reshape(b4 - b1.*b4 - b2.*b4 + b3.*b4 + b1.*b2,2,[]));
isequal(LHS,RHSa,RHSb,RHSc);

%% example eq. (79)
LHS = min(reshape(b1.*b2.*b3 + b1.*b3 - b2,2,[]));
RHS = min(reshape((b4 - b1.*b4 - b3.*b4 + b2.*b4 + 2*b1.*b3)-b2,2,[]));
isequal(LHS,RHS);

%% PTR-YXKK: b1b2b3 + b1b3 − b2 -> b1b2 + ba1 + 2ba2 + 2(1 − b1)(1 − ba1) + 2(1 − b2)(1 − ba1) + 2(1 − b3)(1 − ba2) + 2(1 − ba2)(1 − ba1) − 2(1 − b3) − 1 + b1b3 − b2
b=dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1); b2=b(:,2); b3=b(:,3); ba1=b(:,4); ba2=b(:,5);
LHS=min(reshape(b1.*b2.*b3 + b1.*b3 - b2,4,[]));
RHS=min(reshape(b1.*b2 + ba1 + 2*ba2 + 2*(1-b1).*(1-ba1) + 2*(1-b2).*(1-ba1) + 2*(1-b3).*(1-ba2) + 2*(1-ba2).*(1-ba1) - 2*(1-b3) - 1 + b1.*b3 - b2,4,[]));
isequal(LHS,RHS);

%% SFR-BCR-1: (b1b2 + b1b3 + b1b4 + b2b3 + b2b4 + b3b4) − 3(b1b2b3 + b1b2b4 + b1b3b4 + b2b3b4) + 6b1b2b3b4 -> (−3 + b1 + b2 + b3 + b4 − ba1 + 3ba2)^2
b= dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1); b2=b(:,2); b3=b(:,3); b4=b(:,4); ba1=b(:,5); ba2=b(:,6);
LHS= min(reshape( (b1.*b2 + b1.*b3 + b1.*b4 + b2.*b3 + b2.*b4 + b3.*b4) - 3*(b1.*b2.*b3 + b1.*b2.*b4 + b1.*b3.*b4 + b2.*b3.*b4) + 6*b1.*b2.*b3.*b4 ,4,[]));
RHS= min(reshape( (-3 + b1 + b2 + b3 + b4 - ba1 + 3*ba2).^2 ,4,[]));
isequal(LHS, RHS);

%% SFR-BCR-2: (b1b2 + b1b3 + b1b4 + b2b3 + b2b4 + b3b4)−3(b1b2b3 + b1b2b4 + b1b3b4 + b2b3b4) + 6b1b2b3b4 -> (1 − b1 − b2 − b3 − b4 − ba1 + 3ba2)^2
b= dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1); b2=b(:,2); b3=b(:,3); b4=b(:,4); ba1=b(:,5); ba2=b(:,6);
LHS= min(reshape( (b1.*b2 + b1.*b3 + b1.*b4 + b2.*b3 + b2.*b4 + b3.*b4) - 3*(b1.*b2.*b3 + b1.*b2.*b4 + b1.*b3.*b4 + b2.*b3.*b4) + 6*b1.*b2.*b3.*b4 ,4,[]));
RHS= min(reshape( (1 - b1 - b2 - b3 - b4 - ba1 + 3*ba2).^2 ,4,[]));
isequal(LHS, RHS);

%% ZZZ-TI-CBBK

z = [1 0; 0 -1];
x = [0 1;1 0];
zi = kron(z,eye(2^3));
zj = kron(kron(eye(2),z),eye(2^2));
zk = kron(kron(eye(2^2),z),eye(2));
za = kron(eye(2^3),z);
xa = kron(eye(2^3),x);

alpha = 1;

for delta=1:100:1000;
alpha_I = (delta + ((alpha/6)^(2/5))*(delta^(3/5)))/2;
alpha_zi = -(((alpha*7/6)+((alpha/6)^(3/5))*(delta^(2/5)))-((alpha*delta^4)/6)^(1/5))/2;
alpha_zj = alpha_zi;
alpha_zk = alpha_zi;
alpha_za = (delta - ((alpha/6)^(2/5))*(delta^(3/5)))/2;
alpha_xa = ((alpha*delta^4)/6)^(1/5);
alpha_zzia = -(((alpha*7/6)+((alpha/6)^(3/5))*(delta^(2/5)))+((alpha*delta^4)/6)^(1/5))/2;
alpha_zzja = alpha_zzia;
alpha_zzka = alpha_zzja;

LHS = alpha*zi*zj*zk;
RHS = alpha_I + alpha_zi*zi + alpha_zj*zj + alpha_zk*zk + alpha_za*za + alpha_xa*xa + alpha_zzia*zi*za + alpha_zzja*zj*za + alpha_zzka*zk*za;
min(eig(LHS))-min(eig(RHS));
end

%% NP - SJ (4z5 - 3x1 + 2*z1*y2*x5 + 9*x1*x2*x3*x4 - x1*y2*z3*x5 -> 9*xa1 + 4*za2*z5 - 3*za3*x1 - za3*xa2 + 2*xa3*x5)
x = [0 1;1 0];
y = [0 -1i;1i 0];
z = [1 0;0 -1];

x1 = kron(x,eye(128));x2 = kron(kron(eye(2),x),eye(64));x3 = kron(kron(eye(4),x),eye(32));x4 = kron(kron(eye(8),x),eye(16));x5 = kron(kron(eye(16),x),eye(8));
xa1 = kron(kron(eye(32),x),eye(4));xa2 = kron(kron(eye(64),x),eye(2));xa3 = kron(eye(128),x);
y2 = kron(kron(eye(2),y),eye(64));
z1 = kron(z,eye(128));z3 = kron(kron(eye(4),z),eye(32));z5 = kron(kron(eye(16),z),eye(8));
za2 = kron(kron(eye(64),z),eye(2));
za3 = kron(eye(128),z);

LHS = 4*z5 - 3*x1 + 2*z1*y2*x5 + 9*x1*x2*x3*x4 - x1*y2*z3*x5;
RHS = 9*xa1 + 4*za2*z5 - 3*za3*x1 - za3*xa2 + 2*xa3*x5;

max(eig(LHS)-eig(RHS))<1e-13; % gives 1.

%% P(3->2)-DC1

x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
x1 = kron(x,eye(32)); x2 = kron(kron(eye(2),x),eye(16));
xa1 = kron(kron(eye(8),x),eye(4)); xa2 = kron(kron(eye(16),x),eye(2)); xa3 = kron(eye(32),x);
y3 = kron(kron(eye(4),y),eye(8));
z1 = kron(z,eye(32)); z2 = kron(kron(eye(2),z),eye(16));
za1 = kron(kron(eye(8),z),eye(4)); za2 = kron(kron(eye(16),z),eye(2)); za3 = kron(eye(32),z);

for delta = 1:1e2:1e3
alpha = 1/(8*delta);
alpha_ss = 1/(6*(delta)^(1/3));
alpha_sx = -1/(6*(delta)^(2/3));
alpha_zz = -1/(24*delta);

LHS = (x1 + 3*x2)*(z1 + z2)*y3;
RHS = alpha*eye(64) + alpha_ss*(6*x1*x2 + 2*z1*z2) + 12*alpha_ss*eye(64) + alpha_sx*(x1 + 3*x2)*xa1 + alpha_sx*(z1 + z2)*xa2 + alpha_sx*y3*xa3 + alpha_zz*(za1*za2 + za1*za3 + za2*za3);

min(eig(LHS))-min(eig(RHS))
end

%% P(3->2)-DC2

x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];
x1 = kron(x,eye(32)); x2 = kron(kron(eye(2),x),eye(16)); x3 = kron(kron(eye(4),x),eye(8)); x4 = kron(kron(eye(8),x),eye(4));
xa1 = kron(kron(eye(16),x),eye(2)); xa2 = kron(eye(32),x);
y3 = kron(kron(eye(4),y),eye(8)); y4 = kron(kron(eye(8),y),eye(4));
z1 = kron(z,eye(32)); z2 = kron(kron(eye(2),z),eye(16));
za1 = kron(kron(eye(16),z),eye(2)); za2 = kron(eye(32),z);

for delta = 1:1e2:1e3
alpha = -1/(2*delta);
alpha_z = (1/(4*(delta)^(2/3))) - 1;
alpha_y = (1/(4*(delta)^(2/3))) - 1;
alpha_12_zx = 1/(delta^(2/3));
alpha_zx = 1/(delta^(2/3));
alpha_11_xx = 1/(delta^(1/3));
alpha_xx = 1/(delta^(2/3));
alpha_yz = 1/(4*delta^(2/3));

LHS = x1*z2*y3 - 3*x1*x2*y4 + z1*x2;
RHS = alpha + alpha_z*(za1 + za2) + alpha_y*(y3 + y4) + alpha_12_zx*z1*x2 + alpha_zx*z2*xa1 + alpha_11_xx*x1*x2 + alpha_xx*(x1*xa1 + x1*xa2 + x2*xa2) + alpha_yz*(y3*za1 + y4*za2) + 0.482582936542672*eye(64);

min(eig(LHS))-min(eig(RHS))
end

%% P(3->2)-KKR, example.

x = [0 1;1 0]; y = [0 -1i;1i 0]; z = [1 0;0 -1];
z1 = kron(z,eye(512));z2 = kron(kron(eye(2),z),eye(256));
x1 = kron(x,eye(512));x2 = kron(kron(eye(2),x),eye(256));
y3 = kron(kron(eye(4),y),eye(128));y4 = kron(kron(eye(8),y),eye(64));
xa_11 = kron(kron(eye(16),x),eye(32)); xa_12 = kron(kron(eye(32),x),eye(16)); xa_13 = kron(kron(eye(64),x),eye(8));
xa_21 = kron(kron(eye(128),x),eye(4)); xa_22 = kron(kron(eye(256),x),eye(2)); xa_23 = kron(eye(512),x);
za_11 = kron(kron(eye(16),z),eye(32)); za_12 = kron(kron(eye(32),z),eye(16)); za_13 = kron(kron(eye(64),z),eye(8));
za_21 = kron(kron(eye(128),z),eye(4)); za_22 = kron(kron(eye(256),z),eye(2)); za_23 = kron(eye(512),z);

for delta = 1:1e8:1e9

alpha = -(1/8)*(delta);
alpha_ss = -(1/6)*(delta)^(1/3);
alpha_sx = (1/6)*(delta)^(2/3);
alpha_zz = (1/24)*(delta);

LHS = z1*x2 - x1*z2*y3 - 3*x1*x2*y4;
RHS = z1*x2 - 4*alpha*eye(1024) - 12*alpha_ss*eye(1024) - alpha_sx*(x1*xa_11 + z2*xa_12 + y3*xa_13) - 3*alpha_sx*(x1*xa_21 + x2*xa_22 + y4*xa_23) - alpha_zz*(za_11*za_12 + za_11*za_13 + za_12*za_13) - 3*alpha_zz*(za_21*za_22 + za_21*za_23 + za_22*za_23);

abs(min(eig(LHS))-min(eig(RHS)))
end

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
