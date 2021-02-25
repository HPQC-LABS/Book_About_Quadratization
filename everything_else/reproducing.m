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

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:8:end); EL = EL(:,1:8:end); VR = VR(:,1:8:end); ER = ER(:,1:8:end); % reduce repeats due to auxiliary qubits

for col = 1:1:size(VL,2) % compare eigenvectors and eigenvalues
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

max(E_diff); % 1.1546e-14
min(V_diff); % 1.3632

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

for delta = 1:1e8:1e9

alpha = -(1/8)*(delta);
alpha_ss = -(1/6)*(delta)^(1/3);
alpha_sx = (1/6)*(delta)^(2/3);
alpha_zz = (1/24)*(delta);

LHS = z1*x2 - x1*z2*y3 - 3*x1*x2*y4;
RHS = z1*x2 - 4*alpha*eye(1024) - 12*alpha_ss*eye(1024) - alpha_sx*(x1*xa_11 + z2*xa_12 + y3*xa_13) - 3*alpha_sx*(x1*xa_21 + x2*xa_22 + y4*xa_23) - alpha_zz*(za_11*za_12 + za_11*za_13 + za_12*za_13) - 3*alpha_zz*(za_21*za_22 + za_21*za_23 + za_22*za_23);

abs(min(eig(LHS))-min(eig(RHS)));
end

[VL, EL] = eig(LHS); [VR, ER] = eig(RHS);
[DL, indL] = sort(diag(EL)); [DR, indR] = sort(diag(ER));
VL = VL(:,indL); EL = EL(indL,indL); VR = VR(:,indR); ER = ER(indR,indR);
VL = VL(:,1:64:end); EL = EL(:,1:64:end); VR = VR(:,1:64:end); ER = ER(:,1:64:end);

for col = 1:1:size(VL,2)
V_diff(col,1) = sqrt(dot((VL(:,col)-VR(:,col)),(VL(:,col)-VR(:,col))));
E_diff(col,1) = sqrt(dot((EL(:,col)-ER(:,col)),(EL(:,col)-ER(:,col))));
end

min(V_diff); % 1.3954

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

%% P1B1-OT (Example, Incomplete)
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

alpha = 1;
x1 = kron(x,eye(32));
z2 = kron(kron(eye(2),z),eye(16));
y3 = kron(kron(eye(4),y),eye(8));
y4 = kron(kron(eye(8),y),eye(4));
z5 = kron(kron(eye(16),z),eye(2));


s1 = x1; s2 = z2; s3 = y3; s4 = y4; s5 = z5;

za = kron(eye(32),z);
xa = kron(eye(32),x);
r = 2/3;

delta_array = [];
ans_array = [];


for delta = 1:1e6:1e7
   LHS = alpha*x1*z2*y3;
   RHS = (delta*eye(64) - ((alpha/2)^(1/3))*(delta^(2-2*r))*s3)*((1*eye(64) - za)/2) ...
       + ((alpha/2)^(1/3))*((delta^r)/(2^(1/2)))*(-s1 + s2)*xa ... 
       + (1/2*delta)*(alpha/2)^(2/3)*((delta^r)*(-s2 + s1))^2 + (alpha/4)*delta^(-r)*(s1^2 + s2^2)*s3;
   
   delta_array = [delta_array, delta];
   ans_array = [ans_array, abs(min(eig(LHS))-min(eig(RHS)))];
end
[delta_array ; ans_array];

%% P1B1-CBBK (Example, Incomplete)
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

alpha = 3;
x1 = kron(x,eye(16));
y1 = kron(y,eye(16));

x2 = kron(kron(eye(2),x),eye(8));
z2 = kron(kron(eye(2),z),eye(8));

y3 = kron(kron(eye(4),y),eye(4));
y4 = kron(kron(eye(8),y),eye(2));


s1 = x1; s2 = z2; s3 = y3; s4 = y4;

za = kron(eye(16),z);
xa = kron(eye(16),x);

delta_array = [];
ans_array = [];

for delta = 1:1e5:1e6
   LHS = alpha*s1*s2*s3*s4 + y1*x2;
   RHS = (delta*eye(32) + ((alpha/2)^(3/2))*(delta^(1/2))*s4)*((1*eye(32) - za)/2) ...
       - ((alpha^(2/3))/2)*(1 + (sign(alpha)^2))*( ((2*alpha)^(2/3))*(sign(alpha)^2)*eye(32) + (alpha^(1/3))*s4 - (2^(1/3))*(delta^(1/2))*eye(32) )*((1*eye(32) + za)/2) ...
       + ((alpha/2)^(1/3))*(delta^(3/4))*(s1*s2 - sign(alpha)*s3)*xa + sign(alpha)*(2^(1/3))*(alpha^(2/3))*(delta^(1/2) + delta^(3/2))*(s1*s2*s3) + y1*x2;
   
   delta_array = [delta_array, delta];
   ans_array = [ans_array, min(eig(LHS))-min(eig(RHS))];
end
[delta_array ; ans_array];

%% PSD-CN (Incomplete Example)
x = [0 1 ; 1 0]; y = [0 -1i ; 1i 0]; z = [1 0 ; 0 -1];

a_1 = 2;

x1 = kron(eye(64),x);
y2 = kron(kron(eye(2),y),eye(32));
z3 = kron(kron(eye(4),z),eye(16));
y4 = kron(kron(eye(8),y),eye(8));

H_11 = x1*y2;
H_21 = z3*y4;

a_2 = 3;
y1 = kron(eye(64),y);
x2 = kron(kron(eye(2),x),eye(32));
x3 = kron(kron(eye(4),x),eye(16));
z4 = kron(kron(eye(8),z),eye(8));

H_12 = y1*x2;
H_22 = z3*y4;

za_11 = kron(kron(eye(16),z),eye(4));
xa_11 = kron(kron(eye(16),x),eye(4));
za_1 = kron(kron(eye(32),z),eye(2));

array = [];

for delta = 1:10:1000
alpha = delta/2;
alpha_1 = sqrt((a_1*delta)/(2*1));

LHS = a_1*H_11*H_21;
RHS_book = alpha*(eye(128) - za_11*za_1 + alpha_1*xa_11*(H_11 - H_21)) + alpha*(eye(128) - za_1 + (1 - za_1*za_1));
array = [array , abs(min(eig(LHS))-min(eig(RHS_book)))];
end

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
