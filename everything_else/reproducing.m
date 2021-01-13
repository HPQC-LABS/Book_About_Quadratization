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
[eig(LHS) eig(RHS)];

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
