n=6;
allCombos=dec2bin(0:2^n-1) -'0';

b1=allCombos(:,1);
b2=allCombos(:,2);
b3=allCombos(:,3);
b4=allCombos(:,4);
b5=allCombos(:,5);
ba=allCombos(:,6);
c=ones(2^n,1);

LHS=b1.*b2.*b3.*b4.*b5 + b1.*b2.*b3.*b4 + b1.*b3.*b4.*b5 + b1.*b2.*b4.*b5 + b1.*b2.*b3.*b5 + b2.*b3.*b4.*b5;

%%
%b1.*b2+ b1.*b3+ b1.*b4+ b1.*b5+ b2.*b3+ b2.*b4+ b2.*b5+ b3.*b4+ b3.*b5+ b4.*b5+ b1+ b2+ b3+ b4+ b5;
coeffsQ=ones(22,1);
[b1.*b2 b1.*b3 b1.*b4 b1.*b5 b1.*ba b2.*b3 b2.*b4 b2.*b5 b2.*ba b3.*b4 b3.*b5 b3.*ba b4.*b5 b4.*ba b5.*ba b1 b2 b3 b4 b5 ba c]*coeffsQ;

for i=1:1e9
coeffsQ=randi([-3 3],22,1);
RHS=[b1.*b2 b1.*b3 b1.*b4 b1.*b5 b1.*ba b2.*b3 b2.*b4 b2.*b5 b2.*ba b3.*b4 b3.*b5 b3.*ba b4.*b5 b4.*ba b5.*ba b1 b2 b3 b4 b5 ba c]*coeffsQ;
RHS(1:2:2^n-1) = min(RHS(1:2:2^n-1), RHS(2:2:2^n));
RHS(2:2:2^n) = RHS(1:2:2^n-1);
if LHS==RHS; disp('found_a_quadratization'); end
end

%%
RHS=[b1.*b2 b1.*b3 b1.*b4 b1.*b5 b1.*ba b2.*b3 b2.*b4 b2.*b5 b2.*ba b3.*b4 b3.*b5 b3.*ba b4.*b5 b4.*ba b5.*ba b1 b2 b3 b4 b5 ba c];

