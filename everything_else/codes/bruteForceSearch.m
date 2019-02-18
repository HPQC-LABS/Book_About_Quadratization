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
LHS = LHS(2:2:2^n);
%%
%b1.*b2+ b1.*b3+ b1.*b4+ b1.*b5+ b2.*b3+ b2.*b4+ b2.*b5+ b3.*b4+ b3.*b5+ b4.*b5+ b1+ b2+ b3+ b4+ b5;
coeffsQ=ones(22,1);
[b1.*b2 b1.*b3 b1.*b4 b1.*b5 b1.*ba b2.*b3 b2.*b4 b2.*b5 b2.*ba b3.*b4 b3.*b5 b3.*ba b4.*b5 b4.*ba b5.*ba b1 b2 b3 b4 b5 ba c]*coeffsQ;

allbits = [b1.*b2 b1.*b3 b1.*b4 b1.*b5 b1.*ba b2.*b3 b2.*b4 b2.*b5 b2.*ba b3.*b4 b3.*b5 b3.*ba b4.*b5 b4.*ba b5.*ba b1 b2 b3 b4 b5 ba c];

warning('off', 'MATLAB:rankDeficientMatrix');

t = cputime;
inittime = t;
count = 0;
oddnumbers = int64(1:2:2^n-1);
for i = int64(0): int64(2^(2^(n-1)-1))
    indexlist = bitget(i, 2^(n-1) : -1: 1);
    indexlist = indexlist + oddnumbers;
    
    if sum(allbits(indexlist, 5)) >= sum(allbits(indexlist, 9)) && sum(allbits(indexlist,9))>=sum(allbits(indexlist,12)) && sum(allbits(indexlist,12))>=sum(allbits(indexlist,14)) && sum(allbits(indexlist,14))>=sum(allbits(indexlist,15))
        A = allbits(indexlist, :);
        count = count +1;
        coeffsQ = linsolve(allbits(indexlist, :), LHS);
        RHS = allbits*coeffsQ;
        RHS = min(RHS(1:2:2^n -1), RHS(2:2:2^n));
        
        if max(abs(LHS-RHS))<1e-5 disp('found_a_quadratization'); ans = i; end
    end
    if mod(i,1000000)==0
        fprintf("up to %d, tried %d, time per million i %f, total time %f\n", i, count, cputime-t, cputime - inittime);
        t = cputime;
    end
end

%%
RHS=[b1.*b2 b1.*b3 b1.*b4 b1.*b5 b1.*ba b2.*b3 b2.*b4 b2.*b5 b2.*ba b3.*b4 b3.*b5 b3.*ba b4.*b5 b4.*ba b5.*ba b1 b2 b3 b4 b5 ba c];
warning('on', 'all');
