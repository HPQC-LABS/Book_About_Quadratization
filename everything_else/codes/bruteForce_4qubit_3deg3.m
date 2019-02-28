
n=5;
allCombos=dec2bin(0:2^n-1) -'0';
b=mat2cell(allCombos,2^n,ones(1,n));
c=ones(2^n,1);

allbits=c;
for i=1:n
    for j=i+1:n
        allbits=[allbits b{i}.*b{j}];
    end
end
for i = 1 : n
    allbits=[allbits b{i}];
end

allbits=[allbits c];
allbits(:,1)=[];

warning('off');
N = 4;
goodtrinomial = zeros(3, 0);
quad = zeros(2*N+1, 2*N+1, 2*N+1, n*(n+1)/2+1);
for alpha123 = [-N:-1 1:N]
    for alpha234 = [-N:-1 1:N]
        for alpha341 = [-N:-1 1:N]
            LHS=alpha123.*b{1}.*b{2}.*b{3} + alpha234.*b{2}.*b{3}.*b{4} + alpha341.*b{3}.*b{4}.*b{1};
            coeffsQ=ones(n*(n+1)/2+1,1);
            
            oddnumbers = int64(1:2:2^n-1);
            for i = int64(0): int64(2^(2^(n-1)-1))
                indexlist = bitget(i, 2^(n-1) : -1: 1);
                indexlist = indexlist + oddnumbers;
                
                % may use symmetry to skip some cases
                A = allbits(indexlist, :);
                coeffsQ = A\LHS(1:2:2^n-1);
                RHS = allbits*coeffsQ;
                RHS = min(RHS(1:2:2^n -1), RHS(2:2:2^n));
                
                if max(abs(LHS(1:2:2^n-1)-RHS))<1e-5
                    quad(1+N+alpha123,1+N+alpha234,1+N+alpha341, :) = coeffsQ; %quadratisation coefficients are stored here
                    goodtrinomial = [goodtrinomial, [alpha123; alpha234; alpha341]]; 
 
                end
                
                if max(abs(quad(1+N+alpha123,1+N+alpha234,1+N+alpha341, :))) > 1e-5
                    break;
                end
            end
            
            %%
        end
    end
end
warning('on', 'all');
