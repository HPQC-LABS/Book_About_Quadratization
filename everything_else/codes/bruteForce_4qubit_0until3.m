tic
warning('off');

N = 5;
quadratisations = cell((2*N+1)^5, 1);
hasQuad = zeros((2*N+1)^5, 1);
restartId = 0;
perCheck = 3200;

n=5;
allCombos=dec2bin(0:2^n-1) -'0';
b=mat2cell(allCombos,2^n,ones(1,n));
c=ones(2^n,1);

allbits=c;
for i=1:n
    for j=i+1:n
        allbits=[allbits b{i}.*b{j}]; %[b1b2, b1b3, b1b4, b1ba, b2b3, b2b4, ..., b4, ba, 1]
    end
end
for i = 1 : n
    allbits=[allbits b{i}];
end

allbits=[allbits c];
allbits(:,1)=[];
fileid1 = fopen('allQuadratisations.txt', 'wt');
fileid2 = fopen('failedFunctions.txt', 'wt');
for checkpoint = restartId : floor(((2*N+1)^5-1)/perCheck)
    parfor k = checkpoint*perCheck : min((2*N+1)^5 - 1, (checkpoint+1)*perCheck - 1)
        warning('off');
        alpha123 = floor(k/(2*N+1)^4) - N;
        alpha234 = mod(floor(k/(2*N+1)^3), 2*N+1) - N;
        alpha341 = mod(floor(k/(2*N+1)^2), 2*N+1) - N;
        alpha412 = mod(floor(k/(2*N+1)), 2*N+1) - N;
        alpha1234 = mod(k, 2*N+1) - N;
        temp = zeros(n*(n+1)/2+1i, 0);
        LHS=alpha123.*b{1}.*b{2}.*b{3} + alpha234.*b{2}.*b{3}.*b{4} + alpha341.*b{3}.*b{4}.*b{1} ...
            + alpha412 .* b{4} .* b{1} .* b{2} + alpha1234 .* b{1} .* b{2} .* b{3} .* b{4};
        
        if max(abs(LHS)) < 1e-5 % need some coeff non-zero
            continue;
        end
        
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
                hasQuad(k+1) = 1;
                temp = [temp, coeffsQ];
            end
        end
        quadratisations{k+1}{1} = [alpha123; alpha234; alpha341; alpha412; alpha1234];
        quadratisations{k+1}{2} = temp;
        
    end
    
    for k = checkpoint*perCheck : min((2*N+1)^5 - 1, (checkpoint+1)*perCheck - 1)
        alpha123 = floor(k/(2*N+1)^4) - N;
        alpha234 = mod(floor(k/(2*N+1)^3), 2*N+1) - N;
        alpha341 = mod(floor(k/(2*N+1)^2), 2*N+1) - N;
        alpha412 = mod(floor(k/(2*N+1)), 2*N+1) - N;
        alpha1234 = mod(k, 2*N+1) - N;
        if max(abs([alpha123, alpha234, alpha341, alpha412, alpha1234])) < 1e-5 % need some coeff non-zero
            continue;
        end
        temp = quadratisations{k+1}{2};
        if hasQuad(k+1) == 0
            fprintf(fileid2, "no quadratisation for alpha123 = %d, alpha234 = %d, alpha341 = %d, alpha412 = %d, alpha1234 = %d\n", ...
                alpha123, alpha234, alpha341, alpha412, alpha1234);
        end
        fprintf(fileid1, '\n%db1b2b3 + %db2b3b4 + %db3b4b1 + %db4b1b2 + %db1b2b3b4 has quadratisations:\n', alpha123, alpha234, alpha341, alpha412, alpha1234);
        fprintf(fileid1, "%+.1fb1b2 %+.1fb1b3 %+.1fb1b4 %+.1fb1ba %+.1fb2b3 %+.1fb2b4 %+.1fb2ba %+.1fb3b4 %+.1fb3ba %+.1fb4ba %+.1fb1 %+.1fb2 %+.1fb3 %+.1fb4 %+.1fba %+.1f\n", temp(1,:), temp(2,:), temp(3,:), temp(4,:), temp(5,:), temp(6,:), temp(7,:), temp(8,:), temp(9,:), temp(10,:), temp(11,:), temp(12,:), temp(13,:), temp(14,:), temp(15,:), temp(16,:));
    end
    fprintf('progress %.3f%%, restart id = %d\n', min((checkpoint+1)*perCheck/((2*N+1)^5 - 1), 1) * 100, checkpoint);
end

fclose(fileid1);
fclose(fileid2);
toc
warning('on', 'all');
