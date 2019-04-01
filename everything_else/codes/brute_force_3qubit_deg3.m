tic
warning('off');

n = 3; % number of qubits
N = 2; % coefficients of terms in H are in [-N,N]
quadratisations = cell((2*N+1)^2, 1);
hasQuad = zeros((2*N+1)^2, 1);
sigma = cell(4,1);
sigma{1} = [0,1;1,0];
sigma{2} = [0,-sqrt(-1);sqrt(-1),0];
sigma{3} = [1,0;0,-1];
sigma{4} = eye(2);

allbits = zeros(2^n,0);
%All 3-qubit combinations of X,Y,Z,I that are up to quadratic
for i=1:4
    for j=1:4
        for k=1:4
            if (i==4||j==4||k==4)
                allbits = [allbits kron(sigma{i},kron(sigma{j},sigma{k}))]; %[x1x2,x1y2,x1z2,x1x3,x1y3,...,x3,y3,z3,1];
            end
        end
    end
end
allbits_size = size(allbits,2)/2^n;

fileid1 = fopen('allQuadratisations.txt', 'wt');
fileid2 = fopen('failedFunctions.txt', 'wt');

for k = 0:(2*N+1)^2-1
    alphazzx = mod(floor(k/(2*N+1)), 2*N+1) - N;
    alphazxx = mod(k, 2*N+1) - N;
    temp = zeros(allbits_size,0);
    
    LHS = alphazzx*kron(sigma{3},kron(sigma{3},sigma{1}))...
            + alphazxx*kron(sigma{3},kron(sigma{1},sigma{1}));
	
    if max(max(abs(LHS))) < 1e-5 % need some coeff non-zero
            continue;
    end
	
    coeffsQ = ones(allbits_size,1);
    numbers = int64(1:2^n);
    
    for i = int64(1): int64(2^(2^n)-1)
        indexlist = bitget(i, 2^n : -1: 1);
        indexlist = indexlist.*numbers;
        indexlist = indexlist(indexlist~=0);
        
        A = allbits(indexlist,:)';
        A = reshape(A,2^n,[]);
        B = A(:,1:allbits_size);
        for m = 1:size(indexlist,2)-1
            B = [B;A(:,1+m*allbits_size:(m+1)*allbits_size)];
        end
        
        coeffsQ = B\reshape(LHS(indexlist,:)',[],1);
        RHS = allbits*kron(coeffsQ,eye(2^n));
        
        if ( (max(max(abs(LHS-RHS)))<1e-5) && (max(coeffsQ)>1e-5) )
            hasQuad(k+1) = 1;
            temp = [temp, coeffsQ];
        end
    end
    
    quadratisations{k+1}{1} = [alphazzx; alphazxx];
    quadratisations{k+1}{2} = temp;
end

for k=0:(2*N+1)^2-1
    alphazzx = mod(floor(k/(2*N+1)), 2*N+1) - N;
    alphazxx = mod(k, 2*N+1) - N;
    if max(abs([alphazzx, alphazxx])) < 1e-5 % need some coeff non-zero
        continue;
    end
    temp = quadratisations{k+1}{2};
    if hasQuad(k+1) == 0
        fprintf(fileid2, "no quadratisation for %+dz1z2x3 %+dz1x2x3\n",alphazzx, alphazxx);
    else
    fprintf(fileid1, '\n%+dz1z2x3 %+dz1x2x3 has quadratisations:\n', alphazzx, alphazxx);
    fprintf(fileid1, "%+.1fx1x2 %+.1fx1y2 %+.1fx1z2 %+.1fx1x3 %+.1fx1y3 %+.1fx1z3 %+.1fx1 %+.1fy1x2 %+.1fy1y2 "+...
        "%+.1fy1z2 %+.1fy1x3 %+.1fy1y3 %+.1fy1z3 %+.1fy1 %+.1fz1x2 %+.1fz1y2 %+.1fz1z2 %+.1fz1x3 %+.1fz1y3 "+...
        "%+.1fz1z3 %+.1fz1 %+.1fx2x3 %+.1fx2y3 %+.1fx2z3 %+.1fx2 %+.1fy2x3 %+.1fy2y3 %+.1fy2z3 %+.1fy2 %+.1fz2x3 "+...
        "%+.1fz2y3 %+.1fz2z3 %+.1fz2 %+.1fx3 %+.1fy3 %+.1fz3 %+.1f\n",temp(1,:),temp(2,:),temp(3,:),temp(4,:),...
        temp(5,:),temp(6,:),temp(7,:),temp(8,:),temp(9,:),temp(10,:),temp(11,:),temp(12,:),temp(13,:),temp(14,:),...
        temp(15,:),temp(16,:),temp(17,:),temp(18,:),temp(19,:),temp(20,:),temp(21,:),temp(22,:),temp(23,:),...
        temp(24,:),temp(25,:),temp(26,:),temp(27,:),temp(28,:),temp(29,:),temp(30,:),temp(31,:),temp(32,:),...
        temp(33,:),temp(34,:),temp(35,:),temp(36,:),temp(37,:));
    end
end

fclose(fileid1);
fclose(fileid2);
toc
warning('on', 'all');