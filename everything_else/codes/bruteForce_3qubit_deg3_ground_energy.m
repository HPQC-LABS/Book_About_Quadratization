tic
warning('off');

H = { ["ZZX" "ZXX" "ZZY" "ZYY"] }; %Types of Hamiltonians
n = 3; % number of qubits
N = 3; % coefficients of terms in H are in [-N,N]

N_of_Types = size(H,2);
for T=1:N_of_Types
    N_of_terms(T) = size(H{T},2);
end
restartId = 0;
perCheck = 100;
sigma = cell(4,1);
sigma{1} = [0 1 ; 1 0];
sigma{2} = [0 -1i ; 1i 0];
sigma{3} = [1 0 ; 0 -1];
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

for T=1:N_of_Types
    quadratisations = cell((2*N+1)^N_of_terms(T), 1);
    hasQuad = zeros((2*N+1)^N_of_terms(T), 1);
    
    for checkpoint = restartId : floor(((2*N+1)^N_of_terms(T)-1)/perCheck)
        parfor k = checkpoint*perCheck : min((2*N+1)^N_of_terms(T) - 1, (checkpoint+1)*perCheck - 1)
            warning('off');
            alpha = zeros(N_of_terms(T));
            for t = 1:N_of_terms(T)
                alpha(t) = mod(floor(k/(2*N+1)^(t-1)), 2*N+1) - N;
            end
            
            temp = zeros(allbits_size,0);
            
            LHS = zeros(2^n,2^n);
            for t = 1:N_of_terms(T)
                LHS = LHS + alpha(t)*kron(sigma{H{T}{t}(1)-87},...
                    kron(sigma{H{T}{t}(2)-87},sigma{H{T}{t}(3)-87}));
            end
            
            if max(max(abs(LHS))) < 1e-5 % need some coeff non-zero
                    continue;
            end
            [V,d] = eig(LHS);
            d_gs = min(diag(d));
            
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
                [V,d] = eig(RHS);
                dmin = min(diag(d));
                
                if ( (abs(d_gs-dmin) <1e-5) && (max(coeffsQ)>1e-5) )
                    hasQuad(k+1) = 1;
                    temp = [temp, coeffsQ];
                end
                
            end
            quadratisations{k+1}{1} = alpha;
            quadratisations{k+1}{2} = temp;
        end

        for k=0:(2*N+1)^N_of_terms(T)-1
            alpha = zeros(N_of_terms(T));
            for m = 1:N_of_terms(T)
                alpha(m) = mod(floor(k/(2*N+1)^(m-1)), 2*N+1) - N;
            end
            if max(abs(alpha)) < 1e-5 % need some coeff non-zero
                continue;
            end
            
            if hasQuad(k+1) == 0
                fprintf(fileid2, "no quadratisation for");
                for m = 1:N_of_terms(T)
                    fprintf(fileid2, " %+d%c1%c2%c3", alpha(m), H{T}{m}(1), H{T}{m}(2), H{T}{m}(3));
                end
                fprintf(fileid2, "\n");
            else
                temp = quadratisations{k+1}{2};
                fprintf(fileid1, "\n");
                for m = 1:N_of_terms(T)
                    fprintf(fileid1, "%+d%c1%c2%c3 ", alpha(m), H{T}{m}(1), H{T}{m}(2), H{T}{m}(3));
                end
                fprintf(fileid1, 'has quadratisations:\n');
                %print possible quadratisations
                for count=1:size(temp,2)
                    m = 0;
                    for i=1:4
                        for j=1:4
                            for l=1:4
                                if (i==4||j==4||l==4)
                                    m = m+1;
                                    if(abs(temp(m,count))>1e-5)
                                        fprintf(fileid1, "% +.1f", temp(m,count));
                                        if(i~=4)
                                            fprintf(fileid1, "%c1", char(87+i));
                                        end
                                        if(j~=4)
                                            fprintf(fileid1, "%c2", char(87+j));
                                        end
                                        if(l~=4)
                                            fprintf(fileid1, "%c3", char(87+l));
                                        end
                                    end
                                end
                            end
                        end
                    end
                    fprintf(fileid1, "\n");
                end
            end
        end
        fprintf('progress %.3f%%, restart id = %d\n', min((checkpoint+1)*perCheck/((2*N+1)^N_of_terms(T) - 1), 1) * 100, checkpoint);
    end
end

fclose(fileid1);
fclose(fileid2);
toc
warning('on', 'all');
