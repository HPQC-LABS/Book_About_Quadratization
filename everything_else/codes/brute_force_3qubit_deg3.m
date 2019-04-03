tic
warning('off');

H = { ["ZZX ZXX"] ["ZZX" "ZXX" "ZZY" "ZYY"] }; %Types of Hamiltonians
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
            V_gs = V(:,diag(d)==min(diag(d))); %ground state of H
            R_gs = rref(V_gs');
            gs_num_of_vec = size(V_gs,2);
            
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
                V = V(:,diag(d)==min(diag(d))); %ground state of RHS
                
                if ( size(V,2)==gs_num_of_vec ) %i.e. same number of eigenvectors
                    R = rref(V');
                    if ( (max(max(abs(R-R_gs)))<1e-5) && (max(coeffsQ)>1e-5) )
                        hasQuad(k+1) = 1;
                        temp = [temp, coeffsQ];
                    end
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
        fprintf('progress %.3f%%, restart id = %d\n', min((checkpoint+1)*perCheck/((2*N+1)^N_of_terms(T) - 1), 1) * 100, checkpoint);
    end
end

fclose(fileid1);
fclose(fileid2);
toc
warning('on', 'all');
