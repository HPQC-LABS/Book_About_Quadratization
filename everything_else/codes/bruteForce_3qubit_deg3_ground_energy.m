tic
warning('off');

H = { ["ZZX" "ZXX" "ZZY" "ZYY"] }; %Types of Hamiltonians  
N1 = 2; N2 = 2; % coefficients of terms in H are (a + bi) with {a,b} = [-N1,N1] x [-N2,N2] 

n = 3; % number of qubits
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

fileid = cell(4,1);
fileid{1} = fopen('Quads_ground_energy.txt', 'wt');
fileid{2} = fopen('Quads_ground_energy_some_states.txt', 'wt');
fileid{3} = fopen('Quads_ground_energy_all_degens.txt', 'wt');
fileid{4} = fopen('Quads_ground_energy_plus_more.txt', 'wt');
fileid{5} = fopen('failedFunctions.txt', 'wt');

for T=1:N_of_Types
    quadratisations = cell( ((2*N1+1)*(2*N2+1))^N_of_terms(T), 1);
    hasQuad = zeros( ((2*N1+1)*(2*N2+1))^N_of_terms(T), 1);
    
    for checkpoint = restartId : floor( (((2*N1+1)*(2*N2+1))^N_of_terms(T)-1)/perCheck )
        parfor k = checkpoint*perCheck : min( ((2*N1+1)*(2*N2+1))^N_of_terms(T) - 1, (checkpoint+1)*perCheck - 1)
            warning('off');
            alpha = zeros(N_of_terms(T));
            for t = 1:N_of_terms(T)
                temp = mod(floor( k / ((2*N1+1)*(2*N2+1))^(t-1) ), (2*N1+1)*(2*N2+1) );
                a = floor(temp/(2*N2+1))-N1;
                b = mod(temp,2*N2+1) - N2;
                alpha(t) = a + 1i*b;
            end
            
            temp1 = zeros(allbits_size,0);
            temp2 = zeros(allbits_size,0);
            temp3 = zeros(allbits_size,0);
            temp4 = zeros(allbits_size,0);
            
            LHS = zeros(2^n,2^n);
            for t = 1:N_of_terms(T)
                LHS = LHS + alpha(t)*kron(sigma{H{T}{t}(1)-87},...
                    kron(sigma{H{T}{t}(2)-87},sigma{H{T}{t}(3)-87}));
            end
            
            if (max(max(abs(LHS))) < 1e-5) || ( ~ishermitian(LHS) ) % need some coeff non-zero
                    continue;
            end
            
            [V, d] = eig(LHS);
            LHS_spectrum = uniquetol( diag(d) , 1e-5 );
            LHS_gs = V(:, diag(d)==min(diag(d)) );
            
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
                if ~ishermitian(RHS)
                    continue;
                end
                
                [V, d] = eig(RHS);
                RHS_spectrum = uniquetol( diag(d) , 1e-5 );
                RHS_gs = V(:, diag(d)==min(diag(d)) );
                Delta_E = abs( LHS_spectrum(1) - RHS_spectrum(1) );
                
                if ( ( Delta_E < 1e-5 ) && ( max(coeffsQ) > 1e-5 ) )
                    r = rank( [LHS_gs, RHS_gs] );
                    hasQuad(k+1) = 1;
                    if ( size(LHS_gs,2) + size(RHS_gs,2) - r ) == 0
                        temp1 = [temp1, coeffsQ];
                    elseif  ( ( r == size(RHS_gs,2) ) && ( size(LHS_gs,2) <= r ) )
                        temp3 = [temp3, coeffsQ];
                    else
                        temp2 = [temp2, coeffsQ];
                    end
                    for m=2:min( size(LHS_spectrum,2),size(RHS_spectrum,2) )
                        if abs( LHS_spectrum(m) - RHS_spectrum(m) ) > 1e-5
                            break;
                        end
                    end
                    if m > 2
                        temp4 = [temp4, coeffsQ];
                    end
                end
            end
            quadratisations{k+1}{1} = alpha;
            %remove multiple entries
            if ~sum(sum(real(temp1)~=temp1))
                temp1 = uniquetol(temp1',1e-5,'ByRows',true)';
            end
            if ~sum(sum(real(temp2)~=temp2))
                temp2 = uniquetol(temp2',1e-5,'ByRows',true)';
            end
            if ~sum(sum(real(temp3)~=temp3))
                temp3 = uniquetol(temp3',1e-5,'ByRows',true)';
            end
            if ~sum(sum(real(temp4)~=temp4))
                temp4 = uniquetol(temp4',1e-5,'ByRows',true)';
            end
            quadratisations{k+1}{2} = temp1;
            quadratisations{k+1}{3} = temp2;
            quadratisations{k+1}{4} = temp3;
            quadratisations{k+1}{5} = temp4;
        end
        
        for k = checkpoint*perCheck : min( ((2*N1+1)*(2*N2+1))^N_of_terms(T) - 1, (checkpoint+1)*perCheck - 1)
            alpha = zeros(N_of_terms(T));
            for t = 1:N_of_terms(T)
                temp = mod(floor( k / ((2*N1+1)*(2*N2+1))^(t-1) ), (2*N1+1)*(2*N2+1) );
                a = floor(temp/(2*N2+1))-N1;
                b = mod(temp,2*N2+1) - N2;
                alpha(t) = a + 1i*b;
            end
            if norm(alpha) < 1e-5 % need some coeff non-zero
                continue;
            end
            
            if hasQuad(k+1) == 0
                fprintf(fileid{5}, "No quadratisation for (%d%+di)%c1%c2%c3", real(alpha(1)), imag(alpha(1)), H{T}{1}(1), H{T}{1}(2), H{T}{1}(3));
                for m = 2:N_of_terms(T)
                    fprintf(fileid{5}, " + (%d%+di)%c1%c2%c3", real(alpha(m)), imag(alpha(m)), H{T}{m}(1), H{T}{m}(2), H{T}{m}(3));
                end
                fprintf(fileid{5}, "\n");
            else
                for tempm=1:4
                    temp = quadratisations{k+1}{tempm+1};
                    if size(temp,2)
                        plus_flag = false;
                        for m = 1:N_of_terms(T)
                                if plus_flag
                                    fprintf(fileid{tempm}, " + ");
                                end
                                fprintf(fileid{tempm}, "(%d%+di)%c1%c2%c3", real(alpha(m)), imag(alpha(m)), H{T}{m}(1), H{T}{m}(2), H{T}{m}(3));
                                plus_flag = true;
                        end
                        fprintf(fileid{tempm}, " has quadratisations:\n");
                        %print possible quadratisations
                        for count=1:size(temp,2)
                            m = 0;
                            plus_flag = false;
                            for i=1:4
                                for j=1:4
                                    for l=1:4
                                        if (i==4||j==4||l==4)
                                            m = m+1;
                                            if norm(temp(m,count)) > 1e-5
                                                if imag(temp(m,count)) > 1e-5
                                                    if plus_flag
                                                        fprintf(fileid{tempm}, " + ");
                                                    end
                                                    fprintf(fileid{tempm}, "(%.1f%+.1fi)", real(temp(m,count)), imag(temp(m,count)));
                                                else
                                                    fprintf(fileid{tempm}, " %+ .1f", real(temp(m,count)));
                                                end
                                                plus_flag = true;
                                                if(i~=4)
                                                    fprintf(fileid{tempm}, "%c1", char(87+i));
                                                end
                                                if(j~=4)
                                                	fprintf(fileid{tempm}, "%c2", char(87+j));
                                                end
                                                if(l~=4)
                                                	fprintf(fileid{tempm}, "%c3", char(87+l));
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            fprintf(fileid{tempm}, "\n");
                        end
                        fprintf(fileid{tempm}, "\n");
                    end
                end
            end
        end
        fprintf('progress %.3f%%, restart id = %d\n', min((checkpoint+1)*perCheck/(((2*N1+1)*(2*N2+1))^N_of_terms(T) - 1), 1) * 100, checkpoint);
    end
end

fclose(fileid{1});
fclose(fileid{2});
fclose(fileid{3});
fclose(fileid{4});
fclose(fileid{5});
toc
warning('on', 'all');
