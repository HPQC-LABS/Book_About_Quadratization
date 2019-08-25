n = 4;
coeffs_range = [0 1 -1 2 -2];
allCombos = dec2bin(0:2^n-1) -'0';

b = cell(n,1);
for i = 1:n
    b{i} = allCombos(:,i);
end

LHS = b{1}.*b{2}.*b{3} + b{2}.*b{3}.*b{4};
%LHS = LHS(2:2:2^n); %when using aux
mask = (LHS==0); % positions of Ground states

file_good    = fopen('good.txt', 'w');
fileRHS_good = fopen('RHS__good.txt', 'w');
file_bad     = fopen('bad.txt', 'w');
fileRHS_bad  = fopen('RHS_bad.txt', 'w');

allbits = [];
for i = 1:n
    for j = i+1:n
        allbits = [allbits b{i}.*b{j}];
    end
end

for i = 1:n
    allbits = [allbits b{i}];
end

coeffs_size = n*(n+1)/2;
allcoeffs = repmat(coeffs_range,1,coeffs_size);
base = size(coeffs_range,2);
unit_mod_base = 1:base:base*coeffs_size;

count1 = 0; count2 = 0;
restartId = 0;
perCheck = 100000;
wanted_quad = zeros(1,base^coeffs_size);
t = cputime;
t_init = t;
for checkpoint = restartId : floor( int64(base^coeffs_size-1)/perCheck )
    parfor k = int64(checkpoint*perCheck) : min( int64(base^coeffs_size-1), int64((checkpoint+1)*perCheck) - 1 )
        indexlist = dec2base(k,base,coeffs_size)-'0';
        indexlist = indexlist + unit_mod_base;
        coeffs = allcoeffs(indexlist);
        conflicts = [];

        RHS = allbits * coeffs';
        %RHS = min(RHS(1:2:2^n-1),RHS(2:2:2^n));%when using aux
        RHS = RHS - min(RHS); % make sure Ground Energy is zero

        mask_RHS = (RHS==0);
        conflicts(1) = sum( mask~=mask_RHS ); %number of ground state conflicts
        %conflicts(2) = sum( LHS~=RHS ); %overall number of conflicts
        
        if conflicts(1) <= 5
            wanted_quad(k+1) = 1;
            count1 = count1 + 1;
        end
        if conflicts(1) >= 6 && conflicts(1) <=8
            wanted_quad(k+1) = 1;
            count2 = count2 + 1;
        end
    end
    for k = int64(checkpoint*perCheck) : min( int64(base^coeffs_size-1), int64((checkpoint+1)*perCheck) - 1 )
        if wanted_quad(k+1)
            indexlist = dec2base(k,base,coeffs_size)-'0';
            indexlist = indexlist + unit_mod_base;
            coeffs = allcoeffs(indexlist);

            RHS = allbits * coeffs';
            %RHS = min(RHS(1:2:2^n-1),RHS(2:2:2^n));%when using aux
            const_term = -min(RHS);
            RHS = RHS + const_term;
            
            mask_RHS = (RHS==0);
            conflicts(1) = sum( mask~=mask_RHS ); %number of ground state conflicts
            %conflicts(2) = sum( LHS~=RHS ); %overall number of conflicts
            
            if conflicts(1) <= 5
                %fprintf(fileRHS_good,'%2d,%2d: ',conflicts(1),conflicts(2));
                %fprintf(fileRHS_good,    '%2d  ', RHS);    fprintf(fileRHS_good,'\n');
                fprintf(file_good, '%2d  ', coeffs); fprintf(file_good,'\n');
            end
            if conflicts(1) >= 6 && conflicts(1) <= 8
                %fprintf(fileRHS_bad,'%2d,%2d: ',conflicts(1),conflicts(2));
                %fprintf(fileRHS_bad,    '%2d  ', RHS);    fprintf(fileRHS_bad,'\n');
                fprintf(file_bad, '%2d  ', coeffs); fprintf(file_bad,'\n');
            end
        end
    end
    fprintf('progress %.4f%%, restart id = %d, step time = %.3f, total time = %.3f, count1 = %d, count2 = %d\n',...
      min(double((checkpoint+1)*perCheck)/(base^coeffs_size-1), 1) * 100, checkpoint,cputime - t,cputime - t_init,count1,count2);
    t = cputime;
end

fclose(fileRHS_good);
fclose(file_good);
fclose(fileRHS_bad);
fclose(file_bad);
