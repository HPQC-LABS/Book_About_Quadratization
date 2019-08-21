n = 7;
query = [1 4;
         1 1;
         2 5;
         2 2];
% each row is a query (conflict type,threshold) where:
%       conflict type is 1 for ground state conflicts
%                        2 for overall conflicts
%       threshold is the maximum number of conflicts
%       that are tolerated to print out the result
% for example: 2 5 means that 59/64 (~92%) states are preserved
%              1 1 means that 56/57 (~98%) of the ground states are
%              preserved

number_queries = size(query,1);
allCombos = dec2bin(0:2^n-1) -'0';

b = cell(n,1);
for i = 1:n
    b{i} = allCombos(:,i); % b{n} is the auxiliary 'ba'
end

%choose LHS to be quadratised with one aux
LHS = b{1}.*b{2}.*b{3}.*b{4} + b{3}.*b{4}.*b{5}.*b{6};
LHS = LHS(2:2:2^n);
mask = (LHS==0); % positions of Ground states

for i = 1:number_queries
    fileCoeffs{i} = fopen(strcat('coeffs',int2str(i),'.txt'), 'w');
    fileRHS{i} = fopen(strcat('RHS',int2str(i),'.txt'), 'w');
end

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
allcoeffs = repmat([0 1 -1],1,coeffs_size)';

%use this to exploit symmetry of this particular LHS
index_temp = [0 0; 0 1; 0 2; 1 1; 1 2; 2 2];
index_init = [];
for i = 1:6
    for j = 1:6
        for k = i:6
            for l = 1:3
                index_init = [index_init; [index_temp(i,:) index_temp(j,:) index_temp(k,:) index_temp(l,2)] ];
            end
        end
    end
end

count = 0;
unit_mod3 = 1:3:3*28;

restartId = 0;
perCheck = 1000;
t_init = t;
t = cputime;
for checkpoint = restartId : floor( int64(3^21-1)/perCheck )
	parfor i = int64(checkpoint*perCheck) : min( int64(3^21-1), int64((checkpoint+1)*perCheck) - 1 )
        index_last = dec2base(i,3,21)-'0';
        for j = 1:378
            indexlist = [index_last index_init(j,:)];
            indexlist = indexlist + unit_mod3;
            coeffs = allcoeffs(indexlist);
            conflicts = [];

            RHS = allbits * coeffs;
            RHS = min(RHS(1:2:2^n -1), RHS(2:2:2^n)); %minimizing over 'ba'
            RHS = RHS - min(RHS); % make sure Ground Energy is zero

            mask_RHS = (RHS==0);
            conflicts(1) = sum( mask~=mask_RHS ); %number of ground state conflicts

            dif = LHS - RHS;
            conflicts(2) = sum( dif~=0); %overall number of conflicts

            for k = 1:number_queries
                if conflicts(query(k,1)) <= query(k,2)
                    fprintf(fileRHS{k},    '%2d  ', RHS);    fprintf(fileRHS{k},    '\n');
                    fprintf(fileCoeffs{k}, '%2d  ', coeffs); fprintf(fileCoeffs{k}, '\n');
                    %count(k) = count(k) + 1; % not used for parfor implementation
                    count = count + 1;
                end
            end
            if conflicts(1) == 0
                break;
            end
        end
        %{
        % not used for parfor implementation
        if mod(i,1000) == 0
            fprintf('step(x378000) = %4d,  step time: %8.4f,  total time: %8.4f', i/1000,cputime - t,cputime - t_init);
            for k = 1:number_queries
                fprintf('  count%d = %d,',k,count(k));
            end
            fprintf('\n');
            t = cputime;
        end
        %}
	end
    fprintf('progress %.4f%%, restart id = %d, step time = %.3f, total time = %.3f, found %d quads\n', min((checkpoint+1)*perCheck/int64(3^21-1), 1) * 100, checkpoint,cputime - t,cputime - t_init,count);
    t = cputime;
end

for i = 1:number_queries
    fclose(fileRHS{i});
    fclose(fileCoeffs{i});
end

