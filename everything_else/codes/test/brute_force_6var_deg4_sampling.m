n = 7;
query = [1 5 10; %bad quads to be used for training
         1 2 4;
         1 0 1;
         2 6 12; %bad quads to be used for training
         2 3 5;
         2 0 2];
% each row is a query (conflict type,range) where:
%       conflict type is 1 for ground state conflicts
%                        2 for overall conflicts
%       range gives the range of the number of conflicts
%       required to print out the result
% for example: 2 0 5 means that at least 59/64 (~92%) states are preserved
%              1 0 1 means that at least 56/57 (~98%) of the ground states are
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
index_init = index_init - 1;

count = 0;
restartId = 0;
perCheck = 100;
percentage = 0.0001;
t_init = cputime;
for checkpoint = restartId : floor( percentage * int64(3^21-1)/perCheck )
	for i = 1:perCheck
        sampling_number = randi([1,10]); %think again about this
        index_last = zeros(1,21);
        index_temp = 2*randi([0,1],1,sampling_number)-1;
        index = randi([1,21],1,sampling_number);
        for j = 1:sampling_number
            index_last(index(j)) = index_temp(j);
        end
        
        for j = 1:378
            coeffs = [index_last index_init(j,:)];
            conflicts = [];
            
            RHS = allbits * coeffs';
            RHS = min(RHS(1:2:2^n -1), RHS(2:2:2^n)); %minimizing over 'ba'
            RHS = RHS - min(RHS); % make sure Ground Energy is zero
            
            mask_RHS = (RHS==0);
            conflicts(1) = sum( mask~=mask_RHS ); %number of ground state conflicts

            dif = LHS - RHS;
            conflicts(2) = sum( dif~=0); %overall number of conflicts

            for k = 1:number_queries
                if ( query(k,2) <= conflicts(query(k,1)) ) && ( conflicts(query(k,1)) <= query(k,3) )
                    fprintf(fileRHS{k},    '%2d  ', RHS);    fprintf(fileRHS{k},    '\n');
                    fprintf(fileCoeffs{k}, '%2d  ', coeffs); fprintf(fileCoeffs{k}, '\n');
                    count = count + 1;
                end
            end
        end
	end
    fprintf('progress %.3f%%, restart id = %d, found %d quads\n', min( double((checkpoint+1)*perCheck)/(percentage*(3^21-1)), 1) * 100, checkpoint,count);
end

for i = 1:number_queries
    fclose(fileRHS{i});
    fclose(fileCoeffs{i});
end

