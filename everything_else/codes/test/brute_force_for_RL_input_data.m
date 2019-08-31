n = 5;
coeffs_range = -1:1;
allCombos = dec2bin(0:2^n-1) -'0';

b = cell(n,1);
for i = 1:n
    b{i} = allCombos(:,i);
end

LHS = b{1}.*b{2}.*b{3} + b{2}.*b{3}.*b{4};
LHS = LHS';
LHS = LHS(2:2:2^n); %when using aux

allbits = [];
for i = 1:n
    for j = i+1:n
        allbits = [allbits b{i}.*b{j}];
    end
end

for i = 1:n
    allbits = [allbits b{i}];
end
allbits = allbits';

coeffs_size = n*(n+1)/2;
base = size(coeffs_range,2);
init = int2str( (base-1)/2 );

data_percentage = 1;
conflicts_threshold = 40;

restartId = 0;
perCheck = 100000;

data_size = floor(data_percentage*base^coeffs_size);
progress_const = 100 * perCheck / data_size;
good_coeffs = [];
good_prcntg = [];
good_count = 0;
t = cputime;
t_init = t;
for checkpoint = restartId : floor( (data_size-1)/perCheck )
    k = int64(checkpoint*perCheck) : min( floor(data_size-1), int64((checkpoint+1)*perCheck) - 1 );
    %k = randperm( data_size-1 , perCheck); % use this for random sample
    coeffs = ndec2base(k,base,coeffs_size) - init;
    
    RHS = coeffs*allbits;
    RHS = min(RHS(:,1:2:2^n-1),RHS(:,2:2:2^n)); % when using aux
    RHS = RHS - min(RHS,[],2);
    
    conflicts_percent = mean( RHS ~= LHS , 2 ) * 100; % percentage of overall conflicts
    index_good = (conflicts_percent <= conflicts_threshold);
	
    %input(checkpoint*perCheck+1:min( floor(data_size), int64((checkpoint+1)*perCheck)),:) = coeffs; %for sampling
    %target(checkpoint*perCheck+1:min( floor(data_size), int64((checkpoint+1)*perCheck))) = index_good; %for sampling
    if any(index_good)
        good_coeffs = [good_coeffs; coeffs(index_good,:)];
        good_prcntg = [good_prcntg; conflicts_percent(index_good)];
        good_count = good_count + sum(index_good);
    end
    if mod(checkpoint,10) == 0
        fprintf('progress %.4f%%, restart id = %d, step time = %.3f, total time = %.3f, good = %d\n',...
            min((checkpoint+1)*progress_const, 100), checkpoint,cputime - t,cputime - t_init,good_count);
        t = cputime;
    end
end

fprintf('progress %.3f%%, restart id = %d, step time = %.3f, total time = %.3f, good = %d\n',...
    min((checkpoint+1)*progress_const, 100), checkpoint,cputime - t,cputime - t_init,good_count);

input = good_coeffs(randperm(size(good_coeffs,1),3000),:);

RHS = good_coeffs*allbits;
RHS = min(RHS(:,1:2:2^n-1),RHS(:,2:2:2^n)); % when using aux
const_term = -min(RHS,[],2);
RHS = RHS + const_term;

accuracy_percent = mean( RHS == LHS , 2 ) * 100; % percentage of overall conflicts
index_good = (accuracy_percent <=65);

temp = good_coeffs(index_good,:);
reset_state = temp(randperm(size(temp,1),1),:);
save('data.mat','input','LHS','allbits','reset_state');
