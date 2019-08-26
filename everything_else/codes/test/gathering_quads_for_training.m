n = 5;
coeffs_range = -2:2;
allCombos = dec2bin(0:2^n-1) -'0';

b = cell(n,1);
for i = 1:n
    b{i} = allCombos(:,i);
end

LHS = b{1}.*b{2}.*b{3}.*b{4} + b{2}.*b{3}.*b{4}.*b{5};
LHS = LHS';
%LHS = LHS(2:2:2^n); %when using aux

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
base = size(coeffs_range,2);
init = int2str( (base-1)/2 );

data_percentage = 0.001;
conflicts_threshold = 50;
restartId = 0;
perCheck = 100000;
input  = zeros(floor(data_percentage*base^coeffs_size),coeffs_size);
target = zeros(floor(data_percentage*base^coeffs_size),1);
t = cputime;
t_init = t;
for checkpoint = restartId : floor( (data_percentage*base^coeffs_size-1)/perCheck )
    k = int64(checkpoint*perCheck) : min( int64(data_percentage*base^coeffs_size-1), int64((checkpoint+1)*perCheck) - 1 );
    %k = randperm( base^coeffs_size-1 , perCheck); % use this for random sample
    coeffs = ndec2base(k,base,coeffs_size) - init;
    
    RHS = coeffs*allbits';
    %RHS = min(RHS(:,1:2:2^n-1),RHS(:,2:2:2^n)); % when using aux
    const_term = -min(RHS,[],2);
    RHS = RHS + const_term;
    
    conflicts_percent = mean( RHS ~= LHS , 2 ) * 100; % percentage of overall conflicts
    index_good = (conflicts_percent <= conflicts_threshold);
	
    % training data
    input(k+1,:) = coeffs;
    target(k+1)  = index_good;
    
    if mod(checkpoint,100) == 0
    fprintf('progress %.4f%%, restart id = %d, step time = %.3f, total time = %.3f, good = %d\n',...
      min(double((checkpoint+1)*perCheck)/(data_percentage*base^coeffs_size-1), 1) * 100, checkpoint,cputime - t,cputime - t_init,sum(target));
    t = cputime;
    end
end
