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
good_RHS = [];
good_coeffs = [];
good_prcntg = [];
t = cputime;
t_init = t;
for checkpoint = restartId : floor( (data_percentage*base^coeffs_size-1)/perCheck )
    k = int64(checkpoint*perCheck) : min( int64(data_percentage*base^coeffs_size-1), int64((checkpoint+1)*perCheck) - 1 );
	coeffs = ndec2base(k,base,coeffs_size) - init;
    
    index_pred = predict(coeffs);
    coeffs_pred = coeffs(index_pred >= 0.5,:);
    
    if ~sum(size(coeffs_pred) == 0)
        RHS_pred = coeffs_pred*allbits';
        %RHS_pred = min(RHS_pred(:,1:2:2^n-1),RHS_pred(:,2:2:2^n)); % when using aux
        const_term = -min(RHS_pred,[],2);
        RHS_pred = RHS_pred + const_term;
        
        conflicts_percent = mean( RHS_pred ~= LHS , 2 ) * 100; % percentage of overall conflicts
        index_pred_good = (conflicts_percent <= conflicts_threshold);
        
        conflicts_percent_good = conflicts_percent(index_pred_good);
        RHS_pred_good = RHS_pred(index_pred_good,:);
        coeffs_pred_good = horzcat( coeffs_pred(index_pred_good,:), const_term(index_pred_good) );
        
        good_RHS = [good_RHS; RHS_pred_good];
        good_coeffs = [good_coeffs; coeffs_pred_good];
        good_prcntg = [good_prcntg; conflicts_percent_good];
    end
    
    if mod(checkpoint,10) == 0
    fprintf('progress %.4f%%, restart id = %d, step time = %.3f, total time = %.3f, good = %d\n',...
      min(double((checkpoint+1)*perCheck)/(data_percentage*base^coeffs_size-1), 1) * 100, checkpoint,cputime - t,cputime - t_init,size(good_prcntg,1));
    t = cputime;
    end
end
