function [verified,const_terms] = verify(n,aux,LHS,allbits,coef)
    coeffs_size = n*(n+1)/2;
    
    temp = [coef{:}];
    coef = cell2mat(temp);
    coef = reshape(coef,[],coeffs_size);
    
	temp = [allbits{:}];
    allbits = cell2mat(temp);
    allbits = reshape(allbits,2^n,coeffs_size);
    
    temp = [LHS{:}];
    LHS = cell2mat(temp);
    LHS = reshape(LHS,1,[]);
    
    RHS = coef*allbits';
    if aux
        RHS = min(RHS(:,1:2:2^n-1),RHS(:,2:2:2^n)); % when using aux
    end

    const_terms = -min(RHS,[],2);
    RHS = RHS + const_terms + min(LHS);
    
    % checking if every single state is equal to the min
    % of the corresponding states given by the COMIGs
    if min(RHS) == LHS
        verified = 1;
        fprintf('Hooray!!!\n');
    else
        verified = 0;
        fprintf('Try Again\n');
    end
end