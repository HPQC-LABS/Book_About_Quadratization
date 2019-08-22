n = 10;
allCombos = dec2bin(0:2^n-1) -'0';

b = cell(n,1);
for i = 1:n
    b{i} = allCombos(:,i);
end

%Input LHS and RHS of the gadget
%example: NTR-KZFD(Kolmogorov & Zabin,2004;Freedman&Drineas,2005)
k = n-1;
LHS = -1;
RHS = 0;
for i = 1:k
    LHS = LHS.*b{i};
    RHS = RHS + b{i};
end
RHS = ( (k-1) - RHS ).*b{n};

LHS = LHS(2:2:2^n);  % when minimizing over the aux.
RHS = min(RHS(1:2:2^n -1), RHS(2:2:2^n)); %minimizing over the aux.

if min(LHS) ~= min(RHS)
    fprintf('Ground Energy not preserved.\nLHS:  %d,  RHS:  %d\n',min(LHS),min(RHS));
    LHS = LHS - min(LHS);
    RHS = RHS - min(RHS); %drop both ground energies to zero to investigate the gadget up to a constant
end

% positions of Ground states
LHS_GS = (LHS == 0);
RHS_GS = (RHS == 0);

category = zeros(1,4);
category_name = {'GS'  'FS'  'GSD'  'FSD'};
if any( RHS(LHS_GS) == 0 )
    category(1) = 1;
    if RHS == LHS
        category = ones(1,4);
    else
        if LHS_GS == RHS_GS
            category(3) = 1;
        end
        flag = true;
        LHS_spectrum = uniquetol(LHS, 1e-5 );
        for i=1:size(LHS_spectrum,1)
            LHS_mask = ( LHS == LHS_spectrum(i) );
            if ~any( RHS(LHS_mask) == LHS_spectrum(i) )
                flag = false;
                break;
            end
        end
        if flag
            category(2) = 1;
        end
    end
end

if category == zeros(1,4)
    fprintf('The Gadget preserves no states');
else
    fprintf('The Gadget is in: ');
    for i = 1:4
        if category(i) == 1
            fprintf('%s, ',category_name{i});
        end
    end
end
fprintf('\n');

%example output: "The Gadget is in: GS, FS, GSD, FSD,"

