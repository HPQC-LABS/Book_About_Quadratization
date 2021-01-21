%% PTR_BG Pg. 17, Eqns 50

b = dec2bin(2^6-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);ba2=b(:,6);
LHS = b1.*b2.*b3.*b4;
RHS = ba1.*(2+ b1 -b2 - b3 - b4) + ba2.*(1 + b2 - b3 - b4) + b3.*b4;
%{
%RESULTS%
%LHS shares all of its states with the RHS

%COMMENTS%
Original seems to be the correct derivation, need to determine where additional
term at the end comes from

%ADDITIONAL TESTS%
k = 4
%%%Derivation Test
%b = dec2bin(2^6-1:-1:0)-'0';
%b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);ba2=b(:,6);
%LHS = b1.*b2.*b3.*b4;
%RHS = ba1.*(2+ b1 -b2 - b3 - b4) + ba2.*(1 + b2 - b3 - b4);

%k = 5
%%%Original Test
%b = dec2bin(2^8-1:-1:0)-'0';
%b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);b5=b(:,5);ba1=b(:,6);ba2=b(:,7);ba3=b(:,8);
%LHS = b1.*b2.*b3.*b4.*b5;
%RHS = ba1.*(3 + b1 - b2 - b3 - b4 - b5) + ba2.*(2 + b2 - b3 - b4 - b5) + ba3.*(1 + b3 - b4 - b5) + b4.*b5;
%}

%% PTR-BCR-1: Pg. 20, Eqns 59

b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = 1/2*(b1 + b2 + b3 + b4 - 2*ba1).*(b1 + b2 + b3 + b4 - 2*ba1 - 1);

%{
%RESULTS%
%Ground states are the same, but other states differ

%COMMENTS%
Both derivations use the same formula for the example, need to somehow test
the expanded version to find mistakes

%ADDITIONAL TESTS%
%n = 6 with m = 2   
%%%Original
%b = dec2bin(2^8-1:-1:0)-'0';
%b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);b5=b(:,5);b6=b(:,6);ba1=b(:,7);ba2=b(:,8);
%LHS = b1.*b2.*b3.*b4.*b5.*b6;
%RHS = 1/2*(b1 + b2 + b3 + b4 + b5 + b6 - 2*(ba1+ba2)).*(b1 + b2 + b3 + b4 + b5 + b6 - 2*(ba1+ba2) - 1);
%}

%% PTR-BCR-2 Pg. 21, Eqns 62

%STILL YET TO WORK%

%Original
%b = dec2bin(2^6-1:-1:0)-'0';
%b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);ba2=b(:,6);
%LHS = b1.*b2.*b3.*b4;
%RHS = 1/2*(4 + b1 + b2 + b3 + b4 - ba1 - 2*ba2).*(3 + b1 + b2 + b3 + b4 - ba1 - 2*ba2);

%{
%RESULTS%
%!!!LHS and RHS Share no common states!!! 

%COMMENTS%
Not sure how the original example was determined with the formula. Still
need to continue testing with the base equation as for some reason it
didn't work

%ADDITIONAL TESTS%
%New
b = dec2bin(2^8-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);ba2=b(:,6);ba3=b(:,7);ba4=b(:,8);
LHS = b1.*b2.*b3.*b4;
RHS = (28 + b1 + b2 + b3 + b4 - (ba1 + 2*ba2 + 4*ba3 + 8*ba4)).^2;
%}

%% PTR_Ishikawa Pg. 18, Eqn 52

b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = ba.*(3 - 2*b1 - 2*b2 - 2*b3 - 2*b4) + b1.*b2 + b1.*b3 + b1.*b4 + b2.*b3 + b2.*b4 +  b3.*b4;

%{
%RESULTS%
%LHS shares all of its states with the RHS

%COMMENTS%
Original seems to be the correct equation, need to determine how they
achieve the example equation

%ADDITIONAL TESTS%
%Expanding the formula test
b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = ba.*(3 - 2*b1 - 2*b2 - 2*b3 - 2*b4) + b1.*b2 + b1.*b3 + b1.*b4 + b2.*b3 + b2.*b4 +  b3.*b4 + b1 + b2 + b3 + b4;
%}

%% SFR-ABCG-4 Pg. --, Eqns --

b = dec2bin(2^5-1:-1:0)-'0';
b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);ba1=b(:,5);
LHS = b1.*b2.*b3.*b4;
RHS = + b4.*b2 + ba1.*(b3-1);
%{
%RESULTS%
%LHS shares all of its states with the RHS

%COMMENTS%
This method seems to be working for me. There are certain conditions I've
put on what is alpha_i in the equation which is that it is equal to the
number of terms in that iteration. If someone wants to look this over that
would be great.

%ADDITIONAL TESTS%
%k = 4, m = 1
%RHS = + b4.*b2 + ba1.*(b4+b3-2);

%k = 5, m = 1 
%b = dec2bin(2^6-1:-1:0)-'0';
%b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);b5=b(:,5);ba1=b(:,6);
%LHS = b1.*b2.*b3.*b4.*b5;
%RHS = + b4.*b2 + ba1.*(b5-1);
%RHS = + b4.*b2 + ba1.*(b5+b1-2);

%k = 8, m = 2 
%b = dec2bin(2^10-1:-1:0)-'0';
%b1=b(:,1);b2=b(:,2);b3=b(:,3);b4=b(:,4);b5=b(:,5);b6=b(:,6);b7=b(:,7);b8=b(:,8);ba1=b(:,9);ba2=b(:,10);
%LHS = b1.*b2.*b3.*b4.*b5.*b6.*b7.*b8;
%RHS = + b4.*b2 + ba1.*(b5-1)+ ba2.*(b8+b7-2); 
%}

