function S = ndec2base(D,B,N)
%NDEC2BASE Convert decimal integer to base B string.
%   NDEC2BASE(D,B) returns the representation of D as a string in
%   base B.  D must be a non-negative integer array smaller than 2^52
%   and B must be an integer between 2 and 36.
%
%   FDEC2BASE(D,B,N) produces a representation with at least N digits.
%
%   NDEC2BASE is designed to be a direct replacement for dec2base. It is
%   optimized for speed and conservation of memory.
% 
%   Examples
%       fdec2base(23,3) returns '212'
%       fdec2base(23,3,5) returns '00212'
%

narginchk(2,3);

D = D(:);
if ~(isnumeric(D) || ischar(D))
    error('First argument must be an array of integers, 0 <= D <= 2^52.');
end
if ~isscalar(B) || ~(isnumeric(B) || ischar(B)) || B ~= floor(B) || B < 2 || B > 36
    error('Second argument must be an integer, 2 <= B <= 36.');
end
if nargin == 3
    if ~isscalar(N) || ~(isnumeric(N) || ischar(N)) || N ~= floor(N) || N < 0
        error('Third argument must be a positive integer.');
    end
end

D = double(D);
B = double(B);
l = length(D);

maxd = max(D);
no = max(1,round(log2(maxd+1)/log2(B)));
while B^no <= maxd, no=no+1; end

if nargin == 3 && N>no, nmax = N;else nmax = no; end
S = zeros(l,no,'uint8');

% for large arrays a lot faster than matrix multiplication
while no > 1
    no  = no - 1;
    base = B^no;
    col  = nmax-no;
    S(:,col) = uint8(D/base-0.5);
    D  = D - double(S(:,col))*base;
end

S(:,nmax) = uint8(D);
if B>10
    symbols = uint8('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ');
    % conserve memory + faster than s = symbols(s+1);
    S=S+1;
    maxn = size(S,1);
    for a = 1:10000:maxn
        ind = a:min(a+9999,maxn);
        S(ind,:) = symbols(S(ind,:));
    end
    S = char(S);
else
    S = S+uint8('0');
    S = char(S);
end
end