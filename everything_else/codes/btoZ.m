function latexstr = btoZ(fb, quadb)
% btoZ('- b1*b2*b3*b4 - b2*b3*b4*b5 - b3*b4*b5*b1 - b4*b5*b1*b2 - b5*b1*b2*b3 + 3*b1*b2*b3*b4*b5', '3*ba - b1*ba - b2*ba - b3*ba - b4*ba - b5*ba') 

fb = str2sym(fb);
quadb = str2sym(quadb);
n = length(symvar(fb + quadb)) - 1;

B = sym('b', [1 n]);
syms ba
Zs = sym('Z', [1 n]);
syms Z_a;

degf = polynomialDegree(fb);

for k = 1:n
    fb = subs(fb, B(k), (1-Zs(k))/2);
    quadb = subs(quadb, B(k), (1-Zs(k))/2);
end
fb = subs(fb, ba, (1-Z_a)/2);
quadb = subs(quadb, ba, (1-Z_a)/2);
fz = expand(fb)*(2^degf);
quadz = expand(quadb)*(2^degf);
[Cf, Tf] = coeffs(fz);

l = length(Tf);
for k = 1:l 
    if polynomialDegree(Tf(k)) <= 2
        quadz = quadz - Cf(k) * Tf(k);
        fz = fz - Cf(k) * Tf(k);
    end
end

latexstr = sprintf('%s \\\\\\mapsto %s', latex(fz), latex(quadz));
end