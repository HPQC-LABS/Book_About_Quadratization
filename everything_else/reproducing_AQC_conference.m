z=[1 0 ; 0 -1];
z1 = kron(z,eye(32));
z2 = kron(kron(eye(2),z),eye(16));
z3 = kron(kron(eye(4),z),eye(8));
z4 = kron(kron(eye(8),z),eye(4));
z5 = kron(kron(eye(16),z),eye(2));
za = kron(eye(32),z);


H_5local=z1*z2*z3 + z1*z2*z4 + z1*z2*z5 + z2*z3*z4 + z1*z4*z5 + z2*z3*z5 + z2*z4*z5 + z1*z2*z3*z4+z1*z2*z3*z5+z1*z2*z4*z5+z1*z3*z4*z5+z2*z3*z4*z5-3*z1*z2*z3*z4*z5;
H_2local=3*(z1+z2+z3+z4+z5)-8*za+3*(z1*z2+z1*z3+z1*z4+z2*z3+z1*z5+z2*z4+z2*z5+z3*z4+z3*z5+z4*z5) -8*(z1*za+z2*za+z3*za+z4*za+z5*za)+15;

[states_H5 energies_H5]=eig(H_5local);
[states_H2 energies_H2]=eig(H_2local);


