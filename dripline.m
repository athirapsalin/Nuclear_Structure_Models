
av=15.5;
as=16.8;
ac=0.72;
aa=23;

     A=2:250;
     Z1=A./(2+(0.0156.*A.^(2/3)));
     Z=round(Z1); 
     N=A-Z;
     plot(N,Z)
     xlabel('Neutron Number (N)')
     ylabel('Proton Number (Z)')
     title('Nuclear Dripline')
     hold on


syms E1(Z_p) E4(N_n)
%%proton dripline
for N1=1:155
 E1=(av-(as*(N1+Z_p)^(-1/3))-(ac*Z_p*(Z_p-1)*(N1+Z_p)^(-4/3))-((aa.*((N1-Z_p)^2))/(N1+Z_p)^2))*(N1+Z_p);
 E2=diff(E1,Z_p);
 E3=E2==0;
 Z_p_1(N1)=vpasolve(E3,Z_p);
 N_p(N1)=N1;
end

 plot(N_p,floor(Z_p_1))
 hold on

for Z1=1:108
 E4=(av-(as*(N_n+Z1)^(-1/3))-(ac*Z1*(Z1-1)*(N_n+Z1)^(-4/3))-((aa.*((N_n-Z1)^2))/(N_n+Z1)^2))*(N_n+Z1);
 E5=diff(E4,N_n);
 E6=E5==0;
 N_n_1(Z1)=vpasolve(E6,N_n);
 Z_n(Z1)=Z1;
end
plot(floor(N_n_1),Z_n)
hold on
% not necessary to define as a function

