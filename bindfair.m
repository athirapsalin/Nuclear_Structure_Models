function bindfair()
av=15.5;
as=16.8;
ac=0.72;
aa=23;
ap=34;
B=0;
 A=2:250;
     Z1=A./(2+(0.0156.*A.^(2/3)));
     Z=round(Z1); 
     if mod(Z,2)==0 & mod(A,2)==0
         delta=ap.*A.^(-7/4);
     elseif mod(Z,2)==1 & mod(A,2)==0
         delta=-ap.*A.^(-7/4);
     else
         delta=0;
     end
     B=av.*A./A-(as.*A.^(-1/3))-(ac.*Z.*(Z-1).*A.^(-4/3))-((aa.*((A-2.*Z).^2))./A.^2)+delta;
     
     plot(A,B)
     xlabel('Mass number (A)')
     ylabel('Binding Energy per Nucleon (BE/A) in MeV')
     title('Binding Energy Curve')
     hold on

end

% not necessary to define as a function

