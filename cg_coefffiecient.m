function cgc1=cgc(j1,j2,j,m1,m2,m)
% prompt='enter j1 ';
% j1=input(prompt)
% prompt='enter j2 ';
% j2=input(prompt)
% prompt='enter m1 ';
% m1=input(prompt)
% prompt='enter m2 ';
% m2=input(prompt)
% prompt='enter j ';
% j=input(prompt)
% prompt='enter m';
% m=input(prompt)


if (j1-m1)>=0&&(j2-m2)>=0&&(j-m)>=0&&(j1+j2-j)>=0&&(j+j1-j2)>=0&&(j+j2-j1)>=0
c1=(((2.*j+1).*(factorial(j1+j2-j)*factorial(j+j1-j2)*factorial(j+j2-j1))./factorial(j1+j2+j+1)).^(1/2))   
c2=(factorial(j1+m1).*factorial(j1-m1).*factorial(j2+m2).*factorial(j2-m2).*factorial(j+m).*factorial(j-m)).^(1/2)

c3=0;
for n=(0:10);
    
    t1=(j1+j2-j-n);
    t2=(j1-m1-n);
    t3=(j2+m2-n);
    t4=(j-j2+m1+n);
    t5=(j-j1-m2+n);
    v=m1+m2;
    if t1>=0&&t2>=0&&t3>=0&&t4>=0&&t5>=0
      
        c3=c3+(((-1)^(n))./(factorial(n)))*((factorial(n)*factorial(t1)*factorial(t2)*factorial(t3)*factorial(t4)*factorial(t5)))^(-1)
    end
        cgc1=eq(m,v)*c1*c2*c3
   
end
else
    cgc1=0

end
end
