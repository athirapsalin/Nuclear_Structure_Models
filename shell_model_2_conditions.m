x1=0:3;
plot(x1,0.*ones(size(x1)))
 C=-0.1;
D=-0.0225;
x2=(5:8);  
 hold on
ylabel('ENERGY')

for N=0:1:6
   E1=N+1.5;
    plot(x1,E1.*ones(size(x1)))
    hold on
     end
  

for i=0:6

    if mod(N,2)==0
       for l=i:-2:0
           if l==0
           E7=(i+1.5);
           plot(x2,E7.*ones(size(x2)))
           break
           end
        E2=(i+1.5)+D.*l.*(l+1);
        E3=E2+C.*(0.5.*l);
        E4=E2-C.*(l+1).*0.5;
        
        plot(x2,E3.*ones(size(x2)))
        hold on
        plot(x2,E4.*ones(size(x2)))
        hold on
       end
    
        else 
        for l=i:-2:0
        
        E5=E2+C.*(0.5.*l);
        E6=E5-C.*(l+1).*0.5;
        plot(x2,E5.*ones(size(x2)))
        hold on
        plot(x2,E6.*ones(size(x2)))
       hold on
       end
        
    end
end
x3=8:10;
    
plot(x3,0.*ones(size(x3)))