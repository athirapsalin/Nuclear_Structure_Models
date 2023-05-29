%%Nilsson Model%% 
%% done on 23-november-2019%%
clear all
c=-0.1;
d=-0.0225;
na=[];
la=[];
mla=[];
msa=[];
mja=[];
t1a=[];
t2a=[];
s=0.5;
delta=[];
energy=[];
A=108;
x=-0.3:0.01:0.3;
%%%%%%Basis Creation %%%%%%%%
for n=0:11
    for l=n:-2:0
        for ml=-l:l
           for ms=-0.5:0.5
              for mj=-(l+s):(l+s)                                          %jz defined using (N or l) with s bcoz defining it using   
                 if (ml+ms)==mj                                            %j will give basis states with same mj values with different j
                     if mj>=0                                              %inorder to avoid degeneracy of jz+and jz- states
                           na=[na,n];
                           la=[la,l];            
                           mla=[mla,ml];
                           msa=[msa,ms];
                           mja=[mja,mj];
                
                     end
                 end
              end
           end
        end
    end
end

ba1=[na;la;mla;msa;mja];
ba2=[na;la;mla;msa;mja]; 
%%%%Calculating Hamiltonian%%%%%%
for xi=1:length(x)
 for r1=1:length(na)   
    for r2=1:length(na)      
    
        if ba1(1,r1)==ba2(1,r2)                                                     %(N1=N2)
        if ba1(2,r1)==ba2(2,r2) && ba1(3,r1)==ba2(3,r2) && ba1(4,r1)==ba2(4,r2) %(l1=l2)&(ml1=ml2)&(ms1=ms2)
            ba1(5,66)
            t1(r1,r2)=(ba2(1,r2)+1.5);                                             %(H0)
             t2(r1,r2)=(ba2(2,r2)*(ba2(2,r2)+1));                                   %(l^2)
        else
             t1(r1,r2)=0;
             t2(r1,r2)=0;
        end
        if ba2(2,r2)==ba1(2,r1)                                                      %(l1=l2)
            t_r2(r1,r2)=(ba2(1,r2)+1.5);                                             %(r^2 for (l1=l2))
                                                                                     %(l.s terms)                                                                        
            if ba1(3,r1)==(ba2(3,r2)+1) && ba1(4,r1)==(ba2(4,r2)-1)                  %(ml1=ml2+1)&(ms1=ms2-1)
                 t3(r1,r2)=0.5*sqrt((ba2(2,r2)-ba2(3,r2))*(ba2(2,r2)+ba2(3,r2)+1));  %          
            elseif ba1(3,r1)==(ba2(3,r2)-1) && ba1(4,r1)==(ba2(4,r2)+1)              %(ml1=ml2-1)&(ms1=ms2+1)
                 t3(r1,r2)=0.5*sqrt((ba2(2,r2)+ba2(3,r2))*(ba2(2,r2)-ba2(3,r2)+1));  %
            elseif ba1(3,r1)==ba2(3,r2) && ba1(4,r1)==ba2(4,r2);                     %(ml1=ml2)&(ms1=ms2)
                 t3(r1,r2)=ba2(3,r2)*ba2(4,r2);                                      %
            else                                                                     %
                 t3(r1,r2)=0;                                                        %
            end                                                                      %
                                                                                     %r^2terms)
        elseif ba1(2,r1)==(ba2(2,r2)-2)                                              %((l1=l2-2)
           t_r2(r1,r2)=sqrt((ba2(1,r2)-ba2(2,r2)+2));                                %
        elseif ba1(2,r1)==(ba2(2,r2)+2)                                              %(l1==l2+2)
           t_r2(r1,r2)=sqrt((ba2(1,r2)-ba2(2,r2))*(ba2(1,r2)+ba2(2,r2)+3));          %
        else                                                                         %
           t3(r1,r2)=0;                                                              %                    
           t_r2(r1,r2)=0;                                                            %
        end
                                                                                     %(y_(20) terms
         if ba1(3,r1)==ba2(3,r2)                                                     %(Ml1=ml2) triangle condition for cgc
            cgc1(r1,r2)=cgc(ba2(2,r2),2,ba1(2,r1),ba2(3,r2),0,ba1(3,r1));            %
            cgc2(r1,r2)=cgc(ba2(2,r2),2,ba1(2,r1),0,0,0);                            %
            t_20_1(r1,r2)=sqrt((5*(2*ba2(2,r2))+1)/(4*3.14*((2*ba1(2,r1))+1)));      % const in Y20
            t_20(r1,r2)=t_20_1(r1,r2).*cgc1(r1,r2).*cgc2(r1,r2);                     %
         else                                                                        %
             t_20(r1,r2)=0;                                                          %
         end
                                                                                     %terms inorder to match the dimension  
        else
             t1(r1,r2)=0;
             t2(r1,r2)=0;
             t3(r1,r2)=0;
             t_r2(r1,r2)=0;
             t_20(r1,r2)=0;
            
         
        end
        t4=(4/3)*sqrt(3.14/3)*t_r2.*t_20;                                      %(r^2)*(y_20) term
        g=(((1+((2/3)*x(xi)))^2)*(1-((4/3)*x(xi))))^(-1/6);                    %(f(delta)) term
        hw0=41*A^(-1/3);%*real(g);                                               %taken real g bcoz g can have 6 values including complex
        hw=(hw0*real(g));                                                      % 
        E=(t1+(d*t2)+(c*t3)-(x(xi)*t4))*hw;                                    %Energy eigen values
    end

 end
 delta=[delta,x(xi)];                                                          %array of delta value
 energy=[energy,sort(eig(E))];                                                 %array of energy eigen values 
end

 figure(1)
 plot(delta,energy,'linewidth',1.5);
 hold on
 title('Single Particle energies as a function of deformation (\delta)')
 xlabel('\delta')
 ylabel('Energy levels')
                                                                                %cg coefficient%
function cgc=cgc(j1,j2,jt,m1,m2,m)
if (j1-m1)>=0&&(j2-m2)>=0&&(jt-m)>=0&&(j1+j2-jt)>=0&&(jt+j1-j2)>=0&&(jt+j2-j1)>=0&&(j1+m1)>=0&&(j2+m2)>=0&&(jt+m)>=0
    c1=real((abs(((2.*jt+1).*(factorial(j1+j2-jt)*factorial(jt+j1-j2)*factorial(jt+j2-j1))./factorial(j1+j2+jt+1))).^(1/2)));   
    c2=real((abs(factorial(j1+m1).*factorial(j1-m1).*factorial(j2+m2).*factorial(j2-m2).*factorial(jt+m).*factorial(jt-m))).^(1/2));
    c3=0;
    for n=(0:10)  
        t1=(j1+j2-jt-n);
        t2=(j1-m1-n);
        t3=(j2+m2-n);
        t4=(jt-j2+m1+n);
        t5=(jt-j1-m2+n);
        v=m1+m2;
        if t1>=0 && t2>=0 && t3>=0 && t4>=0 && t5>=0     
            c3=c3+(((-1)^(n))./(factorial(n)))*((factorial(n)*factorial(t1)*factorial(t2)*factorial(t3)*factorial(t4)*factorial(t5)))^(-1);
        end
        cgc=eq(m,v)*c1*c2*c3;
   
     end
else
    cgc=0;
end
end