%harmonic oscillator
clc

% for n=0:7
% l=n:-2:0;
% end
x=-3:0.1:3

V1=0.5*x.^2;

V2(n,i)=V1-0.225* 1;
plot(x,V1,'r',x,V2,'g')

title('Harmonic Oscillator Potential')
xlabel('X')
ylabel('Potential (V)')