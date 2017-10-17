function y = trapfilter(input_signal,remove_f)
fs=125;%����Ƶ��
Ts=1/fs;
x=input_signal; 
f0=remove_f/60; 
NLen=length(x); 
n=0:NLen-1;
apha=-2*cos(2*pi*f0*Ts);
beta=0.96;
b=[1 apha 1];
a=[1 apha*beta beta^2];
y=dlsim(b,a,x); 