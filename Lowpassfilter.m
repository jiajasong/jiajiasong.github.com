
function y = Lowpassfilter(x,M)
N = length(x);
M=M-mod(M,2);
for i=M+1:N
    y(i)=sum(x(i-M+1:i))/M;
end
ytmp=y(M+1:end);
ytmp1=ytmp(1:M/2)+ytmp(1)-ytmp(M/2+1);
ytmp2=ytmp(N-M-M/2+1:N-M)+ytmp(N-M)-ytmp(N-M-M/2);
y=[ytmp1 ytmp ytmp2];
