function p = gm_1_1_using(X,k) 
%gray model: GM(1,1)
%Example p = gm_1_1([200 250 300 350],2)
if nargout>3,error('Too many output argument.');end 
if nargin==1,k=1;x_orig=X; 
elseif nargin==0|nargin>2 
   error('Wrong number of input arguments.'); 
end 
if length(unique(X))==1
    p = unique(X);
else
    x_orig=X; 
    predict=k;

    %AGO process 
    x=cumsum(x_orig);


    %compute the coefficient(a and u)------------------------ 
    n=length(x_orig); 
    %first generate the matrix B 
    for i=1:(n-1); 
    B(i)=-(x(i)+x(i+1))/2; 
    end 
    B=[B' ones(n-1,1)]; 
    %then generate the matrix Y 
    for i=1:(n-1); 
    y(i)=x_orig(i+1); 
    end 
    Y=y'; 
    %get the coefficient. a=au(1) u=au(2) 
    au=(inv(B'*B))*(B'*Y); 
    %-------------------------------------------------------- 
    %change the grey model to symbolic expression 
    coef1=au(2)/au(1); 
    coef2=x_orig(1)-coef1; 
    coef3=0-au(1); 
    costr1=num2str(coef1); 
    costr2=num2str(abs(coef2)); 
    costr3=num2str(coef3); 
    eq=strcat(costr1,'+',costr2,'e^(',costr3,'*(t-1))');

    %comparison of calculated and observed value 
    for t=1:n+predict 
       mcv(t)=coef1+coef2*exp(coef3*(t-1)); 
    end 
    x_mcv0=diff(mcv); 
    x_mcve=[x_orig(1) x_mcv0]; 
    x_mcv=diff(mcv(1:end-predict)); 
    x_orig_n=x_orig(2:end); 
    x_c_error=x_orig_n-x_mcv; 
    x_error=mean(abs(x_c_error./x_orig_n));

    if x_error>0.2 
       disp('model disqualification!'); 
    elseif x_error>0.1 
       disp('model check out'); 
    else 
       disp('model is perfect!'); 
    end 
    %predicting model and plot gragh 
    % plot(1:n,x_orig,'diamond',1:n+predict,x_mcve);
    p=x_mcve(end-predict+1:end); 
    % xlabel('CURVE OF GREY MODEL ANALYSIS'); 
    % title('GM(1,1)'); 
    % grid on 
    y=eq; 
    e=x_error; 
    p=x_mcve(end-predict+1:end);
end

