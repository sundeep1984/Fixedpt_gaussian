function [error,c2,c1,c0] = minimax(f,x1,x2,m)

%if m<188
if m<47
    step=1/2^15;
   % num_div = 2^15;
    temp=15;
else
    temp = m-46;
    step=1/2^(15-temp);
        %temp = floor((m-188)/4)+1;
   % num_div = 2^(15-temp);
end
%x=x1:((x2-x1)/num_div):x2; % the vector of X values
x= 0:step:(x2-x1-step); % the vector of X values
x=x+x1;


if length(x)>1
X = [ ones(length(x),1) x' (x.^2)'];

if f==1
    Y=norminv(x,0,1);
    
end


    p=polyfit(x,Y,2);
    

error= max(Y'-(X*flipud(p')));
%p=p/(2^temp);
a=p(1);%LLS_c(3);
b=p(2);%LLS_c(2);
c=p(3);%LLS_c(1);
%c2=a;
%c1=b;
%c0=c;
c2 = a*(x2-x1)^2;
c1 = 2*a*x1*(x2-x1) + b*(x2-x1);
c0 = a*x1^2+b*x1+c;

% if temp<15 
%    c2 =  c2/(2^(2*temp));
%    c1 = c1/(2^temp);
% end

else 
    error=0;
    c2=0;
    c1=0;
    c0=norminv(x1,0,1)';
end


