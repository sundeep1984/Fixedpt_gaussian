%script to segment the positive portion of standard gaussian(N(0,1)) ICDF 
%into segments which can be used to generate gaussian samples as described
%in R. Gutierrez, V. Torres, J.Valls, �Hardware Architecture of a Gaussian Noise Generator Based on
%Inversion Method, � IEEE Transactions on Circuits & Systems II, Vol.59, No.8, pp.501-505, 2012

% Since we focus only on region where the gaussian support is positive, the
% ICDF range is from [0.5,1]

% we use norminv to get ICDF values. norminv(0.5)=0 and norminv(1)=inf

% Copyright (C) 2017  Sundeep Venkatraman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% starting inputs
a = 0.5; % a and b are the ends of the interval where the function is to be approximated
b = 1;
f = 1; % f=1 corresponds to ICDF of gaussian

ulp = 1/2^25;
e_max = ulp; % 

%temp variable assignments
coeffvec=[];% should store c0, c1,c2 as well as the x1 and x2 bounds
u=[];
%a=[0.5 0.8853];
%b=1;

%for i=1:2
x1 = a;
x2 = b;
m = 1;
done = 0;
check_x2 = 0;
prev_x2 = a;
oscillating = 0;


while ~done
    [error,c2,c1,c0]=minimax(f,x1,x2,m); % minimax does a least squares fit of a degree 2 polynomial to f over interval [x1,x2]
    
    if error<=e_max
        if x2==b
           done = 1;           
           coeffvec = [coeffvec; [c0, c1, c2] ];
           u = [u; x2];
           
        else
            if oscillating
                coeffvec = [coeffvec; [c0, c1, c2] ];
                if c0==inf | c0==-inf
                    m
                end
                u = [u; x2];
                prev_x2 = x2;
                x1 = x2;
                x2 = b;
                m = m+1;
                oscillating = 0;
            else
                change_x2=abs(x2-prev_x2)/2;
                prev_x2 = x2;
                
                if change_x2 > ulp
                    x2 = x2 + change_x2;
                else
                    x2 = x2 + ulp;
                end
            end
                    
        end
        
        
    else
        change_x2 = abs(x2-prev_x2)/2;
        prev_x2 = x2;
        if change_x2 > ulp
            x2 = x2 - change_x2;
        else
            x2 = x2 - ulp;
            if check_x2 == x2
                oscillating = 1;
            else
                check_x2 = x2;
            end
        end
        
    end
    
    
    
end
done
%u(end)=[];
%coeffvec(end,:)=[];
%ulp = 1/2^10;
%e_max = ulp; % 

%end
%l=length(u);

%div=max(61-l,0);

% if div>0
% interval = (u(end)-u(end-1))/(div+1);
% 
% 
% 
% int=u(end-1):interval:u(end);
% u(end)=[];
% coeffvec(end,:)=[];
% 
% for i=1:div
%     [error,c2,c1,c0]=minimax(f,int(i),int(i+1),m);
%     error
%     coeffvec = [coeffvec; [c0, c1, c2] ];
%     u = [u; int(i+1)];
%     m=m+1;
% end
% 
% end

s = size(coeffvec);

u'

c0binvec = zeros(s(1),21);
c1binvec = zeros(s(1),18);
c2binvec = zeros(s(1),18);
%cvec=dec2bin(zeros(s(1),s(2)),18);
cvec=[];
cnt=1;
% in the code below, n1,n2 are selected based on the range of values of the
% coefficients c0,c1,c2
for i=1:s(1)
    for j=1:s(2)
        
        if j==1
            n1=1;% encoding the coeff as Q1.17
            n2= 17;
            if coeffvec(i,j)<0
                temp = 2^n1+coeffvec(i,j); % 2's complement for -ve number
            else 
                temp = coeffvec(i,j);
            end
            c2binvec(cnt,:) = fliplr(de2bi(round(temp*2^n2),n1+n2));% round it out and represent it in 21 bits
            
        elseif j==2
            n1=1;% encoding the coeff as Q1.17
            n2=17;
            if coeffvec(i,j)<0
                temp = 2^n1+coeffvec(i,j);
            else 
                temp = coeffvec(i,j);
            end
            c1binvec(cnt,:) = fliplr(de2bi(round(temp*2^n2),n1+n2));
            
        elseif j==3
            n1=2;% encoding the coeff as Q2.19
            n2=19;
            if coeffvec(i,j)<0
                temp = 2^n1+coeffvec(i,j);
            else 
                temp = coeffvec(i,j);
            end
            c0binvec(cnt,:) = fliplr(de2bi(round(temp*2^n2),n1+n2));
            
        end
                
    end
    cnt=cnt+1;

end

save coeffs.mat c0binvec c1binvec c2binvec