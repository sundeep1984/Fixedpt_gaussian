% add two fixed pt binary numbers 

function c=addbin(a,na1,na2,b,nb1,nb2,nc1,nc2)% a is Qna1.na2, b is Qnb1.nb2 and c is Qnc1.nc2

if na1>nb1
    b = [zeros(1,na1-nb1) b];
    len1 = na1;
else
    a = [zeros(1,nb1-na1) a];
    len1 = nb1;
end


if na2>nb2
    b = [b zeros(1,na2-nb2)];
    len2 = na2;
else
    a = [a zeros(1,nb2-na2)];
    len2 = nb2;
end

ctemp = zeros(1,len1+len2);

carry=0;
for i=(len1+len2):-1:1
    
    c(i) = xor(xor(a(i),b(i)),carry);
    carry = (a(i) & b(i)) | (carry & b(i)) | (a(i) & carry);

    
end

if nc1>=len1
    c=[repmat(c(1),1,nc1-len1) c];
else
    disp('nc1 needs to be large enough to hold the sum of the MSBs of a and b without truncation or overflow');
end

if nc2>=len2
    c=[c repmat(c(1),1,nc2-len2)];
else
    round_bit=1;% used for rounding out to nc2-th place after radix point
    carry=0;
    for i=(nc1+nc2):-1:1
        c(i) = xor(xor(c(i),round_bit),carry);
        carry = ((c(i) & round_bit)|(round_bit & carry)|(carry | c(i)));
        round_bit=0;
    end
    
    c = c(1:(nc1+nc2));
end

