function c=multbin(a,na1,na2,b,nb1,nb2,nc1,nc2)

na=2*(na1+na2);
nb=2*(nb1+nb2);
clen = 2*max(nb,na)-1;
if na>=nb
   term1 = [repmat(a(1),1,length(a)) a];
   term2 = [repmat(b(1),1,2*length(a)-length(b)) b];
else
   term1 = [repmat(b(1),1,length(b)) b];
   term2 = [repmat(a(1),1,2*length(b)-length(a)) a];
end

ctemp = zeros(1,clen);
ctemp2 = zeros(1,clen);
ctemp(end-length(term1)+1:end)=(term1 & term2(end));

count=1;
for i=length(term2)-1:-1:1

    if term2(i) ~= 0
        temp = [zeros(1,length(term2)-1-count) term1 zeros(1,count)];
        carry=0;
        
        for j=length(ctemp):-1:1
            ctemp2(j) = xor(xor(temp(j),ctemp(j)),carry);
            carry = (temp(j) & ctemp(j)) | (carry & ctemp(j)) | (temp(j) & carry);
        end
        ctemp=ctemp2;
        ctemp2=0;
        
    end
    count=count+1;
    
end
%ctemp
ctemp = ctemp(end-(na+nb)/2+1:end);

if nc1>=na1+nb1
    ctemp=[repmat(ctemp(1),1,nc1-na1-nb1) ctemp];
else
    disp('nc1 needs to be large enough to hold the sum of the MSBs of a and b without truncation or overflow');
    %return ctemp;
end

if nc2>=na2+nb2
    c=[ctemp repmat(ctemp(1),1,nc2-na2-nb2)];
else
    round_bit=1;% used for rounding out to nc2-th place after radix point
    carry=0;
    ctemp2 = zeros(1,length(ctemp));
    for i=(nc1+nc2):-1:1
        ctemp2(i) = xor(xor(ctemp(i),round_bit),carry);
        carry = ((ctemp(i) & round_bit)|(round_bit & carry)|(carry & ctemp(i)));
        round_bit=0;
    end
    
    c = ctemp2(1:(nc1+nc2));
end
