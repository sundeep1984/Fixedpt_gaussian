% script to implement fixed point Gaussian Noise generator based on
% inversion method.
% The  goal is to take 64 bit uniform random inputs and produce 16 bit
% gaussian outputs, with 12 fractional bits

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

function [g,gdec]=Gauss_fixedpt(n,write_testvec)% function will return a length n vector of gaussian binary vectors
num_bits=64;
g=[];
gdec=[];

load coeffs.mat
%write_testvec=0;
testfile = 'vectorfile.txt';

if write_testvec==1
   fptr = fopen(testfile,'w'); 
   
   if fptr==-1
       disp('Could not open file to write vectors');
       exit;
   end
end

segments = [repmat(1:4:56,1,4) 53 53 53 53];

for i=1:n
    %urng=floor(rand(1)*2^num_bits); % generating a uniform random number in [0,1] and scaling it up to represent it as a 64 bitvector
    %urng_bin=de2bi(urng,num_bits);
    urng_bin=ceil(rand(1,num_bits)-0.5);
    lead_zero_inp=urng_bin(1:61);
    zerolocs = find(lead_zero_inp);% this is the leading zero detector which is used to select the segment and coefficients
    
    if length(zerolocs)==0
        zerolocs(1)=57;
    end

    if write_testvec==1
        urng_str=num2str(urng_bin);
        urng_str(strfind(urng_str,' '))='';
        fprintf(fptr,'%s ',urng_str);
        fprintf(fptr,'%s ',dec2bin(zerolocs(1),6));
    end
    
    offset = urng_bin(62:63);
    offsetdec = bi2de(offset);
    
    coeff=segments(zerolocs(1))+offsetdec; % picking out the segment based on the leading zero location
    
    
    if write_testvec==1
        off_str=num2str(offset);
        off_str(strfind(off_str,' '))='';
        fprintf(fptr,'%s ',off_str);
        fprintf(fptr,'%s ',dec2bin(coeff,6));
    end
    
    c0_b = c0binvec(coeff,:);%Q2.19
    
    c1_b = c1binvec(coeff,:);%Q1.17
    
    c2_b = c2binvec(coeff,:);%Q1.17
    
    
    databits = urng_bin(47:61);
    mask = ones(1,15);
    
    if zerolocs(1)>46 & zerolocs(1)<61
        
        mask(1:zerolocs(1)-46) = 0;
        databits(1:zerolocs(1)-46)=0; % masking bits 
        databits = fliplr(databits);
    end
    
    if write_testvec==1
        c0b_str=num2str(c0_b);
        c0b_str(strfind(c0b_str,' '))='';
        fprintf(fptr,'%s ',c0b_str);
        
        c1b_str=num2str(c1_b);
        c1b_str(strfind(c1b_str,' '))='';
        fprintf(fptr,'%s ',c1b_str);
        
        c2b_str=num2str(c2_b);
        c2b_str(strfind(c2b_str,' '))='';
        fprintf(fptr,'%s ',c2b_str);
        
        datab_str=num2str(databits);
        datab_str(strfind(datab_str,' '))='';
        fprintf(fptr,'%s ',datab_str);
        
        mask_str=num2str(mask);
        mask_str(strfind(mask_str,' '))='';
        fprintf(fptr,'%s ',mask_str);
        
    end
    stage1out=multbin(c2_b,1,17,databits,0,15,1,17);
    
    stage2in = addbin(stage1out,1,17,c1_b,1,17,2,15);
    stage2out = multbin(stage2in,2,15,databits,0,15,3,15);
    
    stage3out = addbin(stage2out,3,15,c0_b,2,19,4,12);
    
    signbit = urng_bin(64);
    
    if write_testvec==1
        s1o_str=num2str(stage1out);
        s1o_str(strfind(s1o_str,' '))='';
        fprintf(fptr,'%s ',s1o_str);
        
        s2i_str=num2str(stage2in);
        s2i_str(strfind(s2i_str,' '))='';
        fprintf(fptr,'%s ',s2i_str);
        
        s2o_str=num2str(stage2out);
        s2o_str(strfind(s2o_str,' '))='';
        fprintf(fptr,'%s ',s2o_str);
        
        s3o_str=num2str(stage3out);
        s3o_str(strfind(s3o_str,' '))='';
        fprintf(fptr,'%s ',s3o_str);
        
        sgn_str=num2str(signbit);
        %sgn_str(strfind(sgn_str,' '))='';
        fprintf(fptr,'%s ',sgn_str);
        
    end
    
    if signbit==1
        ftemp=-(bi2de(fliplr(stage3out))/2^12);
        if(ftemp>0)
            ftemp=-ftemp;
        end
        
        stage3out = ~stage3out;
        finalout=addbin(stage3out,4,12,[zeros(1,length(stage3out)-1) 1],4,12,4,12);
        
    else
        finalout = stage3out;
        ftemp = (bi2de(fliplr(stage3out))/2^12);
    end

    if write_testvec==1
        fin_str=num2str(finalout);
        fin_str(strfind(fin_str,' '))='';
        fprintf(fptr,'%s \n',fin_str);
    end

    g = [g; finalout];
    gdec = [gdec; ftemp];
end




