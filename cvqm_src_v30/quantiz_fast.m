function y = quantiz_fast(x,t)
% QUANTIZ_FAST
%   This function quantizes data.  It is used as an alternative to the
%   MATLAB routine, becase at the time of writing it ran significantly
%   faster.
% SYNTAX
%   y = quantiz_fast(x,t)
% DESCRIPTION
%   'x' is the data to be quantized into discreet bins.  
%   't' is the array of quantization thresholds.  
%   Return varaible 'y' contains the quantized data. 
%
%  This code performs n=2^b level quantization of the signal samples in the
%  row vector x to produce y.  y will have the same size as x.  t is a
%  length n-1 row or column vector of  quantization thresholds.
%  This code works for all b=1,2,3,...
%  Quantization rules are:
%           x <= t(1)   gives y = 0
%    t(1) < x <= t(2)   gives y = 1
%                ...
%    t(i) < x <= t(i+1) gives y = i
%              ...
%    t(n-1) < x         gives y = n-1


n = length(x);
y = zeros(1,n);           %output variable starts with all zeros
b = log2(length(t)+1);    %number of bits
for i = 1:b               %loop over bits, each pass determines one bit
    temp = find(t(y+(2^(b-i))) < x);
    %compare each sample of x with corresponding current threshold    
    y(temp) = y(temp)+(2^(b-i));
    %where x is above current threshold, update y to refect this
    %this update is equivalent to setting the current bit under
    %consideration (b-i) to 1.    
end
