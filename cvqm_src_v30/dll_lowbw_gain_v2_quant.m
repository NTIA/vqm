function [out_y, out_cb, out_cr] = dll_lowbw_gain_v2_quant(is_quantize, in_y, in_cb, in_cr);
% DLL_LOWBW_GAIN_V2_QUANT
%   Quantizer & reconstruct original features for lowbw gain & offset.
% SYNTAX
%   [index_y, index_cb, index_cr] = dll_lowbw_gain_v2_quant(1, y, cb, cr); % quantize
%   [y, cb, cr] = dll_lowbw_gain_v2_quant(0, index_y, index_cb, index_cr); % reconstruct
% DESCRIPTION
%   First argument is '1' to quantize, and '0' to reconstruct.
%   'y', 'cb', and 'cr are the original features from
%                      lowbw_gain_v2_original.m (i.e., uncompressed)
%   'index_y', 'index_cb', and 'index_cr' are the quantize indicies to be transmitted.

start = 0.0;  % first code
last = 255.0;  % 255 is the maximum observed in the training data 
high_codes = 1024;  % number of codes for 10-bit quantizer, must make epsilon=1.0 
code_lgo = start:(last-start)/(high_codes-1):last;

%  Generate the partitions, halfway between codes
code_lgo_size = size(code_lgo,2);
part_lgo = (code_lgo(2:code_lgo_size)+code_lgo(1:code_lgo_size-1))/2;

if is_quantize,

    %  Quantize orig_y feature like this
    [out_y] = quantiz_fast(in_y',part_lgo);

    %  Quantize org_cb feature like this
    [out_cb] = quantiz_fast(in_cb'+128,part_lgo);

    %  Quantize org_cr feature like this
    [out_cr] = quantiz_fast(in_cr'+128,part_lgo);

else

    %  Look-up the quantized value like this
    orig2 = code_lgo(1+in_y);
    out_y = orig2';

    %  Look-up the quantized value like this
    orig2 = code_lgo(1+in_cb);
    out_cb = orig2'-128;

    %  Look-up the quantized value like this
    orig2 = code_lgo(1+in_cr);
    out_cr = orig2'-128;
end
