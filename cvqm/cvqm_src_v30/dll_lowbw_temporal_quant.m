function [ti2, ti10, y] = dll_lowbw_temporal_quant(is_quantize, orig_ti2, orig_ti10, orig_y);
% DLL_LOWBW_TEMPORAL_QUANT
%   Quantize & reconstruct original features for low bandwidth temporal
%   registration.
% SYNTAX
%  [index_ti2, index_ti10, index_y] = dll_lowbw_temporal_quant(is_quantize, orig_ti2, orig_ti10, orig_y);
%  [orig_ti2, orig_ti10, orig_y] = dll_lowbw_temporal_quant(is_quantize, index_ti2, index_ti10, index_y);
% DESCRIPTION
%   'is_quantize' is 1 for quantize, 0 to reconstruct.
%   When quantizing, the other three paramters are the features returned by 
%   function dll_lowbw_temporal_original (orig_ti2, orig_ti10, orig_y) 
%   and the return values are the indexes to be transmitted.  When reconstructing,
%   the other three parameters are the indexes, and the return values are
%   the re-constructed features.
%
% Example call to quantize:
%   [ti2_index, ti10_index, y_index] =
%       dll_lowbw_temporal_quant(1, orig_ti2, orig_ti10, orig_y);
% Example call to reconstruct:
%   [orig_ti2, orig_ti10, orig_y] =
%       dll_lowbw_temporal_quant(0, ti2_index, ti10_index, y_index);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Define the quantizers for the three temporal registration features.
%  These designs are for 10 bit linear quantizers.

%  cont feature
start = 0.0;  % first code
last = 255.0;  % 252 is the maximum observed in the training data
high_codes = 4096;  % number of codes for 12-bit quantizer
code_cont = start:(last-start)/(high_codes-1):last;
%  Generate the partitions, halfway between codes
code_cont_size = size(code_cont,2);
part_cont = (code_cont(2:code_cont_size)+code_cont(1:code_cont_size-1))/2;

%  ti2 feature
start = 0.0;  % first code
last = 210.0;  % 209 is the maximum observed in the training data
high_codes = 4096;  % number of codes for 12-bit quantizer
code_ti2 = start:(last-start)/(high_codes-1):last;
%  Generate the partitions, halfway between codes
code_ti2_size = size(code_ti2,2);
part_ti2 = (code_ti2(2:code_ti2_size)+code_ti2(1:code_ti2_size-1))/2;

%  ti10 feature
start = 0.0;  % first code
last = 210.0;  % 209 is the maximum observed in the training data
high_codes = 4096;  % number of codes for 12-bit quantizer
code_ti10 = start:(last-start)/(high_codes-1):last;
%  Generate the partitions, halfway between codes
code_ti10_size = size(code_ti10,2);
part_ti10 = (code_ti10(2:code_ti10_size)+code_ti10(1:code_ti10_size-1))/2;


if is_quantize,
    %  Quantize original features
    [ti2] = quantiz_fast((orig_ti2),part_ti2);

    [ti10] = quantiz_fast((orig_ti10),part_ti10);

    [y] = quantiz_fast((orig_y),part_cont);
else
    %  Undo the quantization
    ti2 = code_ti2(1+orig_ti2);
    ti2 = reshape(ti2,1,length(ti2));

    ti10 = code_ti10(1+orig_ti10);
    ti10 = reshape(ti10,1,length(ti10));

    y = code_cont(1+orig_y);
    y = reshape(y,1,length(y));
end
