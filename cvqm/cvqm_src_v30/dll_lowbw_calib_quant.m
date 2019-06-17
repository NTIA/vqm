function [orig_horiz_profile_out, orig_vert_profile_out] = ...
    dll_lowbw_calib_quant(is_quantize, orig_horiz_profile_in, orig_vert_profile_in)
% DLL_LOWBW_CALIB_QUANT
%   Quantize & reconstruct original features for low bandwidth spatial
%   registration.
% SYNTAX
%   [orig_horiz_profile, orig_vert_profile] = ...
%       dll_lowbw_calib_quant(is_quantize, orig_horiz_profile, orig_vert_profile);
% DESCRIPTION
%   'is_quantize' is 1 for quantize, 0 to reconstruct.
%   Other two paramters are profile sfrom dll_lowbw_calib_original.m (when quantizing) 
%   and indexes from previous call (when reconstructing).
%
% Note: original pixels do not need to be quantized.
%
% Example call to quantize:
%   [index_horiz_profile, index_vert_profile] =
%       dll_lowbw_calib_quant(1, orig_horiz_profile, orig_vert_profile);
% Example call to reconstruct:
%   [orig_horiz_profile, orig_vert_profile] =
%       dll_lowbw_calib_quant(0, index_horiz_profile, index_vert_profile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Define the quantizers for the two spatial registration features.
%  These designs are for 10 bit linear quantizers.

%  profile feature
start = 0.0;  % first code
last = 255.0;  % 252 is the maximum observed in the training data
high_codes = 2^16;  % number of codes for 16-bit quantizer
code_profile = start:(last-start)/(high_codes-1):last;
%  Generate the partitions, halfway between codes
code_profile_size = size(code_profile,2);
part_profile = (code_profile(2:code_profile_size)+code_profile(1:code_profile_size-1))/2;


if is_quantize,
    %  Quantize original features
    [row,col] = size(orig_horiz_profile_in);
    [orig_horiz_profile_out] = quantiz_fast(reshape(orig_horiz_profile_in,1,row*col),part_profile);
    orig_horiz_profile_out = reshape(orig_horiz_profile_out,row,col);

    [row,col] = size(orig_vert_profile_in);
    [orig_vert_profile_out] = quantiz_fast(reshape(orig_vert_profile_in,1,row*col),part_profile);
    orig_vert_profile_out = reshape(orig_vert_profile_out,row,col);
else
    %  Undo the quantization
    [row,col] = size(orig_horiz_profile_in);
    orig_horiz_profile_out = code_profile(1+orig_horiz_profile_in);
    orig_horiz_profile_out = reshape(orig_horiz_profile_out,row,col);

    [row,col] = size(orig_vert_profile_in);
    orig_vert_profile_out = code_profile(1+orig_vert_profile_in);
    orig_vert_profile_out = reshape(orig_vert_profile_out,row,col);
end
