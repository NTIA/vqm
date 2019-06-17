function [si_orig, part_si_min, part_si_max, ...
        hv_feat_orig, part_hv_min, part_hv_max, y_orig, ...
        cb_orig, cr_orig, part_c_min, part_c_max, part_c, ...
        ati_orig, part_ati_min, part_ati_max, part_ati, code_ati ] ...
 = model_lowbw_compression (mode, transmit_file, si_orig, hv_feat_orig, luma_orig, cb_orig, cr_orig, ati_orig);
% MODEL_LOWBW_COMPRESSION
%   Compress and write original features from the LowBW or FastLowBW model.
%   Alternatively reverse that process: read and uncompress features.
% SYNTAX
%   model_lowbw_compression ('compress', transmit_file, si_orig, hv_feat_orig, ...
%       luma_orig, cb_orig, cr_orig, ati_orig );
%   [si_orig, part_si_min, part_si_max, ...
%          hv_feat_orig, part_hv_min, part_hv_max, y_orig, ...
%          cb_orig, cr_orig, part_c_min, part_c_max, part_c, ...
%          ati_orig, part_ati_min, part_ati_max, part_ati, code_ati ] = ...
%       model_lowbw_compression('uncompress', transmit_file);
% DESCRIPTION
%   When the first argument is 'compress', this function compresses the 
%   features associated with the lowbw model or fastlowbw model, and write
%   them to a file for transmission.  Operates on 255 seconds or less.  Longer  
%   sequencesmust be split and written to two different files.  To
%   understand the input arguments, see function model_lowbw_features or
%   model_fastlowbw_features. 
%
%   When the first argument is 'uncompress', this function reverses the
%   process: reads compressed original features from file 'transmit_file'
%   and returns those values on the command line.  This function also
%   returns information about the feature compression, which will be needed
%   by function model_lowbw_parameters.m and model_fastlowbw_parameters.m
%   To properly understand those parameters (i.e., everything returned that
%   isn't mentioned in the 'compress' argument list), examine the code
%   below.

if strcmpi(mode,'compress'),
    [row_si,col_si,time_si] = size(si_orig);
    [row_cb,col_cb,time_cb] = size(cb_orig);
    [si_orig]      = model_lowbw_compression_internal (si_orig,      'si', 'compress');
    [hv_feat_orig] = model_lowbw_compression_internal (hv_feat_orig, 'hv', 'compress');
    [luma_orig]    = model_lowbw_compression_internal (luma_orig,    'y', 'compress');
    [cb_orig]      = model_lowbw_compression_internal (cb_orig,      'cb', 'compress');
    [cr_orig]      = model_lowbw_compression_internal (cr_orig,      'cr', 'compress');
    [ati_orig]     = model_lowbw_compression_internal (ati_orig,     'ati', 'compress');
    
    if time_si > 255 | time_cb > 255,
        error('Time sequence too long.  Split into shorter segments');
    end

    % save data to file
    if exist(transmit_file,'file'),
        delete(transmit_file);
    end
    fid = fopen(transmit_file, 'w');

    % write overhead for SI & HV & Luma
    row_si = uint8(row_si); col_si = uint8(col_si); time_si= uint8(time_si);
    fwrite(fid,row_si,'uint8');
    fwrite(fid,col_si,'uint8');
    fwrite(fid,time_si,'uint8');

    % write out SI & HV & Luma
    fwrite(fid,si_orig, 'ubit9');
    fwrite(fid,hv_feat_orig,'ubit9');
    fwrite(fid,luma_orig, 'ubit8');

    % write Cb & Cr
    fwrite(fid,cb_orig, 'ubit9');
    fwrite(fid,cr_orig, 'ubit9');

    % write overhead for ATI 
    frames = length(ati_orig);
    fwrite(fid,frames,'ushort');

    % write ATI
    fwrite(fid,ati_orig, 'ubit10');

%     fprintf('Total bytes written:  %d  Sequence Size:  (%d,%d,%d) ==> %f kbits/s\n', ftell(fid), ...
%         row_si,col_si,time_si, (ftell(fid) * 8.0) / (double(time_si) * 1000.0));

    fclose(fid);
elseif strcmpi(mode,'uncompress'),
    % load original features
    if ~exist(transmit_file,'file'),
        error(sprintf('Cannot open file ''%s'', file does not exist', transmit_file));
    end
    fid = fopen(transmit_file, 'r');

    % read overhead for SI & HV & Luma
    row = fread(fid,1, 'uint8');
    col = fread(fid,1,'uint8');
    time = fread(fid,1,'uint8');

    % read Si & HV & Luma
    si_orig = fread(fid,row*col*time, 'ubit9');
    hv_feat_orig = fread(fid,row*col*time,'ubit9');
    y_orig = fread(fid,row*col*time, 'ubit8');

    % decompress original features
    [si_orig, part_si_min, part_si_max]      = ...
            model_lowbw_compression_internal (si_orig, 'si', 'decompress', row,col,time);
    [hv_feat_orig, part_hv_min, part_hv_max] = ...
            model_lowbw_compression_internal (hv_feat_orig, 'hv', 'decompress',row,col,time);
    [y_orig]    = model_lowbw_compression_internal (y_orig, 'y', 'decompress',row,col,time);

    % read Cb & Cr
    cb_orig = fread(fid,row*col*time, 'ubit9');
    cr_orig = fread(fid,row*col*time, 'ubit9');

    % decompress original features
    [cb_orig]    = model_lowbw_compression_internal (cb_orig, 'cb', 'decompress',row,col,time);
    [cr_orig, part_c_min, part_c_max, part_c]    = ...
            model_lowbw_compression_internal (cr_orig, 'cr', 'decompress',row,col,time);

    % read overhead for ATI & read ATI 
    frames = fread(fid,1,'ushort');
    ati_orig = fread(fid,frames, 'ubit10');

    % decompress original features
    [ati_orig, part_ati_min, part_ati_max, part_ati, code_ati] = model_lowbw_compression_internal (ati_orig, 'ati', 'decompress',1,1,frames);

    fclose(fid);
else
    error('mode not recognized');
end




function [data, part_min, part_max, part, code] = ...
    model_lowbw_compression_internal (feature, type, flag, row, col, time);
% model_lowbw_compression_internal
%   Compress or decompress the lowbw/sec model features.
% SYNTAX
%   [data] = model_lowbw_compression_internal (feature, type, 'compress');
%   [data, part_min, part_max, part, code] = ...
%       model_lowbw_compression_internal (feature, type, 'decompress', row, col, time);
% DESCRIPTION
%   Compress or decompress the features in matrix 'feature'.  Variable
%   'type' controls the compression (e.g., 'hv', 'si','y','cb', 'cr', or 'ati).
%   'flag' should be set to either 'compress' or 'decompress' or 'none'.
%       ('none' is no compression; return data as-is for saving.)
%
%   Return the compressed or decompressed data in [data].  On
%   decompression, also return the minimum & maximum partition boundaries
%   in 'part_min' and 'part_max'.  On compression, [data] will be a one
%   dimensional array.  On decompression, [data] will be of size
%   (row,col,time).  
%
%   return value 'part' is the partitions; 'code' is the code values

% get quantizer codebook
if strcmpi(type,'hv'),
    [part,code] = model_lowbw_compression_internal_hv;
elseif strcmpi(type,'si'),
    [part,code] = model_lowbw_compression_internal_si;
elseif strcmpi(type,'y'),
    [part,code] = model_lowbw_compression_internal_y;
elseif strcmpi(type,'cb') | strcmp(type,'cr'),
    [part,code] = model_lowbw_compression_internal_cb_cr;
elseif strcmpi(type,'ati'),
    [part,code] = model_lowbw_compression_internal_ati;
end

if strcmpi(type,'none'),
    data = feature;
    return;
elseif strcmpi(flag,'compress'),

    %  convert from 3D to 1D
    [frows,fcols,ftime] = size(feature);
    feature = reshape(feature,1,frows*fcols*ftime);  % row vector for quantizer

    %  Quantize feature
    [index] = quantiz_fast(feature,part);
    data = index;
    return;
    
elseif strcmpi(flag,'decompress'),
    part_min = part(1);
    part_max = part( length(part) );
    
    % convert back from indicies to quantized values
    feature = code(feature+1);
    
    data = reshape(feature,row,col,time);
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [part_hv,code_hv] = model_lowbw_compression_internal_hv;

    %  hv_feat_orig quantizer, 9 bit design
    %  Use a quantizer whose distance between bins (dx) increases as a constant
    %  (err) times the code value.  Only preserve this relationship until the
    %  code value decreases to some low value, then uniformly quantize the
    %  remainder.  Start at 1.0 so the code is exact there.
    %  First bin will be values less than part_hv(1) and last bin will be 
    %  values greater than part_hv(code_hv_size-1).  hv_loss and hv_gain parameter values
    %  whose hv_feat_orig values fall outside of the quantizer range will be 
    %  zeroed.
    start = 1.0;
    err = 0.00709;  % fraction of error in the quantized value
    high_codes = 228;  % number of nonlinear codes >= 1.0
    low_codes = 202;  % number of nonlinear codes < 1.0
    vlow_codes = 82;  % number of uniform codes for very low values
    code_hv = [start];
    %  Generate high codes
    for i = 1:high_codes-1
        code_hv = [code_hv code_hv(i)*(1+err)];
    end
    %  Generate low codes
    for i = 1:low_codes
        code_hv = [code_hv(1)*(1-err) code_hv];
    end
    %  Generate very low codes, uniformly distributed to lowest_code
    lowest_code = 0.0991;
    temp = code_hv(1):(lowest_code-1*code_hv(1))/vlow_codes:lowest_code;
    code_hv = [fliplr(temp(2:vlow_codes+1)) code_hv];
    %  Generate the partitions, halfway between codes
    code_hv_size = size(code_hv,2);
    part_hv = (code_hv(2:code_hv_size)+code_hv(1:code_hv_size-1))/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [part_si,code_si] = model_lowbw_compression_internal_si;

    %  si_orig quantizer, 9 bit design
    %  Use a quantizer whose distance between bins (dx) increases as a constant
    %  (err) times the code value.
    %  First bin will be values less than part_si(1) and last bin will be 
    %  values greater than part_si(code_si_size-1).  si_loss and si_gain parameter values
    %  whose si_orig values fall outside of the quantizer range will be 
    %  zeroed. 
    start = 2.99;
    err = 0.00728;  % fraction of error in the quantized value
    high_codes = 512;  % number of nonlinear codes
    code_si = [start];
    %  Generate high codes
    for i = 1:high_codes-1
        code_si = [code_si code_si(i)*(1+err)];
    end
    %  Generate the partitions, halfway between codes
    code_si_size = size(code_si,2);
    part_si = (code_si(2:code_si_size)+code_si(1:code_si_size-1))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [part_y, code_y] = model_lowbw_compression_internal_y;

    %  For y quantizer, just use uniform 8 bit quantizer where the value
    %  represents the y level (0 to 255).  Only one byte per
    %  macroblock(3,3,2).
    start = 0.0;
    err = 1.00;  % distance between codes
    high_codes = 256;  % number of codes
    code_y = [start];
    %  Generate high codes
    for i = 1:high_codes-1
        code_y = [code_y code_y(i)+err];
    end
    %  Generate the partitions, halfway between codes
    code_y_size = size(code_y,2);
    part_y = (code_y(2:code_y_size)+code_y(1:code_y_size-1))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [part_c, code_c] = model_lowbw_compression_internal_cb_cr;

    %  Non-linear quantizer for cb_orig and cr_orig, similar to hv but bipolar (two sided),
    %  since the cb and cr features are bipolar.
    start = 1.0;
    err = 0.0216;  % fraction of error in the quantized value
    high_codes = 217;  % number of nonlinear codes >= start
    low_codes = 40;  % number of uniform codes for very low values (< start)
    code_c = [start];
    %  Generate high codes
    for i = 1:high_codes-1
        code_c = [code_c code_c(i)*(1+err)];
    end
    %  Generate low codes, uniformly distributed from lowest_code to start
    lowest_code = 0.136;
    temp = code_c(1):(lowest_code-1*code_c(1))/low_codes:lowest_code;
    code_c = [fliplr(temp(2:low_codes+1)) code_c];

    %  zero first code and generate sym neg codes
    code_c(1) = 0.0;
    code_c = [-1.0*fliplr(code_c(2:(high_codes+low_codes-1))) code_c];

    %  Generate the partitions, halfway between codes
    code_c_size = size(code_c,2);
    part_c = (code_c(2:code_c_size)+code_c(1:code_c_size-1))/2;

    %  Correct partition around zero
    low_part_spacing = (start-lowest_code)/low_codes;
    part_c(code_c_size/2) = part_c(code_c_size/2 + 1) - low_part_spacing;
    part_c(code_c_size/2 - 1) = part_c(code_c_size/2 - 2) + low_part_spacing;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [part_ati, code_ati] = model_lowbw_compression_internal_ati;
%  ati_orig quantizer, 10 bit design, uniform quantizer design 
%  First bin will be values less than part_ati(1) and last bin will be 
%  values greater than part_ati(code_ati_size-1).
    start = 0.0;  % first code
    last = 220.0;  % 210 is the maximum observed in the training data 
    high_codes = 1024;  % number of codes for 10-bit quantizer
    code_ati = start:(last-start)/(high_codes-1):last;

    %  Generate the partitions, halfway between codes
    code_ati_size = size(code_ati,2);
    part_ati = (code_ati(2:code_ati_size)+code_ati(1:code_ati_size-1))/2;
