function [new_clip_structs, status, unfiltered_clip_structs] = gain_offset_lowbw (test_structs, clip_structs,varargin)
% LUMINANCE_GAIN_OFFSET_LOWBW
%   Automatic computation of gain/offset; 2007 low-bandwidth, RR approach.
%   Improved accuracy over ITU & 2006 algorithms.  Also computes color.
% changes:  1) block_size changed for QCIF (20->10) and CIF (30->22)
%           2) Use 1/2 of blocks, choosing blocks with small original Y stdev,
%              Thus, avoid edges in original luminance plane.
%           3) return (nan) for plane's gain&offset if range of means is 
%              less than 10 (max to min) (e.g., Cr=nan, Y & Cb valid)
%           4) Compute Cb & Cr gain/offset, in addition to Y
%           5) discard blocks with clipping: Y < 2 & Y > 253, Cb/Cr < -126
%              & Cb/Cr > 126. 
%           6) change epsilon from 0.1 to 1.0 (for gain/offset fit
%           iteration)
% SYNTAX
%  [new_clip_structs, status] = gain_offset_lowbw(test_structs, clip_structs);
%  [...] = luminance_gain_offset_lowbw(...,'PropertyName',...)
%  [new_clip_structs, status, unfiltered_clip_structs] =
%                       luminance_gain_offset_lowbw(...,'HRCFilter',...);
% DESCRIPTION
%  Estimate gain/offset for each processed clip.  All clips must 
%  be temporally and spatially registered first; functions
%  'temporal_registration_lowbw' and 'spatial_registration_lowbw'
%  recommended.  
%
%  Optional properties:
%  'HRCFilter'     Median filter luminance gain & offset by HRC, to
%                  improve accuracy.  Return un-filtered in unfiltered_clip_structs.
%  'NoQuantize'    Do not quantize original features; default behavior.
%  'Quantize'      Quantize original features, as when running in-service.
%  'verbose',      Print status information to screen
%  'quiet',        Print no information to screen.  Default behavior.
%
%    Return the new, calibrated clip structure.  Also return status,
%    which has the following fields defined:
%     'status'                    0 if okay; 1 if a fatal error was encountered.
%                                 See also function "lasterr".
% 

status = 0;
try

    % figure out what optional properties have been requested.
    is_hrc_filter = 0;
    is_quantize = 0;
    verbose = 0;
    cnt = 1;
    while cnt <= nargin-2,
        if strcmpi(varargin{cnt},'hrcfilter'),
            is_hrc_filter = 1;
            cnt = cnt + 1;
        elseif strcmpi(varargin{cnt},'quantize'),
            is_quantize = 1;
            cnt = cnt + 1;
        elseif strcmpi(varargin{cnt},'noquantize'),
            is_quantize = 0;
            cnt = cnt + 1;
        elseif strcmpi(varargin{cnt},'verbose'),
            verbose = 1;
            cnt = cnt + 1;
        elseif strcmpi(varargin{cnt},'quiet'),
            verbose = 0;
            cnt = cnt + 1;
        else
            error('Property Name not recognized.  Aborting.');
        end
    end

    % initialize new clip structure to be exactly like old.
    new_clip_structs = clip_structs;

    % sort clips by HRC.
    offsets = sort_clips_by('hrc', clip_structs,test_structs);

    % loop through each hrc-worth of clips.
    for index = 1:length(offsets),
        % original HRC is a special case.
        % Set luminance gain/offset of the original to required values.
        curr_offsets = offsets{index};
        if strcmp(clip_structs(curr_offsets(1)).hrc,'original'),
            for cnt = 1:length(curr_offsets),
                new_clip_structs(curr_offsets(cnt)).luminance_gain = 1.0;
                new_clip_structs(curr_offsets(cnt)).luminance_offset = 0.0;
                new_clip_structs(curr_offsets(cnt)).cb_gain = 1.0;
                new_clip_structs(curr_offsets(cnt)).cb_offset = 0.0;
                new_clip_structs(curr_offsets(cnt)).cr_gain = 1.0;
                new_clip_structs(curr_offsets(cnt)).cr_offset = 0.0;
            end; 
            continue;
        end

        % loop through each clip of this HRC.  Find its spatial registration.
        for cnt = 1:length(curr_offsets),
            orig_offset = find_original(clip_structs,curr_offsets(cnt));
            [ygain, yoffset, cbgain, cboffset, crgain, croffset, valid] = local_one_clip(test_structs,clip_structs(orig_offset), ...
                clip_structs(curr_offsets(cnt)), verbose, is_quantize); 
            new_clip_structs(curr_offsets(cnt)).luminance_gain = ygain;
            new_clip_structs(curr_offsets(cnt)).luminance_offset = yoffset;
            new_clip_structs(curr_offsets(cnt)).cb_gain = cbgain;
            new_clip_structs(curr_offsets(cnt)).cb_offset = cboffset;
            new_clip_structs(curr_offsets(cnt)).cr_gain = crgain;
            new_clip_structs(curr_offsets(cnt)).cr_offset = croffset;

            if verbose && ~valid,
                fprintf('%s:%s(%s):  Luminance gain & offset could not be computed.\n', ...
                    new_clip_structs(curr_offsets(cnt)).test{1}, ...
                    new_clip_structs(curr_offsets(cnt)).scene{1}, ...
                    new_clip_structs(curr_offsets(cnt)).hrc{1});
            end
            
            if verbose && ~is_hrc_filter,
                fprintf('%s:%s(%s):  ', ...
                    new_clip_structs(curr_offsets(cnt)).test{1}, ...
                    new_clip_structs(curr_offsets(cnt)).scene{1}, ...
                    new_clip_structs(curr_offsets(cnt)).hrc{1});
                fprintf('Y-gain=%f, Y-offset=%f\n', ...
                    new_clip_structs(curr_offsets(cnt)).luminance_gain, ...
                    new_clip_structs(curr_offsets(cnt)).luminance_offset);
                fprintf('Cb-gain=%f, Cb-offset=%f\n', ...
                    new_clip_structs(curr_offsets(cnt)).cb_gain, ...
                    new_clip_structs(curr_offsets(cnt)).cb_offset);
                fprintf('Cr-gain=%f, Cr-offset=%f\n', ...
                    new_clip_structs(curr_offsets(cnt)).cr_gain, ...
                    new_clip_structs(curr_offsets(cnt)).cr_offset);
            end
        end

        % compute overall registration, for greater accuracy.
        if is_hrc_filter,
            unfiltered_clip_structs = new_clip_structs;

            % pick out spatial registration results  for this HRC.
            try
                all = [new_clip_structs(curr_offsets).luminance_gain];
                all = all( logical(isfinite(all))  );
                new_ygain = st_collapse('50%', all');
                all = [new_clip_structs(curr_offsets).luminance_offset];
                all = all( logical(isfinite(all))  );
                new_yoffset = st_collapse('50%', all');
            catch
                new_ygain = 1.0;
                new_yoffset = 0.0;
                status = 1;
                if verbose,
                    fprintf('Luminance gain & offset undefined for %s:*(%s)\n', ...
                        new_clip_structs(curr_offsets(1)).test{1}, ...
                        new_clip_structs(curr_offsets(1)).hrc{1}); 
                end
            end
            try
                all = [new_clip_structs(curr_offsets).cb_gain];
                all = all( logical(isfinite(all))  );
                new_cbgain = st_collapse('50%', all');
                all = [new_clip_structs(curr_offsets).cb_offset];
                all = all( logical(isfinite(all))  );
                new_cboffset = st_collapse('50%', all');
            catch
                new_cbgain = nan;
                new_cboffset = nan;
                if verbose,
                    fprintf('Cb gain & offset undefined for %s:*(%s)\n', ...
                        new_clip_structs(curr_offsets(1)).test{1}, ...
                        new_clip_structs(curr_offsets(1)).hrc{1}); 
                end
            end
            
            try
                all = [new_clip_structs(curr_offsets).cr_gain];
                all = all( logical(isfinite(all))  );
                new_crgain = st_collapse('50%', all');
                all = [new_clip_structs(curr_offsets).cr_offset];
                all = all( logical(isfinite(all))  );
                new_croffset = st_collapse('50%', all');
            catch
                new_crgain = nan;
                new_croffset = nan;
                if verbose,
                    fprintf('Cr gain & offset undefined for %s:*(%s)\n', ...
                        new_clip_structs(curr_offsets(1)).test{1}, ...
                        new_clip_structs(curr_offsets(1)).hrc{1}); 
                end
            end

            if verbose,
                fprintf('\n%s:*(%s):\n\t\t', ...
                    new_clip_structs(curr_offsets(1)).test{1}, ...
                    new_clip_structs(curr_offsets(1)).hrc{1});
                fprintf('Y gain=%f, offset=%f\n', ...
                    new_ygain, new_yoffset);
                fprintf('\t\tCb gain=%f, offset=%f\n', ...
                    new_cbgain, new_cboffset);
                fprintf('\t\tCr gain=%f, offset=%f\n\n', ...
                    new_crgain, new_croffset);
            end

            % copy results back to each of these clips.
            for cnt = 1:length(curr_offsets),
                new_clip_structs(curr_offsets(cnt)).luminance_gain = new_ygain;
                new_clip_structs(curr_offsets(cnt)).luminance_offset = new_yoffset;
                new_clip_structs(curr_offsets(cnt)).cb_gain = new_cbgain;
                new_clip_structs(curr_offsets(cnt)).cb_offset = new_cboffset;
                new_clip_structs(curr_offsets(cnt)).cr_gain = new_crgain;
                new_clip_structs(curr_offsets(cnt)).cr_offset = new_croffset;
            end
        end
    end
    
    0;


    
catch
    new_clip_structs = [];
    status = 1;
    if verbose,
        fprintf('\n%s\n', lasterr);
        fprintf('\tAborting.\n'); 
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_gain, y_offset, cb_gain, cb_offset, cr_gain, cr_offset, status] = ...
    local_one_clip(test_structs,orig_clip, proc_clip, verbose, is_quantize);
% Compute luminance gain & offset (low bandwidth) for one clip.
%
% Algorithm:
%
%   Use one frame every second (approx).  
%   Search over ALL frames simultaneously.
%   Sub-sample into blocks.
%

status = 1;

% if strcmp(orig_clip.video_standard,'progressive'),
%     is_field = 0;
% else
%     is_field = 1;
% end

% clear out any spatial registration & scaling.
proc_clip.luminance_gain = 1.0;
proc_clip.luminance_offset = 0.0;
proc_clip.cb_gain = 1.0;
proc_clip.cb_offset = 0.0;
proc_clip.cr_gain = 1.0;
proc_clip.cr_offset = 0.0;

% check for an alignment.  if none exists, set it to the full clip & re-set
% after spatial registration.

clear_proc_align = 0;
clear_both_align = 0;
if isnan(proc_clip.align_start) && isnan(proc_clip.align_stop) && isnan(orig_clip.align_start) && isnan(orig_clip.align_stop),
    clear_both_align = 1;
    proc_clip.align_start = proc_clip.loc_start;
    orig_clip.align_start = orig_clip.loc_start;
    length = min(proc_clip.loc_stop - proc_clip.loc_start, orig_clip.loc_stop - orig_clip.loc_start); 
    proc_clip.align_stop = proc_clip.loc_start + length;
    orig_clip.align_stop = orig_clip.loc_start + length;
elseif isnan(proc_clip.align_start) && isnan(proc_clip.align_stop) && ~isnan(orig_clip.align_start) && ~isnan(orig_clip.align_stop),
    clear_proc_align = 1;
    proc_clip.align_start = (orig_clip.align_start-orig_clip.loc_start) + proc_clip.loc_start;
    length = orig_clip.align_stop - orig_clip.align_start; 
    proc_clip.align_stop = proc_clip.align_start + length;
    while proc_clip.align_stop > proc_clip.loc_stop,
        proc_clip.align_stop = proc_clip.align_stop - 1;
        proc_clip.align_start = proc_clip.align_start - 1;
    end
    if proc_clip.align_start < 1,
        proc_clip.align_start = 1;
    end
end

[y_gain, y_offset, cb_gain, cb_offset, cr_gain, cr_offset, status] = local_core_algorithm(test_structs,orig_clip, ...
    proc_clip, verbose, is_quantize);


if clear_both_align,
    proc_clip.align_start = NaN;
    orig_clip.align_start = NaN;
    proc_clip.align_stop = NaN;
    orig_clip.align_stop = NaN;
end
if clear_proc_align,
    proc_clip.align_start = NaN;
    proc_clip.align_stop = NaN;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_gain, y_offset, cb_gain, cb_offset, cr_gain, cr_offset, valid] = local_core_algorithm(test_structs,orig_clip, proc_clip, verbose, is_quantize);
% Compute spatial registration for one clip of fields..

if proc_clip.image_size.rows < 216,
    block_size = 10;
elseif proc_clip.image_size.rows < 384,
    block_size = 22;
else
    block_size = 46;
end

tic;

% is_field = 1;
% loop over all frames of this clip. 

total = total_tslices(proc_clip,1/proc_clip.fps);
% Choose video frames.
for loop = 1:floor(total/proc_clip.fps),
    % read in this original time-slice. 
    [y, cb, cr] = read_tslice(test_structs, orig_clip, 1/orig_clip.fps, ...
        (loop-1)*round(orig_clip.fps)+1, 'hsize', block_size, 'vsize', block_size);

    % subsample this original frame
    [orig_y(:,:,loop),orig_y_std(:,:,loop)] = block_statistic(y, block_size, block_size, 'mean','std');
    [orig_cb(:,:,loop)] = block_statistic(cb, block_size, block_size, 'mean');
    [orig_cr(:,:,loop)] = block_statistic(cr, block_size, block_size, 'mean');

    % read this processed time-slice
    [y, cb, cr] = read_tslice(test_structs, proc_clip, 1/proc_clip.fps, ...
        (loop-1)*round(orig_clip.fps)+1, 'hsize', block_size, 'vsize', block_size);
    
    % subsample this processed frame
    [proc_y(:,:,loop)] = block_statistic(y, block_size, block_size, 'mean');
    [proc_cb(:,:,loop)] = block_statistic(cb, block_size, block_size, 'mean');
    [proc_cr(:,:,loop)] = block_statistic(cr, block_size, block_size, 'mean');
    
    0;
end
clear y cb cr;

[r1,c1,t1] = size(orig_y);
[r2,c2,t2] = size(proc_y);

if r1 ~= r2 || c1 ~= c2,
    error('CVR must be identical for all original & processed clips associated with one source');
end

orig_y = reshape(orig_y, r1*c1*t1, 1);
proc_y = reshape(proc_y, r1*c1*t1, 1);

orig_y_std = reshape(orig_y_std, r1*c1*t1, 1);

orig_cb = reshape(orig_cb, r1*c1*t1, 1);
orig_cr = reshape(orig_cr, r1*c1*t1, 1);
proc_cb = reshape(proc_cb, r1*c1*t1, 1);
proc_cr = reshape(proc_cr, r1*c1*t1, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use only the half of blocks where Y plane 
% has low stdev.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,b]=sort(orig_y_std);
b = b(1:floor(length(b) / 2));
orig_y = orig_y(b);
proc_y = proc_y(b);
orig_cb = orig_cb(b);
proc_cb = proc_cb(b);
orig_cr = orig_cr(b);
proc_cr = proc_cr(b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eliminate blocks with clipping
% MUST be done second< so that above indicies
% are identical for y, cb, & cr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = find(orig_y >= 2 & orig_y <= 253 & proc_y >= 2 & proc_y <= 253 );
if length(b) ~= length(orig_y),
    orig_y = orig_y(b);
    proc_y = proc_y(b);
end

b = find(orig_cb >= -126 & orig_cb <= 126 & proc_cb >= -126 & proc_cb <= 126 );
if length(b) ~= length(orig_cb),
    orig_cb = orig_cb(b);
    proc_cb = proc_cb(b);
end

b = find(orig_cr >= -126 & orig_cr <= 126 & proc_cr >= -126 & proc_cr <= 126 );
if length(b) ~= length(orig_cr),
    orig_cr = orig_cr(b);
    proc_cr = proc_cr(b);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quantize, if requested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_quantize,
    start = 0.0;  % first code
    last = 255.0;  % 255 is the maximum observed in the training data 
    high_codes = 1024;  % number of codes for 10-bit quantizer, must make epsilon=1.0 
    code_lgo = start:(last-start)/(high_codes-1):last;

    %  Generate the partitions, halfway between codes
    code_lgo_size = size(code_lgo,2);
    part_lgo = (code_lgo(2:code_lgo_size)+code_lgo(1:code_lgo_size-1))/2;

    %  Quantize orig_y feature like this
    [index_lgo] = quantiz_fast(orig_y',part_lgo);
    %  Look-up the quantized value like this
    orig2 = code_lgo(1+index_lgo);
    orig_y = orig2';

    %  Quantize org_cb feature like this
    [index_lgo] = quantiz_fast(org_cb'+128,part_lgo);
    %  Look-up the quantized value like this
    orig2 = code_lgo(1+index_lgo);
    org_cb = orig2'-128;

    %  Quantize org_cr feature like this
    [index_lgo] = quantiz_fast(org_cr'+128,part_lgo);
    %  Look-up the quantized value like this
    orig2 = code_lgo(1+index_lgo);
    org_cr = orig2'-128;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gain & offset with final pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (max(max(max(orig_y)))-min(min(min(orig_y)))) < 10,
    y_gain = 1.0;
    y_offset = 0.0;
    valid = 0;
elseif (max(max(max(proc_y)))-min(min(min(proc_y)))) <= 0
    y_gain = 1.0;
    y_offset = 0.0;
    valid = 0;
else

    % compute initial gain via linear regression
    y = proc_y;
    x = [ones(length(y),1) orig_y];
    b = x\y;
    r = y - x*b;

    done = 0;
    prev_b = b;
    counter = 0;
    while ~done && counter < 10000,
        counter = counter + 1;

        epsilon = 1.0;
        cost = 1.0 ./ (abs(r) + epsilon); % cost vector, reciprocal of errors
        cost = cost ./ sqrt(sum(cost)); % normalize for unity norm
        cost = (cost.^2);

        xp = x' .* repmat(cost,1,2)';
        b = inv(xp*x)*xp*y;
        r = y - x*b;

        if abs(prev_b(2) - b(2)) < 0.0001,
            done = 1;
        else
            prev_b = b;
        end
    end

    if counter < 10000,
        y_gain = b(2);
        y_offset = b(1);
        valid = 1;
    else
        y_gain = 1.0;
        y_offset = 0.0;
        valid = 0;
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gain & offset with final pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (max(max(max(orig_cb)))-min(min(min(orig_cb)))) < 10,
    cb_gain = nan;
    cb_offset = nan;
elseif (max(max(max(proc_cb)))-min(min(min(proc_cb)))) <= 0
    cb_gain = nan;
    cb_offset = nan;
else

    % compute initial gain via linear regression
    y = proc_cb;
    x = [ones(length(y),1) orig_cb];
    b = x\y;
    r = y - x*b;

    done = 0;
    prev_b = b;
    counter = 0;
    while ~done && counter < 10000,
        counter = counter + 1;
        
        epsilon = 1.0;
        cost = 1.0 ./ (abs(r) + epsilon); % cost vector, reciprocal of errors
        cost = cost ./ sqrt(sum(cost)); % normalize for unity norm
        cost = (cost.^2);

        xp = x' .* repmat(cost,1,2)';
        b = inv(xp*x)*xp*y;
        r = y - x*b;

        if abs(prev_b(2) - b(2)) < 0.0001,
            done = 1;
        else
            prev_b = b;
        end
    end

    if counter < 10000,
        cb_gain = b(2);
        cb_offset = b(1);
    else
        cb_gain = nan;
        cb_offset = nan;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gain & offset with final pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (max(max(max(orig_cr)))-min(min(min(orig_cr)))) < 10,
    cr_gain = nan;
    cr_offset = nan;
elseif (max(max(max(proc_cr)))-min(min(min(proc_cr)))) <= 0
    cr_gain = nan;
    cr_offset = nan;
else

    % compute initial gain via linear regression
    y = proc_cr;
    x = [ones(length(y),1) orig_cr];
    b = x\y;
    r = y - x*b;

    done = 0;
    prev_b = b;
    counter =0;
    while ~done && counter < 10000,
        counter = counter + 1;
        
        epsilon = 1.0;
        cost = 1.0 ./ (abs(r) + epsilon); % cost vector, reciprocal of errors
        cost = cost ./ sqrt(sum(cost)); % normalize for unity norm
        cost = (cost.^2);

        xp = x' .* repmat(cost,1,2)';
        b = inv(xp*x)*xp*y;
        r = y - x*b;

        if abs(prev_b(2) - b(2)) < 0.0001,
            done = 1;
        else
            prev_b = b;
        end
    end

    if counter < 10000,
        cr_gain = b(2);
        cr_offset = b(1);
    else
        cr_gain = nan;
        cr_offset = nan;
    end
end

if verbose,
    fprintf('%20s: y = %8.5f,%8.4f\tCb = %8.5f,%8.4f\tCr = %8.5f,%8.4f\n', ...
        proc_clip.scene{1}, y_gain, y_offset, cb_gain, cb_offset, cr_gain, cr_offset);
    0;
end
