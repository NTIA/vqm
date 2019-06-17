function [new_clip_structs, status, unfiltered_clip_structs] = luminance_gain_offset_lowbw (test_structs, clip_structs,varargin);
% LUMINANCE_GAIN_OFFSET_LOWBW
%   Automatic computation of luminance gain/offset; 2006 low-bandwidth, RR approach.
% SYNTAX
%  [new_clip_structs, status] = luminance_gain_offset_lowbw(test_structs, clip_structs);
%  [...] = luminance_gain_offset_lowbw(...,'PropertyName',...)
%  [new_clip_structs, status, unfiltered_clip_structs] =
%                       luminance_gain_offset_lowbw(...,'HRCFilter',...);
% DESCRIPTION
%  Estimate luminance gain/offset for each processed clip.  All clips must 
%  be temporally and spatially registered first; functions
%  'temporal_registration_lowbw' and 'spatial_registration_lowbw'
%  recommended.  
%
%  Optional properties:
%  'HRCFilter'     Median filter luminance gain & offset by HRC, to
%                  improve accuracy.  Return un-filtered in unfiltered_clip_structs.
%  'Quantize'      Quantize original features, as when running in-service.
%  'verbose',       Print status information to screen
%  'quiet',         Print no information to screen.  Default behavior.
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
        if strcmp(lower(varargin{cnt}),'hrcfilter'),
            is_hrc_filter = 1;
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'quantize'),
            is_quantize = 1;
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'verbose'),
            verbose = 1;
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'quiet'),
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
            end; 
            continue;
        end

        % loop through each clip of this HRC.  Find its spatial registration.
        for cnt = 1:length(curr_offsets),
            orig_offset = find_original(clip_structs,curr_offsets(cnt));
            [gain, offset] = lgol_one_clip(test_structs,clip_structs(orig_offset), ...
                clip_structs(curr_offsets(cnt)), verbose, is_quantize); 
            new_clip_structs(curr_offsets(cnt)).luminance_gain = gain;
            new_clip_structs(curr_offsets(cnt)).luminance_offset = offset;

            if verbose & ~is_hrc_filter,
                fprintf('%s:%s(%s):  ', ...
                    new_clip_structs(curr_offsets(cnt)).test{1}, ...
                    new_clip_structs(curr_offsets(cnt)).scene{1}, ...
                    new_clip_structs(curr_offsets(cnt)).hrc{1});
                fprintf('gain=%f, offset=%f\n\n', ...
                    new_clip_structs(curr_offsets(cnt)).luminance_gain, ...
                    new_clip_structs(curr_offsets(cnt)).luminance_offset);
            end
        end

        % compute overall registration, for greater accuracy.
        if is_hrc_filter,
            unfiltered_clip_structs = new_clip_structs;

            % pick out spatial registration results  for this HRC.
            all = [new_clip_structs(curr_offsets).luminance_gain];
            new_gain = st_collapse('50%', all');
            all = [new_clip_structs(curr_offsets).luminance_offset];
            new_offset = st_collapse('50%', all');

            if verbose,
                fprintf('%s:*(%s):  ', ...
                    new_clip_structs(curr_offsets(1)).test{1}, ...
                    new_clip_structs(curr_offsets(1)).hrc{1});
                fprintf('gain=%f, offset=%f\n\n', ...
                    new_gain, new_offset);
            end

            % copy results back to each of these clips.
            for cnt = 1:length(curr_offsets),
                new_clip_structs(curr_offsets(cnt)).luminance_gain = new_gain;
                new_clip_structs(curr_offsets(cnt)).luminance_offset = new_offset;
            end
        end
    end


    
catch
    new_clip_structs = [];
    status = 1;
    if verbose,
        fprintf('\n%s\n', lasterr);
        fprintf('\tAborting.\n'); 
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_gain, y_offset, status] = ...
    lgol_one_clip(test_structs,orig_clip, proc_clip, verbose, is_quantize);
% Compute luminance gain & offset (low bandwidth) for one clip.
%
% Algorithm:
%
%   Use one frame every second (approx).  
%   Search over ALL frames simultaneously.
%   Sub-sample into blocks.
%


if strcmp(orig_clip.video_standard,'progressive'),
    is_field = 0;
else
    is_field = 1;
end

% clear out any spatial registration & scaling.
proc_clip.luminance_gain = 1.0;
proc_clip.luminance_offset = 0.0;

% check for an alignment.  if none exists, set it to the full clip & re-set
% after spatial registration.

clear_proc_align = 0;
clear_both_align = 0;
if isnan(proc_clip.align_start) & isnan(proc_clip.align_stop) & isnan(orig_clip.align_start) & isnan(orig_clip.align_stop),,
    clear_both_align = 1;
    proc_clip.align_start = proc_clip.loc_start;
    orig_clip.align_start = orig_clip.loc_start;
    length = min(proc_clip.loc_stop - proc_clip.loc_start, orig_clip.loc_stop - orig_clip.loc_start); 
    proc_clip.align_stop = proc_clip.loc_start + length;
    orig_clip.align_stop = orig_clip.loc_start + length;
elseif isnan(proc_clip.align_start) & isnan(proc_clip.align_stop) & ~isnan(orig_clip.align_start) & ~isnan(orig_clip.align_stop),,
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

[y_gain, y_offset] = lgol_core_algorithm(test_structs,orig_clip, ...
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
function [y_gain, y_offset] = lgol_core_algorithm(test_structs,orig_clip, proc_clip, verbose, is_quantize);
% Compute spatial registration for one clip of fields..

if proc_clip.image_size.rows < 216,
    block_size = 20;
elseif proc_clip.image_size.rows < 384,
    block_size = 30;
else
    block_size = 46;
end

tic;

is_field = 1;
% loop over all frames of this clip. 

total = total_tslices(proc_clip,1/proc_clip.fps);
% Choose video frames.
for loop = 1:floor(total/proc_clip.fps),
    % read in this original time-slice. 
    orig_y = read_tslice(test_structs, orig_clip, 1/orig_clip.fps, ...
        (loop-1)*round(orig_clip.fps)+1, 'hsize', block_size, 'vsize', block_size);

    % subsample this original frame
    orig(:,:,loop) = block_statistic(orig_y, block_size, block_size, 'mean');

    % read this processed time-slice
    proc_y = read_tslice(test_structs, proc_clip, 1/proc_clip.fps, ...
        (loop-1)*round(orig_clip.fps)+1, 'hsize', block_size, 'vsize', block_size);
    
    % subsample this processed frame
    proc(:,:,loop) = block_statistic(proc_y, block_size, block_size, 'mean');
end
clear orig_y proc_y;

[r1,c1,t1] = size(orig);
[r2,c2,t2] = size(proc);

if r1 ~= r2 | c1 ~= c2,
    error('CVR must be identical for all original & processed clips associated with one source');
end

orig = reshape(orig, r1*c1*t1, 1);
proc = reshape(proc, r1*c1*t1, 1);


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

    %  Quantize orig feature like this
    [index_lgo] = quantiz_fast(orig',part_lgo);

    %  Look-up the quantized value like this
    orig2 = code_lgo(1+index_lgo);
    orig = orig2';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gain & offset with final pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute initial gain via linear regression
y = proc;
x = [ones(length(y),1) orig];
b = x\y;
r = y - x*b;

done = 0;
prev_b = b;
while ~done,
	epsilon = 0.1;
	cost = 1.0 ./ (abs(r) + epsilon); % cost vector, reciprocal of errors
	cost = cost ./ sqrt(sum(cost)); % normalize for unity norm
	cost = (cost.^2);
%     temp = eye(length(y),length(y));
%     for i=1:length(y),
%         temp(i,i) = temp(i,i) * cost(i);
%     end
%     cost = temp;
% 	
% 	b = inv(x'*cost*x)*x'*cost*y;
    xp = x' .* repmat(cost,1,2)';
	b = inv(xp*x)*xp*y;
    r = y - x*b;
    
    if abs(prev_b(2) - b(2)) < 0.0001,
        done = 1;
    else
        prev_b = b;
    end
end

valid = 1;
y_gain = b(2);
y_offset = b(1);

