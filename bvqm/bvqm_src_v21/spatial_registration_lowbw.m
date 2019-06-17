function [new_clip_structs, status, unfiltered_clip_structs] = spatial_registration_lowbw(test_structs, clip_structs,varargin)
% SPATIAL_REGISTRATION_LOWBW
%   Automatic computation of spatial registration, scaling, and luminance
%   gain/offset; 2005 low-bandwidth, RR approach.
% SYNTAX
%  [new_clip_structs] = spatial_registration_lowbw(test_structs, clip_structs);
%  [...] = spatial_registration_lowbw(...,'PropertyName',...)
%  [new_clip_structs, status, alternate_clip_structs] = spatial_registration_lowbw(...)
% DESCRIPTION
%  Estimate spatial registration, scaling registration, and luminance 
%  gain/offset for each processed clip.  All clips must be temporally
%  registered first; function 'temporal_registration_sequence' recommended.
%
%  Optional properties:
%  'HRCFilter'     Median filter spatial registration by HRC, to
%                  improve spatial registration accuracy.  Return
%                  un-filtered clips in 'unfiltered_clip_structs'
%  'NoScaling'     Presume no scaling for all clips.  Use with care.
%  'Uncert',<max_shift_horiz>,<max_shift_vert>,<max_scale_horiz>,<max_scale_vert>
%                   Override default search range.  Defaults by image size:
%              NTSC max_shift=20, vert max_scale=60, horiz max_scale=100
%               CIF max_shift=8,  vert max_scale=40, horiz max_scale=60
%              QCIF max_shift=4,  vert max_scale=40, horiz max_scale=60
%  'gain'           Compute luminance Gain/offset using the pixels &
%                   profiles.  This option is not recommended.  See
%                   algorithm luminance_gain_offset_lowbw.
%  'quantize'       Apply quantization to original features
%  'noquantize'     Do not quantize features (default)
%  'verbose',       Print status information to screen
%  'quiet',         Print no information to screen.  Default behavior.
%
%    Return the new, calibrated clip structure.  Also return status,
%    which has the following fields defined:
%     'status.error'              0 if okay; 1 if a fatal error was encountered.
%                                 If status.error is not 0, the other structure fields 
%                                 will be unavailable.  See also function "lasterr".
%     'status.scale'            Number of video clips that should be re-run
%                                 with a larger scaling uncertainty.
%     'status.shift'              Number of video clips that should be re-run
%                                 with a larger shift uncertainty.
%     'status.luminance'          Number of video clips for which luminance
%                                 calculation failed.
%
%   NOTE: this function can take align_start values with a 0.5 fraction, as
%   produced by temporal_registration_lowbw(..., 'field_select');  In this
%   case, reframe timing will be assumed.  The output clip-set will
%   contain only integer align_start values.
% 

status.error = 0;
status.scale = 0;
status.shift = 0;
status.luminance = 0;
do_quantize = 0;
try

    % set up constants.
    max_shift_horiz = 20; % maximum search 20 pixels in any direction.
    max_shift_vert = 20;
    max_scale_horiz = 100; % maximum scaling search, 10%
    max_scale_vert = 60;  % maximum scaling search, 6%
    default_max = 1;
    do_scaling = 1;
    
    % figure out what optional properties have been requested.
    is_hrc_filter = 0;
    do_gain = 0;
    verbose = 0;
    cnt = 1;
    while cnt <= nargin-2,
        if strcmp(lower(varargin{cnt}),'hrcfilter'),
            is_hrc_filter = 1;
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'noscaling'),
            do_scaling = 0;
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'gain'),
            do_gain = 1;
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'quantize'),
            do_quantize = 1;
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'noquantize'),
            do_quantize = 0;
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'uncert') && cnt <= nargin - 6,
            default_max = 0;
            max_shift_horiz = varargin{cnt+1};
            max_shift_vert = varargin{cnt+2};
            max_scale_horiz = varargin{cnt+3};
            max_scale_vert = varargin{cnt+4};
            cnt = cnt + 5;
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

    % seed random number generator.
    rand('state', sum(100*clock));

    % initialize new clip structure to be exactly like old.
    new_clip_structs = clip_structs;

    % sort clips by HRC.
    offsets = sort_clips_by('hrc', clip_structs,test_structs);

    % loop through each hrc-worth of clips.
    for index = 1:length(offsets),
        % original HRC is a special case.
        % Set spatial registration of the original to (0,0), by definition.
        curr_offsets = offsets{index};
        if strcmp(clip_structs(curr_offsets(1)).hrc,'original'),
            for cnt = 1:length(curr_offsets),
                new_clip_structs(curr_offsets(cnt)).spatial.horizontal = 0;
                new_clip_structs(curr_offsets(cnt)).spatial.vertical = 0;
                new_clip_structs(curr_offsets(cnt)).scale.horizontal = 1000;
                new_clip_structs(curr_offsets(cnt)).scale.vertical = 1000;
            end; 
            continue;
        end

        % loop through each clip of this HRC.  Find its spatial registration.
        for cnt = 1:length(curr_offsets),
            orig_offset = find_original(clip_structs,curr_offsets(cnt));
            
            % Set maximum shift & scale search limits, if defaults chosen
            % these depend upon image size.
            if default_max,
                if clip_structs(orig_offset).image_size.rows <= 216,
                    max_shift_horiz = 4; 
                    max_shift_vert = 4;
                    max_scale_horiz = 60; 
                    max_scale_vert = 40;  
                elseif clip_structs(orig_offset).image_size.rows <= 384,
                    max_shift_horiz = 8; 
                    max_shift_vert = 8;
                    max_scale_horiz = 60;
                    max_scale_vert = 40; 
                else
                    max_shift_horiz = 20; % maximum search 20 pixels in any direction.
                    max_shift_vert = 20;
                    max_scale_horiz = 100; % maximum scaling search, 10%
                    max_scale_vert = 60;  % maximum scaling search, 6%
                end
            end
            if ~do_scaling,
                max_scale_horiz = 0;
                max_scale_vert = 0;
            end
            
            [hold1, hold2, hold3, hold4, status] = ...
                sas_one_clip(test_structs,clip_structs(orig_offset), ...
                clip_structs(curr_offsets(cnt)), max_shift_horiz, max_shift_vert, ...
                max_scale_horiz,max_scale_vert,...
                verbose, status, do_gain, do_quantize); 
            new_clip_structs(curr_offsets(cnt)).spatial = hold1;
            new_clip_structs(curr_offsets(cnt)).scale = hold2;
            if do_gain,
                new_clip_structs(curr_offsets(cnt)).luminance_gain = hold3;
                new_clip_structs(curr_offsets(cnt)).luminance_offset = hold4;
            end
            % discard fractional alignment indicating reframing
            new_clip_structs(curr_offsets(cnt)).align_start =  ...
                floor(new_clip_structs(curr_offsets(cnt)).align_start);
            new_clip_structs(curr_offsets(cnt)).align_stop =  ...
                ceil(new_clip_structs(curr_offsets(cnt)).align_start);
        end
    end
    
    % copy results into "unfiltered" variable
    unfiltered_clip_structs = new_clip_structs;
    
    % loop through each hrc-worth of clips.
    for index = 1:length(offsets),
        % original doesn't need any more work.
        curr_offsets = offsets{index};
        if strcmp(clip_structs(curr_offsets(1)).hrc,'original'),
            continue;
        end

        % compute overall registration, for greater accuracy.
        if is_hrc_filter,

            % pick out spatial registration results  for this HRC.
            all = [new_clip_structs(curr_offsets).spatial];
            new_spatial.vertical = st_collapse('50%', [all.vertical]');
            new_spatial.horizontal = st_collapse('50%', [all.horizontal]');
            all = [new_clip_structs(curr_offsets).scale];
            new_scale.vertical = st_collapse('50%', [all.vertical]');
            new_scale.horizontal = st_collapse('50%', [all.horizontal]');
            if do_gain,
                all = [new_clip_structs(curr_offsets).luminance_gain];
                new_gain = st_collapse('50%', all');
                all = [new_clip_structs(curr_offsets).luminance_offset];
                new_offset = st_collapse('50%', all');
            end

            if verbose,
                fprintf('==> %s:*(%s):  ', ...
                    new_clip_structs(curr_offsets(1)).test{1}, ...
                    new_clip_structs(curr_offsets(1)).hrc{1});
                fprintf('Scale: h=%d, v=%d,  Shift: h=%d, v=%d', ...
                    new_scale.horizontal, new_scale.vertical, ...
                    new_spatial.horizontal, new_spatial.vertical);
                if do_gain,
                    fprintf(' gain=%f, offset=%f', new_gain, new_offset);
                end
                fprintf('\n\n');
            end

            % copy results back to each of these clips.
            for cnt = 1:length(curr_offsets),
                new_clip_structs(curr_offsets(cnt)).spatial.horizontal = new_spatial.horizontal;
                new_clip_structs(curr_offsets(cnt)).spatial.vertical = new_spatial.vertical;
                new_clip_structs(curr_offsets(cnt)).scale.horizontal = new_scale.horizontal;
                new_clip_structs(curr_offsets(cnt)).scale.vertical = new_scale.vertical;
                if do_gain,
                    new_clip_structs(curr_offsets(cnt)).luminance_gain = new_gain;
                    new_clip_structs(curr_offsets(cnt)).luminance_offset = new_offset;
                end
            end
        end
        
        % check for shift or scale values at the maximum.
        for cnt = 1:length(curr_offsets),
            if new_clip_structs(curr_offsets(cnt)).spatial.horizontal == max_shift_horiz || ...
                    new_clip_structs(curr_offsets(cnt)).spatial.vertical == max_shift_vert,
                status.shift = status.shift + 1;
            end
            if abs(new_clip_structs(curr_offsets(cnt)).scale.horizontal - 1000) == max_scale_horiz || ...
                    abs(new_clip_structs(curr_offsets(cnt)).scale.vertical - 1000) == max_scale_vert,
                status.scale = status.scale + 1;
            end
        end
    end
    
    % fix endpoint of these clip structures
    new_clip_structs = fix_temporal(test_structs, new_clip_structs,'EndPoint');
    unfiltered_clip_structs = fix_temporal(test_structs, unfiltered_clip_structs,'EndPoint');


    
catch
    new_clip_structs = [];
    status.error = 1;
    if verbose,
        fprintf('\n%s\n', lasterr);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_spatial, new_scale, y_gain, y_offset, status] = ...
    sas_one_clip(test_structs,orig_clip, proc_clip, max_shift_horiz, max_shift_vert, ...
    max_scale_horiz, max_scale_vert, verbose, status, do_gain, do_quantize);
% Compute spatial registration for one clip.
%
% Algorithm:
%
%   Use one frame every second (approx).  
%       horizontal and vertical profiles AND
%       80% as many randomly subsampled pixels as there are profile pixels
%           - randomly distributed over all frames
%   Search over ALL frames simultaneously.
%   Search original +- 0 second (yes! ZERO); 
%   shift processed +- 20 pixels/lines.
%   scale processed by +/- 6% (nearest neighbor)
%
%   When have final scale & shift, compute luminance gain & offset with
%   those pixels & profiles, too.


if strcmp(orig_clip.video_standard,'progressive'),
    is_field = 0;
else
    is_field = 1;
end

% clear out any spatial registration & scaling.
orig_clip = clear_spatial_registration(orig_clip);
proc_clip = clear_spatial_registration(proc_clip);


% error checks & corrections
if mod(max_shift_horiz,2),
    max_shift_horiz = max_shift_horiz + 1;
end
if mod(max_shift_vert,2),
    max_shift_vert = max_shift_vert + 1;
end

% compute PVR and OROI given the above.
max_pixels_horiz = max_shift_horiz + (max_scale_horiz / 1000) * ...
    proc_clip.image_size.cols;
max_pixels_horiz = ceil(max_pixels_horiz);
max_pixels_horiz = max_pixels_horiz + mod(max_pixels_horiz,2);

max_pixels_vert = max_shift_vert + (max_scale_vert / 1000) * ...
    proc_clip.image_size.rows;
max_pixels_vert = ceil(max_pixels_vert);
max_pixels_vert = max_pixels_vert + mod(max_pixels_vert,2);

[pvr, oroi] = find_pvr_oroi_guess(proc_clip, max_pixels_horiz, max_pixels_vert, is_field);


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

[new_spatial, new_scale, y_gain, y_offset, status] = sas_core_algorithm(test_structs,orig_clip, ...
    proc_clip, max_shift_horiz, max_shift_vert, max_scale_horiz, max_scale_vert, max_pixels_horiz, ...
    max_pixels_vert, pvr, oroi, verbose, status, do_gain, do_quantize);

% fprintf('\t%s:%s(%s) new shift (v=%d, h=%d)\n', ...
%     proc_clip.test{1}, proc_clip.scene{1}, proc_clip.hrc{1}, ...
%     new_spatial.vertical, new_spatial.horizontal);

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
function [new_spatial, new_scale, y_gain, y_offset, status] = sas_core_algorithm(test_structs,orig_clip, ...
    proc_clip, max_shift_horiz, max_shift_vert, max_scale_horiz, max_scale_vert,...
    max_pixels_horiz, max_pixels_vert, pvr, oroi, verbose, status, do_gain, do_quantize);
% Compute spatial registration for one clip of fields..

tic;

new_spatial.horizontal = 0;
new_spatial.vertical = 0;
new_scale.horizontal = 1000;
new_scale.vertical = 1000;


is_field = 1;
% loop over all frames of this clip. 
if verbose,
    fprintf('%s:%s(%s)\n', proc_clip.test{1}, proc_clip.scene{1}, proc_clip.hrc{1});
end

total = total_tslices(proc_clip,1/proc_clip.fps);
if mod(proc_clip.align_start,1) ~= 0,
    total = total - 1;
end
% Choose video frames.
for loop = 1:floor(total/round(proc_clip.fps)),
    % read in this original time-slice.
    orig(:,:,loop) = read_tslice(test_structs, orig_clip, 1/orig_clip.fps, ...
        (loop-1)*round(orig_clip.fps)+1, 'sroi', oroi.top, oroi.left, oroi.bottom, oroi.right);

    [row,col,time] = size(orig);
    
    % read this processed time-slice
    if mod(proc_clip.align_start,1) == 0,
        proc(:,:,loop) = read_tslice(test_structs, proc_clip, 1/proc_clip.fps, ...
            (loop-1)*round(orig_clip.fps)+1, 'sroi', pvr.top, pvr.left, pvr.bottom, pvr.right);
    else
        % change time order as per reframing, but don't swap fields
        proc(:,:,loop) = read_tslice(test_structs, proc_clip, 1/proc_clip.fps, ...
            (loop-1)*round(orig_clip.fps)+1, 'sroi', pvr.top, pvr.left, pvr.bottom, pvr.right, ...
            'align_start',floor(proc_clip.align_start));
        hold = read_tslice(test_structs, proc_clip, 1/proc_clip.fps, ...
            (loop-1)*round(orig_clip.fps)+2, 'sroi', pvr.top, pvr.left, pvr.bottom, pvr.right, ...
            'align_start',floor(proc_clip.align_start));
        [rowp2,colp2,time]=size(proc);
        if strcmp(proc_clip.video_standard,'interlace_lower_field_first'),
            proc(2:2:rowp2,:,loop) = hold(2:2:rowp2,:);
        else
            proc(1:2:rowp2,:,loop) = hold(1:2:rowp2,:);
        end
    end
    
    clear orig_y proc_y row col time cnt;
end

% fprintf('Frames selected and loaded.  '); toc;

[orig_horiz_profile, orig_vert_profile] = sas_profile_images(orig);
[proc_horiz_profile, proc_vert_profile] = sas_profile_images(proc);

% quantize
if do_quantize,
    %  profile feature
    start = 0.0;  % first code
    last = 255.0;  % 252 is the maximum observed in the training data
    high_codes = 2^16;  % number of codes for 16-bit quantizer
    code_profile = start:(last-start)/(high_codes-1):last;
    %  Generate the partitions, halfway between codes
    code_profile_size = size(code_profile,2);
    part_profile = (code_profile(2:code_profile_size)+code_profile(1:code_profile_size-1))/2;


    %  Quantize original features
    [index_temp] = quantiz_fast((orig_horiz_profile),part_profile);
    orig_horiz_profile = code_profile(1+index_temp);
    orig_horiz_profile = reshape(orig_horiz_profile,1,length(orig_horiz_profile));

    [index_temp] = quantiz_fast((orig_vert_profile),part_profile);
    orig_vert_profile = code_profile(1+index_temp);
    orig_vert_profile = reshape(orig_vert_profile,1,length(orig_vert_profile));

end

% now, we have orig with original frames, and proc with processed frames.
% reshape rows & columns into one dimension.
[rows,cols,time] = size(orig);
orig_y = reshape(orig, rows*cols*time,1);

[rowp,colp,timep] = size(proc);
proc_y = reshape(proc, rowp*colp*timep,1);

clear orig proc;

%list of coordinates for profiles
list_horiz_profile = (1:cols) + max_pixels_horiz;
list_vert_profile = (1:rows) + max_pixels_vert;

% Choose the % of pixels to be used.
[seed_state, list_row, list_col, list_time, list_o] = ...
    sas_choose_pixels (orig_y, rows, cols, time, max_pixels_horiz, max_pixels_vert);

% pick our original pixels
orig_pixels = orig_y(list_o);

length_of_list = length(list_o);

    
% random search.
loop = 1;
best_scale_horiz = 1000;
best_scale_vert = 1000;
best_shift_horiz = 0;
best_shift_vert = 0;
best_value = inf;
best_loop = -1;

values = zeros(2*max_scale_horiz+1,2*max_scale_vert+1,2*max_shift_horiz+1,2*max_shift_vert+1);
values(:,:,:,:) = NaN;

random_tries = 15000; % hope it needn't be that large!
loop = 0;
cnt_used = 0;
while loop <= random_tries,

    % randomly choose stretch, shift, and time.

    % for the first 10% of tries, do a flat random search over all possibilities.
    if loop < random_tries / 10,
        scale_horiz = round( -max_scale_horiz + (2 * max_scale_horiz + 1) * rand );
        scale_vert =  round( -max_scale_vert + (2 * max_scale_vert + 1) * rand );
        shift_horiz = round( -max_shift_horiz + (2 * max_shift_horiz + 1) * rand );
        shift_vert =  round( -max_shift_vert + (2 * max_shift_vert + 1) * rand );
    % weight more near best stretch/shift/time found so far.  
    else
        scale_horiz = best_scale_horiz + round( 2 * randn );
        scale_vert = best_scale_vert + round( 2 * randn );
        shift_horiz = best_shift_horiz + round( 2 * randn );
        shift_vert = best_shift_vert + round( 2 * randn );
    end

    if max_scale_horiz == 0,
        scale_horiz = 0;
    end
    if max_scale_vert == 0,
        scale_vert = 0;
    end

    % If this point is out of the legal range, choose again.
    if abs(scale_horiz) > max_scale_horiz | abs(scale_vert) > max_scale_vert | ...
            abs(shift_horiz) > max_shift_horiz | abs(shift_vert) > max_shift_vert,
        continue;
    end
    % check whether this stretch/shift/time has been computed already. 
    if ~isnan( values(scale_horiz + max_scale_horiz + 1, scale_vert + max_scale_vert + 1, ...
        shift_horiz + max_shift_horiz + 1, shift_vert + max_shift_vert + 1) ),
        loop = loop + 1;
        continue;
    end
    cnt_used = cnt_used + 1;

    % convert from chosen original coordinates (in smaller image) into
    % processed coordinates.  
    % First, scale.  The "-0.4" is because the resampling routine we are
    % training for his this factor.
    curr_list_row = list_row * 1000 / (scale_vert+1000) + (rowp/2 - (1000/(scale_vert+1000)) * (rowp/2));
    curr_list_col = list_col * 1000 / (scale_horiz+1000) + (colp/2 - (1000/(scale_horiz+1000)) * (colp/2));
    
    % Second, shift.  Round to nearest pixel (use floor of +0.5 for speed)
    curr_list_row = floor(curr_list_row + shift_vert + 0.5);
    curr_list_col = floor(curr_list_col + shift_horiz + 0.5);
    list_p = (list_time-1)*rowp*colp + (curr_list_col-1)*rowp + curr_list_row;

    % compute the difference value of the random pixels.
    pixel_list = orig_pixels - proc_y(list_p);
    
    % scale the profiles.
    curr_list_row = list_vert_profile * 1000 / (scale_vert+1000) + (rowp/2 - (1000/(scale_vert+1000)) * (rowp/2));
    curr_list_col = list_horiz_profile * 1000 / (scale_horiz+1000) + (colp/2 - (1000/(scale_horiz+1000)) * (colp/2));

    % Second, shift the profiles.  Round to nearest pixel (use floor of +0.5 for speed)
    curr_list_row = floor(curr_list_row + shift_vert + 0.5);
    curr_list_col = floor(curr_list_col + shift_horiz + 0.5);

    % compute the difference value of the profiles.  Append that to random
    % pixels' differneces.
    vert_diff = orig_vert_profile - proc_vert_profile(curr_list_row,:);
    horiz_diff = orig_horiz_profile - proc_horiz_profile(curr_list_col,:);
    pixel_list = [ pixel_list; reshape(vert_diff,rows*time,1); reshape(horiz_diff,cols*time,1)];

    % Compute and record standard deviation of difference.
    curr_value = std(pixel_list);
    values(scale_horiz + max_scale_horiz + 1, scale_vert + max_scale_vert + 1, ...
        shift_horiz + max_shift_horiz + 1, shift_vert + max_shift_vert + 1) = curr_value;

    % keep track of the best one!
    % if find a tie, change to it if the scaling factor is closer to no
    % scaling & no shifting
    want = 0;
    if curr_value < best_value,
        want = 1;
    elseif curr_value == best_value,
        if abs(scale_horiz) <= abs(best_scale_horiz) & abs(scale_vert) <= abs(best_scale_vert) & ...
                abs(shift_horiz) <= abs(best_shift_horiz) & abs(shift_vert) <= abs(best_shift_vert),
            want = 1;
        end
    end
    if want,
        best_scale_horiz = scale_horiz;
        best_scale_vert = scale_vert;
        best_shift_horiz = shift_horiz;
        best_shift_vert = shift_vert;
        best_value = curr_value;
        best_loop = loop;
    end

%         fprintf('  (%d) -- Scale: H%d V%d, Shift: H%d V%d\n', loop, ...
%             best_scale_vert, best_scale_horiz, best_shift_vert, best_shift_horiz);

    loop = loop + 1;
end

scale_horiz = best_scale_horiz;
scale_vert = best_scale_vert;
shift_horiz = best_shift_horiz;
shift_vert = best_shift_vert;

if do_gain,

    % First scale coordinates for random pixels.  
    % First, scale. 
    curr_list_row = list_row * 1000 / (scale_vert+1000) + (rowp/2 - (1000/(scale_vert+1000)) * (rowp/2));
    curr_list_col = list_col * 1000 / (scale_horiz+1000) + (colp/2 - (1000/(scale_horiz+1000)) * (colp/2));
    
    % Second, shift.  Round to nearest pixel (use floor of +0.5 for speed)
    curr_list_row = floor(curr_list_row + shift_vert + 0.5);
    curr_list_col = floor(curr_list_col + shift_horiz + 0.5);
    list_p = (list_time-1)*rowp*colp + (curr_list_col-1)*rowp + curr_list_row;

    % Third, scale the profiles.
    curr_list_row = list_vert_profile * 1000 / (scale_vert+1000) + (rowp/2 - (1000/(scale_vert+1000)) * (rowp/2) );
    curr_list_col = list_horiz_profile * 1000 / (scale_horiz+1000) + (colp/2 -(1000/(scale_horiz+1000)) * (colp/2));

    % Forth, shift the profiles.  Round to nearest pixel (use floor of +0.5 for speed)
    curr_list_row = floor(curr_list_row + shift_vert + 0.5);
    curr_list_col = floor(curr_list_col + shift_horiz + 0.5);

    % pick off scaled & shifted random pixels and profiles 
    orig_pixel_list = [orig_pixels ;reshape(orig_vert_profile,rows*time,1) ;reshape(orig_horiz_profile,cols*time,1) ];
    proc_pixel_list = [proc_y(list_p) ;reshape(proc_vert_profile(curr_list_row,:), rows*time,1); reshape(proc_horiz_profile(curr_list_col,:), cols*time,1) ];
    
    clear orig_y proc_y;
    pack

    % compute the luminance gain & offset for those pixels.
    [valid, y_gain, y_offset] = sas_compute_gain_offset(orig_pixel_list, proc_pixel_list);
    if ~valid,
        if verbose,
            fprintf('Gain/offset failed\n');
        end
        status.luminance = status.luminance + 1;
        y_gain = 1.0;
        y_offset = 0.0;
    end
else
    y_gain = 1.0;
    y_offset = 0.0;
end

% send back all results.
new_spatial.horizontal = best_shift_horiz;
new_spatial.vertical = best_shift_vert;
new_scale.horizontal = best_scale_horiz+1000;
new_scale.vertical = best_scale_vert+1000;

if verbose,
    if do_gain,
        fprintf('Scale: h=%d, v=%d,  Shift: h=%d, v=%d, gain=%f, offset=%f\n', ...
            new_scale.horizontal, new_scale.vertical, ...
            new_spatial.horizontal, new_spatial.vertical, ...
            y_gain, y_offset);
    else
        fprintf('Scale: h=%d, v=%d,  Shift: h=%d, v=%d\n', ...
            new_scale.horizontal, new_scale.vertical, ...
            new_spatial.horizontal, new_spatial.vertical);
    end
end

 
% temp = find(values < best_value + 1);
%     fprintf('%d within 1 of best std.  seed=%d\n', length(temp), seed_state);
;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clip] = clear_spatial_registration(old_clip);
% Erase all spatial registration information in a clip, except for temporal
% alignment (which is presumed to be a reasonable guess) and stretch
% registration.

clip = old_clip;

clip.scale.vertical = 1000;
clip.scale.horizontal = 1000;
clip.spatial.vertical = 0;
clip.spatial.horizontal = 0;
clip.luminance_offset = 0;
clip.luminance_gain = 1.0;
clip.cvr.top = 1;
clip.cvr.left = 1;
clip.cvr.bottom = clip.image_size.rows;
clip.cvr.right = clip.image_size.cols;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pvr] = find_pvr_guess(clip);
% return the best guess for pvr based only on image size.
% i.e., discard overscan.
% 'fchoice' is 1 for 'field' or 0 for 'frame' depending on which is desired.

if (clip.image_size.rows == 486 | clip.image_size.rows == 480 ) ...
        & clip.image_size.cols == 720,
    pvr.top = 19;
    pvr.bottom = 486 - 18;
    pvr.left = 23;
    pvr.right = 720 - 22;
elseif clip.image_size.rows == 576 & clip.image_size.cols == 720,
    pvr.top = 15;
    pvr.bottom = 576 - 14;
    pvr.left = 23;
    pvr.right = 720 - 22;
else
    pvr.top = 1;
    pvr.bottom = clip.image_size.rows;
    pvr.left = 1;
    pvr.right = clip.image_size.cols;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pvr, oroi] = find_pvr_oroi_guess(clip, max_pixels_horiz, max_pixels_vert, is_field);
% Find best guess at PVR and OROI, given image size (within clip) and
% size of maximum search ('max_pixels_horiz'), in pixels.
% 'is_field' is true for 'field', 1 for 'frame', depending on which is desired.

pvr = find_pvr_guess(clip);
if is_field,
    oroi.top = pvr.top + max_pixels_vert;
    oroi.bottom = pvr.bottom - max_pixels_vert;
else
    oroi.top = pvr.top + max_pixels_vert;
    oroi.bottom = pvr.bottom - max_pixels_vert;
end
oroi.left = pvr.left + max_pixels_horiz;
oroi.right = pvr.right - max_pixels_horiz;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [seed_state, list_row, list_col, list_time, list_o] = ...
    sas_choose_pixels (orig_y, rows, cols, time, max_pixels_horiz, max_pixels_vert);


rand('seed',sum(100*clock));
seed_state = round(rand * 255);

% Choose 80% as many random points as profile points.
need = ceil(0.80 * (rows+cols)*time );

rand('state', seed_state);
list_row = round(rand(1,need) * (rows) + 0.5);
list_col = round(rand(1,need) * (cols) + 0.5);
list_time = round(rand(1,need) * (time) + 0.5);

% limit to range available.  'rand' is unlikely to exceed that range, but
% it is possible.
list_row = max( min(list_row,rows), 1);
list_col = max( min(list_col,cols), 1);
list_time = max( min(list_time,time), 1);

list_o = (list_time-1)*rows*cols + (list_col-1)*rows + list_row;

% change coordinates to be for processed, where there are more rows &
% columns.
list_row = list_row + max_pixels_vert;
list_col = list_col + max_pixels_horiz;



%%%%%%%%%%%
function [horiz_profile, vert_profile] = sas_profile_images(y);
% Compute the horizontal profile for each frame (averaging each column).
% Put all such together into one giant image.
% Return that image in profile_image.

[rows,cols,time] = size(y);

horiz_profile = zeros(cols, time);
vert_profile = zeros(rows, time);
for cnt = 1:time,
    horiz_profile(:,cnt) = mean(y(:,:,cnt))';
    vert_profile(:,cnt) = mean(y(:,:,cnt), 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gain & offset with final pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [valid, y_gain, y_offset] = sas_compute_gain_offset(orig, proc);
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

