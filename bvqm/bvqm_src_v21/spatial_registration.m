function [new_clip_structs, status, unfiltered_clip_structs] = spatial_registration(test_structs, clip_structs,varargin)
% SPATIAL_REGISTRATION
%   Automatic computation of spatial registration.
%   2004 ITU standard, General Model FRTV Calibration
% SYNTAX
%  [new_clip_structs, unfiltered_clip_structs] = spatial_registration(test_structs, clip_structs);
%  [...] = spatial_registration(...,'PropertyName',...)
%  [new_clip_structs] = spatial_registration(...,'NoHRCFilter',...)
% DESCRIPTION
%  Perform spatial registration on each processed clip in 'clip_structs'.
%
%  Optional properties:
%  'NoHRCFilter'     Don't median filter spatial registration by HRC, to
%                           improve spatial registration accuracy.
%                    Don't return 'unfiltered_clip_structs', as it would be
%                    identical to new_clip_structs;
%  'verbose'         Print progress and results information to the screen.
%  'quiet'           Print no information to screen.  Default behavior.
%
% status returns the number of clips for which spatial registration could
% not be found.  These will be left 'nan' in the unfiltered_clip_struct 
% (if present) and replaced with 0 in new_clip_structs.  

% figure out what optional properties have been requested.
is_hrc_filter = 1;
is_verbose = 0;

cnt = 1;
while cnt <= nargin-2,
    if strcmp(lower(varargin{cnt}),'nohrcfilter'),
        is_hrc_filter = 0;
        cnt = cnt + 1;
    elseif strcmp(lower(varargin{cnt}),'verbose'),
        is_verbose = 1;
        cnt = cnt + 1;
    elseif strcmp(lower(varargin{cnt}),'quiet'),
        is_verbose = 0;
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
    % Set spatial registration of the original to (0,0), by definition.
    curr_offsets = offsets{index};
    if strcmp(clip_structs(curr_offsets(1)).hrc,'original'),
        for cnt = 1:length(curr_offsets),
            new_clip_structs(curr_offsets(cnt)).spatial.horizontal = 0;
            new_clip_structs(curr_offsets(cnt)).spatial.vertical = 0;
        end; 
        continue;
    end
    
    % loop through each clip of this HRC.  Find its spatial registration.
    for cnt = 1:length(curr_offsets),
        orig_offset = find_original(clip_structs,curr_offsets(cnt));
        new_clip_structs(curr_offsets(cnt)).spatial = ...
            spatial_registration_one_clip(test_structs,clip_structs(orig_offset), ...
            clip_structs(curr_offsets(cnt)), is_verbose); 
    end
    
    % compute overall registration, for greater accuracy.
    if is_hrc_filter,
        unfiltered_clip_structs = new_clip_structs;
        
        % pick out spatial registration results  for this HRC.
        all = [new_clip_structs(curr_offsets).spatial];
        new_spatial.vertical = st_collapse('50%', [all.vertical]');
	    new_spatial.horizontal = st_collapse('50%', [all.horizontal]');

        % print result.
        if is_verbose,
            fprintf('Median filter:  v=%d, h=%d\n', new_spatial.vertical, new_spatial.horizontal);
        end
        
        % copy results back to each of these clips.
        for cnt = 1:length(curr_offsets),
            new_clip_structs(curr_offsets(cnt)).spatial.horizontal = new_spatial.horizontal;
            new_clip_structs(curr_offsets(cnt)).spatial.vertical = new_spatial.vertical;
        end
    end
end

% replace 'nan' with 0 shift in returned structure.  
status = 0;
for cnt=1:length(new_clip_structs),
    if isnan(new_clip_structs(cnt).spatial.horizontal) || isnan(new_clip_structs(cnt).spatial.vertical),
        status = status + 1;
        new_clip_structs(cnt).spatial.horizontal = 0;
        new_clip_structs(cnt).spatial.vertical = 0;
        if is_verbose,
            fprintf('Spatial shift calulation failed for %s:%s(%s)\n', ...
                new_clip_structs(cnt).test{1}, ...
                new_clip_structs(cnt).scene{1}, ...
                new_clip_structs(cnt).hrc{1});
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_spatial] = spatial_registration_one_clip(test_structs,orig_clip, proc_clip, is_verbose);
% Compute spatial registration for one clip.

if strcmp(orig_clip.video_standard,'progressive'),
    is_field = 0;
else
    is_field = 1;
end

% clear out any spatial registration (etc. temporal registration guess).
orig_clip = clear_spatial_registration(orig_clip);
proc_clip = clear_spatial_registration(proc_clip);

% set up constants.
max_delta = round(orig_clip.fps); % plus or minus 30 frames, constant for now.
max_delta_shift = max_delta; % consider one processed frame a second.
max_search = 20; % maximum search 20 pixels in any direction.

% compute PVR and OROI given the above.
[pvr, oroi] = find_pvr_oroi_guess(proc_clip, max_search, is_field);

% error checks & corrections
if mod(max_search,2),
    max_search = max_search + 1;
end

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

if strcmp(orig_clip.video_standard,'progressive'),
	new_spatial = spatial_registration_by_frames(test_structs,orig_clip, ...
        proc_clip, max_delta, max_delta_shift, max_search, pvr, oroi, is_verbose);
elseif strcmp(orig_clip.video_standard,'interlace_lower_field_first') | strcmp(orig_clip.video_standard,'interlace_upper_field_first'),
	new_spatial = spatial_registration_by_fields(test_structs,orig_clip, ...
        proc_clip, max_delta, max_delta_shift, max_search, pvr, oroi, is_verbose);
else
    error('Video standard not recognized for clip %s:%s(%s)', proc_clip.test{1}, proc_clip.scene{1}, proc_clip.hrc{1});
end

if is_verbose,
    fprintf('\t%s:%s(%s) new shift (v=%d, h=%d)\n', ...
        proc_clip.test{1}, proc_clip.scene{1}, proc_clip.hrc{1}, ...
        new_spatial.vertical, new_spatial.horizontal);
end

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
function [new_spatial] = spatial_registration_by_fields(test_structs,orig_clip, ...
    proc_clip, max_delta, max_delta_shift, max_search, pvr, oroi, is_verbose);
% Compute spatial registration for one clip of fields..

is_field = 1;
% loop over all frames of this clip.  Find a first good spatial
% registration.
if is_verbose,
    fprintf('%s:%s(%s)\n', proc_clip.test{1}, proc_clip.scene{1}, proc_clip.hrc{1});
end
total = total_tslices(proc_clip,1/proc_clip.fps);

for proc_frame = (max_delta + 1):max_delta_shift:(total-max_delta),
    
    % read in this original time-slice.  
    % If working in fields, split appart now.
    [orig_y, one_y, two_y, len_tslice, orig_frame] = read_in_needed_fields(max_delta, ...
        test_structs, orig_clip, proc_frame, oroi, proc_clip, pvr);
        
    % do broad search for temporal shift.
    [best_time_f1(1), best_vert_f1(1), best_horiz_f1(1), best_std] = ...
        search_spatial_temporal_range( orig_y, one_y, ...
        [1:4:(max_delta*4+2)], [-8 0 0 8], [0 -8 0 0], is_verbose);

    % if best time is against a boundary, move it in 2 frames
    [best_time_f1(1)] = move_in_from_end(best_time_f1(1), max_delta*4+2, 4);

    % do broad search for spatial shift.
    ssh = [-12 -8 -6 -4 2 0 2 4 6 8 12  -8 -4  0  4  8  -8 -4  0  4  8  0  0  -8 8];
    ssv = [  0  0  0  0 0 0 0 0 0 0  0  -4 -4 -4 -4 -4  -8 -8 -8 -8 -8  1 -1   8 8];
    list_of_times = (best_time_f1(1)-4):2:(best_time_f1(1)+4);
    [best_time_f1(1), best_vert_f1_f1(1), best_horiz_f1(1), best_std] = ...
        search_spatial_temporal_range( orig_y, one_y, list_of_times, ssh, ssv, is_verbose);
    
    % if best time is against a boundary, move it in 1 frames
    [best_time_f1(1)] = move_in_from_end(best_time_f1(1), max_delta*4+2, 2);

    % do repeated fine searches.
    [best_time_f1(1), best_vert_f1(1), best_horiz_f1(1), best_std, found] = do_repeat_fine_search_field(...
        best_time_f1(1), best_vert_f1(1), best_horiz_f1(1), max_search, orig_clip, one_y, orig_y, 5, is_verbose);   

    % record results, if there are any.
    if found,
        break;
    end
end

if ~found,
    new_spatial.vertical = NaN;
    new_spatial.horizontal = NaN;
    return;
end

% note that best...(1) contains the initial estimate & won't be used for
% median filtering.  Record above results as f1 and f2 starting points.
best_time_f2(1) = best_time_f1(1);
best_vert_f2(1) = best_vert_f1(1);
best_horiz_f2(1) = best_horiz_f1(1);


% Now, find shift for each frequency'th frame.
cnt1 = 2;
cnt2 = 2;
first_time = 1;
for proc_frame = (proc_frame):max_delta_shift:(total-max_delta),

    % read in this original time-slice.  split appart into fields.
    % Don't re-read first frame, alread read in.
    if first_time,
        first_time = 0;
    else
        [orig_y, one_y, two_y, len_tslice, orig_frame] = read_in_needed_fields(max_delta, ...
            test_structs, orig_clip, proc_frame, oroi, proc_clip, pvr);
    end

    % try 3 fine searches for field one shift.
    [best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1), best_std, found] = do_repeat_fine_search_field(...
        best_time_f1(cnt1-1), best_vert_f1(cnt1-1), best_horiz_f1(cnt1-1), max_search, orig_clip, one_y, orig_y, 3, is_verbose);   
% fprintf('\tDEBUG: 3fine t=%f,v=%d,h=%d\n', best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1));
% if ~found, fprintf('\tfine failed\n');    end;

    if ~found,
        % do broad search for temporal shift, including last valid spatial shift.
        [best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1), best_std] = ...
            search_spatial_temporal_range( orig_y, one_y, ...
            [1:4:(max_delta*4+2)], [-8 0 0 8 best_horiz_f1(cnt1-1)], [0 -8 0 0 best_vert_f1(cnt1-1)], is_verbose);
% fprintf('\tDEBUG: times t=%f,v=%d,h=%d\n', best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1));
% if ~found, fprintf('\time failed\n');   end; 
	
        % do repeated fine searches.
        [best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1), best_std, found] = do_repeat_fine_search_field(...
            best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1), max_search, orig_clip, one_y, orig_y, 5, is_verbose);   
% fprintf('\tDEBUG: 5fine t=%f,v=%d,h=%d\n', best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1));
% if ~found, fprintf('\tfine failed\n');   end; 
    end

    % record field one results, if there are any.
    if found,
%         fprintf('\t\tproc F%d (f1) -> orig F%d (f%d), v=%d h=%d;  ', ...
%             proc_frame, orig_frame+ceil(best_time_f1(cnt1)/2), 2 - mod(best_time_f1(cnt1),2), ...
%             best_vert_f1(cnt1), best_horiz_f1(cnt1));
        cnt1 = cnt1 + 1;
    end
    
    
    % try 3 fine searches for field two.
    [best_time_f2(cnt2), best_vert_f2(cnt2), best_horiz_f2(cnt2), best_std, found] = do_repeat_fine_search_field(...
        best_time_f2(cnt2-1), best_vert_f2(cnt2-1), best_horiz_f2(cnt2-1), max_search, orig_clip, two_y, orig_y, 3, is_verbose);   
    
    if ~found,
        % do broad search for temporal shift, including last valid spatial shift.
        [best_time_f2(cnt2), best_vert_f2(cnt2), best_horiz_f2(cnt2), best_std] = ...
            search_spatial_temporal_range( orig_y, two_y, ...
            [1:4:(max_delta*4+2)], [-8 0 0 8 best_horiz_f2(cnt2-1)], [0 -8 0 0 best_vert_f2(cnt2-1)], is_verbose);
	
        % do repeated fine searches.
        [best_time_f2(cnt2), best_vert_f2(cnt2), best_horiz_f2(cnt2), best_std, found] = do_repeat_fine_search_field(...
            best_time_f2(cnt2), best_vert_f2(cnt2), best_horiz_f2(cnt2), max_search, orig_clip, two_y, orig_y, 5, is_verbose);   
    end

    % record results, if there are any.
    if found,
%         fprintf('proc F%d (f2) -> orig F%d (f%d), v=%d h=%d\n', ...
%             proc_frame, orig_frame+ceil(best_time_f2(cnt2)/2), 2 - mod(best_time_f2(cnt2),2), ...
%             best_vert_f2(cnt2), best_horiz_f2(cnt2));
        cnt2 = cnt2 + 1;
    else
%         fprintf('\n');
    end
end

if cnt2 <= 2 | cnt1 <= 2,
    new_spatial.vertical = NaN;
    new_spatial.horizontal = NaN;
else
    if is_verbose,
        % median filter the results from each frame.
        fprintf('F1 V: ');
        for i=1:cnt1-1,
            fprintf('%d ', best_vert_f1(i));
        end
        fprintf('\nF2 V: ');
        for i=1:cnt2-1,
            fprintf('%d ', best_vert_f2(i));
        end
        fprintf('\nF1 H: ');
        for i=1:cnt1-1,
            fprintf('%d ', best_horiz_f1(i));
        end
    end

    new_spatial.vertical = st_collapse('50%', best_vert_f1(2:cnt1-1)') + st_collapse('50%', best_vert_f2(2:cnt2-1)');
	new_spatial.horizontal = st_collapse('50%', best_horiz_f1(2:cnt1-1)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clip] = clear_spatial_registration(old_clip);
% Erase all spatial registration information in a clip, except for temporal
% alignment (which is presumed to be a reasonable guess) and stretch
% registration.

clip = old_clip;

clip.spatial.vertical = 0;
clip.spatial.horizontal = 0;
clip.luminance_offset = 0;
clip.luminance_gain = 1.0;
clip.cvr.top = 1;
clip.cvr.left = 1;
clip.cvr.bottom = clip.image_size.rows;
clip.cvr.right = clip.image_size.cols;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [std_of_diff] = compare_proi_with_oroi (orig_y, proc_y, shift, is_verbose);
% Compare original (containing just ORIO) with one shift of the processed 
% image ('proc_y') corresponding to one PROI.  Use standard deviation of
% the difference.  Minimize to get the best match.

% find size of OROI
[o_rows, o_cols] = size(orig_y);
[p_rows, p_cols] = size(proc_y);
max_search_rows = (p_rows - o_rows) / 2;
max_search_cols = (p_cols - o_cols) / 2;

% if request is out of range, warn & return NaN.
if  abs(shift.vertical) > max_search_rows | abs(shift.horizontal) > max_search_cols,
    std_of_diff = NaN;
    if is_verbose,
        fprintf('WARNING:  searched outside of legal range.\n');
    end
    return;
end

% Find location of shifted OROI in the processed image.
top_proc = max_search_rows + shift.vertical;
rows = (top_proc+1):(top_proc+o_rows);

left_proc = max_search_cols + shift.horizontal;
cols = (left_proc+1):(left_proc+o_cols);

% Pick out the OROI from the processed image.  Subtract OROI and return the
% standard deviation of the difference.

diff = orig_y - proc_y(rows,cols);
std_of_diff = std(reshape(diff, o_rows*o_cols, 1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gain] = compute_gain_level_offset (orig_y, proc_y, max_search, shift);
% Compute gain and level offset.  The original image, 'orig_y', must
% contain only the OROI.  The processed image, 'proc_y', must contain OROI
% with a border of 'max_search' pixels on all sides.  Variable 'shift' contains
% the current best guess for the processed image's search, 'shift.horizontal' 
% and 'shift.vertical', where a positive number means the processed image
% has been shifted down or to the right.
%
% Despite the name, level offset is not computed, since it is not needed.

% find size of OROI
[o_rows, o_cols] = size(orig_y);

% Find location of shifted OROI in the processed image.
top_proc = max_search + shift.vertical;
rows = (top_proc+1):(top_proc+o_rows);

left_proc = max_search + shift.horizontal;
cols = (left_proc+1):(left_proc+o_cols);

% Pick out the OROI from the processed image, and compute STD and MEAN.
proc_oroi = reshape( proc_y(rows,cols), o_rows*o_cols, 1);
std_proc = std(proc_oroi);

% Compute mean & std of original, too.
orig_oroi = reshape( orig_y, o_rows*o_cols, 1);
std_orig = std(orig_oroi);

% compute gain.
gain = std_proc / std_orig;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_time, best_vert, best_horiz, best_std] = do_one_fine_search_field(...
    old_best_time, old_best_vert, old_best_horiz, max_search, ...
    orig_clip, proc_y, orig_y, is_verbose);
%

is_field = 1;

% do fine search for spatial shift.
[ssv, ssh] = find_list_of_fine_shifts(old_best_vert, old_best_horiz,...
    max_search, is_field);
if strcmp(orig_clip.video_standard,'interlace_lower_field_first'),
    list_of_times = (old_best_time-2):1:(old_best_time+2);
else
    % pal.  Fields are out of order!  Draw a picture.  Each f1&f2 must
    % be swapped, and then do the above.
    if mod(old_best_time,2),
        list_of_times = [old_best_time-2 old_best_time:1:(old_best_time+2)];
    else
        list_of_times = [(old_best_time-2):1:old_best_time (old_best_time+2)];
    end
end
[best_time, best_vert, best_horiz, best_std] = ...
    search_spatial_temporal_range( orig_y, proc_y, ...
    list_of_times, ssh, ssv, is_verbose);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_time, best_vert, best_horiz, best_std, found] = do_repeat_fine_search_field(...
    old_best_time, old_best_vert, old_best_horiz, max_search, ...
    orig_clip, proc_y, orig_y, total, is_verbose);

% do repeated fine searches.  
%

% initialize return variables and previous results.
best_time = NaN;
best_vert = NaN;
best_horiz = NaN;
best_std = NaN;

btime(1) = NaN;
bvert(1) = NaN;
bhoriz(1) = NaN;

btime(2) = old_best_time;
bvert(2) = old_best_vert;
bhoriz(2) = old_best_horiz;

% for each of 'total' tries, do a fine search.  Check end conditions.
for cnt=3:(2+total),
    [btime(cnt), bvert(cnt), bhoriz(cnt), best_std] = do_one_fine_search_field(...
        btime(cnt-1), bvert(cnt-1), bhoriz(cnt-1), max_search, orig_clip, proc_y, orig_y, is_verbose);
    if (btime(cnt) == btime(cnt-1) & bvert(cnt) == bvert(cnt-1) & bhoriz(cnt) == bhoriz(cnt-1)) | ...
        (btime(cnt) == btime(cnt-2) & bvert(cnt) == bvert(cnt-2) & bhoriz(cnt) == bhoriz(cnt-2)),
        found = 1;
        best_time = btime(cnt);
        best_vert = bvert(cnt);
        best_horiz = bhoriz(cnt);
        return;
    end
end
found = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [list_of_vert, list_of_horiz] = find_list_of_fine_shifts(prev_vert, prev_horiz,...
    max_search, is_field);
% Given a previous shift, find the list of vertical & horizontal shifts to
% be considered.  Keep within search boundary.

% list of all fine shifts around the current point.
hshifts = [-2 -2 -2 -1 -1 -1  0  0  0  0  0  1  1  1  2  2  2] + prev_horiz;
vshifts = [-2  0  2 -1  0  1 -2 -1  0  1  2 -1  0  1  -2 0  2] + prev_vert;

% find limits of search range
maxv = max_search;
if is_field,
    maxh = max_search / 2;
else
    maxh = max_search;
end

% Always search 0,0 shift.
list_of_vert =  [0];
list_of_horiz = [0];

% add in h & v shifts, unless out of range.
for cnt = 1:length(hshifts),
    if abs(vshifts(cnt)) <= maxv & abs(hshifts(cnt)) <= maxh,
        list_of_vert =  [list_of_vert  vshifts(cnt)];
        list_of_horiz = [list_of_horiz hshifts(cnt)];
    end
end


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
function [pvr, oroi] = find_pvr_oroi_guess(clip, max_search, is_field);
% Find best guess at PVR and OROI, given image size (within clip) and
% size of maximum search ('max_search'), in pixels.
% 'is_field' is true for 'field', 1 for 'frame', depending on which is desired.

pvr = find_pvr_guess(clip);
oroi.top = pvr.top + max_search;
oroi.left = pvr.left + max_search;
oroi.bottom = pvr.bottom - max_search;
oroi.right = pvr.right - max_search;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [best_delta, best_horiz, best_vert, best_std] = fine_search_spatial_temporal_range( ...
%     test_structs,orig_clip, proc_clip, orig_frame, proc_frame, ...
%     prev_shift, max_search, fchoice, field);
% % Search a fine local regionof spatial-temporal shifts for best h-v
% % spatial shift.   
% %
% % 'test_structs' is as GTests
% % 'orig_clip' is an original clip as in GClips
% % 'proc_clip' is an processed clip as in GClips, whose original clip is
% %       specified in 'orig_clip'.
% % 'orig_frame' is the original clip's frame number.
% % 'proc_frame' is the processed clip's frame number.
% % 'prev_shift' is the previous best shift.
% % 'max_search' is the maximum search both horizontally & vertically, 
% % 'fchoice' is 'field' or 'frame', depending if should use fields or frames
% %   when searching.
% % 'field' is ignored when 'fchoice' is frame.  1 for field one, 2 for field
% %   two.
% [pvr, oroi] = find_pvr_oroi_guess(proc_clip, max_search, fchoice);
% 
% % 
% v = prev_shift.vertical;
% h = prev_shift.horizontal;
% list_of_vert =  [0];
% list_of_horiz = [0];
% for i= -1:1,
%     for j= -1:1,
%         list_of_vert =  [list_of_vert i];
%         list_of_horiz = [list_of_horiz j];
%     end
% end
% for i= -2:2:2,
%     for j= -2:2:2,
%         if i ~= 0 | j~= 0,
%             list_of_vert =  [list_of_vert i];
%             list_of_horiz = [list_of_horiz j];
%         end
%     end
% end
% 
% % initialize variables.
% best_horiz = 0;
% best_vert = 0;
% best_std = inf;
% 
% % read in the processed frame.  Keep only current field.
% if strcmp(fchoice,'field'),
% 	out = read_tslice(test_structs,proc_clip, 1/proc_clip.fps, proc_frame, ...
%         'sroi', pvr.top*2-1, pvr.left, pvr.bottom*2, pvr.right);
%     [one, two] = split_into_fields(out);
%     if field == 1,
%         out = one;
%     else
%         out = two;
%     end
%     
%     best_std = inf;
%     
%     % consider every other input frame's field one.
% 	for loop = -1:1,
%         % read in the input field
% 	    in = read_tslice(test_structs,orig_clip, 1/orig_clip.fps, orig_frame + loop, ...
%            'sroi', oroi.top*2-1, oroi.left, oroi.bottom*2, oroi.right );
% 	    [in_one,in_two] = split_into_fields(in);
%         
%         % search field one.
%         if field == 2 & loop == -1,
%             % skip field one, because only want +- two fields.
%         else
%             [best_std, best_horiz, best_vert, changed] = loop_over_shifts_image_pair( ...
%                 in_one, out, list_of_horiz, list_of_vert, max_search, ...
%                 best_std, best_horiz, best_vert, is_verbose);
%             if changed,
%                 best_delta = loop;
%             end
%         end
%         % search field one.
%         if field == 1 & loop == 1,
%             % skip field two, because only want +- two fields.
%         else
%             [best_std, best_horiz, best_vert, changed] = loop_over_shifts_image_pair( ...
%                 in_two, out, list_of_horiz, list_of_vert, max_search, ...
%                 best_std, best_horiz, best_vert, is_verbose);
%             if changed,
%                 best_delta = loop;
%             end
%         end
% 	end
% 
% else
%     out = read_tslice(test_structs,proc_clip, 1/proc_clip.fps, proc_frame, ...
%         'sroi', pvr.top, pvr.left, pvr.bottom, pvr.right);
% 	for loop = -1:1,
%         % read in the input field
% 	    in = read_tslice(test_structs,orig_clip, 1/orig_clip.fps, orig_frame + loop, ...
%            'sroi', oroi.top, oroi.left, oroi.bottom, oroi.right );
%         [best_std, best_horiz, best_vert, changed] = loop_over_shifts_image_pair( ...
%               in, out, list_of_horiz, list_of_vert, max_search, ...
%               best_std, best_horiz, best_vert, is_verbose);
%         if changed,
%             best_delta = loop;
%         end
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_std, best_horiz, best_vert, changed] = ...
    loop_over_shifts_image_pair(in, out, list_of_horiz, list_of_vert, max_search,...
                    curr_std, curr_horiz, curr_vert, is_verbose);
% Loop over all shifts for one image pair, finding the best shift.
% update from previous best.

% for each shift, compute comparison statistic
best_std = curr_std;
best_horiz = curr_horiz;
best_vert = curr_vert;
changed = 0;
for cnt = 1:length(list_of_vert),
    shift.vertical = list_of_vert(cnt);;
    shift.horizontal = list_of_horiz(cnt);
    
    curr = compare_proi_with_oroi (in, out, max_search, shift, is_verbose);
    %fprintf('Delay %d, horiz %d vert %d --> compares at %f\n', ...
    %    loop, list_of_horiz(cnt), list_of_vert(cnt), curr);
    
    % update best found.
    if curr < best_std,
        best_std = curr;
        best_horiz = list_of_horiz(cnt);
        best_vert = list_of_vert(cnt);
        changed = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_time] = move_in_from_end(real_best_time, array_length, distance);
% move best_time in from ends of array the requested number of units.

best_time = real_best_time;

% if best time is against a boundary, move it in 'distance' frames
while best_time <= distance,
    best_time = best_time + distance;
end

while best_time > array_length - distance,
    best_time = best_time - distance;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [orig_y, one_y, two_y, len_tslice, orig_frame] = read_in_needed_fields(max_delta, ...
    test_structs, orig_clip, proc_frame, oroi, proc_clip, pvr)
%

% read in this original time-slice.  
% If working in fields, split appart now.
for cnt=1:max_delta*2+1,
    hold = read_tslice(test_structs, orig_clip, 1/orig_clip.fps, ...
        (proc_frame-max_delta) + (cnt-1), 'sroi', ...
        oroi.top, oroi.left, oroi.bottom, oroi.right);
    [one,two] = split_into_fields(hold);
    orig_y(:,:,cnt*2-1) = one;
    orig_y(:,:,cnt*2) = two;
    % IGNORE ntsc/pal ordering.  field numbering will be correct,
    % and the few PAL clips will simply have the time-history of
    % each field switched.  This will only be a issue during fine
    % searches.
end
len_tslice = max_delta*4+2;
orig_frame = proc_frame - max_delta - 1;

% read in the processed field (or frame)
proc_y(:,:) = read_tslice(test_structs, proc_clip, 1/proc_clip.fps, ...
    proc_frame, 'sroi', pvr.top, pvr.left, pvr.bottom, pvr.right);
[one_y, two_y] = split_into_fields(proc_y);
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_delta, best_vert, best_horiz, best_std] = search_spatial_temporal_range( ...
    orig_y, proc_y, list_of_times, list_of_horiz, ...
    list_of_vert, is_verbose);
% Search the specified range of spatial-temporal shifts for best h-v
% spatial shift.  If request field search, will only consider fields of the
% same type (1 or 2).
%
% 'orig_y' is as original frames or fields, xyt
% 'proc_y' is processed frames or fields, xyt
% 'list_of_times' is the list of array offsets that should be searched.  
% 'list_of_horiz' is the list of horizontal shifts.  See list_of_vert.
% 'list_of_vert' is the list of vertical shifts.  Paired item-wise with
%   elements from 'list_of_horiz'.
% 'max_search' is the maximum search both horizontally & vertically, 

best_std = inf;

% Make sure list_of_times doesn't contain any out-of-range values.
[rows,cols,time] = size(orig_y);
list_of_times = list_of_times(find(list_of_times >= 1 & list_of_times <= time));

% consider every other input frame's field one.
for loop = list_of_times,
    % read in the input field
    
    % for each shift, compute comparison statistic
    for cnt = 1:length(list_of_vert),
        shift.vertical = list_of_vert(cnt);;
        shift.horizontal = list_of_horiz(cnt);
        
        curr = compare_proi_with_oroi (orig_y(:,:,loop), proc_y, shift, is_verbose);
        %fprintf('offset %d, horiz %d vert %d --> compares at %f\n', ...
        %    loop, list_of_horiz(cnt), list_of_vert(cnt), curr);
        
        % update best found.
        if curr < best_std,
            best_std = curr;
            best_delta = loop;
            best_horiz = list_of_horiz(cnt);
            best_vert = list_of_vert(cnt);
        end
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_spatial] = spatial_registration_by_frames(test_structs,orig_clip, ...
    proc_clip, max_delta, max_delta_shift, max_search, pvr, oroi, is_verbose);
% Compute spatial registration for one clip of frames..

% loop over all frames of this clip.  Find a first good spatial
% registration.
if is_verbose,
    fprintf('%s:%s(%s)\n', proc_clip.test{1}, proc_clip.scene{1}, proc_clip.hrc{1});
end
total = total_tslices(proc_clip,1/proc_clip.fps);

for proc_frame = (max_delta + 1):max_delta_shift:(total-max_delta),
    
    % read in this original time-slice.  
    % If working in fields, split appart now.
    [orig_y, proc_y, len_tslice, orig_frame] = read_in_needed_frames(max_delta, ...
        test_structs, orig_clip, proc_frame, oroi, proc_clip, pvr);
        
    % do broad search for temporal shift.
    [best_time_fr(1), best_vert_fr(1), best_horiz_fr(1), best_std] = ...
        search_spatial_temporal_range( orig_y, proc_y, ...
        [1:2:(max_delta*2+2)], [-8 0 0 8], [0 -8 0 0], is_verbose);

    % if best time is against a boundary, move it in 2 frames
    [best_time_fr(1)] = move_in_from_end(best_time_fr(1), max_delta*2+2, 2);

    % do broad search for spatial shift.
    ssh = [-12 -8 -6 -4 2 0 2 4 6 8 12  -8 -4  0  4  8  -8 -4  0  4  8  0  0  -8 8];
    ssv = [  0  0  0  0 0 0 0 0 0 0  0  -4 -4 -4 -4 -4  -8 -8 -8 -8 -8  1 -1   8 8];
    list_of_times = (best_time_fr(1)-2):2:(best_time_fr(1)+2);
    [best_time_fr(1), best_vert_fr_fr(1), best_horiz_fr(1), best_std] = ...
        search_spatial_temporal_range( orig_y, proc_y, list_of_times, ssh, ssv, is_verbose);
    
    % if best time is against a boundary, move it in 1 frames
    [best_time_fr(1)] = move_in_from_end(best_time_fr(1), max_delta*2+2, 1);

    % do repeated fine searches.
    [best_time_fr(1), best_vert_fr(1), best_horiz_fr(1), best_std, found] = do_repeat_fine_search_field(...
        best_time_fr(1), best_vert_fr(1), best_horiz_fr(1), max_search, orig_clip, proc_y, orig_y, 5, is_verbose);   

    % record results, if there are any.
    if found,
        break;
    end
end

if ~found,
    new_spatial.vertical = NaN;
    new_spatial.horizontal = NaN;
    return;
end

% note that best...(1) contains the initial estimate & won't be used for
% median filtering.  

% Now, find shift for each frequency'th frame.
cnt = 2;
first_time = 1;
for proc_frame = (proc_frame):max_delta_shift:(total-max_delta),

    % read in this original time-slice.  split appart into fields.
    % Don't re-read first frame, alread read in.
    if first_time,
        first_time = 0;
    else
        [orig_y, proc_y, len_tslice, orig_frame] = read_in_needed_frames(max_delta, ...
            test_structs, orig_clip, proc_frame, oroi, proc_clip, pvr);
    end

    % try 3 fine searches for field one shift.
    [best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt), best_std, found] = do_repeat_fine_search_field(...
        best_time_fr(cnt-1), best_vert_fr(cnt-1), best_horiz_fr(cnt-1), max_search, orig_clip, proc_y, orig_y, 3, is_verbose);   
% fprintf('\tDEBUG: 3fine t=%f,v=%d,h=%d\n', best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt));
% if ~found, fprintf('\tfine failed\n');    end;

    if ~found,
        % do broad search for temporal shift, including last valid spatial shift.
        [best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt), best_std] = ...
            search_spatial_temporal_range( orig_y, proc_y, ...
            [1:2:(max_delta*2+1)], [-8 0 0 8 best_horiz_fr(cnt-1)], [0 -8 0 0 best_vert_fr(cnt-1)], is_verbose);
% fprintf('\tDEBUG: times t=%f,v=%d,h=%d\n', best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt));
% if ~found, fprintf('\time failed\n');   end; 
	
        % do repeated fine searches.
        [best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt), best_std, found] = do_repeat_fine_search_field(...
            best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt), max_search, orig_clip, proc_y, orig_y, 5, is_verbose);   
% fprintf('\tDEBUG: 5fine t=%f,v=%d,h=%d\n', best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt));
% if ~found, fprintf('\tfine failed\n');   end; 
    end

    % record field one results, if there are any.
    if found,
%         fprintf('\t\tproc F%d (f1) -> orig F%d (f%d), v=%d h=%d;  ', ...
%             proc_frame, orig_frame+ceil(best_time_fr(cnt)), 1 - best_time_fr(cnt) ...
%             best_vert_fr(cnt), best_horiz_fr(cnt));
        cnt = cnt + 1;
    end
    
end

if cnt <= 2,
    new_spatial.vertical = NaN;
    new_spatial.horizontal = NaN;
else
    if is_verbose,
        % median filter the results from each frame.
        fprintf('F V: ');
        for i=1:cnt-1,
            fprintf('%d ', best_vert_fr(i));
        end
        fprintf('\nF H: ');
        for i=1:cnt-1,
            fprintf('%d ', best_horiz_fr(i));
        end
        fprintf('\n');
    end

    new_spatial.vertical = st_collapse('50%', best_vert_fr(2:cnt-1)');
	new_spatial.horizontal = st_collapse('50%', best_horiz_fr(2:cnt-1)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [orig_y,proc_y, len_tslice, orig_frame] = read_in_needed_frames(max_delta, ...
    test_structs, orig_clip, proc_frame, oroi, proc_clip, pvr)
%

% read in this original time-slice.  
for cnt=1:max_delta*2+1,
    orig_y(:,:,cnt) = read_tslice(test_structs, orig_clip, 1/orig_clip.fps, ...
        (proc_frame-max_delta) + (cnt-1), 'sroi', ...
        oroi.top, oroi.left, oroi.bottom, oroi.right);
end
len_tslice = max_delta*2+1;
orig_frame = proc_frame - max_delta - 1;

% read in the processed frame
proc_y(:,:) = read_tslice(test_structs, proc_clip, 1/proc_clip.fps, ...
    proc_frame, 'sroi', pvr.top, pvr.left, pvr.bottom, pvr.right);
   

