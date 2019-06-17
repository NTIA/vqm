function [horiz, vert, status] = dll_itu_spatial(uncert_sec, frequency); 
% DLL_ITU_SPATIAL
%   Automatic computation of spatial registration.
% SYNTAX
%   [horiz, vert] = dll_itu_spatial; 
%   [horiz, vert] = dll_itu_spatial(uncert_sec, frequency); 
% DESCRIPTION
%  Perform spatial registration on an original and processed video sequence.
%  Use standard spatial registration algorithm from 2003 NTIA report 
%  & ANSI Rec. 801.03.
%
%  Precondition:  dll_video must have been initialized with fn=1 for the
%  original video file, and fn=2 for the processed video file.
%
%  Optional input argument 'uncert_sec' contains the alignment uncertanty, in
%  seconds.  By default 'uncert_sec' = 1.0 (one second).
%  Optional input argument 'frequency' contains the frequency of images used, in
%  seconds.  By default 'frequency' = 1.0 (images 1sec appart are used).
%
%  Return values are defined as follows:
%   'horiz'     Horizontal shift in pixels
%   'vert'      Vertical shift in frame lines
%   'status'    Status, 1 for success & 0 for failure.

% handle optional properties.
if ~exist('frequency','var'),
    frequency = 1.0;
end
if ~exist('uncert_sec','var'),
    uncert_sec = 1.0;
end

fn1=1;

[video_standard] = dll_video('get_video_standard', fn1);

% set up constants.
[rows,cols, fps] = dll_video('size',fn1);
max_delta = round(fps * uncert_sec); % plus or minus 1-sec (or as chosen)
max_delta_shift = round(fps * frequency); % consider one processed frame a second (or as chosen).
max_search = 20; % maximum search 20 pixels in any direction.

% compute PVR and OROI given the above.
pvr = find_pvr_guess(rows,cols);
oroi.top = pvr.top + max_search;
oroi.left = pvr.left + max_search;
oroi.bottom = pvr.bottom - max_search;
oroi.right = pvr.right - max_search;


% error checks & corrections
if mod(max_search,2),
    max_search = max_search + 1;
end

if strcmp(video_standard,'progressive'),
	[horiz, vert, status] = spatial_registration_by_frames( ...
        max_delta, max_delta_shift, max_search, pvr, oroi);
elseif strcmp(video_standard,'interlace_lower_field_first') | strcmp(video_standard,'interlace_upper_field_first'),
	[horiz, vert, status] = spatial_registration_by_fields( ...
        max_delta, max_delta_shift, max_search, pvr, oroi);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [horiz, vert, status] = spatial_registration_by_fields(...
    max_delta, max_delta_shift, max_search, pvr, oroi);
% Compute spatial registration for one clip of fields..

is_field = 1;
% loop over all frames of this clip.  Find a first good spatial
% registration.
total = min( dll_video('total_frames',1), dll_video('total_frames',2));

for proc_frame = (max_delta + 1):max_delta_shift:(total-max_delta),
    
    % read in this original time-slice.  
    % If working in fields, split appart now.
    [orig_y, one_y, two_y, len_tslice, orig_frame] = read_in_needed_fields(max_delta, ...
        proc_frame, oroi, pvr);
        
    % do broad search for temporal shift.
    [best_time_f1(1), best_vert_f1(1), best_horiz_f1(1), best_std] = ...
        search_spatial_temporal_range( orig_y, one_y, ...
        [1:4:(max_delta*4+2)], [-8 0 0 8], [0 -8 0 0]);

    % if best time is against a boundary, move it in 2 frames
    [best_time_f1(1)] = move_in_from_end(best_time_f1(1), max_delta*4+2, 4);

    % do broad search for spatial shift.
    ssh = [-12 -8 -6 -4 2 0 2 4 6 8 12  -8 -4  0  4  8  -8 -4  0  4  8  0  0  -8 8];
    ssv = [  0  0  0  0 0 0 0 0 0 0  0  -4 -4 -4 -4 -4  -8 -8 -8 -8 -8  1 -1   8 8];
    list_of_times = (best_time_f1(1)-4):2:(best_time_f1(1)+4);
    [best_time_f1(1), best_vert_f1_f1(1), best_horiz_f1(1), best_std] = ...
        search_spatial_temporal_range( orig_y, one_y, list_of_times, ssh, ssv);
    
    % if best time is against a boundary, move it in 1 frames
    [best_time_f1(1)] = move_in_from_end(best_time_f1(1), max_delta*4+2, 2);

    % do repeated fine searches.
    [best_time_f1(1), best_vert_f1(1), best_horiz_f1(1), best_std, found] = do_repeat_fine_search_field(...
        best_time_f1(1), best_vert_f1(1), best_horiz_f1(1), max_search, one_y, orig_y, 5);   

    % record results, if there are any.
    if found,
        break;
    end
end

if ~found,
    vert = NaN;
    horiz = NaN;
    status = 0;
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
            proc_frame, oroi, pvr);
    end

    % try 3 fine searches for field one shift.
    [best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1), best_std, found] = do_repeat_fine_search_field(...
        best_time_f1(cnt1-1), best_vert_f1(cnt1-1), best_horiz_f1(cnt1-1), max_search, one_y, orig_y, 3);   
% fprintf('\tDEBUG: 3fine t=%f,v=%d,h=%d\n', best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1));
% if ~found, fprintf('\tfine failed\n');    end;

    if ~found,
        % do broad search for temporal shift, including last valid spatial shift.
        [best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1), best_std] = ...
            search_spatial_temporal_range( orig_y, one_y, ...
            [1:4:(max_delta*4+2)], [-8 0 0 8 best_horiz_f1(cnt1-1)], [0 -8 0 0 best_vert_f1(cnt1-1)]);
% fprintf('\tDEBUG: times t=%f,v=%d,h=%d\n', best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1));
% if ~found, fprintf('\time failed\n');   end; 
	
        % do repeated fine searches.
        [best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1), best_std, found] = do_repeat_fine_search_field(...
            best_time_f1(cnt1), best_vert_f1(cnt1), best_horiz_f1(cnt1), max_search, one_y, orig_y, 5);   
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
        best_time_f2(cnt2-1), best_vert_f2(cnt2-1), best_horiz_f2(cnt2-1), max_search, two_y, orig_y, 3);   
    
    if ~found,
        % do broad search for temporal shift, including last valid spatial shift.
        [best_time_f2(cnt2), best_vert_f2(cnt2), best_horiz_f2(cnt2), best_std] = ...
            search_spatial_temporal_range( orig_y, two_y, ...
            [1:4:(max_delta*4+2)], [-8 0 0 8 best_horiz_f2(cnt2-1)], [0 -8 0 0 best_vert_f2(cnt2-1)]);
	
        % do repeated fine searches.
        [best_time_f2(cnt2), best_vert_f2(cnt2), best_horiz_f2(cnt2), best_std, found] = do_repeat_fine_search_field(...
            best_time_f2(cnt2), best_vert_f2(cnt2), best_horiz_f2(cnt2), max_search, two_y, orig_y, 5);   
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
    vert = 0;
    horiz = 0;
    status = 0;
else
    vert = st_collapse('50%', best_vert_f1(2:cnt1-1)') + st_collapse('50%', best_vert_f2(2:cnt2-1)');
	horiz = st_collapse('50%', best_horiz_f1(2:cnt1-1)');
    status = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [std_of_diff] = compare_proi_with_oroi (orig_y, proc_y, shift);
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
    proc_y, orig_y);
%

is_field = 1;

% do fine search for spatial shift.
[ssv, ssh] = find_list_of_fine_shifts(old_best_vert, old_best_horiz,...
    max_search, is_field);
[video_standard] = dll_video('get_video_standard', 1);
if strcmp(video_standard,'interlace_lower_field_first'),
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
    list_of_times, ssh, ssv);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_time, best_vert, best_horiz, best_std, found] = do_repeat_fine_search_field(...
    old_best_time, old_best_vert, old_best_horiz, max_search, ...
    proc_y, orig_y, total);

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
        btime(cnt-1), bvert(cnt-1), bhoriz(cnt-1), max_search, proc_y, orig_y);
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
function [pvr] = find_pvr_guess(rows, cols);
% return the best guess for pvr based only on image size.
% i.e., discard overscan.
% 'fchoice' is 1 for 'field' or 0 for 'frame' depending on which is desired.

if (rows == 486 | rows == 480 ) & cols == 720,
    pvr.top = 19;
    pvr.bottom = 486 - 18;
    pvr.left = 23;
    pvr.right = 720 - 22;
elseif rows == 576 & cols == 720,
    pvr.top = 15;
    pvr.bottom = 576 - 14;
    pvr.left = 23;
    pvr.right = 720 - 22;
else
    pvr.top = 1;
    pvr.bottom = rows;
    pvr.left = 1;
    pvr.right = cols;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_std, best_horiz, best_vert, changed] = ...
    loop_over_shifts_image_pair(in, out, list_of_horiz, list_of_vert, max_search,...
                    curr_std, curr_horiz, curr_vert);
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
    
    curr = compare_proi_with_oroi (in, out, max_search, shift);
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
    proc_frame, oroi, pvr)
%

% read in this original time-slice.  
% If working in fields, split appart now.
fn1 = 1;
fn2 = 2;
[rows,cols, fps] = dll_video('size',fn1);


dll_video('set_rewind',fn1);
dll_video('set_tslice',fn1, 1/fps);
dll_video('discard', fn1, (proc_frame-max_delta-1)/fps);  
for cnt=1:max_delta*2+1,

    hold = dll_video('tslice', fn1);
    hold = hold(oroi.top:oroi.bottom, oroi.left:oroi.right);

    [one,two] = split_into_fields(hold);
    orig_y(:,:,cnt*2-1) = one;
    orig_y(:,:,cnt*2) = two;
    % IGNORE ntsc/pal ordering.  field numbering will be correct,
    % and the few PAL clips will simply have the time-history of
    % each field switched.  This will only be a issue during fine
    % searches.
end
dll_video('rewind',fn1);

len_tslice = max_delta*4+2;
orig_frame = proc_frame - max_delta - 1;

% read in the processed field (or frame)
dll_video('set_rewind',fn2);
dll_video('set_tslice',fn2, 1/fps);
dll_video('discard', fn2, (proc_frame-1)/fps);  

hold = dll_video('tslice', fn2);
hold = hold(pvr.top:pvr.bottom, pvr.left:pvr.right);

[one_y, two_y] = split_into_fields(hold);
dll_video('rewind',fn2);
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [best_delta, best_vert, best_horiz, best_std] = search_spatial_temporal_range( ...
    orig_y, proc_y, list_of_times, list_of_horiz, ...
    list_of_vert);
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
        
        curr = compare_proi_with_oroi (orig_y(:,:,loop), proc_y, shift);
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
function [horiz, vert, status] = spatial_registration_by_frames( ...
    max_delta, max_delta_shift, max_search, pvr, oroi);
% Compute spatial registration for one clip of frames..

% loop over all frames of this clip.  Find a first good spatial
% registration.
total = min( dll_video('total_frames',1), dll_video('total_frames',2));

for proc_frame = (max_delta + 1):max_delta_shift:(total-max_delta),
    
    % read in this original time-slice.  
    % If working in fields, split appart now.
    [orig_y, proc_y, len_tslice, orig_frame] = read_in_needed_frames(max_delta, ...
        proc_frame, oroi, pvr);
        
    % do broad search for temporal shift.
    [best_time_fr(1), best_vert_fr(1), best_horiz_fr(1), best_std] = ...
        search_spatial_temporal_range( orig_y, proc_y, ...
        [1:2:(max_delta*2+2)], [-8 0 0 8], [0 -8 0 0]);

    % if best time is against a boundary, move it in 2 frames
    [best_time_fr(1)] = move_in_from_end(best_time_fr(1), max_delta*2+2, 2);

    % do broad search for spatial shift.
    ssh = [-12 -8 -6 -4 2 0 2 4 6 8 12  -8 -4  0  4  8  -8 -4  0  4  8  0  0  -8 8];
    ssv = [  0  0  0  0 0 0 0 0 0 0  0  -4 -4 -4 -4 -4  -8 -8 -8 -8 -8  1 -1   8 8];
    list_of_times = (best_time_fr(1)-2):2:(best_time_fr(1)+2);
    [best_time_fr(1), best_vert_fr_fr(1), best_horiz_fr(1), best_std] = ...
        search_spatial_temporal_range( orig_y, proc_y, list_of_times, ssh, ssv);
    
    % if best time is against a boundary, move it in 1 frames
    [best_time_fr(1)] = move_in_from_end(best_time_fr(1), max_delta*2+2, 1);

    % do repeated fine searches.
    [best_time_fr(1), best_vert_fr(1), best_horiz_fr(1), best_std, found] = do_repeat_fine_search_field(...
        best_time_fr(1), best_vert_fr(1), best_horiz_fr(1), max_search, proc_y, orig_y, 5);   

    % record results, if there are any.
    if found,
        break;
    end
end

if ~found,
    vert = NaN;
    horiz = NaN;
    status = 0;
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
            proc_frame, oroi, pvr);
    end

    % try 3 fine searches for field one shift.
    [best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt), best_std, found] = do_repeat_fine_search_field(...
        best_time_fr(cnt-1), best_vert_fr(cnt-1), best_horiz_fr(cnt-1), max_search, proc_y, orig_y, 3);   
% fprintf('\tDEBUG: 3fine t=%f,v=%d,h=%d\n', best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt));
% if ~found, fprintf('\tfine failed\n');    end;

    if ~found,
        % do broad search for temporal shift, including last valid spatial shift.
        [best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt), best_std] = ...
            search_spatial_temporal_range( orig_y, proc_y, ...
            [1:2:(max_delta*2+1)], [-8 0 0 8 best_horiz_fr(cnt-1)], [0 -8 0 0 best_vert_fr(cnt-1)]);
% fprintf('\tDEBUG: times t=%f,v=%d,h=%d\n', best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt));
% if ~found, fprintf('\time failed\n');   end; 
	
        % do repeated fine searches.
        [best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt), best_std, found] = do_repeat_fine_search_field(...
            best_time_fr(cnt), best_vert_fr(cnt), best_horiz_fr(cnt), max_search, proc_y, orig_y, 5);   
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
    vert = 0;
    horiz = 0;
    status = 0;
else
    vert = st_collapse('50%', best_vert_fr(2:cnt-1)');
	horiz = st_collapse('50%', best_horiz_fr(2:cnt-1)');
    status = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [orig_y,proc_y, len_tslice, orig_frame] = read_in_needed_frames(max_delta, ...
    proc_frame, oroi, pvr)
%

% read in this original time-slice.  
% If working in fields, split appart now.
fn1 = 1;
fn2 = 2;
[rows,cols, fps] = dll_video('size',fn1);


dll_video('set_rewind',fn1);
dll_video('set_tslice',fn1, 1/fps);
dll_video('discard', fn1, (proc_frame-max_delta-1)/fps);  

for cnt=1:max_delta*2+1,
    hold = dll_video('tslice', fn1);
    orig_y(:,:,cnt) = hold(oroi.top:oroi.bottom, oroi.left:oroi.right);
end
dll_video('rewind',fn1);

len_tslice = max_delta*4+2;
orig_frame = proc_frame - max_delta - 1;

% read in the processed field (or frame)
dll_video('set_rewind',fn2);
dll_video('set_tslice',fn2, 1/fps);
dll_video('discard', fn2, (proc_frame-1)/fps);  

hold = dll_video('tslice', fn2);
proc_y = hold(pvr.top:pvr.bottom, pvr.left:pvr.right);

dll_video('rewind',fn2);
    
