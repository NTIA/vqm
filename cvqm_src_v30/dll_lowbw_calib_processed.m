function [shift, scale, status] = ...
    dll_lowbw_calib_processed(fn, seed_state, num_sec, orig_pixels, ...
    orig_horiz_profile, orig_vert_profile, no_scaling);
% DLL_LOWBW_CALIB_PROCESSED
%   Calculate processed features and low bandwidth calibration.
% SYNTAX
%   [shift, scale, status] = ...
%       dll_lowbw_calib_processed(fn, seed_state, num_sec, orig_pixels, ...
%       orig_horiz_profile, orig_vert_profile, no_scaling);
% DESCRIPTION
%  Calcualte processed features needed for low bandwidth calibration, then
%  estimate spatial registration and scaling registration for each processed 
%  clip.  Video clip must be temporally registered first.
%  'fn' is the file identifier from dll_video, which should always be fn=2.
%  'seed_state' is as returned by dll_lowbw_initialize, which should be
%  called anew each time this function is called.
%  'num_sec' is the number of seconds of video from file fn=1 and 2 that should
%  be used for calibration.  
%  'orig_pixels','orig_horiz_profile', and 'orig_vert_profile' are results
%  from dll_lowbw_calib_original, for the original video (i.e., fn=1).
%   'no_scaling' is 1 to pre sume no spatial scaling,
%   'no_scaling' is 0 to calculate spatial scaling
% 
%   status.error of 1 indicates error,
%   status.scale of 1 indicates scale returned is equal to search limit,
%   status.shift of 1 indicates shift returned is equal to search limit,

num_sec = floor(num_sec);

% uncompress input arguments.
seed_state = double(seed_state);
% % orig_pixels = single(orig_pixels);
% % orig_horiz_profile = double( orig_horiz_profile ) / 257;
% % orig_vert_profile = double( orig_vert_profile ) / 257;

%
status.error = 1;
status.scale = 0;
status.shift = 0;
status.luminance = 0;
shift.horizontal = 0;
shift.vertical = 0;
scale.horizontal = 1000;
scale.vertical = 1000;


% try
    % set up constants.
    max_scale = 100; % maximum scaling search, 10%


    [rows,cols, fps] = dll_video('size', fn); 
    durration = floor( dll_calib_video('total_sec',2) );
    if durration < num_sec,
        num_sec = durration;
    end

    if rows <= 216,
        max_shift_horiz = 4;  % maximum search in any direction, in # pixels
        max_shift_vert = 4;  
        max_scale_horiz = 60; % maximum scaling search, 
        max_scale_vert = 40;  
    elseif rows <= 384,
        max_shift_horiz = 8;
        max_shift_vert = 8;
        max_scale_horiz = 60;  
        max_scale_vert = 40;  
    else
        max_shift_horiz = 20;
        max_shift_vert = 20;
        max_scale_horiz = 100;  
        max_scale_vert = 60;  
    end

    % error checks & corrections
    if mod(max_shift_horiz,2),
        max_shift_horiz = max_shift_horiz + 1;
    end
    if mod(max_shift_vert,2),
        max_shift_vert = max_shift_vert + 1;
    end

    % compute PVR and OROI given the above.
    max_pixels_horiz = max_shift_horiz + (max_scale_horiz / 1000) * cols;
    max_pixels_horiz = ceil(max_pixels_horiz);
    max_pixels_horiz = max_pixels_horiz + mod(max_pixels_horiz,2);
    
    max_pixels_vert = max_shift_vert + (max_scale_vert / 1000) * rows; 
    max_pixels_vert = ceil(max_pixels_vert);
    max_pixels_vert = max_pixels_vert + mod(max_pixels_vert,2);

    [pvr, oroi] = find_pvr_oroi_guess(rows, cols, max_pixels_horiz, max_pixels_vert);

    % run the calibration algorithm
    [shift, scale, status] = sas_core_algorithm( ...
        num_sec, fps, max_shift_horiz, max_shift_vert, max_scale_horiz, ...
        max_scale_vert, max_pixels_horiz, max_pixels_vert, status, pvr, oroi, ...
        seed_state, orig_pixels, orig_horiz_profile, orig_vert_profile, no_scaling, fn);
       
    % check for results next to the search limit.
    if shift.horizontal == max_shift_horiz | shift.vertical == max_shift_vert,
        status.shift = status.shift + 1;
    end
    if abs(scale.horizontal - 1000) == max_scale_horiz | abs(scale.vertical - 1000) == max_scale_horiz,
        status.scale = status.scale + 1;
    end
    
    status.error = 0;

% catch
%     status.error = 1;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Algorithm:
%
%   Use one frame every second (approx).  
%       horizontal and vertical profiles AND
%       80% as many randomly subsampled pixels as there are profile pixels
%           - randomly distributed over all frames
%   Search over ALL frames simultaneously.
%   Search original +- 0 second (yes! ZERO); 
%   shift processed +- 20 pixels/lines for NTSC (as specified else)
%   scale processed by +/- 10% for NTSC (as specified else)
%       rescale using nearest neighbor
%
%   When have final scale & shift, compute luminance gain & offset with
%   those pixels & profiles, too.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_spatial, new_scale, status] = sas_core_algorithm( ...
    num_sec, fps, max_shift_horiz, max_shift_vert, max_scale_horiz, max_scale_vert, ...
    max_pixels_horiz, max_pixels_vert, status, pvr, oroi, ...
    seed_state, orig_pixels, orig_horiz_profile, orig_vert_profile, no_scaling, fn);
% Compute spatial registration for one clip of fields..



new_spatial.horizontal = 0;
new_spatial.vertical = 0;
new_scale.horizontal = 1000;
new_scale.vertical = 1000;


% Choose the % of pixels to be used.
rows = oroi.bottom - oroi.top + 1;
cols = oroi.right - oroi.left + 1;
[list_row, list_col, list_time, list_o] = ...
    sas_choose_pixels (seed_state, rows, cols, num_sec, max_pixels_horiz, max_pixels_vert);


% load frames
dll_video('set_rewind', fn);
for loop = 1:num_sec,
    y(:,:,loop) = dll_calib_video('sec', fn, 1.0/fps); 
    if loop ~= num_sec,
        dll_video('discard', fn, (fps-1)/fps); 
    end
end
dll_video('rewind', fn);

y = double(y);

% cut out PVR
proc = y(pvr.top:pvr.bottom, pvr.left:pvr.right, :);

% Compute processed profiles.
[proc_horiz_profile, proc_vert_profile] = sas_profile_images(proc);

% reshape rows & columns into one dimension.
[rowp,colp,timep] = size(proc);
proc_y = reshape(proc, rowp*colp*timep,1);
clear proc;

%list of coordinates for profiles
list_horiz_profile = (1:cols) + max_pixels_horiz;
list_vert_profile = (1:rows) + max_pixels_vert;


if no_scaling,
    max_scale_vert = 0;
    max_scale_horiz = 0;
end

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
    curr_list_row = list_row * 1000 / (scale_vert+1000) + (rowp/2 - (1000/(scale_vert+1000)) * rowp/2);
    curr_list_col = list_col * 1000 / (scale_horiz+1000) + (colp/2 - (1000/(scale_horiz+1000)) * colp/2);
    
    % Second, shift.  Round to nearest pixel (use floor of +0.5 for speed)
    curr_list_row = floor(curr_list_row + shift_vert + 0.5);
    curr_list_col = floor(curr_list_col + shift_horiz + 0.5);
    list_p = (list_time-1)*rowp*colp + (curr_list_col-1)*rowp + curr_list_row;

    % compute the difference value of the random pixels.
    pixel_list = orig_pixels - proc_y(list_p);
    
    % scale the profiles.
    curr_list_row = list_vert_profile * 1000 / ...
        (scale_vert+1000) + (rowp/2 - (1000/(scale_vert+1000)) * rowp/2);
    curr_list_col = list_horiz_profile * 1000 / ...
        (scale_horiz+1000) + (colp/2-(1000/(scale_horiz+1000)) * colp/2);

    % Second, shift the profiles.  Round to nearest pixel (use floor of +0.5 for speed)
    curr_list_row = floor(curr_list_row + shift_vert + 0.5);
    curr_list_col = floor(curr_list_col + shift_horiz + 0.5);

    % compute the difference value of the profiles.  Append that to random
    % pixels' differneces.
    vert_diff = orig_vert_profile - proc_vert_profile(curr_list_row,:);
    horiz_diff = orig_horiz_profile - proc_horiz_profile(curr_list_col,:);
    pixel_list = [ pixel_list; reshape(vert_diff,rows*num_sec,1); reshape(horiz_diff,cols*num_sec,1)];

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

    loop = loop + 1;
end

scale_horiz = best_scale_horiz;
scale_vert = best_scale_vert;
shift_horiz = best_shift_horiz;
shift_vert = best_shift_vert;


% send back all results.
new_spatial.horizontal = best_shift_horiz;
new_spatial.vertical = best_shift_vert;
new_scale.horizontal = best_scale_horiz+1000;
new_scale.vertical = best_scale_vert+1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pvr] = find_pvr_guess(rows, cols);
% return the best guess for pvr based only on image size.
% i.e., discard overscan.
% 'fchoice' is 1 for 'field' or 0 for 'frame' depending on which is desired.

if (rows == 486 | rows == 480 ) ...
        & cols == 720,
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
function [pvr, oroi] = find_pvr_oroi_guess(rows, cols, max_pixels_horiz, max_pixels_vert);
% Find best guess at PVR and OROI, given image size and
% size of maximum search ('max_pixels_horiz'), in pixels.

pvr = find_pvr_guess(rows, cols);
oroi.top = pvr.top + max_pixels_vert;
oroi.bottom = pvr.bottom - max_pixels_vert;
oroi.left = pvr.left + max_pixels_horiz;
oroi.right = pvr.right - max_pixels_horiz;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [list_row, list_col, list_time, list_o] = ...
    sas_choose_pixels (seed_state, rows, cols, time, max_pixels_horiz, max_pixels_vert);

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


