function [orig_pixels, orig_horiz_profile, orig_vert_profile] = ...
    dll_lowbw_calib_original(fn, seed_state, num_sec)
% DLL_LOWBW_CALIB_ORIGINAL
%   Calcualte original features needed for low bandwidth calibration.
% SYNTAX
%  [orig_pixels, orig_horiz_profile, orig_vert_profile] = ...
%       fast_calibration(fn, seed_state, num_sec)
% DESCRIPTION
%  Calcualte original features needed for low bandwidth calibration:
%  estimate spatial registration and scaling registration for each processed 
%  clip.  Video clip must be temporally registered first.
%  'fn' is the file identifier from dll_video, which should always be fn=1.
%  'seed_state' is as returned by dll_lowbw_initialize, which should be
%  called anew each time this function is called.
%  'num_sec' is the number of seconds of video from file fn=1 that should
%  be used for calibration.  
%
%  Return values 'orig_pixels', 'orig_horiz_profile', and
%  'orig_vert_profile' are required by function dll_lowbw_calib_processed.


num_sec = floor(num_sec);

% set up constants.
max_scale = 100; % maximum scaling search, 10%

[rows,cols, fps] = dll_video('size', fn);
durration = floor( dll_video('total_frames',1) / fps );
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



% Choose the % of pixels to be used.
rows = oroi.bottom - oroi.top + 1;
cols = oroi.right - oroi.left + 1;
[list_row, list_col, list_time, list_o] = ...
    sas_choose_pixels (seed_state, rows, cols, num_sec, max_pixels_horiz, max_pixels_vert);

% load frames
dll_video('set_rewind', fn);
for loop = 1:num_sec,
    y(:,:,loop) = dll_video('sec', fn, 0, 1.0/fps); % don't reframe
    if loop ~= num_sec,
        dll_video('discard', fn, (fps-1)/fps); 
    end
end
dll_video('rewind', fn);

y = double(y);

% cut out OROI
orig = y(oroi.top:oroi.bottom, oroi.left:oroi.right, :);

% Compute original profiles.
[orig_horiz_profile, orig_vert_profile] = sas_profile_images(orig);

% reshape rows & columns into one dimension.
[rows,cols,time] = size(orig);
orig_y = reshape(orig, rows*cols*time,1);
clear orig;

%list of coordinates for profiles
list_horiz_profile = (1:cols) + max_pixels_horiz;
list_vert_profile = (1:rows) + max_pixels_vert;

% pick our original pixels
orig_pixels = orig_y(list_o);


% % % Limit precision on return variables
% % orig_pixels = char(orig_pixels);
% % orig_horiz_profile = uint16( round(orig_horiz_profile * 257));
% % orig_vert_profile = uint16( round(orig_vert_profile * 257));

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

rand('state', double(seed_state));
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


