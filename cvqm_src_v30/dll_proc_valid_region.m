function [pvr] = dll_proc_valid_region(ovr);
% DLL_PROC_VALID_REGION
%   Stand-alone DLL code.  Calculates PVR - processed valid region
% SYNTAX
%   [pvr] = dll_proc_valid_region(ovr);
% DESCRIPTION
%   Compute processed valid region.  Use video from dll_video(fn=1).
%   'ovr' is returned by function dll_orig_valid_region

% fetch control variables.
[rows,cols, fps] = dll_video('size', 2); 
[frames] = dll_video('total_frames',2);

% compute control variables.
half_sec_frames = floor(round(fps) / 2);
half_sec_skip = (half_sec_frames - 1) / fps;
curr = 1;
image_size.rows = rows;
image_size.cols = cols;

% initialize using OVR and spatial registration results.
[curr_valid_region, junk, standard] = valid_region_initialize(rows,cols);
max_valid_region = ovr;

% set SROI for processed video read
dll_calib_video('max_roi');

% set rewind point
dll_video('set_rewind', 2);
% loop through frames, improving valid region estimate.
for cnt = 1:half_sec_frames:(frames - half_sec_frames),
    y = dll_calib_video('sec', 2, 1/fps);
    dll_video('discard', 2, half_sec_skip);
    curr = curr + 1;
    [curr_valid_region] = vr_search (max_valid_region, curr_valid_region, y, standard, image_size);
end
% rewind
dll_video('rewind', 2);

if standard,
    % go in by safety margin
    curr_valid_region.top = curr_valid_region.top + 1;
    curr_valid_region.bottom = curr_valid_region.bottom - 1;
    curr_valid_region.left = curr_valid_region.left + 5;
    curr_valid_region.right = curr_valid_region.right - 5;
end

% odd top/left coordinates, even bottom/right coordinates
curr_valid_region.top = curr_valid_region.top + (1-mod(curr_valid_region.top,2));
curr_valid_region.left = curr_valid_region.left + (1-mod(curr_valid_region.left,2));
curr_valid_region.bottom = curr_valid_region.bottom - mod(curr_valid_region.bottom,2);
curr_valid_region.right = curr_valid_region.right - mod(curr_valid_region.right,2);


% print result if debugging
%fprintf('VR = (%d,%d) (%d,%d)\n', curr_valid_region.top, curr_valid_region.left, ...
%  curr_valid_region.bottom, curr_valid_region.right);

% error check.  override curr_valid_region if results were too small.
if curr_valid_region.bottom - curr_valid_region.top < ...
        (max_valid_region.bottom - max_valid_region.top)/2 || ...
        curr_valid_region.right - curr_valid_region.left < ...
        (max_valid_region.right - max_valid_region.left)/2
    curr_valid_region = max_valid_region;
end

pvr = curr_valid_region;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [curr_valid_region, max_valid_region, standard] = valid_region_initialize(rows, cols);
% initialize two variables, given the image size.

if rows == 486 & cols == 720, 
	% initialize maximum valid region.
	max_valid_region.top = 7;
	max_valid_region.left = 7;
	max_valid_region.bottom = rows - 4;
	max_valid_region.right = cols - 6;
    standard = 1;
elseif rows == 480 & cols == 720, 
	% initialize maximum valid region.
	max_valid_region.top = 7;
	max_valid_region.left = 7;
	max_valid_region.bottom = rows - 2;
	max_valid_region.right = cols - 6;
    standard = 1;
elseif rows == 576 & cols == 720, 
	% initialize maximum valid region.
	max_valid_region.top = 7;
	max_valid_region.left = 17;
	max_valid_region.bottom = rows - 6;
	max_valid_region.right = cols - 16;
    standard = 1;
elseif rows == 720 & cols == 1280, 
	% initialize maximum valid region.
	max_valid_region.top = 7;
	max_valid_region.left = 17;
	max_valid_region.bottom = rows - 6;
	max_valid_region.right = cols - 16;
    standard = 1;
elseif rows == 1080 & cols == 1920, 
	% initialize maximum valid region.
	max_valid_region.top = 7;
	max_valid_region.left = 17;
	max_valid_region.bottom = rows - 6;
	max_valid_region.right = cols - 16;
    standard = 1;
else
    max_valid_region.top = 1;
	max_valid_region.left = 1;
	max_valid_region.bottom = rows;
	max_valid_region.right = cols;
    standard = 0;
end

% initialize current valid region.
curr_valid_region.top = rows/2-1;
curr_valid_region.left = cols/2-1;
curr_valid_region.bottom = rows/2+1;
curr_valid_region.right = cols/2+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_curr_valid_region] = vr_search_standard (max_valid_region, curr_valid_region, y);
% search bounderies for one image.

% search for left side.
locn = max_valid_region.left + 1;
prev = mean(y(:,locn - 1));
while locn < curr_valid_region.left,
    next = mean(y(:,locn));
    if next < 20 | next - 2 > prev,
        locn = locn + 1;
        prev = next;
    else
        break;
    end
end
curr_valid_region.left = locn;

% search for top side.
locn = max_valid_region.top + 1;
prev = mean(y(locn - 1,:));
while locn < curr_valid_region.top,
    next = mean(y(locn,:));
    if next < 20 | next - 2 > prev,
        locn = locn + 1;
        prev = next;
    else
        break;
    end
end
curr_valid_region.top = locn;

% search for right side.
locn = max_valid_region.right - 1;
prev = mean(y(:,locn + 1));
while locn > curr_valid_region.right,
    next = mean(y(:,locn));
    if next < 20 | next - 2 > prev,
        locn = locn - 1;
        prev = next;
    else
        break;
    end
end
curr_valid_region.right = locn;

% search for bottom side.
locn = max_valid_region.bottom - 1;
prev = mean(y(locn + 1,:));
while locn > curr_valid_region.bottom,
    next = mean(y(locn,:));
    if next < 20 | next - 2 > prev,
        locn = locn - 1;
        prev = next;
    else
        break;
    end
end
curr_valid_region.bottom = locn;

% return updated CVR
new_curr_valid_region = curr_valid_region;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_curr_vr] = vr_search_noborder (max_vr, curr_vr, y, image_size);
% search bounderies for one image.

% search bounderies for one image.
% max_vr MUST BE exactly equal to the image size.  This
% algorihtm is intended for CIF, QCIF, VGA, and other video where the
% entire image is displayed.

% don't discard more than 4% of the rows or columns on any one border.
max_discard_rows = ceil(image_size.rows * 0.04);
max_discard_cols = ceil(image_size.cols * 0.04);

% search for left side.  Allow left side not to move in, even by one.
for locn = max_vr.left:max_discard_cols,
    if mean(y(:,locn)) < 20 | mean(y(:,locn)) + 20 < mean(y(:,locn+1)),
        % is invalid -- still increasing
    else
        break;
    end
end
curr_vr.left = locn;

% search for right side. Allow right side not to move in, even by one.
for locn = max_vr.right:-1:image_size.cols - max_discard_cols + 1,
    if mean(y(:,locn)) < 20 | mean(y(:,locn)) + 20 < mean(y(:,locn-1)),
        % is invalid -- still increasing
    else
        break;
    end
end
curr_vr.right = locn;

% search for top side. Allow top side not to move in, even by one.
for locn = max_vr.top:max_discard_rows,
    if mean(y(locn,:)) < 20 | mean(y(locn,:)) + 20 < mean(y(locn+1,:)),
        % is invalid -- still increasing
    else
        break;
    end
end
curr_vr.top = locn;

% search for bottom side.  Allow the bottom not to move in, even by one.
for locn = max_vr.bottom:-1:image_size.rows - max_discard_rows + 1,
    if mean(y(locn,:)) < 20 | mean(y(locn,:)) + 20 < mean(y(locn-1,:)),
        % is invalid -- still increasing
    else
        break;
    end
end
curr_vr.bottom = locn;

% return updated CVR
new_curr_vr = curr_vr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_curr_vr] = vr_search (max_vr, curr_vr, y, standard, image_size);
% search bounderies for one image.

if standard,
    [new_curr_vr] = vr_search_standard (max_vr, curr_vr, y);
else
    [new_curr_vr] = vr_search_noborder (max_vr, curr_vr, y, image_size);
end
