function [new_clip_structs, status, unfiltered_clip_structs] = valid_region_lowbw(test_structs, clip_structs, varargin)
% valid_region_lowbw
%   Calculate valid region for each video clip.  2005 low-banwdith algorithm.
% SYNTAX
%  [new_clip_structs, status, unfiltered_clip_structs] = valid_region_lowbw(test_structs, clip_structs);
%  [...] = valid_region_lowbw(...,'PropertyName',PropertyValue,...)
% DEFINITION
%   Calculate OVR and PVR for all clips.  Spatial registration must already
%   be correctly set for each processed clip.  Return variable
%   'unfiltered_clip_structs' contains OVR and PVR without any adjustments
%   (no scene filtering, and no coordinate odd/even correction)
%
%   Algorithm improved from 2004 standard, with additional support for QCIF,
%   CIF, and VGA video streams.
%
% Optional properties:
%   'Frequency',value,      Examine frames at a frequency of 'value',
%                           specified in seconds.  By default, this
%                           will be in seconds.
%   'unaligned',            Used temporally unaligned video sequences.
%                           This is the default behavior.
%   'aligned',              Used temporally ALIGNED video sequences.
%                           Skip unaligned clips.
%   'verbose',              Print results
%   'quiet',                Don't print results
%
%  Return clip structure, modified with new Common Valid Regions.  Also
%  return status, 0 if all works correctly and 1 if an error has occured.

unfiltered_clip_structs = clip_structs;

try
    status = 0;
    % handle optional properties.
    frequency = 0.5;
    verbose = 0;
    flag = 'unaligned';
    cnt = 1;
    while cnt <= nargin-2,
        if strcmp(lower(varargin{cnt}),'frequency'),
            frequency = varargin{cnt+1};
            cnt = cnt + 2;
        elseif strcmp(lower(varargin{cnt}),'verbose'),
            verbose = 1;
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'quiet'),
            verbose = 0;
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'unaligned'),
            flag = 'unaligned';
            cnt = cnt + 1;
        elseif strcmp(lower(varargin{cnt}),'aligned'),
            flag = 'aligned';
            cnt = cnt + 1;
        else
            error('Property Name not recognized.  Aborting.');
        end
    end
    % erase previous CVR information
    for i=1:length(clip_structs),
        clip_structs(i).cvr.top = nan;
        clip_structs(i).cvr.left = nan;
        clip_structs(i).cvr.bottom = nan;
        clip_structs(i).cvr.right = nan;
    end
    warning('off');
    new_clip_structs = clip_structs;


    % sort clips by scene.
    offsets = sort_clips_by('scene', clip_structs, test_structs);

    for scene = 1:length(offsets),
        curr = offsets{scene};

        if ~strcmp('original',clip_structs(curr(1)).hrc{1}),
            error('Original clip missing for test %s scene %s', clip_structs(curr(1)).test{1}, clip_structs(curr(1)).scene{1});
        end

        % calculate OVR for the original scene, examining one frame every
        % frequency seconds.
        [curr_vr, max_vr, standard] = vrp_initialize(clip_structs(curr(1)).image_size);
        clip_structs(curr(1)).cvr = max_vr;
        tst = find_test(test_structs, clip_structs(curr(1)).test{1});
        add_frames = round(clip_structs(curr(1)).fps * frequency);

        frames = total_tslices(clip_structs(curr(1)),1.0 / clip_structs(curr(1)).fps, flag);
        for loop = 1:add_frames:frames,
            y = read_tslice(test_structs(tst), clip_structs(curr(1)), 1.0 / clip_structs(curr(1)).fps, loop, 'full',flag);
            curr_vr = vrp_search (max_vr, curr_vr, y, standard, clip_structs(curr(1)).image_size);
        end

        % max valid region for all processec clips is the original's valid
        % region.
        max_vr = curr_vr;
        unfiltered_clip_structs(curr(1)).cvr = curr_vr;

        % fprintf('OVR is (%d,%d) (%d,%d)\n', curr_vr.top, curr_vr.left, curr_vr.bottom, curr_vr.right); 

        % now, search through all processed clips.
        for hrc = 2:length(curr),
            [curr_vr] = vrp_initialize(clip_structs(curr(1)).image_size);
            clip_structs(curr(hrc)).cvr = vrp_cvr (max_vr,...
                clip_structs(curr(hrc)).spatial, clip_structs(curr(hrc)).image_size, ...
                clip_structs(curr(hrc)));

            if tst ~= find_test(test_structs,clip_structs(curr(hrc)).test{1}),
                error('sort_clips_by error, clips from different tests grouped together!');
            end

            frames = total_tslices(clip_structs(curr(hrc)),1.0 / clip_structs(curr(hrc)).fps, flag);
            for loop = 1:add_frames:frames,
                y = read_tslice(test_structs(tst), clip_structs(curr(hrc)), 1.0 / clip_structs(curr(hrc)).fps, loop, 'full', flag);
                curr_vr = vrp_search (max_vr, curr_vr, y, standard, clip_structs(curr(1)).image_size);
            end
            pvr(hrc) = curr_vr;
            
            if standard,
                unfiltered_clip_structs(curr(hrc)).cvr.top = curr_vr.top+1;
                unfiltered_clip_structs(curr(hrc)).cvr.left = curr_vr.left+5;
                unfiltered_clip_structs(curr(hrc)).cvr.bottom = curr_vr.bottom-1;
                unfiltered_clip_structs(curr(hrc)).cvr.right = curr_vr.right-5;
            else
                unfiltered_clip_structs(curr(hrc)).cvr.top = curr_vr.top;
                unfiltered_clip_structs(curr(hrc)).cvr.left = curr_vr.left;
                unfiltered_clip_structs(curr(hrc)).cvr.bottom = curr_vr.bottom;
                unfiltered_clip_structs(curr(hrc)).cvr.right = curr_vr.right;
            end
            % fprintf('PVR is (%d,%d) (%d,%d)\n', curr_vr.top, curr_vr.left, curr_vr.bottom, curr_vr.right); 
        end
        
        if length(curr) ==1,
            if verbose,
            	fprintf('%s:%s   WARNING: No Processed Clips Defined!\n', clip_structs(curr(1)).test{1}, clip_structs(curr(1)).scene{1});
            end
            continue;
        end
        
        % combine HRC results.
        pvr = pvr(2:length(pvr));
        cvr.top = max([pvr.top]);
        cvr.left = max([pvr.left]);
        cvr.bottom = min([pvr.bottom]);
        cvr.right = min([pvr.right]);

        if standard,
            % go in by safety margin
            cvr.top = cvr.top + 1;
            cvr.bottom = cvr.bottom - 1;
            cvr.left = cvr.left + 5;
            cvr.right = cvr.right - 5;
        end

        % odd top/left coordinates, even bottom/right coordinates
        cvr.top = cvr.top + (1-mod(cvr.top,2));
        cvr.left = cvr.left + (1-mod(cvr.left,2));
        cvr.bottom = cvr.bottom - mod(cvr.bottom,2);
        cvr.right = cvr.right - mod(cvr.right,2);

        % record new CVR
        for hrc = 1:length(curr),
            new_clip_structs(curr(hrc)).cvr = cvr;
        end

        % print result if requested
        if verbose,
            fprintf('%s:%s  PVR = (%d,%d) (%d,%d)\n', clip_structs(curr(1)).test{1}, clip_structs(curr(1)).scene{1}, ....
                cvr.top, cvr.left, cvr.bottom, cvr.right);
        end
        clear pvr cvr;

    end


    warning('on');
catch
    if verbose,
        fprintf('ERROR: %s\n', lasterr);
    end
    status = 1;
    new_clip_structs = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [curr_vr, max_vr, standard] = vrp_initialize(image_size);
% initialize two variables, given the image size.

if image_size.rows == 486 & image_size.cols == 720, 
	% initialize maximum valid region.
	max_vr.top = 7;
	max_vr.left = 7;
	max_vr.bottom = image_size.rows - 4;
	max_vr.right = image_size.cols - 6;
    standard = 1;
elseif image_size.rows == 480 & image_size.cols == 720, 
	% initialize maximum valid region.
	max_vr.top = 7;
	max_vr.left = 7;
	max_vr.bottom = image_size.rows - 2;
	max_vr.right = image_size.cols - 6;
    standard = 1;
elseif image_size.rows == 576 & image_size.cols == 720, 
	% initialize maximum valid region.
	max_vr.top = 7;
	max_vr.left = 17;
	max_vr.bottom = image_size.rows - 6;
	max_vr.right = image_size.cols - 16;
    standard = 1;
elseif image_size.rows == 720 & image_size.cols == 1280, 
	% initialize maximum valid region.
	max_vr.top = 7;
	max_vr.left = 17;
	max_vr.bottom = image_size.rows - 6;
	max_vr.right = image_size.cols - 16;
    standard = 1;
elseif image_size.rows == 1080 & image_size.cols == 1920, 
	% initialize maximum valid region.
	max_vr.top = 7;
	max_vr.left = 17;
	max_vr.bottom = image_size.rows - 6;
	max_vr.right = image_size.cols - 16;
    standard = 1;
else
    max_vr.top = 1;
	max_vr.left = 1;
	max_vr.bottom = image_size.rows;
	max_vr.right = image_size.cols;
    standard = 0;
end

% initialize current valid region.
curr_vr.top = image_size.rows/2-1;
curr_vr.left = image_size.cols/2-1;
curr_vr.bottom = image_size.rows/2+1;
curr_vr.right = image_size.cols/2+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_curr_vr] = vrp_search (max_vr, curr_vr, ...
    y, standard, image_size);
% search bounderies for one image.

if standard,
    [new_curr_vr] = vrp_search_standard (max_vr, curr_vr, y);
else
    [new_curr_vr] = vrp_search_noborder (max_vr, curr_vr, y, image_size);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_curr_vr] = vrp_search_noborder (max_vr, curr_vr, y, image_size);
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
function [new_curr_vr] = vrp_search_standard (max_vr, curr_vr, y);
% search bounderies for one image.  Standard algorithm

% search for left side.
locn = max_vr.left + 1;
prev = mean(y(:,locn - 1));
while locn < curr_vr.left,
    next = mean(y(:,locn));
    if next < 20 | next - 2 > prev,
        locn = locn + 1;
        prev = next;
    else
        break;
    end
end
curr_vr.left = locn;

% search for top side.
locn = max_vr.top + 1;
prev = mean(y(locn - 1,:));
while locn < curr_vr.top,
    next = mean(y(locn,:));
    if next < 20 | next - 2 > prev,
        locn = locn + 1;
        prev = next;
    else
        break;
    end
end
curr_vr.top = locn;

% search for right side.
locn = max_vr.right - 1;
prev = mean(y(:,locn + 1));
while locn > curr_vr.right,
    next = mean(y(:,locn));
    if next < 20 | next - 2 > prev,
        locn = locn - 1;
        prev = next;
    else
        break;
    end
end
curr_vr.right = locn;

% search for bottom side.
locn = max_vr.bottom - 1;
prev = mean(y(locn + 1,:));
while locn > curr_vr.bottom,
    next = mean(y(locn,:));
    if next < 20 | next - 2 > prev,
        locn = locn - 1;
        prev = next;
    else
        break;
    end
end
curr_vr.bottom = locn;

% return updated CVR
new_curr_vr = curr_vr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cvr] = vrp_cvr (max_vr,spatial_shift, image_size, clip_struct);
% find cvr closest to max_vr that is valid given the spatial
% shift.

cvr = max_vr;

if spatial_shift.horizontal <= 0,
    cvr.left = max(cvr.left, -spatial_shift.horizontal+1);
else
    cvr.right = min(cvr.right, image_size.cols - spatial_shift.horizontal);
end

if spatial_shift.vertical <= 0,
    cvr.top = max(cvr.top, -spatial_shift.vertical+is_reframing_indicated(clip_struct)+1);
else
    cvr.bottom = min(cvr.bottom, image_size.rows - spatial_shift.vertical - is_reframing_indicated(clip_struct));
end








