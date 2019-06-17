function [ti2, ti10, ymean, is_white_clip, is_black_clip] = ...
    dll_lowbw_temporal_features(fn, durration, pvr);
% DLL_LOWBW_TEMPORAL_FEATURES
%   Step 1: Computer features needed for temporal registration
%   original code, high memory requirements
% SYNTAX
%   [ti2, ti10, ymean, is_white_clip, is_black_clip] = ...
%       dll_lowbw_temporal_features(fn, durration);
%   [...] = dll_lowbw_temporal_features(fn, durration, pvr);
% DESCRIPTION
%   Input arguments are 'fn' from dll_video, to locate the video file;
%   'durration' in seconds, the durration of video to be used (less than
%   file length).  'pvr' from dll_proc_valid_region, if known, should also
%   be specified.  A default PVR will be used if this information is not
%   known.  WARNING:  'pvr' (if present) and 'durration' must be identical
%   for original & processed function calls.
%
%   Return values are three features: ti2, ti10 & ymean.
%   Also whether the video appears to contain white level clipping
%   ('is_white_clip') or black level clipping ('is_black_clip').

if ~exist('pvr','var'),
    [row,col] = dll_video('size', fn);  
    
    % discard 4% if PVR not defined.
    hold = round(0.04*row);
    hold = hold + mod(hold,2); % next even number
    pvr.top = hold + 1;
    pvr.bottom = row - hold;

    hold = round(0.04*col);
    hold = hold + mod(hold,2); % next even number
    pvr.left = hold + 1;
    pvr.right = col - hold;
end

% get images.
[fps] = dll_video('fps', fn);  
% durration = floor(durration * fps);
[y_frames] = dll_calib_video('peek', fn, durration);
[video_standard] = dll_video('get_video_standard', fn);

if strcmp(video_standard,'interlace_upper_field_first'), % e.g., if fps == 25,
    fld_num(1) = 1;
    fld_num(2) = 2;
    progressive = 0;
elseif strcmp(video_standard,'interlace_lower_field_first'),
    fld_num(1) = 2;
    fld_num(2) = 1;
    progressive = 0;
elseif strcmp(video_standard,'progressive'),
    fld_num(1) = 1;
    progressive = 1;
else
    warning('video standard not recognized');
    fld_num(1) = 1;
    progressive = 1;
end

y_frames = y_frames(pvr.top:pvr.bottom, pvr.left:pvr.right, :);
if ~progressive,
    % Reshape frames into fields.
    [rows,cols,time] = size(y_frames);
%     y_fields = reshape(y_frames,rows/2,2,cols,time);
    y_fields = reshape(y_frames,2,rows/2,cols,time);
    y_fields = permute(y_fields,[2 1 3 4]);
    rows = rows / 2;
else
    % Change name; don't reshape frames into fields.
    [rows,cols,time] = size(y_frames);
    y_fields = reshape(y_frames,rows,1,cols,time);
end
clear y_frames;

% Compute TI2 -- TI on fields, two fields / one frame appart.
cnt = 1;
for loop = 2:time,
    for fld = 1:length(fld_num),
        ti2(cnt) = sqrt(mean(reshape(double(y_fields(:,fld_num(fld),:,loop)) - double(y_fields(:,fld_num(fld),:,loop-1)), 1, rows*cols).^2));
        cnt = cnt + 1;
    end
end

% Compute TI10 -- TI on fields, ten fields / five frames appart.
cnt = 1;
for loop = 6:time,
    for fld = 1:length(fld_num),
        ti10(cnt) = sqrt(mean(reshape(double(y_fields(:,fld_num(fld),:,loop)) - double(y_fields(:,fld_num(fld),:,loop-5)), 1, rows*cols).^2));
        cnt = cnt + 1;
    end
end

% compute white & black clipping
[rows2,timeA,cols2,timeB] = size(y_fields);
total_pixels = rows2 * cols2;
for loop = 1:time,
    for fld = 1:length(fld_num),
        list_white(loop,fld) = length( find(y_fields(:,fld_num(fld),:,loop) >= 254));
        list_black(loop,fld) = length( find(y_fields(:,fld_num(fld),:,loop) <= 1));
    end
end

list_white = list_white ./ total_pixels * 100;
is_white_clip = max(max(list_white)) >= 5.0;

list_black = list_black ./ total_pixels * 100;
is_black_clip = max(max(list_black)) >= 5.0;

% Compute ymean -- average fields luminance.
cnt = 1;
for loop = 5:time,
    for fld = 1:length(fld_num),
        ymean(cnt) = mean(mean(double(y_fields(:,fld_num(fld),:,loop))));
        cnt = cnt + 1;
    end
end


ti2 = single(ti2);
ti10 = single(ti10);
ymean = single(ymean);

