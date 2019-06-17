function [ti2, ti10, ymean, is_white_clip, is_black_clip] = ...
    dll_lowbw_temporal_features(fn, durration, pvr);
% DLL_LOWBW_TEMPORAL_FEATURES
%   Step 1: Compute features needed for temporal registration
% SYNTAX
%   [ti2, ti10, ymean, is_white_clip, is_black_clip] = ...
%       dll_lowbw_temporal_features(fn, durration);
%   [...] = dll_lowbw_temporal_features(fn, durration, pvr);
% DESCRIPTION
%   Input arguments are 'fn' from dll_video, to locate the video file;
%   'durration' in seconds, the durration of video to be used (less than
%   file length).  'pvr' computed by dll_proc_valid_region, if known, should
%   also be specified.  A default PVR will be used if this information is not
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

% initialize
[fps] = dll_video('fps', fn); 
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

rows = pvr.bottom-pvr.top+1;
cols = pvr.right-pvr.left+1;

% loop through frames
dll_video('set_rewind', fn);
dll_video('set_tslice', fn, 1.0/fps);
dll_calib_video('sroi',pvr,0);
if ~progressive,
    buffer = zeros(rows/2*cols,6,2);
else
    buffer = zeros(rows*cols,6,2);
end

ti2_cnt = 1;
ti10_cnt = 1;
ymean_cnt = 1;
white_cnt = 1;
curr = 1;
for loop = 1:floor(fps*durration),
    % read frame
    y_frames = dll_calib_video('tslice', fn);
    
    if ~progressive,
        % Reshape frames into fields.
        y_fields = reshape(y_frames,2,rows/2,cols);
        y_fields = permute(y_fields,[2 1 3]);
        frows = rows / 2;
        buffer(:,curr,1) = reshape(y_fields(:,fld_num(1),:),rows/2*cols,1);
        buffer(:,curr,2) = reshape(y_fields(:,fld_num(2),:),rows/2*cols,1);
    else
        % Change name; don't reshape frames into fields.
        [rows,cols,time] = size(y_frames);
        y_fields = reshape(y_frames,rows,1,cols);
        frows = rows;
        buffer(:,curr,1) = reshape(y_fields(:,fld_num(1),:),rows*cols,1);
    end
    clear y_frames;

    % compute feature TI2
    for fld = 1:length(fld_num),
        if loop >= 2,
            other = curr-1;
            if other < 1,
                other = other + 6;
            end
            ti2(ti2_cnt) = sqrt(mean( (buffer(:,curr,fld) - buffer(:,other,fld)).^2 ));
            ti2_cnt = ti2_cnt + 1;
        end
    end

    % compute feature TI10
    for fld = 1:length(fld_num),
        if loop >= 6,
            other = curr-5;
            if other < 1,
                other = other + 6;
            end
            ti10(ti10_cnt) = sqrt(mean( (buffer(:,curr,fld) - buffer(:,other,fld)).^2 ));
            ti10_cnt = ti10_cnt + 1;
        end
    end

    % compute feature Ymean
    for fld = 1:length(fld_num),
        ymean(ymean_cnt) = mean(buffer(:,curr,fld));
        ymean_cnt = ymean_cnt + 1;
    end

    % white & black clipping
    for fld = 1:length(fld_num),
        list_white(white_cnt) = length( find(buffer(:,curr,fld) >= 254) );
        list_black(white_cnt) = length( find(buffer(:,curr,fld) <= 1) );
        white_cnt = white_cnt + 1;
    end
    
    curr = curr + 1;
    if curr > 6,
        curr = 1;
    end

end
dll_video('rewind', fn);



list_white = list_white ./ (frows * cols) * 100;
is_white_clip = max(max(list_white)) >= 5.0;

list_black = list_black ./ (frows * cols) * 100;
is_black_clip = max(max(list_black)) >= 5.0;



ti2 = single(ti2);
ti10 = single(ti10);
ymean = single(ymean);

