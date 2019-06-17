function [clips_out] = fix_temporal (test_struct, clips_in, varargin);
% FIX_TEMPORAL
%   Adjust / fix temporal registration for a clipset.
% SYNTAX
%   [clips_out] = fix_temporal (test_struct, clips_in);
%   [..] = fix_temporal (...,'PropertyName',PropertyValue, ...);
% DESCRIPTION
%   'clips_in' is a clip structure, formatted like GClips.  Modify the
%   temporal registration and requested and return in 'clips_out'.
%   'test_struct' is the test structure, like GTests.
%   This function prints no information to the screen.
%
%   Optional properties (one or more must be chosen):
%  
%   'StartPoint'  Keep start frame alignment.  Replace starting frames,
%                (moving as little as possible) but moving inward as needed
%                so that the first frame always 1+.
%   'EndPoint'   Keep starting frame alignment; replace ending frame
%                alignment.  Set the stop frame alignment to include
%                as much of the video sequence as possible; with an
%                identical length for all HRCs from one SRC.
%   'FirstFrame' Alignment is unknown or incorrect.  Set all clips to
%                first frames align, then make sure all video sequences
%                from each video sequence have the same length.
%   'Spatial'    Spatial registration has changed.  Shift end of aligned 
%                video by 1 frame where needed.
%   'Start', N,  Move all starting alignment point forward 'N' frames.
%                Negative numbers will move alignments backward.
%                Warning:  use with care, no error checks.
%   'Stop', N,   Move all ending alignment point forward by 'N' frames.
%                Negative numbers will move alignment backward.
%                Warning:  use with care, no error checks.
%   'Original', Start, Stop,
%                Keep relative alignments between each original & all
%                processed versions thereof.  Move original's align_start
%                and align_stop to 'Start' and 'Stop' arguments,
%                respectively.  Move all processed versions accordingly.
%                Will throw an error if this cannot be done.
%                When chosen, all other options ignored.
%   'HalfOkay',  If present, align_start values for processed clips that
%                contain one-half (0.5) will be okay (e.g., 1.5, -3.5).
%                These are returned by temporal_registration_lowbw when
%                called with the 'field' option, to indicate reframing is
%                likely.

if length(varargin) == 0,
    error('choose one or more property string for function fix_temporal.');
end

start_point = 0;
end_point = 0;
first_frame = 0;
half_okay = 0;


is_spatial = 0;
move_start = 0;
move_stop = 0;

original_range = 0;


cnt = 1;
while cnt <= length(varargin),
    if strcmpi(varargin{cnt},'StartPoint'),
        start_point = 1;
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'EndPoint'),
        end_point = 1;
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'FirstFrame'),
        first_frame = 1;
        is_spatial = 1;
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'Spatial'),
        is_spatial = 1;
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'Start'),
        move_start = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin{cnt},'Stop'),
        move_stop = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin{cnt},'original'),
        original_range = 1;
        orig_start = varargin{cnt+1};
        orig_stop = varargin{cnt+2};
        cnt = cnt + 3;
    elseif strcmpi(varargin{cnt},'HalfOkay'),
        half_okay = 1;
        cnt = cnt + 1;
    else
        error('Property Name not recognized');
    end
end

clips_out = clips_in;

for cnt=1:length(clips_in),
    if mod(clips_in(cnt).align_start,1) ~= 0,
        if mod(clips_in(cnt).align_start,1) == 0.5 & half_okay & ...
                ~strcmp(clips_in(cnt).hrc{1},'original'),
            ; % 0.5 fraction okay because HalfOkay flag specified
        elseif first_frame,
            ; % don't worry about HalfOkay, alignment to be discarded.
        else
            error('Invalid clip; align_start must be an integer');
        end
    end
end


% user specified range for all originals.  Keep alignments.
if original_range,
    list = sort_clips_by('scene', clips_out, test_struct);
    for i=1:length(list),
        if ~strcmpi('original',clips_out(list{i}(1)).hrc{1}),
            error('No original defined for one or more scenes; operation unavailable');
        end
        move_start = orig_start - clips_out(list{i}(1)).align_start;
        move_stop = orig_stop - clips_out(list{i}(1)).align_stop;
        for j=list{i},
            clips_out(j).align_start = clips_out(j).align_start + move_start; 
            clips_out(j).align_stop = clips_out(j).align_stop + move_stop; 
        end
    end
    % check
    for i=1:length(clips_out),
        if clips_out(i).align_start < clips_out(i).loc_start | clips_out(i).align_stop > clips_out(i).loc_stop,
            error('Insufficient frames available to implement ''fix_temporal'' request');
        end
    end
    
    return;
end


% Set default alignment of "first frames align" for ALL clips.
if first_frame,
    for cnt = 1:length(clips_out),
        clips_out(cnt).align_start = clips_out(cnt).loc_start;
    end
    
    % also, fix end points.
    end_point = 1;
end


if start_point,
    order = sort_clips_by('scene', clips_out, test_struct);
    for cnt = 1:length(order),
        hold = [clips_out(order{cnt})];
        hold = min([hold.align_start]);
        if hold < 1,
            hold = 1-hold;
            for loop = order{cnt},
                clips_out(loop).align_start = clips_out(loop).align_start + hold;
            end
        end
    end
end


% Set end points to last available frame, where all clips from one scene
% have the same length
if end_point,
    order = sort_clips_by('scene', clips_out, test_struct);
    % loop through each video scene
    for cnt = 1:length(order),
        curr = order{cnt};
        % find the minimum scene length
        list = [];
        for loop = 1:length(curr),
            list = [list (clips_out(curr(loop)).loc_stop - clips_out(curr(loop)).align_start)];
        end
        total = min(list);
        % set align_start to first frames; set align_stop for minimum
        % length.
        for loop = 1:length(curr),
            clips_out(curr(loop)).align_stop = clips_out(curr(loop)).align_start + total;
        end
    end
    
    % also, fix spatial registration off-by-one temporal registration.
    is_spatial = 1;
end


% make sure odd vertical shift clips have an extra frame.
if is_spatial,
    order = sort_clips_by('scene', clips_out, test_struct);
    % loop through each video scene
    for cnt = 1:length(order),
        curr = order{cnt};
        
        % error checks.
        if ~strcmpi(clips_out(curr(1)).hrc{1},'original'),
            continue;
        end
        if isnan(clips_out(curr(1)).align_stop) | isnan(clips_out(curr(1)).align_start),
            continue;
        end
        if strcmp(clips_in(curr(1)).video_standard,'progressive'),
            continue;
        end

        % loop through each HRC.
        len_orig = clips_out(curr(1)).align_stop - clips_out(curr(1)).align_start;
        move_all_back = 0;
        for cnt = 2:length(curr),
            % error check
            if isnan(clips_out(curr(cnt)).align_stop) | isnan(clips_out(curr(cnt)).align_start),
                continue;
            end
            % compute length
            len_proc = clips_out(curr(cnt)).align_stop - clips_out(curr(cnt)).align_start;
            if mod( clips_out(curr(cnt)).spatial.vertical,2),
                % need one extra
                if len_proc == len_orig,
                    clips_out(curr(cnt)).align_stop = clips_out(curr(cnt)).align_stop + 1;
                    if clips_out(curr(cnt)).align_stop > clips_out(curr(cnt)).loc_stop,
                        move_all_back = 1;
                    end
                end
            else
                % need same frames
                if len_proc + 1 == len_orig,
                    clips_out(curr(cnt)).align_stop = clips_out(curr(cnt)).align_stop - 1;
                end
            end
        end
        % move all back, because went past the last available frame.
        if move_all_back,
            for cnt1 = 1:length(curr),
                clips_out(curr(cnt1)).align_stop = clips_out(curr(cnt1)).align_stop - 1;
            end
        end
    end
    
end

if move_start ~= 0,
    for cnt = 1:length(clips_out),
        clips_out(cnt).align_start = clips_out(cnt).align_start + move_start;
    end
end

if move_stop ~= 0,
    for cnt = 1:length(clips_out),
        clips_out(cnt).align_stop = clips_out(cnt).align_stop + move_stop;
    end
end


