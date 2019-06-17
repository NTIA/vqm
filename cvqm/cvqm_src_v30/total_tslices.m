function [total] = total_tslices(clip_struct, tslice_length_sec, varargin)
% TOTAL_TSLICES  
%  Compute the total number of time-slices available for a video clip,
%  given a time-slice length in seconds.
% SYNTAX
%  [total] = total_tslices(clip_struct, tslice_length_sec);
%  [total] = total_tslices(..., 'ProptertyName',PropertyValue, ...);
% DESCRIPTION
%  Given a time-slice length in seconds ('tslice_length_sec'), compute the total
%  number of time-slices available for this clip.  Description of clip passed
%  in 'clip_struct' (of the same format as GClips).  
%  
%  Optional properties include:
%  'align_start',value      Specified 'value' overrides clip_struct.align_start
%  'align_stop', value      Specified 'value' overrides clip_struct.align_stop
%  'unaligned'                Count timeslices based on total range of frames
%                                 (loc_start & loc_stop) instead of aligned frames.
%                                 Cannot be used with 'align_start' & 'align_stop'.
%  'aligned'                   Count timeslices based on aligned segment.
%                                 This is the default behavior.
%  'TIS', sec,              See read_tslice for more information.

is_start = clip_struct.align_start;
is_stop = clip_struct.align_stop;
is_tis = 0;
tis_skip = 0;
tis_frames = 0;

cnt = 1;
while cnt <= nargin - 2,
    if strcmp(lower(varargin{cnt}),'align_start'),
        is_start = varargin{cnt+1};
        cnt = cnt + 2;
        if is_start < clip_struct.loc_start,
            error('When override ''align_start'' in ''total_tslices'', must adhere to constraints of clip''s loc_start');
        end
    elseif strcmp(lower(varargin{cnt}),'align_stop'),
        is_stop = varargin{cnt+1};
        cnt = cnt + 2;
        if is_stop > clip_struct.loc_stop,
            error('When override ''align_stop'' in ''total_tslices'', must adhere to constraints of clip''s loc_stop');
        end
    elseif strcmp(lower(varargin{cnt}),'unaligned'),
        is_start = clip_struct.loc_start;
        is_stop = clip_struct.loc_stop;
        cnt = cnt + 1;
    elseif strcmp(lower(varargin{cnt}),'aligned'),
        is_start = clip_struct.align_start;
        is_stop = clip_struct.align_stop;
        cnt = cnt + 1;
    elseif strcmp(upper(varargin(cnt)),'TIS') == 1,
        is_tis = 1;
        tis_sec = varargin{cnt+1};
        cnt = cnt + 2;
        
        % compute extra & skip frames 
        [tis_frames] = tslice_conversion(tis_sec, clip_struct.fps);
        [tslice_frames] = tslice_conversion(tslice_length_sec, clip_struct.fps);
        tis_skip = - mod(tis_frames, tslice_frames);
        if tis_skip < 0,
            tis_skip = tis_skip + tslice_frames;
        end
    else
        error('Property value passed into total_tslices not recognized');
    end
end

% return immediately if given an undefined range.
if isnan(is_start) | isnan(is_stop),
    fprintf('Aligned segment undefined for this clip.\n');
    total = 0;
    return;
end

% modify start by one frame, if need to reframe.
if is_reframing_indicated(clip_struct),
    is_start = is_start + 1;
end

% Convert from time-slice length in seconds to time slice length in frames.
[tslice_frames,over_sec] = tslice_conversion(tslice_length_sec, clip_struct.fps);

% move starting point in by the TIS skip & frames
is_start = is_start + tis_skip + tis_frames;

% Increase total until pass the actual amount.  Then, go back one.
total_frames = is_stop - is_start + 1;
for total = 1:total_frames+1,
    prev_tslice_frames = (total-1) * tslice_frames - floor(over_sec * total);
    %fprintf('%d: %d to %d\n', total, prev_tslice_frames+1, prev_tslice_frames + tslice_frames);
    if prev_tslice_frames + tslice_frames + is_start - 1 > is_stop,
        break;
    end
end
total = total - 1;

% error checking
if total == 0,
    error('Time-slice length in seconds is longer than this clip');
end