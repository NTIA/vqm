function [status] = is_reframing_indicated(clip_struct)
% DOES_CLIP_REFRAME
%  Determine whether spatial registration indicates reframing
% SYNTAX
%  [status] = is_reframing_indicated(clip_struct)
% DESCRIPTION
%  [status] = is_reframing_indicated(clip_struct); takes a variable
%  'clip_struct', as defined within GClips.  Return 'status' of 1 if
%  reframing is indicated, return 0 if reframing is not indicated.

% expect 'progressive' or 'interlace_lower_field_first' or 'interlace_upper_field_first' for
% clip_struct.video_standard


% progressive never reframes.
if strcmp(clip_struct.video_standard,'progressive') == 1,
    status = 0;

% interlace, reframing
elseif mod(clip_struct.spatial.vertical,2) == 1,
    status = 1;
  
% interlace, frame alignment
else
    status = 0;
end