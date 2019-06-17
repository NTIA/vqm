function [new_clip_struct, status] = check_gain_offset(clip_struct);
% CHECK_GAIN_OFFSET
%   Check luminance gain & offset values for obviously invalid values.
% SYNTAX
%   [new_clip_struct, status] = check_gain_offset(clip_struct);
% DESCRIPTION
%   Take a clip structure, e.g., GClips.  Find the clips where luminance
%   gain and/or luminance offset are obviously invalid.  Change those to
%   default values.  Return an array of unique HRC names of the invalid clips.
%   Also return the updated clip structure.

status = [];
jc = 1;
new_clip_struct = clip_struct;
for cnt=1:length(clip_struct),
    % mark as failure also if gain is unreasonable.
    if new_clip_struct(cnt).luminance_gain < 0.6 || new_clip_struct(cnt).luminance_gain > 1.6 || ...
            new_clip_struct(cnt).luminance_offset < -80 || new_clip_struct(cnt).luminance_offset > 80 || ...
            isnan(new_clip_struct(cnt).luminance_gain) || isnan(new_clip_struct(cnt).luminance_gain),
        new_clip_struct(cnt).luminance_gain = 1.0;
        new_clip_struct(cnt).luminance_offset = 0.0;
        status{jc} = clip_struct(cnt).hrc{1};
        jc = jc + 1;
    end
end

if ~isempty(status),
    status = unique(status);
end

