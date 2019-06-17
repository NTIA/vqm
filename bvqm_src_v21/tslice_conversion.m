function [tslice_frames, over_sec] = tslice_conversion (tslice_sec, fps)
% TSLICE_CONVERSION
%   Convert from time-slice length in seconds, to time-slice length in
%   frames at a given frame rate.
% SYNTAX
%  [tslice_frames, over_sec] = tslice_conversion (tslice_sec, fps)
% DESCRIPTION
%
%  [tslice_frames, over_sec] = tslice_conversion (tslice_sec, fps)
%  Given the current frame rate in frames per second ('fps') and the length
%  of the current time-slice in seconds ('tslice_sec'), compute the
%  length of the current time-slice in frames ('tslice_frames').  Also return how
%  much the chosen block length exceeds the requested block length, in
%  seconds ('over_sec').
% 
%  NOTE:  if the specified tslice_sec is within (plus or minus) one
%  thousandth of one millisecond of being exactly tslice_frames 
% (e.g., over_sec <= 0.000001), then this function will assume the 
% user intended that exact number of tslice_frames.
% 
%  NOTE:  Any time-slice smaller than one frame will be implemented as
%  one frame per time-slice. 

% compute length of time-slice, in frames
tslice_frames = ceil(tslice_sec * fps);

% compute amount of seconds that above measurement exceeds that requested.
over_sec = tslice_frames - tslice_sec * fps;

% if within 1 ms of an integer tslice_frames with over_sec = 0, use that. 
if tslice_frames == 1,
    over_sec = 0;
elseif over_sec <= 0.000001,
    over_sec = 0;
elseif over_sec >= 0.999999,
    tslice_frames = tslice_frames - 1;
    over_sec = 0;
end
