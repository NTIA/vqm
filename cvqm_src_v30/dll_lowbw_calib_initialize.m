function [seed_state] = dll_lowbw_calib_initialize;
% DLL_LOWBW_CALIB_INITIALIZE
%   Initialize low bandwidth calibration.
% SYNTAX
%   [seed_state] = dll_lowbw_calib_initialize;
% DESCRIPTIONS
%   Initialize low bandwidth calibration.  The value returned by this
%   function ('seed_state') must be passed IDENTICALLY to 
%   dll_lowbw_calib_original and then dll_lowbw_calib_processed.  The next
%   time these two functions are again required (i.e., the next video 
%   sequence), this initialization function should be called again.

rand('seed',sum(100*clock));
seed_state = round(rand * 255);
seed_state = uint8(seed_state);

