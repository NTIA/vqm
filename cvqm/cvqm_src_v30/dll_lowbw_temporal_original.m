function dll_lowbw_temporal_original(fn, delay);
% DLL_LOWBW_TEMPORAL_ORIGINAL
%   step 3: Apply delay to original video
% SYNTAX
%   dll_lowbw_temporal_original(fn, delay);
% DESCRIPTION
%   'fn' from dll_video, fn=1, description of original video sequence.
%   'delay' as returned by dll_lowbw_temporal

% Adjust the image buffer according to the temporal registration just computed.
if delay < 0,
    [fps] = dll_video('fps', fn);  
    dll_video('discard', fn, floor(-delay)/fps); 
end