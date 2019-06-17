function dll_lowbw_temporal_processed(fn, delay);
% DLL_LOWBW_TEMPORAL_PROCESSED
%   step 4: Apply delay to processed video
% SYNTAX
%   dll_lowbw_temporal_processed(fn, delay);
% DESCRIPTION
%   'fn' from dll_video, fn=2, description of processed video sequence.
%   'delay' as returned by dll_lowbw_temporal


% Adjust the image buffer according to the temporal registration just computed.
if delay > 0,
    [fps] = dll_video('fps', fn);  
    dll_video('discard', fn, floor(delay)/fps); 
end
