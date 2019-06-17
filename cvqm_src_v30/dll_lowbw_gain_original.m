function [orig_blocks] = dll_lowbw_gain_original(fn, num_sec)
% DLL_LOWBW_GAIN_ORIGINAL
%   Calculate original features needed for low bandwidth luminance gain &
%   offset.
% SYNTAX
%  [orig_blocks] = dll_lowbw_gain_original(fn, num_sec)
% DESCRIPTION
%  Calculate original features needed for low bandwidth luminance gain &
%  offset.  Video clip must be temporally registered first.  Spatial
%  registration & valid region should also be calculated.
%  'fn' is the file identifier from dll_video, preferably fn=1.
%  'num_sec' is the number of seconds of video from file fn=1 that should
%  be used for calibration.  
%
%  Return values 'orig_blocks' are required by function dll_lowbw_gain_processed.


num_sec = floor(num_sec);

[rows,cols, fps] = dll_video('size', fn);
durration = floor( dll_video('total_frames',1) / fps );
if durration < num_sec,
    num_sec = durration,
end

if rows <= 216,
    block_size = 20;  
elseif rows <= 384,
    block_size = 30;  
else
    block_size = 46;  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Algorithm:
%
%   Use one frame every second (approx).  
%   Search over ALL frames simultaneously.
%   Search original +- 0 second (yes! ZERO); 
%   Use luminance image only, sub-sampled by block_size.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set SROI given specified block size
[temp.image_size.rows,temp.image_size.cols] = dll_video('size',fn);
[temp.cvr] = dll_calib_video('pvr');
extra = 0;
[sroi,vert,horiz] = adjust_requested_sroi (temp, ...
    'vsize',block_size, 'hsize',block_size, 'extra',extra);
dll_calib_video('sroi', sroi, extra);

% allocate space for results
orig_blocks = zeros(vert, horiz, num_sec);

% loop through frames
dll_video('set_rewind', fn);
dll_video('set_tslice', fn, 1.0/fps);
for loop = 1:num_sec,
    % compute mean of frame
    y = dll_calib_video('tslice', fn);
    orig_blocks(:,:,loop) = block_statistic(y, block_size, block_size, 'mean');
    
    % skip over the rest of the frames in this second of video.
    if loop ~= num_sec,
        dll_video('discard', fn, (fps-1)/fps); 
    end
end
dll_video('rewind', fn);


% set SROI to PVR again
dll_calib_video('sroi', temp.cvr, 0);



