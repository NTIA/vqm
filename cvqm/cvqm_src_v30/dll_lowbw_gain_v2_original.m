function [orig_y_blocks, orig_cb_blocks, orig_cr_blocks, yesno] = ...
    dll_lowbw_gain_v2_original(fn, num_sec)
% DLL_LOWBW_GAIN_V2_ORIGINAL
%   Calculate original features needed for low bandwidth YCbCr gain &
%   offset (rrcal version 2).
% SYNTAX
%  [orig_y_blocks, orig_cb_blocks, orig_cr_blocks, yesno] = ...
%       dll_lowbw_gain_v2_original(fn, num_sec)
% DESCRIPTION
%  Calculate original features needed for low bandwidth  gain &
%  offset.  Video clip must be temporally registered first.  Spatial
%  registration & valid region should also be calculated.
%  'fn' is the file identifier from dll_video, preferably fn=1.
%  'num_sec' is the number of seconds of video from file fn=1 that should
%  be used for calibration.  
%
%  Return values are required by function dll_lowbw_gain_processed.
%   'orig_y_blocks' averaged Y blocks selected.
%   'orig_cb_blocks' averaged Cb blocks selected
%   'orig_cr_blocks' averaged Cr blocks selected
%   'yesno' logicals (i.e., booleans) indicating for all blocks, which 1/2
%       of blocks were selected.


num_sec = floor(num_sec);

[rows,cols, fps] = dll_video('size', fn);
durration = floor( dll_video('total_frames',1) / fps );
if durration < num_sec,
    num_sec = durration;
end

if rows <= 216,
    block_size = 10;  
elseif rows <= 384,
    block_size = 22;  
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
orig_y_blocks = zeros(vert, horiz, num_sec);
orig_y_std = zeros(vert, horiz, num_sec);
orig_cb_blocks = zeros(vert, horiz, num_sec);
orig_cr_blocks = zeros(vert, horiz, num_sec);

% loop through frames
dll_video('set_rewind', fn);
dll_video('set_tslice', fn, 1.0/fps);
for loop = 1:num_sec,
    % compute mean of frame
    [y, cb, cr] = dll_calib_video('tslice', fn);
    orig_y_blocks(:,:,loop) = block_statistic(y, block_size, block_size, 'mean');
    orig_y_std(:,:,loop) = block_statistic(y, block_size, block_size, 'std');
    orig_cb_blocks(:,:,loop) = block_statistic(cb, block_size, block_size, 'mean');
    orig_cr_blocks(:,:,loop) = block_statistic(cr, block_size, block_size, 'mean');
   
    % skip over the rest of the frames in this second of video.
    if loop ~= num_sec,
        dll_video('discard', fn, (fps-1)/fps); 
    end
end
dll_video('rewind', fn);


% set SROI to PVR again
dll_calib_video('sroi', temp.cvr, 0);

% pick off 1/2 of blocks with lowest Y stdev

[r1,c1,t1] = size(orig_y_std);

orig_y_std = reshape(orig_y_std, r1*c1*t1, 1);
orig_y_blocks = reshape(orig_y_blocks, r1*c1*t1, 1);
orig_cb_blocks = reshape(orig_cb_blocks, r1*c1*t1, 1);
orig_cr_blocks = reshape(orig_cr_blocks, r1*c1*t1, 1);


[a,b]=sort(orig_y_std);
b = b(1:floor(length(b) / 2));
yesno = logical(orig_y_std <= orig_y_std(b(length(b))));

orig_y_blocks = orig_y_blocks(yesno);
orig_cb_blocks = orig_cb_blocks(yesno);
orig_cr_blocks = orig_cr_blocks(yesno);


