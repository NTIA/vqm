function [y_gain, y_offset, sucess] = ...
    dll_lowbw_gain_processed(fn, num_sec, orig_blocks);
% DLL_LOWBW_GAIN_PROCESSED
%   Calculate processed features and low bandwidth gain/offset.
% SYNTAX
%   [y_gain, y_offset, status] = ...
%       dll_lowbw_gain_processed(fn, seed_state, num_sec, orig_blocks);
% DESCRIPTION
%  Calculate processed features needed for low bandwidth luminance gain &
%  offset.  Video clip must be temporally registered first.  Spatial
%  registration & valid region should also be calculated.
%  'fn' is the file identifier from dll_video, preferably fn=2.
%  'num_sec' is the number of seconds of video from file fn=2 that should
%  be used for calibration.  
%
%  Input values 'orig_blocks' are computed by function dll_lowbw_gain_original.
%  Returned values are 'y_gain' the luminance gain, and 'y_offset', the
%  luminance offset.
%
%  Return value of 'sucess' is 1 if algorithm succeeds, and 0 if algorithm
%  may have failed, and -1 if a catestrophic failure results in a return of
%  gain=1, offset=0.

num_sec = floor(num_sec);

y_gain = 1;
y_offset = 0;


%figure out block size
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
proc_blocks = zeros(vert, horiz, num_sec);

% loop through frames
dll_video('set_rewind', fn);
dll_video('set_tslice', fn, 1.0/fps);
for loop = 1:num_sec,
    % compute mean of frame
    y = dll_calib_video('tslice', fn);
    proc_blocks(:,:,loop) = block_statistic(y, block_size, block_size, 'mean');
    
    % skip over the rest of the frames in this second of video.
    if loop ~= num_sec,
        dll_video('discard', fn, (fps-1)/fps); 
    end
end
dll_video('rewind', fn);

% set SROI to PVR again
dll_calib_video('sroi', temp.cvr, 0);


[r1,c1,t1] = size(orig_blocks);
[r2,c2,t2] = size(proc_blocks);

orig = reshape(orig_blocks, r1*c1*t1, 1);
proc = reshape(proc_blocks, r1*c1*t1, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gain & offset with final pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute initial gain via linear regression
y = proc;
x = [ones(length(y),1) orig];
b = x\y;
r = y - x*b;

done = 0;
prev_b = b;
while ~done,
	epsilon = 0.1;
	cost = 1.0 ./ (abs(r) + epsilon); % cost vector, reciprocal of errors
	cost = cost ./ sqrt(sum(cost)); % normalize for unity norm
	cost = (cost.^2);

    xp = x' .* repmat(cost,1,2)';
	b = inv(xp*x)*xp*y;
    r = y - x*b;
    
    if abs(prev_b(2) - b(2)) < 0.0001,
        done = 1;
    else
        prev_b = b;
    end
end

sucess = 1;
y_gain = b(2);
y_offset = b(1);

% failure causes large offset or small gain.
if y_offset > 20 || y_gain < 0.70,
    sucess = 0;
end
if y_gain < 0.6 || y_gain > 1.6 || y_offset < -80 || y_offset > 80,
    y_gain = 1.0;
    y_offset = 0.0;
    sucess = -1;
end

