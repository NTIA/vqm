function [y_gain, y_offset, cb_gain, cb_offset, cr_gain, cr_offset, sucess] = ...
    dll_lowbw_gain_v2_processed(fn, num_sec, orig_y, orig_cb, orig_cr, yesno);
% DLL_LOWBW_GAIN_V2_PROCESSED
%   Calculate processed features and low bandwidth YCbCr gain/offset (rrcal
%   version 2).
% SYNTAX
%   [y_gain, y_offset, cb_gain, cb_offset, cr_gain, cr_offset, status] = ...
%       dll_lowbw_gain_v2_processed(fn, seed_state, num_sec, ...
%           orig_y_blocks, orig_cb_blocks, orig_cr_blocks, yesno);
% DESCRIPTION
%  Calculate processed features needed for low bandwidth gain &
%  offset.  Video clip must be temporally registered first.  Spatial
%  registration & valid region should also be calculated.
%  'fn' is the file identifier from dll_video, preferably fn=2.
%  'num_sec' is the number of seconds of video from file fn=2 that should
%  be used for calibration.  
%
%  Other input values are computed by function dll_lowbw_gain_original.
%  Returned values are 'y_gain' the luminance gain, and 'y_offset', the
%  luminance offset; and likewise for Cb and Cr.
%
%  Return value of 'sucess' is 1 if algorithm succeeds, and 0 if algorithm
%  may have failed, and -1 if a catestrophic failure results in a return of
%  gain=1, offset=0.  Cb & Cr values set to 'nan' if algorithm failed,
%  otherwise not checked. 


sucess = 1;


num_sec = floor(num_sec);

y_gain = 1;
y_offset = 0;


%figure out block size
num_sec = floor(num_sec);

[rows,cols, fps] = dll_video('size', fn);
durration = floor( dll_video('total_frames',2) / fps );
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
proc_y = zeros(vert, horiz, num_sec);
proc_cb = zeros(vert, horiz, num_sec);
proc_cr = zeros(vert, horiz, num_sec);

% loop through frames
dll_video('set_rewind', fn);
dll_video('set_tslice', fn, 1.0/fps);
for loop = 1:num_sec,
    % compute mean of frame
    [y, cb, cr] = dll_calib_video('tslice', fn);
    proc_y(:,:,loop) = block_statistic(y, block_size, block_size, 'mean');
    proc_cb(:,:,loop) = block_statistic(cb, block_size, block_size, 'mean');
    proc_cr(:,:,loop) = block_statistic(cr, block_size, block_size, 'mean');
    
    % skip over the rest of the frames in this second of video.
    if loop ~= num_sec,
        dll_video('discard', fn, (fps-1)/fps); 
    end
end
dll_video('rewind', fn);

% set SROI to PVR again
dll_calib_video('sroi', temp.cvr, 0);

[r1,c1,t1] = size(proc_y);
proc_y = reshape(proc_y, r1*c1*t1, 1);
proc_y = proc_y(yesno);
proc_cb = reshape(proc_cb, r1*c1*t1, 1);
proc_cb = proc_cb(yesno);
proc_cr = reshape(proc_cr, r1*c1*t1, 1);
proc_cr = proc_cr(yesno);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eliminate blocks with clipping
% MUST be done second< so that above indicies
% are identical for y, cb, & cr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = find(orig_y >= 2 & orig_y <= 253 & proc_y >= 2 & proc_y <= 253 );
if length(b) ~= length(orig_y),
    orig_y = orig_y(b);
    proc_y = proc_y(b);
end

b = find(orig_cb >= -126 & orig_cb <= 126 & proc_cb >= -126 & proc_cb <= 126 );
if length(b) ~= length(orig_cb),
    orig_cb = orig_cb(b);
    proc_cb = proc_cb(b);
end

b = find(orig_cr >= -126 & orig_cr <= 126 & proc_cr >= -126 & proc_cr <= 126 );
if length(b) ~= length(orig_cr),
    orig_cr = orig_cr(b);
    proc_cr = proc_cr(b);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gain & offset with final pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (max(max(max(orig_y)))-min(min(min(orig_y)))) < 10,
    y_gain = 1.0;
    y_offset = 0.0;
    sucess = -1;
elseif (max(max(max(proc_y)))-min(min(min(proc_y)))) <= 0
    y_gain = 1.0;
    y_offset = 0.0;
    sucess = -1;
else

    % compute initial gain via linear regression
    y = proc_y;
    x = [ones(length(y),1) orig_y];
    b = x\y;
    r = y - x*b;

    done = 0;
    prev_b = b;
    counter = 0;
    while ~done && counter < 10000,
        counter = counter + 1;
        
        epsilon = 1.0;
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

    if counter < 10000,
        y_gain = b(2);
        y_offset = b(1);
    else
        y_gain = 1.0;
        y_offset = 0.0;
        sucess = -1;
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gain & offset with final pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (max(max(max(orig_cb)))-min(min(min(orig_cb)))) < 10,
    cb_gain = nan;
    cb_offset = nan;
elseif (max(max(max(proc_cb)))-min(min(min(proc_cb)))) <= 0
    cb_gain = nan;
    cb_offset = nan;
else

    % compute initial gain via linear regression
    y = proc_cb;
    x = [ones(length(y),1) orig_cb];
    b = x\y;
    r = y - x*b;

    done = 0;
    prev_b = b;
    counter = 0;
    while ~done && counter < 10000,
        counter = counter + 1;
        
        epsilon = 1.0;
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

    if counter < 10000,
        cb_gain = b(2);
        cb_offset = b(1);
    else
        cb_gain = nan;
        cb_offset = nan;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gain & offset with final pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (max(max(max(orig_cr)))-min(min(min(orig_cr)))) < 10,
    cr_gain = nan;
    cr_offset = nan;
elseif (max(max(max(proc_cr)))-min(min(min(proc_cr)))) <= 0
    cr_gain = nan;
    cr_offset = nan;
else

    % compute initial gain via linear regression
    y = proc_cr;
    x = [ones(length(y),1) orig_cr];
    b = x\y;
    r = y - x*b;

    done = 0;
    prev_b = b;
    counter =0;
    while ~done && counter < 10000,
        counter = counter + 1;
        
        epsilon = 1.0;
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

    if counter < 10000,
        cr_gain = b(2);
        cr_offset = b(1);
    else
        cr_gain = nan;
        cr_offset = nan;
    end
end

% failure causes large offset or small gain.
if y_offset > 20 || y_gain < 0.70,
    sucess = 0;
end
if y_gain < 0.6 || y_gain > 1.6 || y_offset < -80 || y_offset > 80,
    y_gain = 1.0;
    y_offset = 0.0;
    sucess = -1;
end

