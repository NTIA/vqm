function [y_gain, y_offset, status] = dll_itu_gain_offset(uncert_sec, frequency)
% DLL_ITU_GAIN_OFFSET
%   Compute frame-based luminance gain and offset.
% SYNTAX
%  [y_gain, y_offset, status] = dll_itu_gain_offset;
%  [...] = dll_itu_gain_offset(uncert, freq);
% DESCRIPTION
%  Calculate luminance gain and offset.
%  Use standard temporal registration algorithm from 2003 NTIA report 
%  & ANSI Rec. 801.03, frame based algorithm.
%
%  Precondition:  dll_video must have been initialized with fn=1 for the
%  original video file, and fn=2 for the processed video file.
%  Spatial Shift must have been calculated and set within dll_calib_video.
%  Valid region must have been calculated and set within dll_calib_video.
%  Additionally, dll_calib_video must have been initialized and set to
%  calibration values for fn=2, particularly processed spatial registration.
%
%  Optional input argument 'uncert' contains the alignment uncertanty, in
%  seconds.  By default 'uncert' = 1.0 (one second).
%  Optional input argument 'freq' contains the frequency of images used, in
%  seconds.  By default 'freq' = 0.5 (one-half second).
%
%  Return values are defined as follows:
%   'y_gain'  
%   'y_offset'
%   'status'        0 if failed, 1 if succeed.


% handle optional properties.
if ~exist('frequency','var'),
    frequency = 0.5;
end
if ~exist('uncert_sec','var'),
    uncert_sec = 1.0;
end

fn1 = 1;
fn2 = 2;

% set subsampling size.
hsize = 16;
vsize = 16;

% clear any luminance gain/offset settings.
dll_calib_video('luma',1.0, 0.0);

% set SROI for processed video read
dll_calib_video('max_roi');

% set default SROI
[temp.image_size.rows,temp.image_size.cols, fps] = dll_video('size',fn1);
[temp.cvr] = dll_calib_video('pvr');
extra = 0;
[sroi,vert,horiz] = adjust_requested_sroi (temp, ...
    'vsize',vsize, 'hsize',hsize, 'extra',extra, 'evenodd');
dll_calib_video('sroi', sroi, extra);

% Compute the number of blocks available
number_tslices = min( dll_video('total_frames',fn1), dll_video('total_frames',fn2));

% Check if this is an interlace or progressive sequence
% Allocate memory to hold all features for this clip.
[video_standard] = dll_video('get_video_standard', fn1);
if strcmp(video_standard,'progressive'),
    is_progressive = 1;
    src = zeros(vert*horiz, number_tslices);
elseif strcmp(video_standard,'interlace_lower_field_first') || strcmp(video_standard,'interlace_upper_field_first'),
    is_progressive = 0;
    src = zeros(vert*horiz, number_tslices*2);
end

% set rewind point
dll_video('set_rewind', fn1);
dll_video('set_rewind', fn2);

% set read size
dll_video('set_tslice', fn1, 1.0/fps);
dll_video('set_tslice', fn2, 1.0/fps);

% read in each original field (or frame), and sub-sample it.
num = 0;
for frame=1:number_tslices,
    [y] = dll_calib_video('tslice', fn1);
    if is_progressive,
        num = num + 1;
        src(:,num) = reshape( block_statistic(y,hsize,vsize,'mean'), vert*horiz, 1);
    else
        [one,two] = split_into_fields(y);

        num = num + 1;
        src(:,num) = reshape( block_statistic(one,hsize/2,vsize,'mean'), vert*horiz, 1);

        num = num + 1;
        src(:,num) = reshape( block_statistic(two,hsize/2,vsize,'mean'), vert*horiz, 1);
    end
end

% display_xyt(src,'zoom4');
% pause;

% read in processed fields (or frames), at the requested interval, and sub-sample them.
all_gain = [];
all_offset = [];
add_frames = round(fps * frequency);
uncert = round(fps * uncert_sec);

for frame=(1+add_frames):add_frames:(number_tslices-add_frames),
    dll_video('rewind', fn2);
    dll_video('discard', fn2, (frame-1)/fps);  
    y = dll_calib_video('tslice', fn2);
    if is_progressive,
        one = reshape( block_statistic(y,hsize,vsize,'mean'), vert*horiz, 1);
        num = frame;
    else
        [one,two] = split_into_fields(y);
        one = reshape( block_statistic(one,hsize/2,vsize,'mean'), vert*horiz, 1);
        two = reshape( block_statistic(two,hsize/2,vsize,'mean'), vert*horiz, 1);
        num = frame*2-1;
    end

    [valid, y_gain, y_offset] = luminance_gain_offset_search(src, one, num, uncert);
    if valid,
        all_gain = [all_gain y_gain];
        all_offset = [all_offset y_offset];
    end

    if ~is_progressive,
        [valid, y_gain, y_offset] = luminance_gain_offset_search(src, two, num+1, uncert);
        if valid,
            all_gain = [all_gain y_gain];
            all_offset = [all_offset y_offset];
        end
    end
end

% median filter & record result for this clip.
if length(all_gain) == 0,
    y_gain = 1.0;
    y_offset = 0.0;
    status = 0;
else
    y_gain = median(all_gain);
    y_offset = median(all_offset);
    status = 1;
end

% mark as failure also if gain is unreasonable.
if y_gain < 0.6 || y_gain > 1.6 || y_offset < -80 || y_offset > 80,
    y_gain = 1.0;
    y_offset = 0.0;
    status = 0;
end


% rewind video files
dll_video('rewind', fn1);
dll_video('rewind', fn2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [valid, y_gain, y_offset] = luminance_gain_offset_search(src, deg, num, uncert)
%

[a, src_length] = size(src);
std_of_diff = zeros(1, src_length);
std_of_diff(:) = NaN;

done = 0;
while ~done,
	for i=max(1,num-uncert) : min(src_length,num+uncert),
        % if haven't filled in this delay's temporal registration, do it
        % now.
        if isnan(std_of_diff(i)),
            std_of_diff(i) = std(src(:,i) - deg);
        end
    end
    
    [a,where] = min(std_of_diff);
    
    % if at end point, return failure.
    if where == 1 || where == src_length 
        valid = 0;
        y_gain = 1;
        y_offset = 0;
        return;
    elseif where == num-uncert || where == num+uncert,
        %fprintf('\tat edge of uncertainty.  extending\n');
        uncert = uncert + 15;
    else
        done = 1;
    end
        
end

%  Check added 4/26/10 to make sure that there is some variance in both the
%  src and deg images to prevent failure of the linear fit
if (std(deg)==0 || std(src(:,where))==0)
    valid = 0;
    y_gain = 1;
    y_offset = 0;
    return;
end

% compute initial gain via linear regression
y = deg;
x = [ones(length(y),1) src(:,where)];
b = x\y;
r = y - x*b;

done = 0;
prev_b = b;
while ~done,
	epsilon = 0.1;
	cost = 1.0 ./ (abs(r) + epsilon); % cost vector, reciprocal of errors
	cost = cost ./ sqrt(sum(cost)); % normalize for unity norm
	cost = (cost.^2);
    temp = eye(length(y),length(y));
    for i=1:length(y),
        temp(i,i) = temp(i,i) * cost(i);
    end
    cost = temp;
	
	b = inv(x'*cost*x)*x'*cost*y;
    r = y - x*b;
    
    if abs(prev_b(2) - b(2)) < 0.0001,
        done = 1;
    else
        prev_b = b;
    end
end

valid = 1;
y_gain = b(2);
y_offset = b(1);


