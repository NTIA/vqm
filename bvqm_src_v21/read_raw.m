function [y,cb,cr] = read_raw(test_struct,clip_struct,flag, one, two, varargin);
% SSCQE_READ_RAW
%   Read raw, uncalibrated images from bigyuv-files associated with a
%   clip (GClip).
% SYNTAX
%   [y] = read_raw(test_struct, clip_struct,'tslice', tslice_sec, number);
%   [y] = read_raw(test_struct,clip_struct,'frames', begin, end);
%   [y,cb,cr] = read_raw(...);
%   [...] = read_raw(...,'PropertyName',PropertyValue,...);
% DESCRIPTION
%   Read in raw, uncalibrated images from bigyuv files associated with
%   a video clip.  When the 'tslice' flag is used, the user specifies the
%   length of the time-slice in seconds ('tslice_sec') and the number of
%   the time-slice ('number').  When the 'frames' flag is used, the user
%   specifies the first and last frames, inclusive, to be read ('begin' and
%   'end').  Frame & time-slice numbers are relative to the entire,
%   UNaligned video sequence, starting at clip_struct.loc_start (i.e.,
%   clip_struct.loc_start is frame number 1).
%
%   The luminance plane is returned in 'Y'; the color planes are
%   returned in 'cb' and 'cr' upon request.  Color planes are pixel
%   replicated.  The following optional properties may be requested:
%
%   'sroi'           Spatial region of intrest to be returned.  By default,
%                     the entirity of each image is returned.
%
%   Color image pixels will be pixel replicated, so that Cb and Cr images
%   are not subsampled by 2 horizontally.

%
if length(test_struct) > 1,
    num = search_test_list(test_struct, clip_struct);
    test_struct = test_struct(num);
end
% read values from clip_struct that can be over written by variable argument
% list.
is_whole_image = 1;

is_start = clip_struct.loc_start;
is_stop = clip_struct.loc_stop;

% parse variable argument list (property values)
cnt = 1;
while cnt <= nargin - 5,
    if strcmp(lower(varargin(cnt)),'sroi') == 1,
        sroi.top = varargin{cnt+1};
        sroi.left = varargin{cnt+2};
        sroi.bottom = varargin{cnt+3};
        sroi.right = varargin{cnt+4};
        is_whole_image = 0;
        cnt = cnt + 5;
    else
        error('Property value passed into read_tslice not recognized');
    end
end

% Open image file 
% [test_struct.path{1} clip_struct.file_name{1}]
[fid, message] = fopen([test_struct.path{1} clip_struct.file_name{1}], 'r');
if fid == -1
    fprintf(message);
    error('read_tslice cannot open this clip''s bigyuv file, %s', [test_struct.path{1} clip_struct.file_name{1}]);
end

if strcmp(flag,'tslice'),
    % Find time-slice requested.
    tslice_sec = one;
    number = two;
    
    % Convert from time-slice length in seconds to time-slice length in frames.
    [tslice_frames,over_sec] = tslice_conversion(tslice_sec, clip_struct.fps);

	% Compute number of previous time-slices (to be skipped over) and remove from
	% that the amount of overlap (including this time-slice's overlap).
	prev_tslice_frames = (number-1) * tslice_frames - floor(over_sec * number);
else
    % find range of frames requested.
    prev_tslice_frames = one - 1;
    tslice_frames = two - one + 1;
    number = one;
end

% check that this requested time-slice exists
if prev_tslice_frames + tslice_frames - 1 > is_stop - is_start + 1,
    fclose(fid);
    error('read_tslice, tried to read past the end of this clip');
elseif number < 1,
    fclose(fid);
    error('read_tslice, tried to read before the beginning of this clip');
end;

% go to requested location
if isnan(is_start),
    error('first frame of this clip is undefined (NaN).');
end
num_rows = clip_struct.image_size.rows;
num_cols = clip_struct.image_size.cols;
offset = ((is_start - 1) + prev_tslice_frames); % number of prior images
offset = offset * num_rows * num_cols * 2; %pixels each image
status = fseek(fid, offset, 'bof');

if status == -1,
    fclose(fid);
    error('read_tslice cannot seek requested image location');
end

% initialize memory to hold return images.
y = zeros(clip_struct.image_size.rows,clip_struct.image_size.cols,tslice_frames);

if (nargout == 3),
    cb = y;
    cr = y;
end

% loop through & read in the time-slice of images
this_try = 1;
for cnt = 1:tslice_frames,
    where = ftell(fid);
	[hold_fread,count] = fread(fid, [2*num_cols,num_rows], 'uint8=>uint8');
	if count ~= 2*num_cols*num_rows,
        % try one more time.
        fprintf('Warning: read_tslice could not read entirity of requested image time-slice; re-trying\n');
        %pause(5);
        if where == -1,
            fprintf('Could not determine current location.  Re-try failed.\n');
            error('read_tslice could not read entirity of requested image time-slice');
            fclose(fid);
        end
        fseek(fid, where, 'bof');
    	[hold_fread,count] = fread(fid, [2*num_cols,num_rows], 'uint8=>uint8');
        if count ~= 2*num_cols*num_rows,
            fclose(fid);
            hold = sprintf('time-slice read failed for time-slice %s:%s(%s)\nread_tslice could not read entirity of requested image time-slice\nRead failed for image %d of file %s', ...
                clip_struct.test{1}, clip_struct.scene{1}, clip_struct.hrc{1}, ...
                (offset / (num_rows*num_cols*2)) + cnt, ...
                [test_struct.path{1} clip_struct.file_name{1}]);
            error(hold);
        end
	end
    
    % pick off the Y plane (luminance)
    temp = reshape(hold_fread', num_rows, 2, num_cols);
    uncalib = squeeze(temp(:,2,:));  
    y(:,:,cnt) = double(uncalib);
    
    % If color image planes are requested, pick those off and perform
    % pixel replication to upsample horizontally by 2.
    if nargout == 3,
        temp = reshape(hold_fread,4,num_rows*num_cols/2);

        color = reshape(temp(1,:),num_cols/2,num_rows)';
        color2 = [color ; color];
        uncalib = reshape(color2,num_rows,num_cols);
        cb(:,:,cnt) = double(uncalib) - 128;
        
        color = reshape(temp(3,:),num_cols/2,num_rows)';
        color2 = [color ; color];
        uncalib = reshape(color2,num_rows,num_cols);
        cr(:,:,cnt) = double(uncalib) - 128;
    end
end

fclose(fid);

if ~is_whole_image,
    y = y(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
    if nargout == 3,
        cb = cb(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
        cr = cr(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
    end
end


