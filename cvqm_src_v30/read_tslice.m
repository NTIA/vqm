function [y, cb, cr, is_roi] = read_tslice (test_struct, clip_struct, ...
    tslice_sec, number, varargin)
% READ_TSLICE
%  Read and calibrate a time-slice of images associated with one particular 
%  video clip.  Scope includes spatial & temporal registration, luminance 
%  gain/offset, region of interest, and reframing (when needed).
%  Compatible with Big-YUV and uncompressed AVI files only.
% SYNTAX
%
%  [y] = read_tslice (test_struct, clip_struct, tslice_sec, number);
%  [y] = read_tslice (..., 'PropertyName', PropertyValue, ...);
%  [y] = read_tslice (..., 'Flag', ...);
%  [y, cb, cr] = read_tslice (...);
%  [y, cb, cr, sroi] = read_tslice (...);
%
% DESCRIPTION
% 
%  [y] = read_tslice (test_struct, clip_struct, tslice_sec, number)  
%  reads in one time-slice of images from the test specified in 'test_struct'
%  (of the same format as GTests), from the clip specified in 'clip_struct'
%  (of the same format as GClips).  'tslice_sec' is the length of one 
%  time-slice in seconds; and 'number' indexes through the time-slices available for 
%  this clip, 1 to M.  The images read are calibrated for luminance gain, 
%  luminance offset and spatial registration.
%
%  Return value 'Y' contains the read images as a 3-D matrix:  row, column,
%  time.  Return value PROI contains the processed region of interest,
%  calcualted from this clip's CVR.  Return values 'Cb' and 'Cr' contain
%  the color image planes, also calibrated. (pixel replication, not interp)
%
%  The following optional arguments may be specified: 
%
%  'precision',value,  Return image precision ('double' or 'single').  The
%                      default is 'double', unless the returned video has
%                      more than 650 rows (in which case the default is
%                      'single').
%  'sroi',...       Requested spatial region of interest (SROI), overriding
%                   default values of SROI.  Must be followed by 4 values, 
%                   specifying the region of interest, in the order: top, 
%                   left, bottom, right.  SROI will be adjusted.
%  'align_start',value      Specified 'value' overrides clip_struct.align_start
%  'align_stop', value      Specified 'value' overrides clip_struct.align_stop
%  'unaligned'                Read timeslices based on total range of frames
%                                 (loc_start & loc_stop) instead of aligned frames.
%                                 Cannot be used with 'align_start' & 'align_stop'.
%  'aligned'                    Read timeslices based on aligned segment.
%                                  This is the default behavior.
%  'hsize', value,  Horizontal size of S-T blocks.  SROI must evenly divide
%                   by this value, horizontally.
%  'vsize', value,  Vertically size of S-T blocks.  SROI must evenly divide
%                   by this value, vertizontally.
%  'extra', value,  This many valid pixels are required on all sides of the
%                   SROI, extra pixels for filtering.
%  'yxextra', y, x,  This many valid pixels are required. 'x' on left and
%                   right, 'y' on top and bottom of SROI.
%  'full',          Ignore default SROI.  Return full image, with pixels
%                   outside of CVR black.  Cannot be used with 'extra',
%                   'sroi','hsize', or 'vsize'.
%  'TI',            Need this many extra images before this time-slice.  Extra
%                   image is NOT supplied for the first time-slice.
%  'TIS', sec,      Will be doing a 'sec' wide temporal filter. Thus, need
%                   'sec' seconds extra frames before this
%                   time-slice.  Also, time-slices will be lined up in time
%                   with time-slices that don't have the 'TIS' option.
%                   
%
%  Return the actual SROI used.  If the full image is requested ('full')
%  then return the CVR, instead.  If extra pixels are requested, the
%  returned SROI will NOT include those extra pixels.
%
% EXAMPLE
%  [y] = read_tslice( GTests(2), GClips(57), 1.0/5.0, 1, 'ROI', 21, 21, 468, 700, ...
%                   'align_start', 20, 'align_stop', 280);
% REMARKS
%
% Color image pixels will be pixel replicated, so that Cb and Cr
% images are not subsampled by 2 horizontally.
%
% See function 'default_sroi' for default requested SROI.
% SROI will be modified as specified by function 'adjust_requested_sroi'.
%
% Can read UYVY files (Big-YUV) with 'yuv' suffix, and uncompressed AVI
% files (UYVY or RGB) with 'avi' suffix.

if length(test_struct) > 1,
    num = search_test_list(test_struct, clip_struct);
    test_struct = test_struct(num);
end
% read values from clip_struct that can be over written by variable argument
% list.
is_whole_image = 0;
default_roi = 1;
is_yextra = 0;
is_xextra = 0;
is_hsize = 1;
is_vsize = 1;
is_ti = 0;
is_tis = 0;
tis_skip = 0;

is_start = clip_struct.align_start;
is_stop = clip_struct.align_stop;

precision = 'double';  % Normal default precision unless greater than 650 rows
precision_input = 0;  % boolean will be set to 1 if precision option was requested

% parse variable argument list (property values)
cnt = 1;
while cnt <= nargin - 4,
    if strcmpi(varargin(cnt),'sroi') == 1,
        is_roi.top = varargin{cnt+1};
        is_roi.left = varargin{cnt+2};
        is_roi.bottom = varargin{cnt+3};
        is_roi.right = varargin{cnt+4};
        is_whole_image = 0;
        default_roi = 0;
        if is_within_cvr(is_roi, clip_struct) == 0,
            error('ERROR: read_tslice, ROI invalid, must lie within CVR\n');
        end
        cnt = cnt + 5;
    elseif strcmpi(varargin(cnt),'extra') == 1,
        is_yextra = varargin{cnt+1};
        is_xextra = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt),'yxextra') == 1,
        is_yextra = varargin{cnt+1};
        is_xextra = varargin{cnt+2};
        cnt = cnt + 3;
    elseif strcmpi(varargin(cnt),'TI') == 1,
        is_ti = 1;
        cnt = cnt + 1;
        is_tis = 0; tis_skip = 0;
    elseif strcmpi(varargin(cnt),'TIS') == 1,
        is_tis = 1;
        tis_sec = varargin{cnt+1};
        cnt = cnt + 2;
        is_ti = 0;
        
        % compute extra & skip frames 
        [tis_frames] = tslice_conversion(tis_sec, clip_struct.fps);
        [tslice_frames] = tslice_conversion(tslice_sec, clip_struct.fps);
        tis_skip = - mod(tis_frames, tslice_frames);
        if tis_skip < 0,
            tis_skip = tis_skip + tslice_frames;
        end
    elseif strcmpi(varargin(cnt),'align_start') == 1,
        is_start = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt),'align_stop') == 1,
        is_stop = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt),'unaligned') == 1,
        is_start = clip_struct.loc_start;
        is_stop = clip_struct.loc_stop;
        cnt = cnt + 1;
    elseif strcmpi(varargin(cnt),'aligned') == 1,
        is_start = clip_struct.align_start;
        is_stop = clip_struct.align_stop;
        cnt = cnt + 1;
    elseif strcmpi(varargin(cnt),'full') == 1,
        is_whole_image = 1;
        default_roi = 0;
        is_roi.top = 1;
        is_roi.left = 1;
        is_roi.bottom = clip_struct.image_size.rows;
        is_roi.right = clip_struct.image_size.cols;

        cnt = cnt + 1;
    elseif strcmpi(varargin(cnt),'hsize') == 1,
        is_hsize = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt),'vsize') == 1,
        is_vsize = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt),'precision') == 1,
        precision = varargin{cnt+1};
        precision_input = 1;
        cnt = cnt + 2;
    else
        error('Property value passed into read_tslice not recognized');
    end
end

if default_roi,
    is_roi = adjust_requested_sroi(clip_struct, 'yxextra', is_yextra, is_xextra, ...
        'vsize', is_vsize, 'hsize', is_hsize);
else
    is_roi = adjust_requested_sroi(clip_struct, 'yxextra', is_yextra, is_xextra, ...
        'vsize', is_vsize, 'hsize', is_hsize, 'sroi', is_roi.top, ...
        is_roi.left, is_roi.bottom, is_roi.right);
end
read_roi.top = is_roi.top-is_yextra;
read_roi.left = is_roi.left-is_xextra;
read_roi.bottom = is_roi.bottom+is_yextra;
read_roi.right = is_roi.right+is_xextra;


% Convert from time-slice length in seconds to time-slice length in frames.
[tslice_frames,over_sec] = tslice_conversion(tslice_sec, clip_struct.fps);

% Compute number of previous time-slices (to be skipped over) and remove from
% that the amount of overlap (including this time-slice's overlap).
prev_tslice_frames = (number-1) * tslice_frames - floor(over_sec * number) + tis_skip;

% check that this requested time-slice exists
if prev_tslice_frames + tslice_frames - 1 > is_stop - is_start + 1,
    error('read_tslice, tried to read past the end of this clip');
elseif number < 1,
    error('read_tslice, tried to read before the beginning of this clip');
end;

hold_len = length(clip_struct.file_name{1});
hold_suffix = clip_struct.file_name{1}(hold_len-2:hold_len);

if strcmpi(hold_suffix,'yuv'),
    % Open image file 
    % [test_struct.path{1} clip_struct.file_name{1}]
    [fid, message] = fopen([test_struct.path{1} clip_struct.file_name{1}], 'r');
    if fid == -1
        fprintf(message);
        error('read_tslice cannot open this clip''s bigyuv file, %s', [test_struct.path{1} clip_struct.file_name{1}]);
    end

    % go to requested location
    num_rows = clip_struct.image_size.rows;
    num_cols = clip_struct.image_size.cols;
    if isnan(is_start),
        error('first frame of this clip is undefined (NaN).  Try optional argument "align_start" or "unaligned" ');
    end
    offset = ((is_start - 1) + prev_tslice_frames); % number of prior images
    if is_ti && number > 1, % handle TI read of extra image.
        offset = offset - 1;
        tslice_frames = tslice_frames + 1;
    elseif is_tis, % read extra frames at the END for TIS case.
        tslice_frames = tslice_frames + tis_frames;
    end;
    offset = offset * num_rows * num_cols * 2; %pixels each image
    status = fseek(fid, offset, 'bof');

    if status == -1,
        fclose(fid);
        error('read_tslice cannot seek requested image location. Run "check_clips" on your clip structure.');
    end

    % Set default precision to single for large images
    if ((read_roi.bottom - read_roi.top + 1 > 650) && precision_input==0)
        precision = 'single';
    end
    
    % initialize memory to hold return images.
    if ~is_whole_image,
        y = zeros(read_roi.bottom - read_roi.top + 1,read_roi.right - read_roi.left + 1,tslice_frames, precision);
    else
        y = zeros(clip_struct.image_size.rows,clip_struct.image_size.cols,tslice_frames, precision);
    end

    if (nargout >= 3),
        cb = y;
        cr = y;
    end

    %Read extra frame for reframing.  Do reframing immediately after the fread.
    do_reframing = is_reframing_indicated(clip_struct);
    if do_reframing,
        [prev_fread,count] = fread(fid, [2*num_cols,num_rows], 'uint8=>uint8');
        prev_fread = reshape(prev_fread,num_cols*2,2,num_rows/2);
    end

    % figure out which rows, cols to pick off for ROI & spatial registration.
    if is_whole_image,
        read_roi = clip_struct.cvr;
        is_roi = read_roi;
        if is_sroi_valid(read_roi,clip_struct) == 0 || ...
                isnan(read_roi.top) || isnan(read_roi.left) || isnan(read_roi.bottom) || isnan(read_roi.right),
            warning('Warning: the CVR for this clip is invalid.\nWill use a CVR set to maximum image content.\n');
            read_roi.top = 1;
            read_roi.left = 1;
            read_roi.bottom = clip_struct.image_size.rows;
            read_roi.right = clip_struct.image_size.cols;
            if clip_struct.spatial.horizontal >= 0,
                read_roi.right = read_roi.right - clip_struct.spatial.horizontal;
            else
                read_roi.left = read_roi.left - clip_struct.spatial.horizontal;
            end
            if clip_struct.spatial.vertical >= 0,
                read_roi.bottom = read_roi.bottom - clip_struct.spatial.vertical - is_reframing_indicated(clip_struct);
            else
                read_roi.top = read_roi.top - clip_struct.spatial.vertical + is_reframing_indicated(clip_struct);
            end
            if is_sroi_valid(read_roi,clip_struct) == 0,
                error('ERROR: Code Defect.  The CVR for this clip is STILL invalid.\nWill use a CVR set to maximum image content.\n');
            end            
        end
    end
    rows_raw = [read_roi.top:read_roi.bottom];
    cols_raw = [read_roi.left:read_roi.right];
    if do_reframing,
        if strcmp(clip_struct.video_standard,'interlace_lower_field_first') == 1,
            rows = rows_raw + round(clip_struct.spatial.vertical - 1);
        else
            rows = rows_raw + round(clip_struct.spatial.vertical + 1);
        end
    else 
        rows = rows_raw + round(clip_struct.spatial.vertical);
    end
    cols = cols_raw + clip_struct.spatial.horizontal;

    % loop through & read in the time-slice of images
    for cnt = 1:tslice_frames,
        where = ftell(fid);
        [hold_fread,count] = fread(fid, [2*num_cols,num_rows], 'uint8=>uint8');
        if count ~= 2*num_cols*num_rows,
            % try one more time.
            fprintf('Warning: read_tslice could not read entirety of requested image time-slice; re-trying\n');
            %pause(5);
            if where == -1,
                fprintf('Could not determine current location.  Re-try failed.\n');
                fclose(fid);
                error('read_tslice could not read entirety of requested image time-slice');
            end
            fseek(fid, where, 'bof');
            [hold_fread,count] = fread(fid, [2*num_cols,num_rows], 'uint8=>uint8');
            if count ~= 2*num_cols*num_rows,
                fclose(fid);
                hold = sprintf('time-slice read failed for time-slice %s:%s(%s)\nread_tslice could not read entirety of requested image time-slice\nRead failed for image %d of file %s', ...
                    clip_struct.test{1}, clip_struct.scene{1}, clip_struct.hrc{1}, ...
                    (offset / (num_rows*num_cols*2)) + cnt, ...
                    [test_struct.path{1} clip_struct.file_name{1}]);
                error(hold);
            end
        end

        % reframe, if need be.
        if do_reframing,
            hold_fread = reshape(hold_fread,num_cols*2,2,num_rows/2);
            temp = hold_fread;
            if strcmp(clip_struct.video_standard,'interlace_lower_field_first') == 1,
                hold_fread(:,1,:) = hold_fread(:,2,:);
                hold_fread(:,2,1:num_rows/2-1) = prev_fread(:,1,2:num_rows/2);
            elseif strcmp(clip_struct.video_standard,'interlace_upper_field_first') == 1,
                hold_fread(:,2,:) = hold_fread(:,1,:);
                hold_fread(:,1,2:num_rows/2) = prev_fread(:,2,1:num_rows/2-1);
            else
                fclose(fid);
                error('GClips.video_standard not recognized.  Should be ''interlace_lower_field_first'' or ''interlace_upper_field_first'' since reframing.');
            end
            prev_fread = temp;
            hold_fread = reshape(hold_fread,2*num_cols,num_rows);
        end

        % pick off the Y plane (luminance)
        temp = reshape(hold_fread', num_rows, 2, num_cols);
        uncalib = squeeze(temp(:,2,:));  
        if isfield(clip_struct,'scale'),
            uncalib = resample_image(double(uncalib), clip_struct.scale.vertical, clip_struct.scale.horizontal, 'interlace');
        end
        if is_whole_image,
            if clip_struct.luminance_offset == 0 && clip_struct.luminance_gain == 1.0,
                y(rows_raw,cols_raw,cnt) = double(uncalib(rows,cols));
            else
                y(rows_raw,cols_raw,cnt) = (double(uncalib(rows,cols)) - clip_struct.luminance_offset ) ./ clip_struct.luminance_gain;
            end
        else
            if clip_struct.luminance_offset == 0 && clip_struct.luminance_gain == 1.0,
                y(:,:,cnt) = double(uncalib(rows,cols));
            else
                y(:,:,cnt) = (double(uncalib(rows,cols)) - clip_struct.luminance_offset ) ./ clip_struct.luminance_gain;
            end
        end

        % If color image planes are requested, pick those off and perform
        % pixel replication to upsample horizontally by 2.
        if nargout >= 3,
            temp = reshape(hold_fread,4,num_rows*num_cols/2);

            color = reshape(temp(1,:),num_cols/2,num_rows)';
            color2 = [color ; color];
            uncalib = reshape(color2,num_rows,num_cols);
            if isfield(clip_struct,'scale'),
                uncalib = resample_image(double(uncalib), clip_struct.scale.vertical, clip_struct.scale.horizontal, 'interlace');
            end
            if is_whole_image
                cb(rows_raw,cols_raw,cnt) = double(uncalib(rows,cols)) - 128;
            else
                cb(:,:,cnt) = double(uncalib(rows,cols)) - 128;
            end

            color = reshape(temp(3,:),num_cols/2,num_rows)';
            color2 = [color ; color];
            uncalib = reshape(color2,num_rows,num_cols);
            if isfield(clip_struct,'scale'),
                uncalib = resample_image(double(uncalib), clip_struct.scale.vertical, clip_struct.scale.horizontal, 'interlace');
            end
            if is_whole_image,
                cr(rows_raw,cols_raw,cnt) = double(uncalib(rows,cols)) - 128;
            else
                cr(:,:,cnt) = double(uncalib(rows,cols)) - 128;
            end
        end
    end

    fclose(fid);
    
elseif strcmpi(hold_suffix,'avi'),
    

    % compute number of frames to skip
    offset = (is_start - 1 + prev_tslice_frames); % number of prior images
    if is_ti && number > 1, % handle TI read of extra image.
        offset = offset - 1;
        tslice_frames = tslice_frames + 1;
    elseif is_tis, % read extra frames at the END for TIS case.
        tslice_frames = tslice_frames + tis_frames;
    end;

    % Set default precision to single for large images
    if ((read_roi.bottom - read_roi.top + 1 > 650) && precision_input==0)
        precision = 'single';
    end
    
    % initialize memory to hold return images.
    if ~is_whole_image,
        y = zeros(read_roi.bottom - read_roi.top + 1,read_roi.right - read_roi.left + 1,tslice_frames, precision);
    else
        y = zeros(clip_struct.image_size.rows,clip_struct.image_size.cols,tslice_frames, precision);
    end

    if (nargout >= 3),
        cb = y;
        cr = y;
    end

    %Read extra frame for reframing.  Do reframing immediately after the fread.
    do_reframing = is_reframing_indicated(clip_struct);
    if do_reframing,
        [prev_y,prev_cb,prev_cr] = read_avi('YCbCr',[test_struct.path{1} clip_struct.file_name{1}], ...
            'frames',offset+1, offset+1, '128');  
        offset = offset + 1;
        [num_rows,num_cols] = size(prev_y);
        prev_y = reshape(prev_y,2,num_rows/2,num_cols);
        prev_cb = reshape(prev_cb,2,num_rows/2,num_cols);
        prev_cr = reshape(prev_cr,2,num_rows/2,num_cols);
    end

    % figure out which rows, cols to pick off for ROI & spatial registration.
    if is_whole_image,
        read_roi = clip_struct.cvr;
        is_roi = read_roi;
        if is_sroi_valid(read_roi,clip_struct) == 0 || ...
                isnan(read_roi.top) || isnan(read_roi.left) || isnan(read_roi.bottom) || isnan(read_roi.right),
            warning('Warning: the CVR for this clip is invalid.\nWill use a CVR set to maximum image content.\n');
            read_roi.top = 1;
            read_roi.left = 1;
            read_roi.bottom = clip_struct.image_size.rows;
            read_roi.right = clip_struct.image_size.cols;
            if clip_struct.spatial.horizontal >= 0,
                read_roi.right = read_roi.right - clip_struct.spatial.horizontal;
            else
                read_roi.left = read_roi.left - clip_struct.spatial.horizontal;
            end
            if clip_struct.spatial.vertical >= 0,
                read_roi.bottom = read_roi.bottom - clip_struct.spatial.vertical - is_reframing_indicated(clip_struct);
            else
                read_roi.top = read_roi.top - clip_struct.spatial.vertical + is_reframing_indicated(clip_struct);
            end
            if is_sroi_valid(read_roi,clip_struct) == 0,
                error('ERROR: Code Defect.  The CVR for this clip is STILL invalid.\nWill use a CVR set to maximum image content.\n');
            end            
        end
    end
    rows_raw = [read_roi.top:read_roi.bottom];
    cols_raw = [read_roi.left:read_roi.right];
    if do_reframing,
        if strcmp(clip_struct.video_standard,'interlace_lower_field_first') == 1,
            rows = rows_raw + round(clip_struct.spatial.vertical - 1);
        else
            rows = rows_raw + round(clip_struct.spatial.vertical + 1);
        end
    else 
        rows = rows_raw + round(clip_struct.spatial.vertical);
    end
    cols = cols_raw + clip_struct.spatial.horizontal;

    % loop through & read in the time-slice of images
    if (nargout >= 3)
        [all_y,all_cb,all_cr] = read_avi('YCbCr',[test_struct.path{1} clip_struct.file_name{1}], ...
            'frames',offset+1, offset+tslice_frames, '128');
    else
        [all_y] = read_avi('YCbCr',[test_struct.path{1} clip_struct.file_name{1}], ...
            'frames',offset+1, offset+tslice_frames, '128');
    end
      
    for cnt = 1:tslice_frames,
        if (nargout >= 3)
            curr_y = all_y(:,:,cnt);
            curr_cb = all_cb(:,:,cnt);
            curr_cr = all_cr(:,:,cnt);
        else
            curr_y = all_y(:,:,cnt);
        end

        % reframe, if need be.
        if do_reframing,
            
            hold_y = reshape(curr_y,2,num_rows/2,num_cols);
            temp_y = hold_y;
            if strcmp(clip_struct.video_standard,'interlace_lower_field_first') == 1,
                hold_y(1,:,:) = hold_y(2,:,:);
                hold_y(2,1:num_rows/2-1,:) = prev_y(1,2:num_rows/2,:);
            elseif strcmp(clip_struct.video_standard,'interlace_upper_field_first') == 1,
                hold_y(2,:,:) = hold_y(1,:,:);
                hold_y(1,2:num_rows/2,:) = prev_y(2,1:num_rows/2-1,:);
            else
                error('GClips.video_standard not recognized.  Should be ''interlace_lower_field_first'' or ''interlace_upper_field_first'' since reframing.');
            end
            prev_y = temp_y;
            curr_y = reshape(hold_y,num_rows,num_cols);
            if nargout >= 3,
                hold_cb = reshape(curr_cb,2,num_rows/2,num_cols);
                temp_cb = hold_cb;
                if strcmp(clip_struct.video_standard,'interlace_lower_field_first') == 1,
                    hold_cb(1,:,:) = hold_cb(2,:,:);
                    hold_cb(2,1:num_rows/2-1,:) = prev_cb(1,2:num_rows/2,:);
                elseif strcmp(clip_struct.video_standard,'interlace_upper_field_first') == 1,
                    hold_cb(2,:,:) = hold_cb(1,:,:);
                    hold_cb(1,2:num_rows/2,:) = prev_cb(2,1:num_rows/2-1,:);
                else
                    error('GClips.video_standard not recognized.  Should be ''interlace_lower_field_first'' or ''interlace_upper_field_first'' since reframing.');
                end
                prev_cb = temp_cb;
                curr_cb = reshape(hold_cb,num_rows,num_cols);
                
                hold_cr = reshape(curr_cr,2,num_rows/2,num_cols);
                temp_cr = hold_cr;
                if strcmp(clip_struct.video_standard,'interlace_lower_field_first') == 1,
                    hold_cr(1,:,:) = hold_cr(2,:,:);
                    hold_cr(2,1:num_rows/2-1,:) = prev_cr(1,2:num_rows/2,:);
                elseif strcmp(clip_struct.video_standard,'interlace_upper_field_first') == 1,
                    hold_cr(2,:,:) = hold_cr(1,:,:);
                    hold_cr(1,2:num_rows/2,:) = prev_cr(2,1:num_rows/2-1,:);
                else
                    error('GClips.video_standard not recognized.  Should be ''interlace_lower_field_first'' or ''interlace_upper_field_first'' since reframing.');
                end
                prev_cr = temp_cr;
                curr_cr = reshape(hold_cr,num_rows,num_cols);
            end
            
        end

        % Handle the Y plane (luminance)
        if isfield(clip_struct,'scale'),
            curr_y = resample_image(double(curr_y), clip_struct.scale.vertical, clip_struct.scale.horizontal, 'interlace');
        end
        if is_whole_image,
            if clip_struct.luminance_offset == 0 && clip_struct.luminance_gain == 1.0,
                y(rows_raw,cols_raw,cnt) = curr_y(rows,cols);
            else
                y(rows_raw,cols_raw,cnt) = (double(curr_y(rows,cols)) - clip_struct.luminance_offset ) ./ clip_struct.luminance_gain;
            end
        else
            if clip_struct.luminance_offset == 0 && clip_struct.luminance_gain == 1.0,
                y(:,:,cnt) = double(curr_y(rows,cols));
            else
                y(:,:,cnt) = (double(curr_y(rows,cols)) - clip_struct.luminance_offset ) ./ clip_struct.luminance_gain;
            end
        end

        % If color image planes are requested, pick those off and perform
        % pixel replication to upsample horizontally by 2.
        if nargout >= 3,
            if isfield(clip_struct,'scale'),
                curr_cb = resample_image(curr_cb, clip_struct.scale.vertical, clip_struct.scale.horizontal, 'interlace');
                curr_cr = resample_image(curr_cr, clip_struct.scale.vertical, clip_struct.scale.horizontal, 'interlace');
            end
            if is_whole_image
                cb(rows_raw,cols_raw,cnt) = double(curr_cb(rows,cols));
                cr(rows_raw,cols_raw,cnt) = double(curr_cr(rows,cols));
            else
                cb(:,:,cnt) = double(curr_cb(rows,cols));
                cr(:,:,cnt) = double(curr_cr(rows,cols));
            end
        end
        
    end

else
    error(sprintf('Video file type not recognized.\n''yuv'' suffix required for Big-YUV files;\n''avi'' suffix required for uncompressed AVI files'));
end

