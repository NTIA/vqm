function [one,two,three,four] = dll_calib_video(control, fn, varargin);
% DLL_CALIB_VIDEO
%   This function implements calibrated video file read.  Understands
%   model's SROI.
% SYNTAX
%   [...] = dll_calib_video(control, ...);
%   [...] = dll_calib_video(control, fn, ...);
% DESCRIPTION
%   'fn' is either 1 for original, or 2 for processed (when required)
%   'control' is one of the following strings.  Additional parameters may
%   be required, as specified below:
%
% dll_calib_video('initialize', fn);
%           % initialize calibration.  dll_video(fn) must have been
%           % initialized on this computer, for either original or processed.
%
% [pvr] = dll_calib_video('pvr');
%           % get PVR.
%
% dll_calib_video('sroi', roi, extra);
%           % set Spatial Region of Interest required by model, 'extra' is
%           % the extra pixels needed on all sides for spatial filtering.
%           % roi should contain a new region of interest structure.
%
% dll_calib_video('max_roi');
%           % maximize SROI and PVR, given shifts.
%
% dll_calib_video('calibration', horiz, vert, pvr, gain, offset, horiz_stretch, vert_stretch);
%           % Set calibration values for processed video.  horiz & vert are 
%           % spatial registration; PVR is Destination valid region. 
%           % WARNING: if fn=1 and fn=2 calculated on two different
%           % computers, this call must be made on BOTH computers.
%
% [y,cb,cr] = dll_calib_video('sec', fn, durration);
%           % read [Y,Cb,Cr] images, 'durration' seconds, and calibrate.
%           % Return entire image, pixels outside of PVR replace with black.
%
% [y,cb,cr] = dll_calib_video('tslice', fn);
%           % read the next tslice of [Y,Cb,Cr] images, and calibrate!
%           % Return pixels within SROI, only.
%
% dll_calib_video('clear');  
%           % Clear calibration values.
%
% [y] = dll_calib_video('peek', fn, durration);  Get the Y
%               images without removing them from the buffer.  So, the next
%               call with 'sec' or 'tslice' will get these same frames.
%               Perform shift & scaling & valid region but NOT gain &
%               offset!
%
% dll_calib_video('luma', gain, offset);
%           % Set luminance gain & offset values for processed video.  
%           % WARNING: if fn=1 and fn=2 calculated on two different
%           % computers, this call must be made on BOTH computers.
%
% dll_calib_video('set_reframe', value);
%           % set reframe to yes (value=1) or no (values=0), ignoring
%           % spatial shift.  Next set of spatial shift will over-ride.
%
% value = dll_calib_video('get_reframe');
%           % get whether reframe.  yes (value=1) or no (values=0).
%
% value = dll_calib_video('total_sec', fn)
%           % return the total number of seconds of CALIBRATED video left
%           % in the file, after the current "read" point.
%
% dll_calib_video('image_mode', flag);
%           % set image mode to either 'field' or 'frame'.
%           % This is useful when the definition of upper and lower field
%           % must remain correct.
%
% dll_calib_video('set_vfd', indices, indices_fuzzy);
%           % set the new indicies for the original video to match the
%           % processed video after performing variable frame delay
%           % calibrations.
%
% dll_calib_video('get_vfd');
%           % get the values of the vfd indices pairs in a vector.


persistent CALIB;

if strcmp(control, 'initialize'),
    CALIB.do_reframe = 0;
    CALIB.do_calibration = 0;
    CALIB.horizontal = 0;
    CALIB.vertical = 0;
    CALIB.pvr = dll_default_vr(fn);
    CALIB.gain = 1.0;
    CALIB.offset = 0.0;
    CALIB.horiz_stretch = 1000;
    CALIB.vert_stretch = 1000;
    CALIB.sroi = [];
    CALIB.image_mode = 'frame';
    CALIB.vfd_flag = 0;


elseif strcmp(control,'set_reframe'),
    CALIB.do_reframe = fn;
    CALIB.do_calibration = 1;
    
elseif strcmp(control,'get_reframe'),
    one = CALIB.do_reframe;

elseif strcmp(control,'pvr'),
    if isfield(CALIB,'pvr'),
        [one] = CALIB.pvr;
    else
        error('PVR must be defined (default or actual) prior to model calculation');
    end
    
elseif strcmp(control,'print'),
    CALIB
    CALIB.sroi
    CALIB.pvr

elseif strcmp(control,'max_roi'),
    [rows,cols] = dll_video('size', 2);  
    CALIB.sroi.top = 1 - min(0, CALIB.vertical);
    CALIB.sroi.left = 1 - min(0, CALIB.horizontal);
    CALIB.sroi.bottom = rows - max(0, CALIB.vertical);
    CALIB.sroi.right = cols - max(0, CALIB.horizontal);
    if strcmp(CALIB.image_mode, 'field')
        CALIB.sroi = gciC_adjust_roi(CALIB.sroi);
    end
    CALIB.pvr = CALIB.sroi;

elseif strcmp(control,'sroi'),
    CALIB.sroi = fn; % a new roi variable
    CALIB.sroi.top = CALIB.sroi.top - varargin{1};
    CALIB.sroi.left = CALIB.sroi.left - varargin{1};
    CALIB.sroi.bottom = CALIB.sroi.bottom + varargin{1};
    CALIB.sroi.right = CALIB.sroi.right + varargin{1};
    if strcmp(CALIB.image_mode, 'field')
        CALIB.sroi = gciC_adjust_roi(CALIB.sroi);
    end

elseif strcmp(control,'sec'),
    [one,two,three] = gciC_tslice(fn, CALIB, varargin{1});

elseif strcmp(control, 'tslice'),
    [one,two,three] = gciC_tslice(fn, CALIB);
    one = one(CALIB.sroi.top:CALIB.sroi.bottom, CALIB.sroi.left:CALIB.sroi.right,:);
    if length(two) > 0,
        two = two(CALIB.sroi.top:CALIB.sroi.bottom, CALIB.sroi.left:CALIB.sroi.right,:);
        three = three(CALIB.sroi.top:CALIB.sroi.bottom, CALIB.sroi.left:CALIB.sroi.right,:);
    end

elseif strcmp(control, 'calibration'),
    CALIB.do_calibration = 1;
    CALIB.horizontal = fn;
    CALIB.vertical = varargin{1};
    if strcmp(CALIB.image_mode, 'field')
        CALIB.pvr = gciC_adjust_roi(varargin{2});
    else
        CALIB.pvr = varargin{2};
    end
    CALIB.gain = varargin{3};
    CALIB.offset = varargin{4};
    CALIB.horiz_stretch = varargin{5};
    CALIB.vert_stretch = varargin{6};
    
    % get buffer image if needed for reframing
    CALIB.do_reframe = 0;
    if mod(abs(CALIB.vertical),2),
        if (dll_video('exist',2) && ~strcmp('progressive',dll_video('get_video_standard', 2))) | ...
                (dll_video('exist',2) && ~strcmp('progressive',dll_video('get_video_standard', 1))),
            CALIB.do_reframe = 1;
        end
    end
    
    % set default sroi
    CALIB.sroi = CALIB.pvr;

elseif strcmp(control, 'luma'),
    CALIB.gain = fn;
    CALIB.offset = varargin{1};

elseif strcmp(control, 'clear'),
    CALIB.do_calibration = 0;
    
elseif strcmp(control, 'peek'),
    [one] = gciC_tslice_yonly_nogain(fn, CALIB, varargin{1});

elseif strcmp(control, 'total_sec'),
    [one] = dll_video('total_sec',fn);
    if CALIB.do_reframe && fn == 2,  % processed only.
        one = one - 1 / dll_video('fps', fn); 
    end

elseif strcmp(control, 'image_mode'),
    if strcmp(fn, 'field')
        CALIB.image_mode = 'field';
    elseif strcmp(fn, 'frame')
        CALIB.image_mode = 'frame';
    else
        error('With control image_mode, next argument must be either "field" or "frame".')
    end
    
elseif strcmp(control, 'set_vfd')
    CALIB.vfd = fn;
    CALIB.vfd_fuzzy = varargin{1};
    CALIB.vfd_flag = 1;

elseif strcmp(control, 'get_vfd')
    if CALIB.vfd_flag == 0
        error('You must run vfd calibration before accessing vfd information.')
    end
    one = CALIB.vfd;
    two = CALIB.vfd_fuzzy;
    
end



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,cb,cr] = gciC_tslice(fn, CALIB, durration)
% get requested images.  Do calibration.  Keep & handle reframing buffer.

    % Get time-slice of images & do calibration.
    if fn == 2 && CALIB.do_calibration,
        if exist('durration', 'var'),
            [uy,ucb,ucr] = dll_video('sec', fn, CALIB.do_reframe, durration);
        else
            [uy,ucb,ucr] = dll_video('tslice', fn, CALIB.do_reframe);
        end
        
        y = do_calibration_on_tslice(fn, uy, CALIB.horizontal, CALIB.vertical, ...
            CALIB.pvr, CALIB.gain, CALIB.offset, CALIB.do_reframe, CALIB.horiz_stretch, ...
            CALIB.vert_stretch, 0);
        cb = do_calibration_on_tslice(fn, ucb, CALIB.horizontal, CALIB.vertical, ...
            CALIB.pvr, 1.0, 0.0, CALIB.do_reframe, CALIB.horiz_stretch, CALIB.vert_stretch, 1);
        cr = do_calibration_on_tslice(fn, ucr, CALIB.horizontal, CALIB.vertical, ...
            CALIB.pvr, 1.0, 0.0, CALIB.do_reframe, CALIB.horiz_stretch, CALIB.vert_stretch, 1);
        

    else
        % Get time-slice of images only (no calibration)
        if exist('durration', 'var'),
            [y,cb,cr] = dll_video('sec', fn, 0, durration);
        else
            [y,cb,cr] = dll_video('tslice', fn);
        end
        % If the clips have been processed for vfd
        if CALIB.vfd_flag == 1
            if ~strcmp('progressive', dll_video('get_video_standard', 2)) % interlaced
                [y_one, y_two] = split_into_fields(y); % split the frames into fields
                [~, ~, time] = size(y_two);
                if strcmp('interlace_upper_field_first', dll_video('get_video_standard',1)) % upper field first
                    y_combo(:, :, 1:2:(time*2-1)) = y_two; % put them in a single array
                    y_combo(:, :, 2:2:(time*2)) = y_one;
                    y_combo = y_combo(:, :, CALIB.vfd);    % index using the vfd information
                    y_two = y_combo(:, :, 1:2:end);        % split back up into upper fields and lower fields
                    y_one = y_combo(:, :, 2:2:end);
                    y = join_into_frames(y_one, y_two);    % join the fields into frames
                    [cb_one, cb_two] = split_into_fields(cb);
                    cb_combo(:, :, 1:2:(time*2-1)) = cb_two;
                    cb_combo(:, :, 2:2:(time*2)) = cb_one;
                    cb_combo = cb_combo(:, :, CALIB.vfd);
                    cb_two = cb_combo(:, :, 1:2:end);
                    cb_one = cb_combo(:, :, 2:2:end);
                    cb = join_into_frames(cb_one, cb_two);
                    [cr_one, cr_two] = split_into_fields(cr);
                    cr_combo(:, :, 1:2:(time*2-1)) = cr_two;
                    cr_combo(:, :, 2:2:(time*2)) = cr_one;
                    cr_combo = cr_combo(:, :, CALIB.vfd);
                    cr_two = cr_combo(:, :, 1:2:end);
                    cr_one = cr_combo(:, :, 2:2:end);
                    cr = join_into_frames(cr_one, cr_two);
                else                                                    % lower field first
                    y_combo(:, :, 1:2:(time*2-1)) = y_one;
                    y_combo(:, :, 2:2:(time*2)) = y_two;
                    y_combo = y_combo(:, :, CALIB.vfd);
                    y_one = y_combo(:, :, 1:2:end);
                    y_two = y_combo(:, :, 2:2:end);
                    y = join_into_frames(y_one, y_two);
                    [cb_one, cb_two] = split_into_fields(cb);
                    cb_combo(:, :, 1:2:(time*2-1)) = cb_one;
                    cb_combo(:, :, 2:2:(time*2)) = cb_two;
                    cb_combo = cb_combo(:, :, CALIB.vfd);
                    cb_one = cb_combo(:, :, 1:2:end);
                    cb_two = cb_combo(:, :, 2:2:end);
                    cb = join_into_frames(cb_one, cb_two);
                    [cr_one, cr_two] = split_into_fields(cr);
                    cr_combo(:, :, 1:2:(time*2-1)) = cr_one;
                    cr_combo(:, :, 2:2:(time*2)) = cr_two;
                    cr_combo = cr_combo(:, :, CALIB.vfd);
                    cr_one = cr_combo(:, :, 1:2:end);
                    cr_two = cr_combo(:, :, 2:2:end);
                    cr = join_into_frames(cr_one, cr_two);
                end
            else   % progressive
                y = y(:, :, CALIB.vfd);
                cb = cb(:, :, CALIB.vfd);
                cr = cr(:, :, CALIB.vfd);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c_image = do_calibration_on_tslice(fn, u_image, horizontal, vertical, pvr, ...
    gain, offset, do_reframe, horiz_stretch, vert_stretch, is_color);

if strcmp('interlace_lower_field_first',dll_video('get_video_standard',fn)),
    f1 = 2;
    f2 = 1;
elseif strcmp('interlace_upper_field_first', dll_video('get_video_standard',fn)),
    f1 = 1;
    f2 = 2;
elseif strcmp('progressive', dll_video('get_video_standard',fn)),
end

[rows,cols, frames] = size(u_image);

% do reframing & shift, if required
if do_reframe,
	% reshape the images
    u_image = reshape(u_image, 2,rows/2,cols,frames);

    c_image = u_image(:,:,:,2:frames);
%    c_image = zeros(2,rows/2,cols,frames-1);
    %
    if strcmp('interlace_lower_field_first',dll_video('get_video_standard',fn)),
        c_image(f2,1:rows/2,:,:) = u_image(f1,1:rows/2,:,2:frames);
        c_image(f1,1:rows/2-1,:,:) = u_image(f2,2:rows/2,:,1:(frames-1));
        vertical = vertical - 1;
    else % strcmp('interlace_upper_field_first',dll_video('get_video_standard',fn))
        c_image(f1,2:rows/2,:,:) = u_image(f2,1:((rows/2)-1),:,1:(frames-1));
        c_image(f2,1:rows/2,:,:) = u_image(f1,1:rows/2,:,2:frames);
        vertical = vertical + 1;
    end
    
    u_image = reshape(c_image, rows, cols, frames-1);
end

if horiz_stretch ~= 1000 || vert_stretch ~= 1000,
    [row,col,time] = size(u_image);
    if strcmp('progressive', dll_video('get_video_standard',fn)),
        for cnt = 1:time,
            if is_color,
                u_image2(:,:,cnt) = resample_image(double(u_image(:,:,cnt)), ...
                    vert_stretch, horiz_stretch, 'Fast');
            else
                u_image2(:,:,cnt) = resample_image(double(u_image(:,:,cnt)), ...
                    vert_stretch, horiz_stretch);
            end
        end
    else
        for cnt = 1:time,
            if is_color,
                u_image2(:,:,cnt) = resample_image(double(u_image(:,:,cnt)), ...
                    vert_stretch, horiz_stretch, 'Fast', 'Interlace');
            else
                u_image2(:,:,cnt) = resample_image(double(u_image(:,:,cnt)), ...
                    vert_stretch, horiz_stretch, 'Interlace');
            end
        end
    end
    u_image = u_image2;
    clear u_image2;
end

% undo shift.
c_image = circshift(u_image,[-vertical, -horizontal, 0]);
% undo gain & offset
if gain ~= 1.0 && offset ~= 0.0,
    c_image = double(c_image);
    c_image =  (c_image - offset) / gain;
end

% zero area voided by PVR
c_image(1:pvr.top-1,:) = 0;
c_image(pvr.bottom+1:rows,:) = 0;
c_image(:,1:pvr.left-1) = 0;
c_image(:,pvr.right+1:cols) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = gciC_tslice_yonly_nogain(fn, CALIB, durration);
% get requested images.  Do calibration.  Keep & handle reframing buffer.

    % Get time-slice of images.
    [y] = dll_video('peek', fn, CALIB.do_reframe, durration);
    
    % Do calibration if appropriate.  Skip gain/offset.
    % perform in groups of 10 images, because this can be memory intensive.
    if fn == 2 && CALIB.do_calibration,
        [r,c,t]=size(y);
        y2 = zeros(r,c,t-CALIB.do_reframe,'single');
        for i=1:10:(t-CALIB.do_reframe),
            j=min(i+9,t-CALIB.do_reframe);
            y2(:,:,i:j) = do_calibration_on_tslice(fn, y(:,:,i:(j+CALIB.do_reframe)), ...
                CALIB.horizontal, CALIB.vertical, ...
                CALIB.pvr, 1.0, 0.0, CALIB.do_reframe, CALIB.horiz_stretch, CALIB.vert_stretch, 0);
        end
        y = y2;
        clear y2;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_roi] = gciC_adjust_roi(old_roi)
% keep field 1 in field 1, keep field 2 in field 2
% keep the top of roi odd and bottom of roi even
    new_roi = old_roi;
    if(mod(old_roi.top, 2) == 0)
        new_roi.top = old_roi.top + 1;
    end
    if(mod(old_roi.bottom, 2) == 1)
        new_roi.bottom = old_roi.bottom - 1;
    end
    