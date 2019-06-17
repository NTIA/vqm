function [one,two,three,four] = dll_video(control, fn, varargin);
% DLL_VIDEO
%   This function implements video file read.
% SYNTAX
%   [...] = dll_video(control, fn, ...);
% DESCRIPTION
%   'fn' is either 1 for original, or 2 for processed
%   'control' is one of the following strings.  Additional parameters may
%   be required, as specified below:
%
% dll_video('initialize', fn, file_name, 'avi', video_standard);
% dll_video('initialize', fn, file_name, 'uyvy', video_standard, rows, cols, fps);
%       'video_standard' is 'progressive',
%       'interlace_lower_field_first', or 'interlace_upper_field_first'
%
% dll_video('set_tslice', fn, tslice_len);
%          % Set durration retrieved by command 'tslice'
%
% [rows,cols,fps,durration] = dll_video('size', fn);  
%          % return image size, frames per second, TOTAL file durration
%
% [fps] = dll_video('fps', fn);  
% [fps] = dll_video('fps');  
%          % return frames per second
%
% [video_standard] = dll_video('get_video_standard', fn);
%           % return video standard:  'progressive',
%           'interlace_lower_field_first' or 'interlace_upper_field_first'
%
% dll_video('set_rewind', fn);
%           % Set 'rewind'to go to current point in the file
%
% dll_video('rewind', fn);
%           % Go to point in file specified by 'set_rewind'
%
% dll_video('discard', fn, durration);  
%           % Discard the next 'durration' seconds of images from the buffer.   
%
% [y,cb,cr] = dll_video('tslice', fn, reframe);
%           % read [Y,Cb,Cr] images, the next tslice in durration, 
%           % previously specified via 'set_tslice'
%
% [y,cb,cr] = dll_video('sec', fn, reframe, durration);
% [y] = dll_video('sec', fn, reframe durration);
%           % read [Y,Cb,Cr] images, 'durration' seconds.
%
% [y] = dll_video('peek', fn, reframe durration);  Get the Y
%               images without removing them from the buffer.  So, the next
%               call with 'sec' or 'tslice' will get these same frames.
%               'reframe' is 1 if ONE extra farme neede for reframing, 0 otherwise 
%
% [total] = dll_video('total_frames',fn);
%           % return the total number of frames left in the file, after the
%           % current "read" point.
%
% [total] = dll_video('total_sec',fn);
%           % return the total number of seconds left in the file, after the
%           % current "read" point.
%
% [boolean] = dll_video('exist',fn);
%           % return 1 (true) if file defined/exists, 0 (false) otherwise.
%
% [code] = dll_set_align('delay_8s', delay);
%           % Works only for files exactly 8seconds long.  
%           % With delay=0, discard first 0.8s start and last 0.2s.
%           % Delays around that allowed, from (-0.2) sec to (0.2sec - 1frame).
%           % For original, delay must be set to 0.
%           % After this option is chosen, 'size' will always return a
%           % durration of '7 seconds'.
%           % "code" will be 0 on success, 1 if file was longer than 8sec
%           % and 2 if file is shorter than 8sec -- in which case, the
%           % request cannot be accomodated!


persistent data1;
persistent data2;

% Don't need 'fn' to get 'fps' value.  Return value from whichever
% structure is defined.  Must always be equal, so doesn't matter.
if strcmp(control,'fps') && ~exist('fn'),
    if exist('data1') && isfield(data1,'fps'),
        one = data1.fps;
    elseif exist('data2') && isfield(data2,'fps'),
        one = data2.fps;
    else
        error('function dll_video: no ''fps'' to retreive');
    end
    return;
end


% otherwise, always need 'fn' value, to say if working from source or
% processed video. 
if fn == 1,
    if exist('data1','var')
        data = data1;
    end
elseif fn == 2,
    if exist('data2','var')
        data = data2;
    end
else
    error('function dll_video: fn must be either 1 or 2');
end

%
if strcmp(control, 'initialize'),
    [data] = dll_video__initialize(varargin);
    data.sec7 = 0;
    
elseif strcmp(control,'exist'),
    if exist('data','var'),
        one = 1;
    else
        one = 0;
    end
    
elseif strcmp(control, 'set_tslice'),
    [data.tslice_frames, data.over_sec] = tslice_conversion(varargin{1}, data.fps);

elseif strcmp(control,'size'),
    one = data.rows;
    two = data.cols;
    three = data.fps;
    four = data.durration;
    % if demanding 7-second length, return that length instead of actual
    % length.
    if data.sec7 == 1,
        four = 7;
    end

    
elseif strcmp(control,'fps'),
    one = data.fps;
    
elseif strcmp(control,'get_video_standard'),
    one = data.video_standard;
    
elseif strcmp(control,'set_rewind'),
    data.rewind_curr_frame = data.curr_frame;
    data.rewind_overby = data.overby;
    
elseif strcmp(control,'rewind'),
    data.curr_frame = data.rewind_curr_frame;
    data.overby = data.rewind_overby;
    
elseif strcmp(control,'discard'),
    [tslice_frames, over_sec] = tslice_conversion(varargin{1}, data.fps);
    data.curr_frame = data.curr_frame + tslice_frames;
    data.overby = data.overby + over_sec;
    if data.overby >= 1.0,
        data.curr_frame = data.curr_frame - 1;
        data.overby = data.overby - 1.0;
    end
    
elseif strcmp(control, 'tslice'),
    if length(varargin) == 0,
        [one, two, three, data] = dll_video__next_tslice(data);
    elseif length(varargin) == 1,
        [one, two, three, data] = dll_video__next_tslice(data, varargin{1});
    else
        error('not sure what the second argument is!?!');
    end
    
elseif strcmp(control,'peek'),
    [one] = dll_video__peek_tslice(data, varargin{1}, varargin{2});

elseif strcmp(control, 'sec'),
    if nargout == 1,
        [one, data] = dll_video__next_tslice(data, varargin{1}, varargin{2}, 0);
    else
        [one, two, three, data] = dll_video__next_tslice(data, varargin{1}, varargin{2}, 1);
    end
    
elseif strcmp(control, 'total_frames'),
    one = data.total_frames - data.curr_frame + 1;
    
elseif strcmp(control, 'total_sec'),
    one = (data.total_frames - data.curr_frame + 1) / data.fps;

    
elseif strcmp(control, 'delay_8s'),
    delay_value = varargin{1};

    % can't have delays beyond +/- 0.2
    % leave one extra at end for re-framing
    delay_value = min(round(0.2 * data.fps)-1, delay_value);
    delay_value = max(round(-0.2 * data.fps), delay_value);

    
    % initialize processed video segment to be used: 7 sec + 1frame for
    % re-framing
    data.curr_frame = round(0.8 * data.fps) + 1 + delay_value;
    data.rewind_curr_frame = data.curr_frame;
    data.overby = 0;
    data.sec7 = 1;

    % set return code.  See help above.
    if data.durration < 8,
        one = 2;
    elseif data.durration > 8,
        one = 1;
    else
        one = 0;
    end

else
    error('function dll_video:  control value not recognized');
end

% copy back to persistent variable
if fn == 1,
    data1 = data;
elseif fn == 2,
    data2 = data;
end

%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = dll_video__initialize(list);
% Organize all of the data that will be needed.

data.file_name = list{1};
data.file_type = list{2};
data.video_standard = list{3};

data.curr_frame = 1;
data.overby = 0;

if strcmp(data.file_type,'avi'),
    [info] = aviinfo(data.file_name);
    data.rows = info.Height;
    data.cols = info.Width;
    data.fps = info.FramesPerSecond;
    data.total_frames = info.NumFrames;
    
elseif strcmp(data.file_type,'uyvy'),
    data.rows = list{4};
    data.cols = list{5};
    data.fps = list{6};
    fid = fopen(data.file_name,'r');
    fseek(fid,0, 'eof');
    data.total_frames = ftell(fid)./(data.rows*data.cols*2);
    fclose(fid);
    
else
    error('file type not recognized.  Use ''avi'' or ''uyvy'' ');
end

data.durration = data.total_frames / data.fps;

%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = dll_video__peek_tslice(data, reframe, durration);

if reframe ~=1 & reframe ~= 0,
    error('CALL WRONG -- go back & add reframing argument');
end

[tslice_frames, over_sec] = tslice_conversion(durration, data.fps);
start = data.curr_frame;
stop = data.curr_frame + tslice_frames + reframe - 1;
stop = min(stop, data.total_frames); % don't try to read past end of file!

if strcmp(data.file_type,'avi'),
    [y] = read_avi('YCbCr',data.file_name, '128', ...
        'frames',start, stop);
elseif strcmp(data.file_type,'uyvy'),
    [y] = read_bigyuv(data.file_name, 'frames', start, stop,  ...
        '128','size',data.rows,data.cols);
end


%%%%%%%%%%%%%%%%%%%%%%%%
function [y, cb, cr, data ] = dll_video__next_tslice( data, reframe, durration, get_color);

% by default, get color & luminance
if nargin < 4,
    get_color = 1;
end
if nargin == 1,
    reframe = 0;
end
if reframe ~=1 & reframe ~= 0,
    error('CALL WRONG -- go back & add reframing argument');
end

% if optional 'durration' argument exists, override the default values for
% data.tslice_frames and data.over_sec.
if exist('durration', 'var'),
    [tslice_frames, over_sec] = tslice_conversion(durration, data.fps);
else
    tslice_frames = data.tslice_frames;
    over_sec = data.over_sec;
end

if tslice_frames == 0,
    error('Request made to read zero (0) frames of video.');
end

if strcmp(data.file_type,'avi'),
    if get_color,
        [y,cb,cr] = read_avi('YCbCr',data.file_name, '128', ...
            'frames',data.curr_frame, data.curr_frame + tslice_frames +reframe - 1);
    else
        [y] = read_avi('YCbCr',data.file_name, '128', ...
            'frames',data.curr_frame, data.curr_frame + tslice_frames + reframe - 1);
        cb=0;
        cr=0;
    end
elseif strcmp(data.file_type,'uyvy'),
	% read video from file.
	if get_color,
        [y,cb,cr] = read_bigyuv(data.file_name, 'frames', data.curr_frame, data.curr_frame + tslice_frames + reframe - 1,  ...
            '128','size',data.rows,data.cols);
	else
        [y] = read_bigyuv(data.file_name, 'frames', data.curr_frame, data.curr_frame + tslice_frames + reframe - 1,  ...
            '128','size',data.rows,data.cols);
        cb=0;
        cr=0;
	end
end

% handle overlap
data.curr_frame = data.curr_frame + tslice_frames;
data.overby = data.overby + over_sec;
if data.overby >= 1.0,
    data.curr_frame = data.curr_frame - 1;
    data.overby = data.overby - 1.0;
end

% convert to double -- everything!
y = double(y);

if get_color,
    cb = double(cb);
    cr = double(cr);
end

if ~get_color,
    cb = data;
end