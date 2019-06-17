function [o1, o2, o3, o4, o5] = read_avi( ret_val, varargin )
% readAvi
%   Reads an uncompressed AVI file of a variety of formats, including
%   the following:
%       10-bit : uyvy : yuy2 : yv12 : rgb
%   If FILENAME does not include an extension, then '.avi' will be used.
% SYNTAX
%   [info] = read_avi('Info',filename);
%   [c1,c2,c3] = read_avi(color_out,filename);
%   [...] = read_avi(...,'flag',...);
% DESCRIPTION
%   [info] = read_avi('Info',filename);
%       returns only a struct, containing information about the file.
%       This struct can then be passed as an argument to this function
%       preceeded by the 'struct' flag, and the file will be read.
%   [c1,c2,c3] = read_avi(color_out,filename);
%       returns the color components of the frames read in from the file.
%       The color components will depend on color_out. If no frames are
%       requested in particular, only the first frame is read.
%   [...] = read_avi(...,'Flag',...);
%       designates a flag to be set inside the function. See below for a
%       complete list of possible flags.
% INPUT ARGUMENTS
%   color_out >
%   'Info'  -- return listing of header information (see aviinfo).
%   'RGB'   -- return image planes in the RGB colorspace.
%   'YCbCr' -- return image planes in the YCbCr colorspace.
%
%   filename >
%   Avi file to be opened. If the 'struct' flag has
%   already been given, do NOT also give a filename
%
%   flag >
%   'struct', avi_struct    A struct returned by aviinfo or by this
%                           function with the 'Info' property.
%   'sroi',top,left,bottom,right,   Spatial region of intrest.  By 
%                                   default, all of each image is 
%                                   returned.
%   'frames',start,stop,    Specify the first and last frames, 
%                           inclusive, to be read ('start' and 'stop').
%                           By default, read first frame.
%   '128'               Subtract 128 from all Cb and Cr values.  By 
%                       default, Cb and Cr values are left in the 
%                       [0..255] range.
%   'interp'            Linearly interpolate Cb and Cr values.  By default, 
%                       color planes are pixel replicated.  Note:  
%                       Interpolation is slow. Only implemented for
%                       YUV colorspaces excepting YV12.
%   'audio',['frames' or 'seconds'], start, stop
%               Request audio be returned if it exists. If 'frames' are
%               requested, the audio for the given frames [1..NumFrames]
%               will be returned. If 'seconds' are requested, audio with
%               the given duration [0..TotalTime) will be returned. Feel
%               free to request more than is in the file; I handle it :)
% OUTPUT ARGUMENTS
%   c1 >
%   Depending on the color_out property, could be Info if 'Info', Y if
%   'YCbCr', or R if 'RGB'.
%
%   c2, c3 >
%   Depending on the color_out property, could be Cb and Cr if 'YCbCr'
%   or G and B if 'RGB'.
%
%   c4 >
%   If audio is requested and exists, this is the raw audio data,
%   separated by channels.
%
%   c5 >
%   if audio is requested and exists, this is the Audio Rate.
% EXAMPLES
%---[info] = read_avi('Info','twocops.avi');
%---[r,g,b] = read_avi('RGB','twocops.avi','frames',1,30);
%---[y,cb,cr] = read_avi('YCbCr','twocops.avi','frames',61,90,'128');
%---info = aviinfo('my.avi');
%   [r, g, b] = read_avi('RGB', 'struct', info);
%---[y,cb,cr,aud,rate] = read_avi('YCbCr','my.avi','audio','seconds',0,5);
% NOTES
%   When reading files with the YV12 fourcc, the cb and cr color
%   components will be extrapolated to fit the Y component matrix size.
%   The current extrapolation algorithm simply copies the cb and cr
%   values. A better implementation might include a bi-linear
%   interpolation.

% SIGNATURE
%   Programmer: Zebulon Fross
%   Version:    08/10/2010
%

% Initialization
is_whole_image  =  1;
frame_start     =  1;
frame_stop      =  1;
audio_start     =  0;
audio_stop      =  0;
used_frames     =  0;
sroi            = [];
is_sub128       =  0;
is_interp       =  0;
ret_info        =  0;
ret_ycbcr       =  1;
% Parse return information type flag
if strcmpi(ret_val,'info')
    ret_info = 1;
elseif strcmpi(ret_val,'RGB')
    ret_ycbcr = 0;
elseif ~strcmpi(ret_val,'YCbCr')
    error('Return type flag not recognized');
end

persistent info;

% Validate input/output.
error(nargoutchk(0,5,nargout));
error(nargchk(2,19,nargin));
try
    cnt=1;
    while cnt <= length(varargin),
        if ~ischar(varargin{cnt}),
            error('parameter not recognized');
        end
        if strcmpi(varargin(cnt),'struct') == 1,
            info = varargin{cnt+1};
            cnt = cnt + 2;
        elseif strcmpi(varargin(cnt),'sroi') == 1,
            sroi.top = varargin{cnt+1};
            sroi.left = varargin{cnt+2};
            sroi.bottom = varargin{cnt+3};
            sroi.right = varargin{cnt+4};
            is_whole_image = 0;
            cnt = cnt + 5;
        elseif strcmpi(varargin(cnt),'frames') == 1,
            frame_start = varargin{cnt+1};
            frame_stop = varargin{cnt+2};
            cnt = cnt + 3;
        elseif strcmp(varargin(cnt),'128') == 1,
            is_sub128 = 1;
            cnt = cnt + 1;
            if ~ret_ycbcr,
                error('RGB and ''128'' flag are incompatible');
            end
        elseif strcmpi(varargin(cnt),'interp') == 1,
            is_interp = 1;
            cnt = cnt + 1;
        elseif strcmpi(varargin(cnt),'audio') == 1,
            cnt = cnt + 1;
            if strcmpi(varargin(cnt),'frames') == 1,
                used_frames = 1;
            end
            audio_start = varargin{cnt+1};
            audio_stop  = varargin{cnt+2};
            cnt = cnt + 3;
        else
            % assume file name is given
            filename = varargin{cnt};
            [~,~,ext] = fileparts(filename);
            if isempty(ext)
                filename = strcat(filename,'.avi');
            end
            new_dir = dir(filename);
            if isempty(info) || ...
                    ~strcmp(info.Filename, filename) || ...
                    ~strcmp(info.FileModDate, new_dir.date) || ...
                    ~eq(info.FileSize, new_dir.bytes)
                info = aviinfo(filename); %#ok<FREMO>
            end
            cnt = cnt + 1;
        end
    end
catch e
    error(e.identifier, ...
        'Unable to parse input arguments. Please check calling syntax.');
end

% if frames were given for audio, convert frames to seconds
if used_frames == 1,
    audio_start = (audio_start-1)/info.FramesPerSecond;
    audio_stop  = audio_stop/info.FramesPerSecond;
end
% check for invalid audio request
if (audio_stop < audio_start)
    error('MATLAB:readavialt', 'invalid audio chunk request');
end

% TODO: Is this correct?
if ispc
    file = fopen(info.Filename, 'r', 'l');
else
    file = fopen(info.Filename, 'r', 'b');
end
assert(file >= 0);

if ret_info == 1
    o1 = info;
    return ;
end

% ensure frame request is valid
if frame_start <= 0 || frame_start > info.NumFrames
    error('frame numbers to read must be valid');
end
if frame_stop < frame_start || frame_stop > info.NumFrames
    error('frame numbers to read must be valid');
end

if is_whole_image,
    o1=zeros(info.Height,...
        info.Width,...
        frame_stop-frame_start+1,'single');
else
    o1=zeros(sroi.bottom-sroi.top+1,...
        sroi.right-sroi.left+1,...
        frame_stop-frame_start+1,'single');
end
if nargout > 1
    o2=o1;
    o3=o1;
    o4=[];
    o5= 0;
end

%% gather video data -------------------------------------------------
out_pos = 1;
% read in the requested frames
for ind = frame_start:frame_stop
    % seek to position of frame
    fseek(file, info.vidFrames(1, ind), -1);
    
    if ( strcmp(info.ColorType, 'UYVY') )
        if (strcmpi(info.Codec, 'UYVY') || ...
                strcmpi(info.Codec, 'ffds') || ...
                strcmpi(info.Codec, 'HDYC') || ...
                strcmpi(info.Codec, 'DIB '))
            [y,cb,cr] = read_uyvy_frame(file, is_whole_image, ...
                sroi, is_interp, ...
                info.Height, info.Width);
            %                 % One fluke AVI file had upsidedown images, and so needed
            %                 % images to be flipped. None of the ITS datasets have this
            %                 % problem, so this appears to have been a mistake.
            %                 if (strcmpi(info.Codec, 'DIB '))
            %                 % I am still not sure how to interpret DIB files
            %                    y  = flipud( y);
            %                    cb = flipud(cb);
            %                    cr = flipud(cr);
            %                 end
        elseif (strcmpi(info.Codec, 'v210'))
            [y, cb, cr] = ...
                read_10bit_uyvy(file, is_whole_image, ...
                sroi, is_interp, ...
                info.Height, info.Width);
        elseif (strcmpi(info.Codec, 'yuy2'))
            [y,cb,cr] = read_yuyv_frame(file, is_whole_image, ...
                sroi, is_interp, ...
                info.Height, info.Width);
        elseif (strcmpi(info.Codec, 'yv12'))
            [y, cb, cr] = read_yv12_frame(file, is_whole_image, ...
                sroi, is_interp, ...
                info.Height, info.Width);
            % if the codec is not one of the previous, we will resort to
            % the compressed file reader
        else
            % this piece of code is experimental for reading
            % compressed avi files. The mex file does not work on
            % 64-bit machines. On 32-bit machines, the program
            % returns slightly-faulty data.
            error('could not interpet avi codec');
            %                 try
            %                     eStr = ['The mex function could not read this AVI ' ...
            %                         'file. This could be because the mex function' ...
            %                         ' is not available or because the file''s ' ...
            %                         'codec is not installed on this system.'];
            %                     % determine the correct mex function to call
            %                     if ( strcmp(computer, 'PCWIN') )
            %                         [r, g, b] = read_compressed_avi_32(info.Filename, ...
            %                             ind, ...
            %                             info.Height, ...
            %                             info.Width);
            %                     else
            %                         eStr = ['Your operating system is not supported by' ...
            %                             ' our mex function. Please run this again ' ...
            %                             'on a 32-bit Windows machine.'];
            %                         error(eStr);
            %                     end
            %                     [y, cb, cr] ...
            %                         = rgb2ycbcr_double(double( r ), ...
            %                         double( g ), ...
            %                         double( b ));
            %                 catch e
            %                     rethrow(e);
            %                 end
        end
        % append y, cb, and cr to output variables
        if ret_ycbcr,
            o1(:,:,out_pos) = y;
            if nargout > 1,
                if is_sub128
                    cb = single(cb)-128;
                    cr = single(cr)-128;
                end
                o2(:,:,out_pos) = cb;
                o3(:,:,out_pos) = cr;
            end
        else
            [o1(:,:,out_pos),o2(:,:,out_pos),o3(:,:,out_pos)] ...
                = ycbcr2rgb_double(single(y),...
                single(cb),...
                single(cr));
        end
    elseif strcmp(info.ColorType, 'RGB'),
        if info.BitDepth == 24,
            [r,g,b] = read_rgb24_frame(file, is_whole_image, ...
                sroi, info.Height, ...
                info.Width);
        elseif info.BitDepth == 32,
            [r,g,b] ...
                = read_rgb32_frame(file, is_whole_image, sroi, ...
                info.Height, ...
                info.Width);
            % if the codec is not one of the previous, we will resort to
            % the compressed file reader
        else
            % this piece of code is experimental for reading
            % compressed avi files. The mex file does not work on
            % 64-bit machines. On 32-bit machines, the program
            % returns slightly-faulty data.
            error('could not interpet avi codec');
            %                 try
            %                     eStr = ['The mex function could not read this AVI ' ...
            %                         'file. This could be because the mex function' ...
            %                         ' is not available or because the file''s ' ...
            %                         'codec is not installed on this system.'];
            %                     % determine the correct mex function to call
            %                     if ( strcmp(computer, 'PCWIN') )
            %                         [r, g, b] = read_compressed_avi_32(info.Filename, ...
            %                             ind, ...
            %                             info.Height, ...
            %                             info.Width);
            %                     else
            %                         eStr = ['Your operating system is not supported by' ...
            %                             ' our mex function. Please run this again ' ...
            %                             'on a 32-bit Windows machine.'];
            %                         error(eStr);
            %                     end
            %                 catch e
            %                     rethrow(e);
            %                 end
        end
        
        if ~ret_ycbcr,
            o1(:,:,out_pos) = r;
            o2(:,:,out_pos) = g;
            o3(:,:,out_pos) = b;
        else
            [y, cb, cr] = rgb2ycbcr_double(r,g,b);
            o1(:,:,out_pos) = y;
            if nargout > 1,
                if is_sub128
                    cb = single(cb)-128;
                    cr = single(cr)-128;
                end
                o2(:,:,out_pos) = cb;
                o3(:,:,out_pos) = cr;
            end
        end
    else
        error('unsupported input format');
    end
    out_pos = out_pos + 1;
end

%% only return audio if requested and if it exists -------------------
if audio_stop > 0 && isfield(info, 'audFrames')
    data = [];
    
    % determine the datatype to read in
    BytesPerSample = info.BitsPerSample/8;
    if (BytesPerSample == 1),
        dtype='uchar'; % unsigned 8-bit
    elseif (BytesPerSample == 2),
        dtype='int16'; % signed 16-bit
    elseif (BytesPerSample == 3)
        dtype='bit24'; % signed 24-bit
    elseif (BytesPerSample == 4),
        % 32-bit 16.8 float (type 1 - 32-bit)
        if (strcmpi(info.AudioFormat, 'Format # 0x1'))
            dtype = 'bit32'; %signed 32-bit
            % 32-bit normalized floating point
        elseif (strcmpi(info.AudioFormat, 'Format # 0x3'))
            dtype = 'float'; % floating point
        elseif (strcmpi(info.AudioFormat, 'Format # 0xFFFE'))
            %32 Bit data with either integer or floating point numbers. Use
            %info.SubFormat to determine whether or not the audio data is
            %32 bit or 32 bit floating point.  Also, a check will still be
            %in place below just incase this information does not exist or
            %is something "weird".
            if(info.SubFormat == 1)
                %PCM data is contained, meaning that the data is 32 bit,
                %not floating point.
                dtype = 'bit32'; %signed 32-bit
            elseif(info.SubFormat == 3)
                %IEEE floating point data is contained, meaning that the
                %data is 32 bit floating point numbers.
                dtype = 'float'; % floating point
            else
                %The SubFormat is formatted in a way that isnt checked by
                %this program, it will be assumed that 32 bit will be used.
                dtype = 'bit32'; %signed 32-bit
                warning('SubFormat not formatted in a way this program understands.  Audio may be incorrect.');
            end
        else 
             dtype = 'bit32';
        end
    end

    skip_samples = uint32(audio_start*info.BytesPerSec/BytesPerSample);
    stop_samples = uint32(audio_stop*info.BytesPerSec/BytesPerSample);
    total_samples = (stop_samples - skip_samples)/info.NumAudioChannels;
    sample_count = 0;
    if sample_count < skip_samples
        ind = 0;
    else
        ind = 1;
    end
    % determine where to start reading the audio by skipping parts
    total_samples_thus_far = 0;
    while sample_count < skip_samples
        ind = ind + 1;
        sz = info.audFrames(2, ind);
        % Total samples in a chunk = size of chunk / Bytes per sample
        TotalSamples = floor(sz/BytesPerSample);
        sample_count = sample_count + TotalSamples;
        total_samples_thus_far = total_samples_thus_far + TotalSamples;
    end
    % determine how much needs to be skipped in the initial chunk
    % seek to initial chunk and skip necessary samples
    
    if(ind == 1)
        fseek(file, info.audFrames(1, ind), -1);
        fseek(file, floor(double(skip_samples)*BytesPerSample), 0);
        
        if(info.NumAudioChannels > 2)
            %Get a Reference Size to know how many samples you need to
            %actually have.
            
            data = [data fread(file, [info.NumAudioChannels, ...
                floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                /(BytesPerSample*info.NumAudioChannels))], dtype)];
            samples_needed = size(data,2);
            %Clear Data
            data = [];
            
            %Now grab all the data in the full chunk and chop off the not
            %needed data in the begining.
            fseek(file, info.audFrames(1, ind), -1);
            data = [data fread(file, [info.NumAudioChannels, ...
                floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                /(BytesPerSample*info.NumAudioChannels))], dtype)];
            data = data(:,(size(data,2)-samples_needed+1):size(data,2));
            if(size(data,2) ~= samples_needed)
                %These two HAVE to be equal for the seeking through the
                %audio to be correct.
                error('Seeking through Audio data has failed!');
            end
            ind = ind + 1;
            try
                fseek(file, info.audFrames(1, ind), -1);
            catch
                %No more audio in the file, it was all read in.  Thats
                %fine, this will be taken care of down below.
            end
        end
    else
        %Re-adjust skip_samples since we are not starting at the beginning
        %of the clip.
        skip_samples_in_file = skip_samples - (total_samples_thus_far - TotalSamples);
        if(skip_samples < 0)
            error('Skip Samples was not calculated correctly')
        end
        fseek(file, info.audFrames(1, ind), -1);
        fseek(file, floor(double(skip_samples_in_file)*BytesPerSample), 0);
        
        %Ok, so when 6 channel audio is present, seeking in the file messes
        %up the order of the audio channels for the first sample that is
        %read in.  Therefore, this program will read the whole chunk in and
        %then only store the number of samples that is needed.
        if(info.NumAudioChannels > 2)
            %Get a Reference Size to know how many samples you need to
            %actually have.
            
            data = [data fread(file, [info.NumAudioChannels, ...
                floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                /(BytesPerSample*info.NumAudioChannels))], dtype)];
            samples_needed = size(data,2);
            %Clear Data
            data = [];
            
            %Now grab all the data in the full chunk and chop off the not
            %needed data in the begining.
            fseek(file, info.audFrames(1, ind), -1);
            data = [data fread(file, [info.NumAudioChannels, ...
                floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                /(BytesPerSample*info.NumAudioChannels))], dtype)];
            data = data(:,(size(data,2)-samples_needed+1):size(data,2));
            if(size(data,2) ~= samples_needed)
                %These two HAVE to be equal for the seeking through the
                %audio to be correct.
                error('Seeking through Audio data has failed!');
            end
            ind = ind + 1;
            try
                fseek(file, info.audFrames(1, ind), -1);
            catch
                %No more audio in the file, it was all read in.  Thats
                %fine, this will be taken care of down below.
            end
        end
    end
    
    % while we haven't read in all that has been requested
    while size(data, 2) < total_samples
        if(isempty(data))
            %Dividing has to do with staying within the audio region!
            %Otherwise, it goes outside of the audio region and noise
            %is added.
            data = [data fread(file, [info.NumAudioChannels, ...
                floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                /(BytesPerSample*info.NumAudioChannels))], dtype)];
        else
            data = [data fread(file, [info.NumAudioChannels, ...
                    floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                    /(BytesPerSample*info.NumAudioChannels))], dtype)];
        end

        ind = ind + 1;
        % break out if there is no more audio to read. We won't
        % penalize for users requesting too much.
        if (ind > size(info.audFrames, 2))
            break;
        end
        fseek(file, info.audFrames(1, ind), -1);
    end
    
    %Chop off extra data before scaling occurs.  This has been moved from
    %the bottom of this function to here.  This is because if the program
    %read to much and left the audio section of the data, noise enters the
    %data.  This noise needs to be removed before the data can be
    %normalized.  Thus, this is what is accomplished here.
    if size(data, 2) < total_samples
        total_samples = size(data, 2);
    end
    % truncate the output just in case we read too much
    data = data(:, 1:total_samples);
    
    data = data';
    
    % Normalize data range: min will hit -1, max will not quite hit +1.
    if BytesPerSample==1,
        data = (data-128)/128;  % [-1,1)
    elseif BytesPerSample==2,
        data = data/32768;      % [-1,1)
    elseif BytesPerSample==3,
        data = data/(2^23);     % [-1,1)
    elseif BytesPerSample==4,
        % Type 3 32-bit is already normalized
        if(~strcmpi(info.AudioFormat, 'Format # 0x3') && info.SubFormat ~= 3)
           data = data/(2^31); % [-1,1)
        end
    end
    
    %Only needed for 32 bit audio - This checks if the 32 bit audio should
    %be int32 (default) or if it should be float 32 bit.  While this can
    %be known from the aviinfo function (info.AudioFormat = 0x3) it is
    %unknown if the audio format is WAVE_FORMAT_EXTENSIBLE(0xfffe).
    %Therefore, this needs to be checked, this function is resposible for
    %checking.  This also provides a "check" for info.SubFormat, just
    %incase it was not formatted or not formatted correctly.
    if(strcmpi(info.AudioFormat, 'Format # 0xFFFE') && (strcmpi(dtype,'bit32') || strcmpi(dtype,'float')))
        if(max(max(isnan(data))) == 1) %The assumed format was not correct.
            %32 bit float will be used
            
            skip_samples = uint32(audio_start*info.BytesPerSec/BytesPerSample);
            stop_samples = uint32(audio_stop*info.BytesPerSec/BytesPerSample);
            total_samples = (stop_samples - skip_samples)/info.NumAudioChannels;
            sample_count = 0;
            if sample_count < skip_samples
                ind = 0;
            else
                ind = 1;
            end
            % determine where to start reading the audio by skipping parts
            total_samples_thus_far = 0;
            while sample_count < skip_samples
                ind = ind + 1;
                sz = info.audFrames(2, ind);
                % Total samples in a chunk = size of chunk / Bytes per sample
                TotalSamples = floor(sz/BytesPerSample);
                sample_count = sample_count + TotalSamples;
                total_samples_thus_far = total_samples_thus_far + TotalSamples;
            end
            % determine how much needs to be skipped in the initial chunk
            % seek to initial chunk and skip necessary samples
            
            if(strcmpi(dtype,'bit32'))
                data = [];
                dtype = 'float';
            else
                data = [];
                dtype = 'bit32';
            end
            
            if(ind == 1)
                fseek(file, info.audFrames(1, ind), -1);
                fseek(file, floor(double(skip_samples)*BytesPerSample), 0);
                
                if(info.NumAudioChannels > 2)
                    %Get a Reference Size to know how many samples you need to
                    %actually have.
                    
                    data = [data fread(file, [info.NumAudioChannels, ...
                        floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                        /(BytesPerSample*info.NumAudioChannels))], dtype)];
                    samples_needed = size(data,2);
                    %Clear Data
                    data = [];
                    
                    %Now grab all the data in the full chunk and chop off the not
                    %needed data in the begining.
                    fseek(file, info.audFrames(1, ind), -1);
                    data = [data fread(file, [info.NumAudioChannels, ...
                        floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                        /(BytesPerSample*info.NumAudioChannels))], dtype)];
                    data = data(:,(size(data,2)-samples_needed+1):size(data,2));
                    if(size(data,2) ~= samples_needed)
                        %These two HAVE to be equal for the seeking through the
                        %audio to be correct.
                        error('Seeking through Audio data has failed!');
                    end
                    ind = ind + 1;
                    try
                        fseek(file, info.audFrames(1, ind), -1);
                    catch
                        %No more audio in the file, it was all read in.  Thats
                        %fine, this will be taken care of down below.
                    end
                end
            else
                %Re-adjust skip_samples since we are not starting at the beginning
                %of the clip.
                skip_samples_in_file = skip_samples - (total_samples_thus_far - TotalSamples);
                if(skip_samples < 0)
                    error('Skip Samples was not calculated correctly')
                end
                fseek(file, info.audFrames(1, ind), -1);
                fseek(file, floor(double(skip_samples_in_file)*BytesPerSample), 0);
                
                %Ok, so when 6 channel audio is present, seeking in the file messes
                %up the order of the audio channels for the first sample that is
                %read in.  Therefore, this program will read the whole chunk in and
                %then only store the number of samples that is needed.
                if(info.NumAudioChannels > 2)
                    %Get a Reference Size to know how many samples you need to
                    %actually have.
                    
                    data = [data fread(file, [info.NumAudioChannels, ...
                        floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                        /(BytesPerSample*info.NumAudioChannels))], dtype)];
                    samples_needed = size(data,2);
                    %Clear Data
                    data = [];
                    
                    %Now grab all the data in the full chunk and chop off the not
                    %needed data in the begining.
                    fseek(file, info.audFrames(1, ind), -1);
                    data = [data fread(file, [info.NumAudioChannels, ...
                        floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                        /(BytesPerSample*info.NumAudioChannels))], dtype)];
                    data = data(:,(size(data,2)-samples_needed+1):size(data,2));
                    if(size(data,2) ~= samples_needed)
                        %These two HAVE to be equal for the seeking through the
                        %audio to be correct.
                        error('Seeking through Audio data has failed!');
                    end
                    ind = ind + 1;
                    try
                        fseek(file, info.audFrames(1, ind), -1);
                    catch
                        %No more audio in the file, it was all read in.  Thats
                        %fine, this will be taken care of down below.
                    end
                end
            end
            
            % while we haven't read in all that has been requested
            while size(data, 2) < total_samples
                if(isempty(data))
                    %Dividing has to do with staying within the audio region!
                    %Otherwise, it goes outside of the audio region and noise
                    %is added.
                    data = [data fread(file, [info.NumAudioChannels, ...
                        floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                        /(BytesPerSample*info.NumAudioChannels))], dtype)];
                else
                    data = [data fread(file, [info.NumAudioChannels, ...
                        floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                        /(BytesPerSample*info.NumAudioChannels))], dtype)];
                end
                
                ind = ind + 1;
                % break out if there is no more audio to read. We won't
                % penalize for users requesting too much.
                if (ind > size(info.audFrames, 2))
                    break;
                end
                fseek(file, info.audFrames(1, ind), -1);
            end
            
            if size(data, 2) < total_samples
                total_samples = size(data, 2);
            end
            % truncate the output just in case we read too much
            data = data(:, 1:total_samples);
            
            data = data';
            
            if(strcmpi(dtype,'bit32'))
                data = data/(2^31); % [-1,1)
            end
        else %Check the histogram of the audio output.  If gaussian in 
             %shape then everything is ok.  If two big peaks, the wrong
             %format was used.  32 bit float will be used.
            two_peaks = 0;
            peak_left = 0;
            peak_right = 0;
           
            %Obtain the historgram
            max_chan = [0,0];
            for i = 1:info.NumAudioChannels
                %Make sure the channel is not just all zeros
                temp = max(max(data(:,i)));
                if(temp > max_chan(1,1))
                    max_chan(1,1) = temp;
                    max_chan(1,2) = i;
                end
            end
             
            clear temp
            temp = histc(data(:,max_chan(1,2)),[min(min(data)):.01:max(max(data))]);
            
            %check if two peaks exist.
            %Split the historgram into two parts (the peaks will be on
            %either end of the spectrum if the wrong format was used.
            temp1 = temp(1:round((size(temp,1))/4));
            temp2 = temp(round(size(temp,1)*(3/4)):size(temp));
            temp3 = temp(round(size(temp,1)*(1/4)):round(size(temp,1)*(3/4)));
            
            %Set middle peak value for comparison
            middle_peak = max(temp3);
            
            %If a side peak is higher than the middle peak or is within 50%
            %of the middle_peak in value, more than one peak exists.
            if(max(temp1) > middle_peak || max(temp1) > (middle_peak - floor(middle_peak*.5)))
                %A peak on the left side exists
                peak_left = 1;
            end
            
            if(max(temp2) > middle_peak || max(temp2) > (middle_peak - floor(middle_peak*.5)))
                %A peak on the right side exists
                peak_right = 1;
            end
            
            if(peak_left == 1 && peak_right == 1)
                %The wrong format was used!
                two_peaks = 1;
                display('32 Bit Floating Point Audio Has Been Detected');
            elseif(peak_left == 1 || peak_right == 1)
                %Only one peak exists which is strange
                warning('Only one peak detected and its not centered!  Audio may not be correct!');
            end
            
            if(two_peaks == 1)
                %Audio format was wrong, use 32 bit floating point format.
                
                skip_samples = uint32(audio_start*info.BytesPerSec/BytesPerSample);
                stop_samples = uint32(audio_stop*info.BytesPerSec/BytesPerSample);
                total_samples = (stop_samples - skip_samples)/info.NumAudioChannels;
                sample_count = 0;
                if sample_count < skip_samples
                    ind = 0;
                else
                    ind = 1;
                end
                % determine where to start reading the audio by skipping parts
                total_samples_thus_far = 0;
                while sample_count < skip_samples
                    ind = ind + 1;
                    sz = info.audFrames(2, ind);
                    % Total samples in a chunk = size of chunk / Bytes per sample
                    TotalSamples = floor(sz/BytesPerSample);
                    sample_count = sample_count + TotalSamples;
                    total_samples_thus_far = total_samples_thus_far + TotalSamples;
                end
                % determine how much needs to be skipped in the initial chunk
                % seek to initial chunk and skip necessary samples
                
                if(strcmpi(dtype,'bit32'))
                    data = [];
                    dtype = 'float';
                else
                    data = [];
                    dtype = 'bit32';
                end
                
                if(ind == 1)
                    fseek(file, info.audFrames(1, ind), -1);
                    fseek(file, floor(double(skip_samples)*BytesPerSample), 0);
                    
                    if(info.NumAudioChannels > 2)
                        %Get a Reference Size to know how many samples you need to
                        %actually have.
                        
                        data = [data fread(file, [info.NumAudioChannels, ...
                            floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                            /(BytesPerSample*info.NumAudioChannels))], dtype)];
                        samples_needed = size(data,2);
                        %Clear Data
                        data = [];
                        
                        %Now grab all the data in the full chunk and chop off the not
                        %needed data in the begining.
                        fseek(file, info.audFrames(1, ind), -1);
                        data = [data fread(file, [info.NumAudioChannels, ...
                            floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                            /(BytesPerSample*info.NumAudioChannels))], dtype)];
                        data = data(:,(size(data,2)-samples_needed+1):size(data,2));
                        if(size(data,2) ~= samples_needed)
                            %These two HAVE to be equal for the seeking through the
                            %audio to be correct.
                            error('Seeking through Audio data has failed!');
                        end
                        ind = ind + 1;
                        try
                            fseek(file, info.audFrames(1, ind), -1);
                        catch
                            %No more audio in the file, it was all read in.  Thats
                            %fine, this will be taken care of down below.
                        end
                    end
                else
                    %Re-adjust skip_samples since we are not starting at the beginning
                    %of the clip.
                    skip_samples_in_file = skip_samples - (total_samples_thus_far - TotalSamples);
                    if(skip_samples < 0)
                        error('Skip Samples was not calculated correctly')
                    end
                    fseek(file, info.audFrames(1, ind), -1);
                    fseek(file, floor(double(skip_samples_in_file)*BytesPerSample), 0);
                    
                    %Ok, so when 6 channel audio is present, seeking in the file messes
                    %up the order of the audio channels for the first sample that is
                    %read in.  Therefore, this program will read the whole chunk in and
                    %then only store the number of samples that is needed.
                    if(info.NumAudioChannels > 2)
                        %Get a Reference Size to know how many samples you need to
                        %actually have.
                        
                        data = [data fread(file, [info.NumAudioChannels, ...
                            floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                            /(BytesPerSample*info.NumAudioChannels))], dtype)];
                        samples_needed = size(data,2);
                        %Clear Data
                        data = [];
                        
                        %Now grab all the data in the full chunk and chop off the not
                        %needed data in the begining.
                        fseek(file, info.audFrames(1, ind), -1);
                        data = [data fread(file, [info.NumAudioChannels, ...
                            floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                            /(BytesPerSample*info.NumAudioChannels))], dtype)];
                        data = data(:,(size(data,2)-samples_needed+1):size(data,2));
                        if(size(data,2) ~= samples_needed)
                            %These two HAVE to be equal for the seeking through the
                            %audio to be correct.
                            error('Seeking through Audio data has failed!');
                        end
                        ind = ind + 1;
                        try
                            fseek(file, info.audFrames(1, ind), -1);
                        catch
                            %No more audio in the file, it was all read in.  Thats
                            %fine, this will be taken care of down below.
                        end
                    end
                end
                
                % while we haven't read in all that has been requested
                while size(data, 2) < total_samples
                    if(isempty(data))
                        %Dividing has to do with staying within the audio region!
                        %Otherwise, it goes outside of the audio region and noise
                        %is added.
                        data = [data fread(file, [info.NumAudioChannels, ...
                            floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                            /(BytesPerSample*info.NumAudioChannels))], dtype)];
                    else
                        data = [data fread(file, [info.NumAudioChannels, ...
                            floor((info.audFrames(1,ind)+info.audFrames(2,ind)-ftell(file))...
                            /(BytesPerSample*info.NumAudioChannels))], dtype)];
                    end
                    
                    ind = ind + 1;
                    % break out if there is no more audio to read. We won't
                    % penalize for users requesting too much.
                    if (ind > size(info.audFrames, 2))
                        break;
                    end
                    fseek(file, info.audFrames(1, ind), -1);
                end
                
                if size(data, 2) < total_samples
                    total_samples = size(data, 2);
                end
                % truncate the output just in case we read too much
                data = data(:, 1:total_samples);
                
                data = data';
                
                if(strcmpi(dtype,'bit32'))
                    data = data/(2^31); % [-1,1)
                end
            end
        end
    end     
            
    % return only what was requested
    o4 = data;
    if size(o4, 1) < total_samples
        total_samples = size(o4, 1);
    end
    % truncate the output just in case we read too much
    o4 = o4(1:total_samples, :);
    % give the user the audio rate
    o5 = info.AudioRate;
    
    if(info.BitsPerSample == 32)
        %Need to let the AVI write function know whether the 32 bit audio
        %is floating point or not.  This information will be placed in o5
        %which displays the audio rate.  This was choosen because it is the
        %only returnable item that doesn't contain video or audio data.
        if(strcmpi(dtype,'bit32'))
            %Its 32 bit (not floating point)
            o5 = [info.AudioRate, 0];
        else
            %Its 32 bit floating point
            o5 = [info.AudioRate, 1];
        end
    end
            
end

fclose(file);
return ;

%% -----------------------------------------------------------------------
function [r,g,b] ...
    = read_rgb24_frame(fid, is_whole_image, sroi, num_rows, num_cols)
% Read one RGB24 frame

% read in image
temp = readAndCheck(fid, [3*num_cols,num_rows], '*uint8');

% flip.
temp = temp(:,num_rows:-1:1);

% pick off the planes
temp = reshape(temp', num_rows, 3, num_cols);

b = single(squeeze(temp(:,1,:)));
g = single(squeeze(temp(:,2,:)));
r = single(squeeze(temp(:,3,:)));

if ~is_whole_image,
    r = r(sroi.top:sroi.bottom, sroi.left:sroi.right);
    g = g(sroi.top:sroi.bottom, sroi.left:sroi.right);
    b = b(sroi.top:sroi.bottom, sroi.left:sroi.right);
end
return ;

%% -----------------------------------------------------------------------
function [r,g,b] = read_rgb32_frame(fid, is_whole_image, ...
    sroi, num_rows, num_cols)
% Read one RGB24 frame

% read in image
temp = readAndCheck(fid, [4*num_cols,num_rows], '*uint8');

% flip.
temp = temp(:,num_rows:-1:1);

% pick off the planes
temp = reshape(temp', num_rows, 4, num_cols);

b = single(squeeze(temp(:,1,:)));
g = single(squeeze(temp(:,2,:)));
r = single(squeeze(temp(:,3,:)));

if ~is_whole_image,
    r = r(sroi.top:sroi.bottom, sroi.left:sroi.right);
    g = g(sroi.top:sroi.bottom, sroi.left:sroi.right);
    b = b(sroi.top:sroi.bottom, sroi.left:sroi.right);
end
return ;

%% -----------------------------------------------------------------------
function [y,cb,cr] = read_uyvy_frame(fid, is_whole_image, sroi, ...
    is_interp, num_rows, num_cols)
% read one YCbCr frame

% read in image
data = readAndCheck(fid, [2*num_cols, num_rows], '*uint8');
% pick off the Y plane (luminance)
temp = reshape(data', num_rows, 2, num_cols);
y = squeeze(temp(:,2,:));

% If color image planes are requested, pick those off and perform
% pixel replication to upsample horizontally by 2.
temp = reshape(data,4,num_rows*num_cols/2);

cb = reshape(temp(1,:),num_cols/2,num_rows)';
cb = [cb ; cb];
cb = reshape(cb, num_rows, num_cols);

cr = reshape(temp(3,:),num_cols/2,num_rows)';
cr = [cr ; cr];
cr = reshape(cr, num_rows, num_cols);

% Interpolate, if requested
if is_interp == 1,
    for i=2:2:num_cols-2,
        % Bug fix 3/16/09
        cb(:,i) = uint8(round( ...
            (double(cb(:,i-1)) + double(cb(:,i+1)))/2));
        % Bug fix 3/16/09
        cr(:,i) = uint8(round( ...
            (double(cr(:,i-1)) + double(cr(:,i+1)))/2));
    end
end

if ~is_whole_image,
    y = y(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
    if nargout == 3,
        cb = cb(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
        cr = cr(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
    end
end
return ;

%% -----------------------------------------------------------------------
function [y,cb,cr] = read_yuyv_frame(fid, is_whole_image, sroi, ...
    is_interp, num_rows, num_cols)
% read one YCbCr frame (yuy2 format)

% read in image
data = readAndCheck(fid, [2*num_cols, num_rows], '*uint8');
% pick off the Y plane (luminance)
temp = reshape(data', num_rows, 2, num_cols);
y = squeeze(temp(:,1,:));

% If color image planes are requested, pick those off and perform
% pixel replication to upsample horizontally by 2.
temp = reshape(data,4,num_rows*num_cols/2);

cb = reshape(temp(2,:),num_cols/2,num_rows)';
cb = [cb ; cb];
cb = reshape(cb, num_rows, num_cols);

cr = reshape(temp(4,:),num_cols/2,num_rows)';
cr = [cr ; cr];
cr = reshape(cr, num_rows, num_cols);

% Interpolate, if requested
if is_interp == 1,
    for i=2:2:num_cols-2,
        % Bug fix 3/16/09
        cb(:,i) = uint8(round( ...
            (double(cb(:,i-1)) + double(cb(:,i+1)))/2));
        % Bug fix 3/16/09
        cr(:,i) = uint8(round( ...
            (double(cr(:,i-1)) + double(cr(:,i+1)))/2));
    end
end

if ~is_whole_image,
    y = y(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
    if nargout == 3,
        cb = cb(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
        cr = cr(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
    end
end
return ;

%% -----------------------------------------------------------------------
function [y,cb,cr] = read_yv12_frame(fid, is_whole_image, sroi, ...
    is_interp, num_rows, num_cols)
% read one YCbCr frame (yv12 format - http://fourcc.org/yuv.php#YV12)

% read in Y plane (luminance)
y = readAndCheck(fid, [num_cols, num_rows], '*uint8')';

% If color image planes are requested, read them in from the
% file and resample them
cr = readAndCheck(fid, num_cols/2 * num_rows/2, '*uint8');
cr = [cr, cr];
cr = reshape(cr', num_cols, num_rows/2)';
cr = [cr, cr];
cr = reshape(cr', num_cols, num_rows)';

cb = readAndCheck(fid, num_cols/2 * num_rows/2, '*uint8');
cb = [cb, cb];
cb = reshape(cb', num_cols, num_rows/2)';
cb = [cb, cb];
cb = reshape(cb', num_cols, num_rows)';

% cannot interpolate yv12
if is_interp == 1,
    error('Cannot interpolate yv12 format');
end

if ~is_whole_image,
    y = y(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
    if nargout == 3,
        cb = cb(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
        cr = cr(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
    end
end
return ;

%% -----------------------------------------------------------------------
function [y,cb,cr] = read_10bit_uyvy(fid, is_whole_image, sroi, ...
    is_interp, num_rows, num_cols)
% read one YCbCr frame from 10-bit data
% NOTE: Look at the following site for explanation:
% http://developer.apple.com/quicktime/icefloe/dispatch019.html#v210

global data;
% Picture of what the bytes look like coming in: 1 char = 1 bit
%   ///-CB-/// ////-Y-/// ///-CR-/// XX     - 1st 32 bits
%   ////-Y-/// ///-CB-/// ////-Y-/// XX     - 2nd 32 bits
%   ///-CR-/// ////-Y-/// ///-CB-/// XX     - 3rd 32 bits
%   ////-Y-/// ///-CR-/// ////-Y-/// XX     - 4th 32 bits

% read in image
sz = num_rows*num_cols*2;
% reads three 10-bit fields and then skips 2 bits. Nifty, right?
data = fread(fid, sz, '3*ubit10', 2);

y = data(2:2:sz);
cb= data(1:4:sz);
cr= data(3:4:sz);

% because this is 10-bit data, the bits need to be shifted right
% 2 places to turn it into 8-bit data.
y  = single(y)/4;
cb = single(cb)/4;
cr = single(cr)/4;

% reshape the extracted information
y = reshape(y, num_cols, num_rows)';

cb = reshape(cb,num_cols/2,num_rows)';
cb = [cb ; cb];
cb = reshape(cb,num_rows,num_cols);

cr = reshape(cr,num_cols/2,num_rows)';
cr = [cr ; cr];
cr = reshape(cr,num_rows,num_cols);

% Interpolate, if requested
if is_interp == 1,
    for i=2:2:num_cols-2,
        % Bug fix 3/16/09
        cb(:,i) = (cb(:,i-1) + cb(:,i+1))/2;
        % Bug fix 3/16/09
        cr(:,i) = (cr(:,i-1) + cr(:,i+1))/2;
    end
end

if ~is_whole_image,
    y = y(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
    if nargout == 3,
        cb = cb(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
        cr = cr(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
    end
end
return ;

%% -----------------------------------------------------------------------
function data = readAndCheck( file, num, datatype )
% Reads data from the specified file
[data, count] = fread(file, num, datatype);
% quick fix for data read into 2D matrices. This needs to be made
% more general though
if (size(num, 2) > 1)
    num = num(1,1)*num(1,2);
end
assert( eq( count, num) );
return ;

%% -----------------------------------------------------------------------
% THIS IS AN ALTERNATE IMPLEMENTATION OF THE 10-BIT READ FUNCTION.
% TO USE IT, SIMPLY UNCOMMENT THE LINES BELOW AND COMMENT-OUT THE
% CURRENT FUNCTION. IT IS SLOWER THAN THE CURRENT VERSION, BUT THERE
% MIGHT BE A WAY TO SPEED IT UP.
% % ----------------------------------------------------------------------
% function [y,cb,cr] = read_10bit_uyvy(fid, is_whole_image, sroi, ...
%                                      is_interp, num_rows, num_cols, sz)
% % read one YCbCr frame from 10-bit data
% % NOTE: Look at the following site for explanation:
% % http://developer.apple.com/quicktime/icefloe/dispatch019.html#v210
%
%     global data;
%     % read in image
%     data = readAndCheck(fid, [4, sz/16], '*uint32');
%     % preallocate - ZEROS and NOT ONES is VERY important here I
%     % think...but maybe it doesn't matter
%      y = zeros([1,num_rows*num_cols], 'uint16');
%     cb = zeros([1,num_cols*num_rows/2], 'uint16');
%     cr = cb;
%
%     % Picture of what the bytes look like coming in: 1 char = 1 bit
%     %   ///-CB-/// ////-Y-/// ///-CR-/// XX     - 1st 32 bits
%     %   ////-Y-/// ///-CB-/// ////-Y-/// XX     - 2nd 32 bits
%     %   ///-CR-/// ////-Y-/// ///-CB-/// XX     - 3rd 32 bits
%     %   ////-Y-/// ///-CR-/// ////-Y-/// XX     - 4th 32 bits
%
%     sz = num_cols*num_rows;
%     % in every 4x4 byte matrix, there are 6 y values that need to
%     % be extracted
%     y  = read10( y, 1, 6, sz, 1, 11:20);
%     y  = read10( y, 2, 6, sz, 2,  1:10);
%     y  = read10( y, 3, 6, sz, 2, 21:30);
%     y  = read10( y, 4, 6, sz, 3, 11:20);
%     y  = read10( y, 5, 6, sz, 4,  1:10);
%     y  = read10( y, 6, 6, sz, 4, 21:30);
%
%     sz = sz/2;
%     % in every 4x4 byte matrix, there are 3 cb and cr values that
%     % need to be extracted
%     cb = read10(cb, 1, 3, sz, 1,  1:10);
%     cb = read10(cb, 2, 3, sz, 2, 11:20);
%     cb = read10(cb, 3, 3, sz, 3, 21:30);
%
%     cr = read10(cr, 1, 3, sz, 1, 21:30);
%     cr = read10(cr, 2, 3, sz, 3,  1:10);
%     cr = read10(cr, 3, 3, sz, 4, 11:20);
%
%     % because this is 10-bit data, the bits need to be shifted right
%     % 2 places to turn it into 8-bit data.
%     y  = single(y)/4;
%     cb = single(cb)/4;
%     cr = single(cr)/4;
%
%     % reshape the extracted information
%     y = reshape(y, num_cols, num_rows)';
%
%     cb = reshape(cb,num_cols/2,num_rows)';
%     cb = [cb ; cb];
%     cb = reshape(cb,num_rows,num_cols);
%
%     cr = reshape(cr,num_cols/2,num_rows)';
%     cr = [cr ; cr];
%     cr = reshape(cr,num_rows,num_cols);
%
%     % Interpolate, if requested
%     if is_interp == 1,
%         for i=2:2:num_cols-2,
%             % Bug fix 3/16/09
%             cb(:,i) = (cb(:,i-1) + cb(:,i+1))/2;
%             % Bug fix 3/16/09
%             cr(:,i) = (cr(:,i-1) + cr(:,i+1))/2;
%         end
%     end
%
%     if ~is_whole_image,
%         y = y(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
%         if nargout == 3,
%             cb = cb(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
%             cr = cr(sroi.top:sroi.bottom, sroi.left:sroi.right, :);
%         end
%     end
% return ;
%
% function ret = read10( V, vp, st, sz, dp, bi )
% % This function serves a VERY specific purpose. It reads the given
% % 10 bits from the global data and sets those bits in the y, cb, or
% % cr matrix that is given by V. Honestly, it is best if you don't
% % try to understand it, ha ha.
%     global data;
%     try
%         V(vp:st:sz) = bitset(V(vp:st:sz), 1, bitget(data(dp,:), bi(1)));
%         V(vp:st:sz) = bitset(V(vp:st:sz), 2, bitget(data(dp,:), bi(2)));
%         V(vp:st:sz) = bitset(V(vp:st:sz), 3, bitget(data(dp,:), bi(3)));
%         V(vp:st:sz) = bitset(V(vp:st:sz), 4, bitget(data(dp,:), bi(4)));
%         V(vp:st:sz) = bitset(V(vp:st:sz), 5, bitget(data(dp,:), bi(5)));
%         V(vp:st:sz) = bitset(V(vp:st:sz), 6, bitget(data(dp,:), bi(6)));
%         V(vp:st:sz) = bitset(V(vp:st:sz), 7, bitget(data(dp,:), bi(7)));
%         V(vp:st:sz) = bitset(V(vp:st:sz), 8, bitget(data(dp,:), bi(8)));
%         V(vp:st:sz) = bitset(V(vp:st:sz), 9, bitget(data(dp,:), bi(9)));
%         V(vp:st:sz) = bitset(V(vp:st:sz),10, bitget(data(dp,:), bi(10)));
%     catch e
%         % there have been some 10-bit videos made by VirtualDub that
%         % are unlike any 10-bit format I've seen.  There should be 6
%         % pixels for every 4 bytes, but these videos have more bytes
%         % than expected.
%         disp([e.message ' - Corrupted 10-bit file']);
%     end
%
%     ret = V;
% return ;
