function info = aviinfo( filename )
%% Reads the header of an AVI file and returns a structure,
% holding all the important information of the header in the following
% format:
%
%   The set of fields for FILEINFO are:
%   
%   Filename           - A string containing the name of the file.
%
%   FileModDate        - A string containing the modification date of the 
%                        file.
%   		      
%   FileSize           - An integer indicating the size of the file in 
%                        bytes.
%
%   Codec              - A four-character code, representing the file's
%                        compression type.
%   		      
%   FramesPerSecond    - An integer indicating the desired frames per 
%                     	 second during playback.
%
%   NumFrames          - An integer indicating the total number of frames 
%                     	 in the movie.
%   		      
%   Height             - An integer indicating the height of AVI movie in
%                     	 pixels.
%   		      
%   Width              - An integer indicating the width of AVI movie in
%                        pixels.
%   		      
%   ImageType          - A string indicating the type of image; either
%                     	 'truecolor' for a truecolor (RGB) image, or
%                     	 'indexed', for an indexed image.
%   		      
%   VideoCompression   - A string containing the compressor used to  
%                     	 compress the AVI file.   If the compressor is 
%                     	 not Microsoft Video 1, Run-Length Encoding, 
%                     	 Cinepak, or Intel Indeo, the four character code 
%                        is returned.
%		      
%   Quality            - A number between 0 and 100 indicating the video
%                     	 quality in the AVI file.  Higher quality numbers
%                     	 indicate higher video quality, where lower
%                     	 quality numbers indicate lower video quality.  
%                     	 This value is not always set in AVI files and 
%                     	 therefore may be inaccurate.
%
%   NumColormapEntries - The number of colors in the colormap. For a
%                        truecolor image this value is zero.
%
%   BitDepth           - A double representing the bit depth of the video.
%
%   ColorType          - A string indicating the video's colorspace.
%
%   vidFrames          - A 2x[NumFrames] matrix containing offsets and
%                        sizes of all frames in the AVI file.
%
% If AUDIO is present in the file, the following fields we will be added:
%
%   AudioFormat        - A string indicating the format in which the
%                        audio is stored.
%
%   AudioRate          - The audio's sampling rate.
%
%   NumAudioChannels   - The number of audio channels in the file.
%
%   BytesPerSec        - The number of bytes per second for audio. This
%                        is mostly for internal use.
%
%   BlockAlign         - A number indicating the block alignment. This is
%                        mostly for internal use to compute BitsPerSample.
%
%   SampleSize         - The size of an audio sample.
%
%   BufferSize         - The size of one audio data chunk.
%
%   BitsPerSample      - The number of bits per sample (per channel).
%
%   audFrames          - A 2x[number of audio chunks] matrix, containing
%                        the offsets and sizes of all audio chunks in the
%                        file.
% EXAMPLES:
% Read in an AVI file to determine it's frame rate
%   % use aviinfo to read file and store in structure
%   inf = aviinfo('my_file.avi');
%   inf.FramesPerSecond
%
% Output:
%   ans = 
%       29.9700
%
% REFERENCES:
%   http://the-labs.com/Video/odmlff2-avidef.pdf
%       Offers a detailed description of the AVI file format.
%
%   http://abcavi.kibi.ru/docs/riff2.pdf
%       An excellent reference for registered AVI audio formats.
%
%   http://www.jmcgowan.com/avi.html#WinProg
%       A not as detailed but easy to understand reference on the AVI
%       file format.
%
%   http://developer.apple.com/quicktime/icefloe/dispatch019.html#v210
%       What you need to know about 10-bit AVI files.
%
    file = fopen(filename, 'r', 'l');
    if file < 0
        error('MATLAB:aviinfoalt', '''%s'' does not exist', filename);
    end
    file_inf            = dir(filename);
    info.Filename       = filename;
    info.FileModDate    = file_inf.date;
    info.FileSize       = file_inf.bytes;
    fpos = 1; % stores the current reading position of data from file

    % format of beginning header in file:
    %   'RIFF' size 'avi  LIST' size 'hdrl avih' size 
    beg_inf = readAndCheck(file, 32, '*uint8');
    fpos = fpos + 4;
%     riff_size = typecast( beg_inf(fpos:fpos+3), 'uint32' );
    fpos = fpos+4;

    if (strcmpi( 'avi ', char(beg_inf(fpos:fpos+3)') ) < 1 )
        error('MATLAB:aviinfoalt', ...
              '''%s'' is not an avi file', info.Filename);
    end
    fpos = fpos+4+4;
    % LIST header contains all header information before the
    % actual movie information
    hdr_end_pos = typecast( beg_inf(fpos:fpos+3), 'uint32' ) ...
                  + fpos-1 + 4;
    fpos = fpos+4+8;

    avi_hdr_size = typecast( beg_inf(fpos:fpos+3), 'uint32' );
%     fpos = fpos+4;
    data = readAndCheck(file, avi_hdr_size, '*uchar');
    
    avi_hdr = readAviHeader( data );
    vids_found = 0;
    auds_found = 0;
% sometimes there might be 2 streams but only one strh LIST, so this
% is a very tricky way to make sure all the strh headers are found
    while vids_found+auds_found < avi_hdr.NumStreams
        stream_size = findlist( file, 'strl' );
        stream_end_pos = ftell(file) + stream_size - 4;
        stream ...
            = readAndCheck(file, stream_end_pos-ftell(file), '*uchar');
        fpos = 1;
        
        % the 'strh' chunk tells us whether this is a header for
        % audio or video and tells us the color encoding if video
        [strh_size, fpos] = findData(stream, fpos, 'strh');
        strh_end_pos = fpos + strh_size;
        stream_type = char( stream(fpos:fpos+3)' );
        fpos = fpos + 4;
%% -----------------------------------------------------------------------
        if ( strcmpi(stream_type, 'vids') )
            vids_found = 1;
            strh = readAudioVideoHeader( stream, fpos );
            info.Codec = strh.Codec;
            info.FramesPerSecond = strh.Rate/strh.Scale;
            fpos = strh_end_pos;
    
            % Some files have a mismatched number of frames in the 
            % MainHeader and VideoStreamHeader.  We should always trust 
            % the VideoStreamHeaders number of frames over the 
            % MainHeader.TotalFrmes. However, there was a bug in 
            % AVIFILE which generates files such that 
            % MainHeader.NumFrames was greater than 
            % VideoStreamHeader's frames. So, to maintain 
            % backward compatibility with AVIFILE generated AVIs only 
            % update the main headers total frames if it is less than 
            % the video stream headers frame count.
            if avi_hdr.NumFrames < strh.Length
               avi_hdr.NumFrames = strh.Length;
            end
            
            % the 'strf' chunk contains information about the
            % planes, bit depth, compression, and color map
            [strf_size, fpos] = findData(stream, fpos, 'strf');
            strf_end_pos = fpos + strf_size;
            strf_hdr = readBitmapHeader(stream, fpos, strf_size);
            fpos = strf_end_pos;
            
            % the 'indx' is either a 'frame index' with offsets and
            % sizes of '00db' chunks or a 'super index' with offsets
            % and sizes to 'frame indexes' throughout the file. This
            % chunk does not always exist.
            if (fpos+3 <= size(stream, 1))
                nextChunk = stream(fpos:fpos+3);
            else
                nextChunk = ' ';
            end
            fpos = fpos + 4;
            if ( strcmpi( char( nextChunk' ), 'indx' ) )
                fpos = fpos + 8;
                num_indexes ...
                    = typecast( stream(fpos:fpos+3), 'uint32' );
                % skip 4 bytes for what was just read and skip '00db' 
                % and empty bytes after it
                fpos = fpos + 4 + 16;
                index_vids = ones(num_indexes, 3); % preallocate
                % This is safe because I believe all 'indx'
                % chunks are super indexes.
                for index_num = 1:num_indexes
                    % assert we aren't overwriting data from other 
                    % streams
                    assert(index_vids(index_num, 1) == 1);
                    % read in offset of index
                    index_vids(index_num, 1) ...
                        = typecast(stream(fpos:fpos+7), 'uint64');
                    fpos = fpos + 8;
                    % read in size of index
                    index_vids(index_num, 2) ...
                        = typecast(stream(fpos:fpos+3), 'uint32');
                    fpos = fpos + 4;
                    % read in number of frames contained in index
                    index_vids(index_num, 3) ...
                        = typecast(stream(fpos:fpos+3), 'uint32');
                    fpos = fpos + 4;
                end
                vidFrames = getvidFrames(file, index_vids, avi_hdr);
            else
            % there is no super index, so we must seek past movi
            % information and find the 'idx1' chunk
                ret = ftell(file);
                fseek(file, hdr_end_pos, -1);
                movi_size = findlist( file, 'movi' );
                header_offset = ftell(file)-4; % -4 for 'movi'
                fseek(file, movi_size-4, 0);
                Itype = readAndCheck(file, 4, '*uchar')';
                if (strcmpi(char(Itype), 'idx1'))
                    index_size = fread(file, 1, 'uint32');
                    index = fread(file, (index_size/4), 'uint32');
                    % idx1 comes in looking like the following:
                    % '00db' or '00dc'
                    % 00 00 00 10 -- I don't know what this is
                    % [offset of frame]
                    % [size of frame]
                    % ... repeat ...
                    index = reshape(index, 4, (index_size/16));

                    % if audio exists, it will be interleaved here
                    % '00db' cast to int == 1650733104
                    % '01wb' cast to int == 1651978544
                    video = index(:, (index(1,:) == 1650733104) | ...
                                     (index(1,:) == 1667510320));
                    audio = index(:, (index(1,:) == 1651978544) | ...
                                     (index(1,:) == 1668755760));
                    vidFrames = video(3:4, :);
                    if (size(audio, 2) > 0)
                        audFrames = audio(3:4, :);
                    end
                    % the offsets in this index may be absolute 
                    % offsets or relative to the header. AVI parsers 
                    % should handle either case
                    fseek(file, vidFrames(1,1), -1);
                    Otype = char(fread(file, 4, 'uchar')');
                    if (~strcmpi(Otype, '00db') && ...
                        ~strcmpi(Otype, '00dc'))
                        vidFrames(1, :) ...
                            = vidFrames(1, :) + header_offset;
                        if (size(audio, 2) > 0)
                            audFrames(1, :) ...
                                = audFrames(1, :) + header_offset;
                        end
                    end
                    vidFrames(1, :) = vidFrames(1, :) + 8;
                    if (size(audio, 2) > 0)
                        audFrames(1, :) = audFrames(1, :) + 8;
                    end
                else
                    % No index was found, so the frames must be
                    % manually located. Issue a warning, so the
                    % user can be aware of the problem.
                    warning('readAvi:index_location', ...
                            'No index was found - reading may be slow');
                    vidFrames = zeros(2, avi_hdr.NumFrames);
                    fseek(file, header_offset+4, -1);
                    frame = 1;
                    while (ftell(file) < movi_size+header_offset);
                        chunk.ckid ...
                            = char( readAndCheck(file, 4, '*uchar')' );
                        chunk.cksize ...
                            = readAndCheck(file, 1, 'uint32');
                        if (strcmpi(chunk.ckid,'00db') || ...
                            strcmpi(chunk.ckid,'00dc'))
                            vidFrames(1, frame) = ftell(file);
                            vidFrames(2, frame) = chunk.cksize;
                            frame = frame + 1;
                        end
                        skipchunk(file,chunk);
                    end
                end
                fseek(file, ret, -1); % return to header positions
            end
%             fpos = fpos + 4;
%% -----------------------------------------------------------------------
        elseif ( strcmpi(stream_type, 'auds') )
            auds_found = 1;
            strh_aud = readAudioVideoHeader(stream, fpos);
            fpos = strh_end_pos;
            [strf_size, fpos] = findData(stream, fpos, 'strf');
            strf_end_pos = fpos + strf_size;
            
            strf_aud = readAudioFormat( stream, fpos );
            fpos = strf_end_pos;
            
            if ~exist('audFrames', 'var')
                % the 'indx' is either a 'frame index' with offsets and
                % sizes of '00db' chunks or a 'super index' with offsets
                % and sizes to 'frame indexes' throughout the file. This
                % chunk does not always exist.
                if (fpos+3 <= size(stream, 1))
                    nextChunk = stream(fpos:fpos+3);
                else
                    nextChunk = ' ';
                end
                fpos = fpos + 4;
                if ( strcmpi( char( nextChunk' ), 'indx' ) )
                    fpos = fpos + 8;
                    num_indexes ...
                        = typecast( stream(fpos:fpos+3), 'uint32' );
                    % skip 4 bytes for what was just read and skip '00db' 
                    % and empty bytes after it
                    fpos = fpos + 4 + 16;
                    index_auds = ones(num_indexes, 3); % preallocate
                    % This is safe because I believe all 'indx'
                    % chunks are super indexes.
                    for index_num = 1:num_indexes
                        % assert we aren't overwriting data from other 
                        % streams
                        assert(index_auds(index_num, 1) == 1);
                        % read in offset of index
                        index_auds(index_num, 1) ...
                            = typecast(stream(fpos:fpos+7), 'uint64');
                        fpos = fpos + 8;
                        % read in size of index
                        index_auds(index_num, 2) ...
                            = typecast(stream(fpos:fpos+3), 'uint32');
                        fpos = fpos + 4;
                        % read in number of frames contained in index
                        index_auds(index_num, 3) ...
                            = typecast(stream(fpos:fpos+3), 'uint32');
                        fpos = fpos + 4;
                    end
                    audFrames = getaudFrames(file, index_auds);
    %                 audFrames = [];
                else
                % there is no super index, so we must seek past movi
                % information and find the 'idx1' chunk
                    ret = ftell(file);
                    fseek(file, hdr_end_pos, -1);
                    movi_size = findlist( file, 'movi' );
                    header_offset = ftell(file)-4; % -4 for 'movi'
                    fseek(file, movi_size-4, 0);
                    Itype = readAndCheck(file, 4, '*uchar')';
                    if (strcmpi(char(Itype), 'idx1'))
                        index_size = fread(file, 1, 'uint32');
                        index = fread(file, (index_size/4), 'uint32');
                        % idx1 comes in looking like the following:
                        % '00db' or '00dc'
                        % 00 00 00 10 -- I don't know what this is
                        % [offset of frame]
                        % [size of frame]
                        % ... repeat ...
                        index = reshape(index, 4, (index_size/16));
                        % just in case the index is interleaved, pick
                        % off all chunks that aren't video chunks.
                        % '01wb' cast to int == 1651978544
                        audio = index(:, (index(1,:) == 1651978544) | ...
                                         (index(1,:) == 1668755760));
                        audFrames = audio(3:4, :);
                        % the offsets in this index may be absolute 
                        % offsets or relative to the header. AVI parsers 
                        % should handle either case
                        fseek(file, audFrames(1,1), -1);
                        Otype = char(fread(file, 4, 'uchar')');
                        if (~strcmpi(Otype, '01wb') && ...
                            ~strcmpi(Otype, '01wc'))
                            audFrames(1, :) ...
                                = audFrames(1, :) + header_offset;
                        end
                        audFrames(1, :) = audFrames(1, :) + 8;
                    else
                        % No index was found, so the frames must be
                        % manually located. Issue a warning, so the
                        % user can be aware of the problem.
                        warning('readAvi:index_location', ...
                                'No index was found - reading may be slow');
                        audFrames = zeros(2, avi_hdr.NumFrames);
                        fseek(file, header_offset+4, -1);
                        frame = 1;
                        while (ftell(file) < movi_size+header_offset);
                            chunk.ckid ...
                                = char( readAndCheck(file, 4, '*uchar')' );
                            chunk.cksize ...
                                = readAndCheck(file, 1, 'uint32');
                            if (strcmpi(chunk.ckid,'01wb') || ...
                                strcmpi(chunk.ckid,'01wc'))
                                audFrames(1, frame) = ftell(file);
                                audFrames(2, frame) = chunk.cksize;
                                frame = frame + 1;
                            end
                            skipchunk(file,chunk);
                        end
                    end
                    fseek(file, ret, -1); % return to header positions
                end
    %             fpos = fpos + 4;
            end
        end
    end
%% prepare output --------------------------------------------------------
    info.NumFrames          = cast(avi_hdr.NumFrames, 'double');
    info.Height             = cast(avi_hdr.Height, 'double');
    info.Width              = cast(avi_hdr.Width, 'double');
    info.ImageType          = strf_hdr.ImageType;
    info.VideoCompression   = strf_hdr.VideoCompression;
    info.Quality            = 0; % TODO
    info.NumColormapEntries = cast(strf_hdr.NumColormapEntries, 'double');
    info.BitDepth           = cast(strf_hdr.BitDepth, 'double');
    % determine the file's color type
    % bytes/pixel possibilities:
    % 3/1  - rgb 24
    % 4/1  - rgb 32
    % 4/2  - uyvy, yuyv
    % 16/6 - uyvy 10-bit
    if strcmpi(info.VideoCompression, 'none')
        if info.BitDepth == 24 || info.BitDepth == 32
            info.ColorType = 'RGB';
        else
            info.ColorType = 'UYVY';
        end
    else
        info.ColorType = 'UYVY';
    end
    
    info.vidFrames    = cast(vidFrames, 'double');
    if exist('strf_aud', 'var')
        info.AudioFormat = strf_aud.Format;
        info.AudioRate = strf_aud.SampleRate;
        info.NumAudioChannels = strf_aud.NumChannels;
        info.BytesPerSec = strf_aud.BytesPerSec;
        info.BlockAlign = strf_aud.BlockAlign;
        info.SampleSize = strh_aud.SampleSize;
        info.BufferSize = strh_aud.BufferSize;
        info.BitsPerSample = 8*ceil(info.BlockAlign/info.NumAudioChannels);
        info.audFrames = audFrames;
        if(isfield(strf_aud,'SubFormat'))
            info.SubFormat = strf_aud.SubFormat;
        end
            
    end
    
    fclose(file);
return ;

%% -----------------------------------------------------------------------
function avi_hdr = readAviHeader( data )
% Reads an 'avih' header
% Input: The header information read in from file

    fpos = 1;
    fpos = fpos+4; % skip microseconds per frame
    avi_hdr.MaxBytePerSec ...
        = double(typecast( data(fpos:fpos+3), 'uint32' ));
    fpos = fpos+8; % skip reserved
    
    flags =  typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos+4;
    flagbits = find( bitget( flags, 1:32 ) );
    for i = 1:length(flagbits)
      switch flagbits(i)
       case 5
        avi_hdr.HasIndex = 'True';
       case 6
        avi_hdr.MustUseIndex = 'True';
       case 9
        avi_hdr.IsInterleaved = 'True';
       case 12
        avi_hdr.TrustCKType = 'True';
       case 17
        avi_hdr.WasCaptureFile = 'True';
       case 18
        avi_hdr.Copywrited = 'True';
      end
    end

    avi_hdr.NumFrames = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos+8; % skip Initial frames
    avi_hdr.NumStreams = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos+8; % skip Suggested Buffer Size
    avi_hdr.Width = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos+4;
    % Height may be negative if AVI is written top down
    avi_hdr.Height ...
        = typecast(abs( ...
          typecast( data(fpos:fpos+3), 'int32' ) ), 'uint32' );
    fpos = fpos+4;
    avi_hdr.Scale = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos+4;
    avi_hdr.Rate = typecast( data(fpos:fpos+3), 'uint32' );
%     fpos = fpos+4;
    % skip start and length
return ;

%% -----------------------------------------------------------------------
function strh = readAudioVideoHeader( stream, fpos )
% reads the STRH chunk information. this can either contain video
% or audio information

    strh.Codec = char( stream(fpos:fpos+3)' );
    fpos = fpos + 4 + 12; % skip flags, reserved, initial frames

    strh.Scale = double(typecast(stream(fpos:fpos+3), 'uint32'));
    fpos = fpos + 4;
    
    strh.Rate = double(typecast(stream(fpos:fpos+3), 'uint32'));
    fpos = fpos + 4 + 4; % skip start
    
    strh.Length = double(typecast(stream(fpos:fpos+3), 'uint32'));
    fpos = fpos + 4;
    
    strh.BufferSize = double(typecast(stream(fpos:fpos+3), 'uint32'));
    fpos = fpos + 4 + 4; % skip quality (it is unreliable)
    
    strh.SampleSize = double(typecast(stream(fpos:fpos+3), 'uint32'));
return ;

%% -----------------------------------------------------------------------
function strf = readBitmapHeader(data, fpos, strf_size)
% Reads the BITMAPINFO header information.

    strf.BitmapHeaderSize = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos + 4;
    strf.Width = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos + 4;
    strf.Height = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos + 4;
    strf.Planes = typecast( data(fpos:fpos+1), 'uint16' );
    fpos = fpos + 2;
    strf.BitDepth = typecast( data(fpos:fpos+1), 'uint16' );
    fpos = fpos + 2;

    % Read Compression.
    compress = typecast( data(fpos:fpos+3), 'uint32' );
    switch compress
    case 0
        compression = 'none';
    case 1
        compression = '8-bit RLE';
    case 2
        compression = '4-bit RLE';
    case 3
        compression = 'bitfields';
    otherwise
        compression = '';
    end
    if isempty(compression)
        code = [char(bitshift(compress,0,8)) ...
                char(bitshift(compress,-8,8)) ...
                char(bitshift(compress,-16,8)) ...
                char(bitshift(compress,-24,8))];
        switch lower(code)
        case 'none'
            compression = 'None';
        case 'rgb '
            compression = 'None';
        case 'raw '
            compression = 'None';  
        case '    '
            compression = 'None';
        case 'rle '
            compression = 'RLE';
        case 'cvid'
            compression = 'Cinepak';
        case 'iv32'
            compression = 'Indeo3';
        case 'iv50'
            compression = 'Indeo5';
        case 'msvc'
            compression = 'MSVC';
        case 'cram'
            compression = 'MSVC';
        otherwise
            compression = code;
        end
    end
    fpos = fpos + 4;
    strf.VideoCompression = compression;
    
    strf.Bitmapsize = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos + 4;
    strf.HorzResoltion = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos + 4;
    strf.VertResolution = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos + 4;
    strf.NumColorsUsed = typecast( data(fpos:fpos+3), 'uint32' );
    fpos = fpos + 4;
    strf.NumImportantColors = typecast( data(fpos:fpos+3), 'uint32' );
%     fpos = fpos + 4;
    strf.NumColormapEntries = ...
        (strf_size - strf.BitmapHeaderSize)/4;
    % 8-bit grayscale
    % 24-bit truecolor
    % 8-bit indexed
    if strf.NumColorsUsed > 0
        strf.ImageType = 'indexed';
    else
        if strf.BitDepth > 8
            strf.ImageType = 'truecolor';
            strf.ImageType = 'truecolor';
        else
            strf.ImageType = 'grayscale';
        end
    end
return ;

%% -----------------------------------------------------------------------
function strf = readAudioFormat( data, fpos )
% Read WAV format chunk information.

    % Read format tag.
    formatTag = typecast( data(fpos:fpos+1), 'uint16' );
    fpos = fpos + 2;

    % Complete list of formats can be found in Microsoft Platform SDK
    % header file "MMReg.h" or in MSDN Library (search for "registered 
    % wave formats").
    switch formatTag
     case  1
      strf.Format = 'PCM';
     case 2
      strf.Format = 'Microsoft ADPCM';
     case 6
      strf.Format = 'CCITT a-law';
     case 7
      strf.Format = 'CCITT mu-law';
     case 17
      strf.Format = 'IMA ADPCM';
     case 34
      strf.Format = 'DSP Group TrueSpeech TM';
     case 49
      strf.Format = 'GSM 6.10';
     case 50
      strf.Format = 'MSN Audio';
     otherwise
      strf.Format = ['Format # 0x' dec2hex(formatTag)];
    end

    % Read number of channels.
    strf.NumChannels = double(typecast( data(fpos:fpos+1), 'uint16'));
    fpos = fpos + 2;

    % Read samples per second.
    strf.SampleRate = double(typecast( data(fpos:fpos+3), 'uint32'));
    fpos = fpos + 4;

    % Read buffer estimation.
    strf.BytesPerSec = double(typecast( data(fpos:fpos+3), 'uint32'));
    fpos = fpos + 4;
    
    % Read block size of data.
    strf.BlockAlign = double(typecast( data(fpos:fpos+1), 'uint16'));
%     fpos = fpos + 2;

%MY EDIT - Information gathered from 
%       http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/WAVE.html
    if(formatTag == 65534)
        %The following section only needs to be added for clips that are of
        %format fffe (WAVE_FORMAT_EXTENSIBLE).  Otherwise, this information
        %does not exist.
        fpos = fpos + 2;
        %Bits per sample
        %     test = double(typecast( data(fpos:fpos+1),'uint16'));
        
        fpos = fpos + 2;
        %Size of extension (cbSize)
        %     test1 = double(typecast( data(fpos:fpos+1),'uint16'));
        
        fpos = fpos + 2;
        %Valid Bits Per Sample
%         strf.ValidBitsPerSample = double(typecast( data(fpos:fpos+1),'uint16'));
        
        fpos = fpos + 2;
        %dwChannelMask (Speaker position mask)
%         strf.SpeakerPos = double(typecast( data(fpos:fpos+3),'uint32'));
        %Quadraphonic = 0x33 = 00110011
        %Positioning is the following:
        %
        % Back Right, Back Left, Low Freq, Front Center, Front Right, Front
        % Left.
        %
        % 5.1 = 63 (decimal => to hex =>  3F => to binary => 00111111
        
        fpos = fpos + 4;
        %SubFormat (GUID (first tow bytes are the data format code))  As of
        %right now, only the first byte is important to use.  This will allow
        %us to tell if audio that is greater than 2 channels has 32 bit audio
        %or 32 bit floating point audio.  If the response is 1, then its just
        %32 bit audio (Baiscally it just is PCM).  If the response is 3, then
        %its 32 bit floating point audio (IEEE floating point).  GUID,
        %SubFormat, is 16 bytes long.  This function is only taking the first
        %bit.  If needed at a later time, this can be modified to include all
        %16 bytes.
        strf.SubFormat = double(typecast( data(fpos),'uint8'));
        %     test2 = double(typecast( data(fpos:fpos+15),'uint64'));
    end
    
return

%% -----------------------------------------------------------------------
function [size, pos] = findData(data, fpos, target)
% Finds a chunk in the data read from file. MAKE SURE the chunk is
% in the data being passed.
% This is quite similar to "findchunk". The only difference is that
% the data is read from the file and then searched.

    chunk.ckid      = '    ';
    chunk.cksize    = 0;

    while( strcmpi( chunk.ckid, target ) == 0 )
        % if the size is odd, it requires a pad byte
        pad = rem(chunk.cksize, 2); 
        % seek past current chunk
        fpos = fpos + chunk.cksize + pad;
        
        chunk.ckid = char( data( fpos:fpos+3 )' );
        fpos = fpos + 4;
        chunk.cksize = typecast( data(fpos:fpos+3), 'uint32' );
        fpos = fpos + 4;
    end
    size = chunk.cksize;
    pos = fpos;
return ;

%% -----------------------------------------------------------------------
function [frame_list] = getvidFrames( file, index_vids, avi_hdr )
% read data from sub indexes throughout file
% This function is safe because if the file has more than 1 RIFF
% chunk, this function won't be called...theoretically.

    ix_data     = ones( size(index_vids, 1), max(index_vids(:,2))/4 );
    frame_list  = ones(2, avi_hdr.NumFrames);
    framesIn    = 1;
    file_pos    = ftell(file);
    % for each video stream, read in positions of frames
    for index_num = 1:size(index_vids, 1)
        % seek to position of sub-index and read data
        fseek(file, index_vids(index_num, 1)+8, -1);
        num_vals = (index_vids(index_num, 2)-8)/4;

        ix_data(index_num, 1:num_vals) ...
            = readAndCheck(file, num_vals, '*uint32');
        fpos = 1 + 3;
        
        % the '00ix' chunk is formatted as such 'offset size offset
        % size ... '. the offsets represent the distance from the 
        % beginning of their parent RIFF chunk which may also have an
        % offset from the beginning of the file which is read here
        temp = uint32(ix_data(index_num, fpos:fpos+1));
        riff_offset = double(typecast(temp, 'uint64')); % ...wow
            
        fpos = fpos + 3;
        frame_list(:,framesIn:(index_vids(index_num, 3)+framesIn-1))...
            = reshape( ix_data(index_num, ...
                       fpos:(index_vids(index_num, 3)*2+fpos-1)), ...
                       2, index_vids(index_num, 3));
        % account for RIFF offset
        frame_list(1,framesIn:(index_vids(index_num, 3)+framesIn-1))...
            = frame_list(1,...
            framesIn:(index_vids(index_num, 3)+framesIn-1))...
            + riff_offset;
        % in case you weren't following. frames now looks like this:
        %   offset  offset  offset  offset  ...
        %   size    size    size    size    ...
        % nifty right? And the offsets are absolute (from beginning
        % of file)
        framesIn = framesIn + index_vids(index_num, 3);
    end
    % seek back to original position before function call
    fseek(file, file_pos, -1);
return ;

%% -----------------------------------------------------------------------
function [frame_list] = getaudFrames( file, index_auds )
% read data from sub indexes throughout file
% This function is safe because if the file has more than 1 RIFF
% chunk, this function won't be called...theoretically.
%
% TODO: This function is dreadfully similar to getvidFrames. You could
%       easily call this function in place of "getvidFrames" because
%       the end result is the same. The only reason I don't rush to do
%       that is because "getvidFrames" is more efficient because it
%       preallocates. Something could be done in the future, however.

    file_pos    = ftell(file);
    frame_list  = [];
    % for each video stream, read in positions of frames
    for index_num = 1:size(index_auds, 1)
        % seek to position of sub-index and read data
        fseek(file, index_auds(index_num, 1)+8, -1);
        num_vals = (index_auds(index_num, 2)-8)/4;
        
        ix_data = readAndCheck(file, num_vals, '*uint32');
        fpos = 1 + 3;
        
        % the '00ix' chunk is formatted as such 'offset size offset
        % size ... '. the offsets represent the distance from the 
        % beginning of their parent RIFF chunk which may also have an
        % offset from the beginning of the file which is read here
        tmp = uint32(ix_data(fpos:fpos+1));
        riff_offset = double(typecast(tmp, 'uint64')); % ...wow
        fpos = fpos + 3;
        
        tmp = double(reshape(ix_data(fpos:size(ix_data)), 2, []));
        tmp(1, :) = tmp(1, :) + riff_offset;
        % we cannot preallocate since we don't know the size
        frame_list = [frame_list, tmp]; %#ok<AGROW>
    end
    % seek back to original position before function call
    fseek(file, file_pos, -1);
return ;

%% -----------------------------------------------------------------------
function size = findlist(fid,listtype)
% Finds a list in the given file

    found = -1;
    size = -1;
    while(found == -1)
        [chunk,~,~] = findchunk(fid,'LIST');
        checktype = fread(fid, 4, '*uchar')';
        if (checktype == listtype) %#ok<BDSCI>
            size = chunk.cksize;
            break;
        else
            fseek(fid,-4,0); %Go back so we can skip the LIST
            skipchunk(fid,chunk);
        end
        if ( feof(fid) )
            return ;
        end
    end
return ;

%% -----------------------------------------------------------------------
function data = readAndCheck( file, num, datatype )
% Reads data from the specified file
    [data, count] = fread(file, num, datatype);
    % count how many elements were actually requested
%    num = sum(sum(sum(num > 0 | num < 1)));
    assert( eq( count, num) );
return ;

%% -----------------------------------------------------------------------
function [chunk,msg,msgID] = findchunk(fid,chunktype)
%FINDCHUNK find chunk in AVI
%   [CHUNK,MSG,msgID] = FINDCHUNK(FID,CHUNKTYPE) finds a chunk of type 
%   CHUNKTYPE in the AVI file represented by FID.  CHUNK is a structure 
%   with fields 'ckid' and 'cksize' representing the chunk ID and chunk 
%   size respectively.  Unknown chunks are ignored (skipped). 

%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2007/07/26 19:30:47 $

    chunk.ckid = '';
    chunk.cksize = 0;
    msg = '';
    msgID='';

    while( strcmp(chunk.ckid,chunktype) == 0 )
      [msg msgID] = skipchunk(fid,chunk);
      if ~isempty(msg)
        fclose(fid);
        error(msgID,msg);
      end
      [id, count] = fread(fid,4,'uchar');
      chunk.ckid = char(id)';
      if (count ~= 4 )
        msg = sprintf('''%s'' did not appear as expected.',chunktype);
        msgID = 'MATLAB:findchunk:unexpectedChunkType';
      end
      [chunk.cksize, count] = fread(fid,1,'uint32');
      if (count ~= 1)
        msg = sprintf('''%s'' did not appear as expected.',chunktype);
        msgID = 'MATLAB:findchunk:unexpectedChunkType';
      end
      if ( ~isempty(msg) ), return; end
    end
return ;

%% -----------------------------------------------------------------------
function [msg msgID] = skipchunk(fid,chunk)
%SKIPCHUNK skip chunk in AVI
%   [MSG MSGID] = SKIPCHUNK(FID,CHUNK) skips CHUNK.cksize bytes in the 
%   AVI file FID.  MSG contains an error message string if the skip 
%   fails, otherwise it is an empty string.

%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2007/07/26 19:31:04 $

    msg = '';
    msgID = '';
    % Determine if pad byte is necessary; % If the chunk size is odd, 
    % there is a pad byte
    if ( rem(chunk.cksize,2) ) 
      pad = 1;
    else 
      pad = 0;
    end

    % Skip the chunk
    status = fseek(fid,chunk.cksize + pad,0);
    if ( status == -1 )
      msg = 'Incorrect chunk size information in AVI file.';
      msgID = 'MATLAB:skipChunk:incorrectChunkSize';
    end
return;

