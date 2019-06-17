function est_hrc_fdf_rr(clip_dir, test, results_file, varargin);
% EST_HRC_FDF_RR 'clip_dir' 'test' 'results_file' options
%
%   Estimate the fraction dropped frames (FDF) of all HRCs in a video test
%   where the video clips are stored in the specified directory.  The video
%   clips must have names that conform to the standard naming convention 
%   (test_scene_hrc.avi or test_scene_hrc.yuv, with no extra '_' or '.' in 
%   the file names).  If AVI files are used, they must be in the 'YCbCr'
%   format, commonly known as UYVY.  'test' is the name of the test,
%   'scene' is the name of the scene, and 'hrc' is the name of the HRC.
%   The original reference clip's HRC name must be 'original'. FDF results 
%   are output to the file specified by the user.
%
% SYNTAX
%   est_hrc_fdf_rr 'clip_dir' 'test' 'results_file' options
%
% DESCRIPTION
%   This function will process all video clips in the user specified
%   clip_dir and test, estimate the average fdf of each clip, and then
%   average these clip results to produce an estimate for each HRC.  The
%   algorithm is a zero reference algorithm that computes the Temporal
%   Information (TI) difference between successive frames, zeros low motion
%   pixels, sums the TI^2 energy of the highest motion pixels remaining in
%   each frame, and then examines these time history waveforms to locate
%   dropped frames.  The algorithm can be run in RR mode provided the
%   original clips are also present in the specified directory and their
%   names have the form test_scene_original.yuv or test_scene_original.avi.
%
%   Any or all of the following optional properties may be requested (the
%   first option is required for yuv files, but not for avi files since
%   this information is read from the avi header).
%
%   'yuv',rows,cols,fps    Specifies the number of rows, cols, and frames
%                          per second (fps) of the yuv files.
%
%   'sroi',top,left,bottom,right,   Only use the specified spatial region 
%                                   of interest (sroi) for the temporal
%                                   difference calculation.  By default,
%                                   all of the image is used.  The sroi is
%                                   inclusive, where top/left start at 1.
%
%   'frames',fstart,fstop  Only use the frames from fstart to fstop
%                          (inclusive) to perform the fdf estimate.
%                          By default, all frames will be used - this may
%                          cause 'out of memory' errors for resolutions
%                          larger than QVGA or CIF on some computers.
%
%   'rr'    Operate in reduced reference mode where the number of dropped
%           frames in the processed clip is reduced by the number of
%           dropped frames in the original clip for the fdf calculation.
%
%   'm_image',m_image   The M_image parameter, default = 30.
%
%   'f_cut',f_cut       The F_cut parameter, default = 0.02.
%
%   'a',a               The a parameter, default = 2.5.
%
%   'b',b               The b parameter, default = 1.25.
%
%   'c',c               The c parameter, default = 0.1.
%
%   'm_drop',m_drop     The M_drop parameter, default = 0.015.
%
%   'm_dip',m_dip       The M_dip parameter, default = 1.0.
%
%   'a_dip',a_dip       The A_dip parameter, default = 3.0.
%
%   'verbose'   Display output and plots during processing.  This includes
%               the average FDF for each HRC that is processed.
%
% RESULTS
%   Each line in the results_file contains the following information for
%   each processed clip, in this order, in Comma-separated values (CSV format:
%       'test'          The test name for the video clip.
%       'scene'         The scene name for the video clip.
%       'hrc'           The hrc name for the video clip.
%       'fdf'           The fraction of dropped frames for the video clip.
%       'fps'           The effective frames per second after correcting for drops. 
%       'start_frame'   The first frame in the file where drops & dips can
%                       be detected by the algorithm, serves as a reference
%                       point for the 'drops' indices.
%       'drops'         The frame indices of the drops and dips, where 
%                       index = 1 is 'start_frame' in the file.
%
% EXAMPLES
%   est_hrc_fdf_rr 'clip_dir' 'test' 'results_file'
%   est_hrc_fdf_rr 'c:\fdf' 'fdf' 'fdf_results.csv' 'yuv' 144 176 30 'sroi' 5 5 140 172 'verbose'
%

if nargin == 0,
    fprintf('EST_HRC_FDF_RR ''clip_dir'' ''test'' ''results_file'' options\n');
    fprintf('\n');
    fprintf('Estimate the fraction dropped frames (FDF) of all HRCs in a video test\n');
    fprintf('where the video clips are stored in the specified directory.  The video\n');
    fprintf('clips must have names that conform to the standard naming convention\n');
    fprintf('(test_scene_hrc.avi or test_scene_hrc.yuv, with no extra ''_'' or ''.'' in\n');
    fprintf('the file names).  If AVI files are used, they must be in the ''YCbCr''\n');
    fprintf('format, commonly known as UYVY.  ''test'' is the name of the test,\n');
    fprintf('''scene'' is the name of the scene, and ''hrc'' is the name of the HRC.\n');
    fprintf('The original reference clip''s HRC name must be ''original''. FDF results\n');
    fprintf('are output to the file specified by the user.\n');
    fprintf('\n');
    fprintf('SYNTAX\n');
    fprintf('est_hrc_fdf_rr ''clip_dir'' ''test'' ''results_file'' options\n');
    fprintf('\n');
    fprintf('DESCRIPTION\n');
    fprintf('This function will process all video clips in the user specified\n');
    fprintf('clip_dir and test, estimate the average fdf of each clip, and then\n');
    fprintf('average these clip results to produce an estimate for each HRC.  The\n');
    fprintf('algorithm is a zero reference algorithm that computes the Temporal\n');
    fprintf('Information (TI) difference between successive frames, zeros low motion\n');
    fprintf('pixels, sums the TI^2 energy of the highest motion pixels remaining in\n');
    fprintf('each frame, and then examines these time history waveforms to locate\n');
    fprintf('dropped frames.  The algorithm can be run in RR mode provided the\n');
    fprintf('original clips are also present in the specified directory and their\n');
    fprintf('names have the form test_scene_original.yuv or test_scene_original.avi.\n');
    fprintf('\n');
    fprintf('Any or all of the following optional properties may be requested (the\n');
    fprintf('first option is required for yuv files, but not for avi files since\n');
    fprintf('this information is read from the avi header).\n');
    fprintf('\n');
    fprintf('''yuv'',rows,cols,fps    Specifies the number of rows, cols, and frames\n');
    fprintf('                       per second (fps) of the yuv files.\n');
    fprintf('\n');
    fprintf('''sroi'',top,left,bottom,right,   Only use the specified spatial region\n');
    fprintf('                                of interest (sroi) for the temporal\n');
    fprintf('                                difference calculation.  By default,\n');
    fprintf('                                all of the image is used.  The sroi is\n');
    fprintf('                                inclusive, where top/left start at 1.\n');
    fprintf('\n');
    fprintf('''frames'',fstart,fstop  Only use the frames from fstart to fstop\n');
    fprintf('                       (inclusive) to perform the fdf estimate.\n');
    fprintf('                       By default, all frames will be used - this may\n');
    fprintf('                       cause ''out of memory'' errors for resolutions\n');
    fprintf('                       larger than QVGA or CIF on some computers.\n');
    fprintf('\n');
    fprintf('''rr''    Operate in reduced reference mode where the number of dropped\n');
    fprintf('        frames in the processed clip is reduced by the number of\n');
    fprintf('        dropped frames in the original clip for the fdf calculation.\n');
    fprintf('\n');
    fprintf('''m_image'',m_image   The M_image parameter, default = 30.\n');
    fprintf('\n');
    fprintf('''f_cut'',f_cut       The F_cut parameter, default = 0.02.\n');
    fprintf('\n');
    fprintf('''a'',a               The a parameter, default = 2.5.\n');
    fprintf('\n');
    fprintf('''b'',b               The b parameter, default = 1.25.\n');
    fprintf('\n');
    fprintf('''c'',c               The c parameter, default = 0.1.\n');
    fprintf('\n');
    fprintf('''m_drop'',m_drop     The M_drop parameter, default = 0.015.\n');
    fprintf('\n');
    fprintf('''m_dip'',m_dip       The M_dip parameter, default = 1.0.\n');
    fprintf('\n');
    fprintf('''a_dip'',a_dip       The A_dip parameter, default = 3.0.\n');
    fprintf('\n');
    fprintf('''verbose''   Display output and plots during processing.  This includes\n');
    fprintf('            the average FDF for each HRC that is processed.\n');
    fprintf('\n');
    fprintf('RESULTS\n');
    fprintf('Each line in the results_file contains the following information for\n');
    fprintf('each processed clip, in this order, in Comma-separated values (CSV format:\n');
    fprintf('    ''test''          The test name for the video clip.\n');
    fprintf('    ''scene''         The scene name for the video clip.\n');
    fprintf('    ''hrc''           The hrc name for the video clip.\n');
    fprintf('    ''fdf''           The fraction of dropped frames for the video clip.\n');
    fprintf('    ''fps''           The effective frames per second after correcting for drops.\n');
    fprintf('    ''start_frame''   The first frame in the file where drops & dips can\n');
    fprintf('                    be detected by the algorithm, serves as a reference\n');
    fprintf('                    point for the ''drops'' indices.\n');
    fprintf('    ''drops''         The frame indices of the drops and dips, where\n');
    fprintf('                    index = 1 is ''start_frame'' in the file.\n');
    fprintf('\n');
    fprintf('EXAMPLES\n');
    fprintf('est_hrc_fdf_rr ''clip_dir'' ''test'' ''results_file''\n');
    fprintf('est_hrc_fdf_rr ''c:\\fdf'' ''fdf'' ''fdf_results.csv'' ''yuv'' 144 176 30 ''sroi'' 5 5 140 172 ''verbose''\n');
    fprintf('\n');
    return;
end

% strip off the extra single quotes ''
clip_dir = eval(clip_dir); 
test = eval(test);
results_file = eval(results_file);

% Parameters that control the dropped frames detection algorithm.

% TI image pixels at or below this magnitude are not considered motion and
% are set equal to zero.  Thus, a pixel must have a TI difference magnitude
% greater than this level to contribute to the motion energy of the frame.
m_image = 30;

% If the average TI^2 value of a frame is less than or equal to a motion
% threshold, m_drop*dfact (a dynamic factor that is a function of the mean 
% TI^2 value computed over all frames in the video sequence), then that
% frame is declared a frame repeat.  
m_drop = 0.015;

% Simple thresholding does not catch all the perceptual dropped frames when
% there is some residual motion (small blocks being updated) since the TI^2
% waveform does not dip below the motion threshold.  If this condition 
% occurs for single frame drops, they can be detected with a dip
% detector.  In order to be declared a valid dip, these negative 1-frame
% width dips must have a dip amplitude of at least a_dip*dfact, and the
% average TI^2 value for the frame must be less than or equal to
% m_dip*dfact.
a_dip = 3;
m_dip = 1;

% Constants that control the computation of the dynamic factor (dfact)
a = 2.5;  % DC shift
b = 1.25;  % Gain
c = 0.1;  % low end limit

% Fraction of scene cuts to eliminate before computing TI2_ave
f_cut = 0.02;

% Disable output warnings
warning off all;

% Add extra \ in clip_dir in case user forgot
clip_dir = strcat(clip_dir,'\');

% Validate input arguments and set their defaults
file_type = 'avi';  % default file type, uncompressed UYVY AVI
is_yuv = 0;
is_whole_image = 1;
is_whole_time = 1;
is_rr = 0;
verbose = 0;
cnt=1;
while cnt <= length(varargin),
    if ~ischar(varargin{cnt}),
        error('Property value passed into est_hrc_fdf_rr is not recognized');
    end
    if strcmpi(eval(char(varargin(cnt))),'yuv') == 1
        rows = str2double(varargin{cnt+1});
        cols = str2double(varargin{cnt+2});
        fps = str2double(varargin{cnt+3});
        is_yuv = 1;
        file_type = 'yuv';
        cnt = cnt + 4;
    elseif strcmpi(eval(char(varargin(cnt))),'sroi') == 1
        top = str2double(varargin{cnt+1});
        left = str2double(varargin{cnt+2});
        bottom = str2double(varargin{cnt+3});
        right = str2double(varargin{cnt+4});
        is_whole_image = 0;
        cnt = cnt + 5;
    elseif strcmpi(eval(char(varargin(cnt))),'frames') == 1
        fstart = str2double(varargin{cnt+1});
        fstop = str2double(varargin{cnt+2});
        nframes = fstop-fstart+1;
        is_whole_time = 0;
        cnt = cnt + 3;
    elseif strcmpi(eval(char(varargin(cnt))),'rr') == 1
        is_rr = 1;
        cnt = cnt + 1;
    elseif strcmpi(eval(char(varargin(cnt))),'m_image') == 1
        m_image = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'f_cut') == 1
        f_cut = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'a') == 1
        a = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'b') == 1
        b = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'c') == 1
        c = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'m_drop') == 1
        m_drop = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'m_dip') == 1
        m_dip = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'a_dip') == 1
        a_dip = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'verbose') == 1
        verbose = 1;
        cnt = cnt + 1;
    else
        error('Property value passed into est_hrc_fdr_rr not recognized');
    end
end

files = dir(clip_dir);  % first two files are '.' and '..'
num_files = size(files,1);

% Find the HRCs and their associated scenes
hrc_list = {};
scene_list = {};
for i=3:num_files
    this_file = files(i).name;
    und = strfind(this_file,'_'); % find underscores and period
    dot = strfind(this_file,'.');
    if(und) % possible standard naming convention file found
        this_test = this_file(1:und(1)-1);  % pick off the test name
        if(~isempty(strmatch(test,this_test,'exact')) && ...
                ~isempty(strmatch(file_type,this_file(dot+1:length(this_file)),'exact')))  % test clip found
            this_scene = this_file(und(1)+1:und(2)-1);
            this_hrc = this_file(und(2)+1:dot(1)-1);
            % See if this HRC already exists and find its list location
            loc = strmatch(this_hrc,hrc_list,'exact');
            if(loc)  % HRC already present, add to scene list for that HRC
                if(size(strmatch(this_scene,scene_list{loc},'exact'),1)==0)
                    scene_list{loc} = [scene_list{loc} this_scene];
                end
            else  % new HRC found
                hrc_list = [hrc_list;{this_hrc}];
                this_loc = size(hrc_list,1);
                scene_list(this_loc) = {{this_scene}};
            end
        end
    end
end

scene_list = scene_list';
num_hrcs = size(hrc_list,1);
if (num_hrcs == 0)
    fprintf('No files with standard naming convention found.\n');
    return
end

%  Reorganize the HRC and scene lists to process the originals first, if
%  they are present.
orig_loc = strmatch('original',hrc_list,'exact');
if (~isempty(orig_loc))
    if (orig_loc > 1)
        top_hrc_list = [hrc_list(orig_loc);hrc_list(1:orig_loc-1)];
        top_scene_list = [scene_list(orig_loc);scene_list(1:orig_loc-1)];
    else
        top_hrc_list = hrc_list(orig_loc);
        top_scene_list = scene_list(orig_loc);
    end
    if (orig_loc < num_hrcs)
        hrc_list = [top_hrc_list; hrc_list(orig_loc+1:num_hrcs)];
        scene_list = [top_scene_list; scene_list(orig_loc+1:num_hrcs)];
    else
        hrc_list = top_hrc_list;
        scene_list = top_scene_list;
    end
end

%Results struct
results = struct('test', {}, 'scene', {}, 'hrc', {}, 'start_frame', {}, 'ti2', {}, ...
    'drops', {}, 'fdf', {}, 'fps', {});

% Process one HRC at a time to compute average fdf
index = 1;  % index used to store results
fid_results = fopen(results_file,'a');  % open results file for appending
fprintf(fid_results,'Test,Scene,HRC,FDF,FPS,Start_Frame,Drops\n');
fclose(fid_results);

% Hold the number of drops and dips in the original for RR mode of operation
orig_scene_list = {};  % contains the original scene name
orig_dropnum = [];  % contains the number of dropped frames in the original
orig_index = 0;  % index used to access original dropped frame information

for i = 1:num_hrcs
    
    fdf_ave = 0;  % initialize the fdf average summer for this HRC
    fdf_scenes = 0;  % initialize number of scenes summer for this HRC
    this_hrc = hrc_list{i};
    num_scenes = size(scene_list{i},2);  % number of scenes in this HRC
    
    for j = 1:num_scenes
        this_scene = scene_list{i}{j};
        results(index).test = test;
        results(index).scene = this_scene;
        results(index).hrc = this_hrc;
        
        % Read video file
        if (~is_yuv)  % AVI file
            
            % Re-generate the processed avi file name
            proc = strcat(clip_dir, test,'_', this_scene, '_', this_hrc, '.avi');
            [avi_info] = read_avi('Info',proc);
            rows = avi_info.Height;
            cols = avi_info.Width;
            
            % Set/Validate the ROI
            if (is_whole_image)
                top = 1;
                left = 1;
                bottom = rows;
                right = cols;
            elseif (top<1 || left<1 || bottom>rows || right>cols)
                display('Requested SROI too large for image size.\n');
                return;
            end
            
            fps = avi_info.FramesPerSecond;  % frames per second of input file
            tframes = avi_info.NumFrames;  % total frames in processed file
            
            if (is_rr)
                orig = strcat(clip_dir, test,'_', this_scene, '_', 'original', '.avi');
                [avi_info_orig] = read_avi('Info',orig);
                tframes_orig = avi_info_orig.NumFrames;
                % Validate that orig and proc have the same number of frames
                if (tframes ~= tframes_orig)
                    display('The orig & proc files must have the same length for the rr option.');
                    return
                end
            end
            
            % Set/Validate the time segment to use
            if (is_whole_time)
                fstart= 1;
                fstop = tframes;
                nframes = tframes;
            elseif (fstart<1 || fstop>tframes)
                display('Requested temporal segment is too large for file size.\');
                return;
            end
            
            % Read in video and clear color planes to free up memory
            [y,cb,cr] = read_avi('YCbCr',proc,'frames',fstart,fstop,...
                'sroi',top,left,bottom,right);
            clear cb cr;
            
        else  % YUV file
            
            % Re-generate the processed YUV file name
            proc = strcat(clip_dir, test,'_', this_scene, '_', this_hrc, '.yuv');
            
            % Set/Validate the ROI
            if (is_whole_image)
                top = 1;
                left = 1;
                bottom = rows;
                right = cols;
            elseif (top<1 || left<1 || bottom>rows || right>cols)
                display('Requested SROI too large for image size.\n');
                return;
            end

            % Find the total frames of the input file
            [fid, message] = fopen(proc, 'r');
            if fid == -1
                fprintf(message);
                error('Cannot open this clip''s bigyuv file, %s', proc);
            end
            % Find last frame.
            fseek(fid,0, 'eof');
            tframes = ftell(fid) / (2 * rows * cols);
            fclose(fid);
            
            if (is_rr)
                orig = strcat(clip_dir, test,'_', this_scene, '_', 'original', '.yuv');
                % Find the total frames of the original file
                [fid, message] = fopen(orig, 'r');
                if fid == -1
                    fprintf(message);
                    error('Cannot open this clip''s bigyuv file, %s', orig);
                end
                % Find last frame.
                fseek(fid,0, 'eof');
                tframes_orig = ftell(fid) / (2 * rows * cols);
                fclose(fid);
                % Validate that orig and proc have the same number of frames
                if (tframes ~= tframes_orig)
                    display('The orig & proc files must have the same length for the rr option.');
                    return
                end
            end
            
            % Set/Validate the time segment to use
            if (is_whole_time)
                fstart= 1;
                fstop = tframes;
                nframes = tframes;
            elseif (fstart<1 || fstop>tframes)
                display('Requested temporal segement is too large for file size.\n');
            end

            % Read in video and clear color planes to free up memory
            [y,cb,cr] = read_bigyuv(proc,'frames',fstart,fstop,...
                'size',rows,cols,'sroi',top,left,bottom,right);
            clear cb cr;
        end
        
        results(index).start_frame = fstart + 1;  % ti(1) at start_frame
        
        % Compute TI over the sroi and time segment
        ti = y(:,:,2:nframes)-y(:,:,1:nframes-1);
        clear y;
        [nrows, ncols, nsamps] = size(ti);
        
        %  Zero pixels below m_image (image motion threshold)
        ti((abs(ti) <= m_image)) = 0.0;
        
        %  Convert to energy and take mean of each frame
        ti2 = mean(reshape(ti,nrows*ncols,nsamps).^2);
        clear ti;
        results(index).ti2 = ti2;
        
        % Calculate the average TI^2 value for the clip, but eliminate 
        % scene cuts and 
        ti2_sort = sort(ti2);
        ti2_ave = mean(ti2_sort(ceil(f_cut*nsamps):floor((1-f_cut)*nsamps)));

        % dynamic thresholding scheme based on ti2_ave
        dfact = a+b*log(ti2_ave);

        % Limit dfact to positive numbers
        if (dfact < c)
            dfact = c;
        end

        % Compute the dynamic thresholds for this video clip
        this_m_drop = dfact*m_drop;
        this_a_dip = dfact*a_dip;
        this_m_dip = dfact*m_dip;

        %  Find dropped frames that fall below motion threshold
        drops = find(ti2 <= this_m_drop);
        
        %  Find additional drops with dip detector, which is the same as finding
        %  the peaks of the negative waveform.
        [dips_mag, dips] = findpeaks2(-1*ti2,'threshold',this_a_dip);
        dips_mag = -1*dips_mag;
        
        % Eliminate dips that are above this_m_dip
        dips = dips(find(dips_mag <= this_m_dip));
        
        % Calculate total of drops and dips, eliminating redundancies
        drops_dips = union(drops,dips);
        results(index).drops = drops_dips;
        
        % Check if original and save number of dropped frames
        if (strmatch('original',this_hrc,'exact'))
            orig_index = orig_index+1;  % index for this original clip
            orig_scene_list{orig_index} = this_scene;
            orig_dropnum(orig_index) = length(drops_dips);
        end

        % calculate the fraction of dropped frames
        this_fdf = length(drops_dips)/(nsamps-2);
        this_fps = fps*(1-this_fdf);  % convert to fps if desired
        if (is_rr)
            % Find the corresponding original dropnum and compensate fdf
            this_match = strmatch(this_scene,orig_scene_list,'exact');
            if(~isempty(this_match))
                orig_fdf = orig_dropnum(this_match)/(nsamps-2);
                if (orig_fdf > 0.9)
                    this_fdf = 'NaN';
                    this_fps = 'NaN';
                else
                    this_fdf = max(0,(this_fdf-orig_fdf)/(1-orig_fdf));
                    this_fps = fps*(1-this_fdf);
                end
            else
                fprintf('Original scene %s is missing.\n', this_scene);
                return
            end
        end

        % Store the fdf information for this clip
        results(index).fdf = this_fdf;
        results(index).fps = this_fps;
        fid_results = fopen(results_file,'a');  % open results file for appending
        fprintf(fid_results,'%s, %s, %s, %5.4f, %5.4f, %5.0f', ...
            results(index).test, results(index).scene, results(index).hrc, ...
            results(index).fdf, results(index).fps, results(index).start_frame);
        this_drops = [results(index).drops];
        if (~isempty(this_drops))
            for k = 1:size(this_drops,2)
                fprintf(fid_results,',%5.0f',this_drops(k));
            end
            fprintf(fid_results,'\n');
        else
            fprintf(fid_results,'\n');
        end
        fclose(fid_results);
        index = index+1;
        
        %  Don't include FDF's that are NAN in the HRC average FDF
        if (~ischar(this_fdf))
            fdf_ave = fdf_ave + this_fdf;
            fdf_scenes = fdf_scenes+1;
        end
        
        pause(0.1);  % to allow user to stop program
        if(verbose)
            plot(ti2, 'k', 'LineWidth', 2);
            grid on;
            set(gca,'LineWidth',1)
            set(gca,'FontName','Ariel')
            set(gca,'fontsize',12)
            hold on
            plot(dips, ti2(dips), '.g', 'markersize', 10);
            plot(drops, ti2(drops), '.r', 'markersize', 10);
            xlabel('Frame');
            ylabel('TI2');
            title('Detected Frame Drops and Dips in Red and Green');
            hold off
            if (~ischar(this_fdf))
                fprintf('Test = %s,   Scene = %s,   HRC = %s,   fdf = %4.2f\n',...
                    test, this_scene, this_hrc, this_fdf);
            else
                fprintf('Test = %s,   Scene = %s,   HRC = %s,   fdf = %s\n',...
                    test, this_scene, this_hrc, this_fdf);
            end
            pause(1);
        end
        
    end
    
    %  Compute average fdf for this HRC if it is available
    if (fdf_scenes >= 1)
        fdf_ave = fdf_ave/fdf_scenes;
        if (verbose)
            fprintf('HRC = %s, fdf_ave = %4.2f\n', this_hrc, fdf_ave);
        end
        fprintf('\n');
    end
    
end

close