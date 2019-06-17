function psnr_search(clip_dir, test, results_file, varargin)
% PSNR_SEARCH 'clip_dir' 'test' 'results_file' options
%   Estimate the Y-channel PSNR (PSNR) of all clips and HRCs (Hypothetical 
%   Reference Circuits) in a video test (input argument test) where the 
%   video clips are stored in the specified directory (clip_dir).  The 
%   video clips must have names that conform to the naming convention 
%   test_scene_hrc.yuv, with no extra '_' or '.' in the file names.  "test"
%   is the name of the test, "scene" is the name of the scene, and "hrc" is
%   the name of the HRC.  The name of the original reference clip for the 
%   PSNR calculation must be "test_scene_original.yuv".  PSNR results are
%   output to the file specified by the user.
%
%   Video files can be stored in "Big YUV format", which is a binary format
%   for storing ITU-R Recommendation BT.601 video sequences.  The format can
%   be used for any image size.  In the Big YUV format, all the frames are
%   stored sequentially in one big binary file. The sampling is 4:2:2 and 
%   image pixels are stored sequentially by video scan line as bytes in the
%   following order: Cb1 Y1 Cr1 Y2 Cb3 Y3 Cr3 Y4…, where Y is the 
%   luminance component, Cb is the blue chrominance component, Cr is the 
%   red chrominance component, and the subscript is the pixel number. The 
%   Y signal is quantized into 220 levels where black = 16 and white = 235,
%   while the Cb and Cr signals are quantized into 225 levels with zero 
%   signal corresponding to 128.  Occasional excursions beyond these levels
%   may occur. For example, Y=15 may be produced by a real system, but in 
%   any case, Y, Cb, and Cr values are always between 0 and 255.  The
%   original and processed Big YUV files must have the same number of frames.
%
%   Alternatively, files can be stored in uncompressed UYVY AVI format.  In
%   this case, the file naming convention must be test_scene_hrc.avi, with 
%   no extra '_' or '.' in the file names.  The name of the original
%   reference clip for the PSNR calculation must be "test_scene_original.avi".
%
%   A peak signal of 255 is used for calculation of PSNR.  Double precision
%   calculations are used everywhere.  A 64-bit operating system with at
%   least 4 GB of free memory is recommended since the entire double
%   precision versions of the original and processed sequences must be held
%   in memory.
%
% SYNTAX
%   psnr_search 'clip_dir' 'test' 'results_file' options
%
% DESCRIPTION
%   This function will process all video clips in the user specified
%   clip_dir and test, estimate the Y-channel PSNR of each clip, and then
%   output these PSNR results to the user specified results_file.  If the
%   verbose option is specified, the program will also average the clip
%   results to produce an average PSNR estimate for each HRC.  
%
%   The algorithm performs an exhaustive search for max PSNR over plus or
%   minus the spatial_uncertainty (in pixels) and plus or minus the 
%   temporal_uncertainty (in frames).  The processed video segment is fixed
%   and the original video segment is shifted over the search range.  For
%   each spatial-temporal shift, a linear fit between the processed pixels 
%   and the original pixels is performed such that the mean square error of
%   [original-gain*processed+offset] is minimized (hence maximizing PSNR).
%
%   Any or all of the following optional properties may be requested (the
%   first option is required for yuv files, not avi files).
%
%   'yuv' rows cols    Specifies the number of rows and cols for the Big
%                      YUV files.
%
%   'sroi' top left bottom right    Only use the specified spatial region 
%                                   of interest (sroi) for the PSNR
%                                   calculation.  This is the sroi of the
%                                   processed sequence, which remains fixed
%                                   over all spatial shifts.  By default,
%                                   sroi is the entire image reduced by the
%                                   spatial uncertainty.  If the user
%                                   inputs a sroi, allowance must be made
%                                   for the spatial search specified by
%                                   'spatial_uncertainty'.
%
%   'frames' fstart fstop  Only use the frames from fstart to fstop
%                          (inclusive) to perform the PSNR estimate.  This
%                          specifies the temporal segment of the processed
%                          sequence, which remains fixed over all temporal
%                          shifts.  By default, the temporal segment is the
%                          entire file reduced by the temporal uncertainty.
%                          If the user inputs an fstart and fstop,
%                          allowance must be made for the temporal search
%                          specified by 'temporal_uncertainty'. 
%
%   'fraction_sampled' p   Specifies the fraction of pixels to be sub
%                        sampled.  This affects only the gain and offset
%                        calculation.  This is a deviation from ITU-T 
%                        Recommendation J.340 because it uses a randomly
%                        sub-sampled version of the video sequence.  There 
%                        is a performance gain in speed and considerable
%                        less memory is used and the results seem to agree
%                        closely with using the full video sequence.  The 
%                        p value ranges from (0:1].  However the lowest
%                        recommended value is .001.  By default, this is
%                        set to 1 for 100 percent sampling.
%                          
%   'spatial_uncertainty' x y   Specifies the spatial uncertainty (plus 
%                               or minus, in pixels) over which to 
%                               search.  The processed remains fixed and
%                               the original is shifted.  By default,
%                               this is set to zero.  These x and y shifts
%                               are always with respect to the video FRAME.
%
%   'temporal_uncertainty' t  Specifies the temporal uncertainty
%                             (plus or minus, in frames) over which
%                             to search.  The processed remains fixed
%                             and the original is shifted.  By
%                             default, this is set to zero.
%
%   'video_standard' v    Specifies the video standard of the clips in the
%                         clip_dir. Values for v are:
%                         'interlace_lower_field_first'
%                         'interlace_upper_field_first' 
%                         'progressive' (default)  
%                         If the processed clips are interlaced and have a
%                         chance of having been 'reframed', the
%                         video_standard option should be used and the
%                         spatial_uncertainty y should be set to at least
%                         1.  If the clips are interlaced but do NOT have
%                         reframing, then the default 'progressive'
%                         video_standard can be used. Specifying an
%                         interlaced video_standard will cause the program
%                         to search reframed conditions in addition to
%                         non-reframed conditions.  When reframing occurs,
%                         a + or - 0.5 frames is added to the temporal
%                         shift (as appropriate).
%       Interlace_upper_field_first reframing is explained as follows:
%       (field_2 of frame_n+1 = field_1 of frame_n, field_1 of frame_n+1 =
%       field_2 of frame_n+1, etc.), where field 1 = lines 2, 4, 6 (late)
%                                          field 2 = lines 1, 3, 5 (early)
%       Interlace_lower_field_first reframing is explained as follows:
%       (field_1 of frame_n+1 = field_2 of frame_n, field_2 of frame_n+1 =
%       field_1 of frame_n+1, etc.), where field 1 = lines 2, 4, 6 (early)
%                                      field 2 = lines 1, 3, 5 (late)
%
%   'verbose'   Display progress during processing.  An update line is
%               printed for spatial and temporal shifts that produce
%               an improved PSNR, the average PSNR for the HRC is printed,
%               and for progressive video, a 3D graph of PSNR (as a
%               function of X and Y shift) is generated for temporal
%               shifts that produce a better PSNR (if both x and y
%               spatial_uncertainty are non-zero).
%
%   'full_results'   Generate an output file for each clip that has every
%                    gain, offset, and PSNR calculations for each spatial
%                    and temporal shift that was examined. These files are
%                    stored in the same directory that is specified by
%                    'results_file'.
%
% 
% RESULTS
%   Each line in the results_file contains the following information for
%   each processed clip, in this order, in Comma-separated values (CSV) format: 
%
%   test    The test name for the video clip.
%   scene   The scene name for the video clip.
%   hrc     The HRC name for the video clip.
%   yshift  The y shift for maximum PSNR.  This is in frames.
%   xshift  The x shift for maximum PSNR.  This is in frames.
%   tshift  The time shift for maximum PSNR.  This is in frames, however
%   when a .5 is present, for interlace, this means that there is a field
%   shift and not a whole frame shift.
%   gain    The gain*processed+offset for maximum PSNR.
%   offset
%   psnr    The maximum PSNR observed over the search range.
%
% EXAMPLES
%   These examples illustrate how to call the routine to process the VQEG
%   MM Phase I test scenes, where test scenes from each subjective
%   experiment are stored in a unique directory.  These three examples
%   illustrate the Big YUV format for QCIF, CIF, and VGA.
%
%   psnr_search 'd:\q01\' 'q01' 'q01_psnr.csv' 'yuv' 144 176 'sroi' 5 5 140 172 'spatial_uncertainty' 1 1 'temporal_uncertainty' 8 'verbose' 
%   psnr_search 'd:\c01\' 'c01' 'c01_psnr.csv' 'yuv' 288 352 'sroi' 8 8 281 345 'spatial_uncertainty' 1 1 'temporal_uncertainty' 8 'verbose' 
%   psnr_search 'd:\v01\' 'v01' 'v01_psnr.csv' 'yuv' 480 640 'sroi' 14 14 467 627 'spatial_uncertainty' 1 1 'temporal_uncertainty' 8 'verbose' 
%
%   This example illustrates the uncompressed UYVY AVI format for QCIF.
%
%   psnr_search 'd:\q01\' 'q01' 'q01_psnr.csv' 'sroi' 5 5 140 172 'spatial_uncertainty' 1 1 'temporal_uncertainty' 8 'verbose' 
%

if nargin == 0,
    fprintf('PSNR_SEARCH ''clip_dir'' ''test'' ''results_file'' options\n');
    fprintf('\n');
    fprintf('Estimate the Y-channel PSNR (PSNR) of all clips and HRCs (Hypothetical\n');
    fprintf('Reference Circuits) in a video test (input argument test) where the\n');
    fprintf('video clips are stored in the specified directory (clip_dir).  The\n');
    fprintf('video clips must have names that conform to the naming convention\n');
    fprintf('test_scene_hrc.yuv, with no extra ''_'' or ''.'' in the file names.  "test"\n');
    fprintf('is the name of the test, "scene" is the name of the scene, and "hrc" is\n');
    fprintf('the name of the HRC.  The name of the original reference clip for the\n');
    fprintf('PSNR calculation must be "test_scene_original.yuv".  PSNR results are\n');
    fprintf('output to the file specified by the user.\n');
    fprintf('\n');
    fprintf('Video files can be stored in "Big YUV format", which is a binary format\n');
    fprintf('for storing ITU-R Recommendation BT.601 video sequences.  The format can\n');
    fprintf('be used for any image size.  In the Big YUV format, all the frames are\n');
    fprintf('stored sequentially in one big binary file. The sampling is 4:2:2 and\n');
    fprintf('image pixels are stored sequentially by video scan line as bytes in the\n');
    fprintf('following order: Cb1 Y1 Cr1 Y2 Cb3 Y3 Cr3 Y4…, where Y is the\n');
    fprintf('luminance component, Cb is the blue chrominance component, Cr is the\n');
    fprintf('red chrominance component, and the subscript is the pixel number. The\n');
    fprintf('Y signal is quantized into 220 levels where black = 16 and white = 235,\n');
    fprintf('while the Cb and Cr signals are quantized into 225 levels with zero\n');
    fprintf('signal corresponding to 128.  Occasional excursions beyond these levels\n');
    fprintf('may occur. For example, Y=15 may be produced by a real system, but in\n');
    fprintf('any case, Y, Cb, and Cr values are always between 0 and 255.  The\n');
    fprintf('original and processed Big YUV files must have the same number of frames.\n');
    fprintf('\n');
    fprintf('Alternatively, files can be stored in uncompressed UYVY AVI format.  In\n');
    fprintf('this case, the file naming convention must be test_scene_hrc.avi, with\n');
    fprintf('no extra ''_'' or ''.'' in the file names.  The name of the original\n');
    fprintf('reference clip for the PSNR calculation must be "test_scene_original.avi".\n');
    fprintf('\n');
    fprintf('A peak signal of 255 is used for calculation of PSNR.  Double precision\n');
    fprintf('calculations are used everywhere.  A 64-bit operating system with at\n');
    fprintf('least 4 GB of free memory is recommended since the entire double\n');
    fprintf('precision versions of the original and processed sequences must be held\n');
    fprintf('in memory.\n');
    fprintf('\n');
    fprintf('SYNTAX\n');
    fprintf('psnr_search ''clip_dir'' ''test'' ''results_file'' options\n');
    fprintf('\n');
    fprintf('DESCRIPTION\n');
    fprintf('This function will process all video clips in the user specified\n');
    fprintf('clip_dir and test, estimate the Y-channel PSNR of each clip, and then\n');
    fprintf('output these PSNR results to the user specified results_file.  If the\n');
    fprintf('verbose option is specified, the program will also average the clip\n');
    fprintf('results to produce an average PSNR estimate for each HRC.\n');
    fprintf('\n');
    fprintf('The algorithm performs an exhaustive search for max PSNR over plus or\n');
    fprintf('minus the spatial_uncertainty (in pixels) and plus or minus the\n');
    fprintf('temporal_uncertainty (in frames).  The processed video segment is fixed\n');
    fprintf('and the original video segment is shifted over the search range.  For\n');
    fprintf('each spatial-temporal shift, a linear fit between the processed pixels\n');
    fprintf('and the original pixels is performed such that the mean square error of\n');
    fprintf('[original-gain*processed+offset] is minimized (hence maximizing PSNR).\n');
    fprintf('\n');
    fprintf('Any or all of the following optional properties may be requested (the\n');
    fprintf('first option is required for yuv files, not avi files).\n');
    fprintf('\n');
    fprintf('''yuv'' rows cols    Specifies the number of rows and cols for the Big\n');
    fprintf('                   YUV files.\n');
    fprintf('\n');
    fprintf('''sroi'' top left bottom right    Only use the specified spatial region\n');
    fprintf('                                of interest (sroi) for the PSNR\n');
    fprintf('                                calculation.  This is the sroi of the\n');
    fprintf('                                processed sequence, which remains fixed\n');
    fprintf('                                over all spatial shifts.  By default,\n');
    fprintf('                                sroi is the entire image reduced by the\n');
    fprintf('                                spatial uncertainty.  If the user\n');
    fprintf('                                inputs a sroi, allowance must be made\n');
    fprintf('                                for the spatial search specified by\n');
    fprintf('                                ''spatial_uncertainty''.\n');
    fprintf('\n');
    fprintf('''frames'' fstart fstop  Only use the frames from fstart to fstop\n');
    fprintf('                       (inclusive) to perform the PSNR estimate.  This\n');
    fprintf('                       specifies the temporal segment of the processed\n');
    fprintf('                       sequence, which remains fixed over all temporal\n');
    fprintf('                       shifts.  By default, the temporal segment is the\n');
    fprintf('                       entire file reduced by the temporal uncertainty.\n');
    fprintf('                       If the user inputs an fstart and fstop,\n');
    fprintf('                       allowance must be made for the temporal search\n');
    fprintf('                       specified by ''temporal_uncertainty''.\n');
    fprintf('\n');
    fprintf('\n');
    fprintf('''fraction_sampled'' p  Specifies the fraction of pixels to be sub\n');
    fprintf('                      sampled.  This affects only the gain and offset\n');
    fprintf('                      calculation.  This is a deviation from ITU-T \n');
    fprintf('                      Recommendation J.340 because it uses a randomly\n');
    fprintf('                      sub-sampled version of the video sequence.  There\n');
    fprintf('                      is a performance gain in speed and considerable\n');
    fprintf('                      less memory is used and the results seem to agree\n');
    fprintf('                      closely with using the full video sequence.  The \n');
    fprintf('                      p value ranges from (0:1].  However the lowest\n');
    fprintf('                      recommended value is .001.  By default, this is\n');
    fprintf('                      set to 1 for 100 percent sampling.\n');
    fprintf('\n');
    fprintf('''spatial_uncertainty'' x y   Specifies the spatial uncertainty (plus\n');
    fprintf('                            or minus, in pixels) over which to\n');
    fprintf('                            search.  The processed remains fixed and\n');
    fprintf('                            the original is shifted.  By default,\n');
    fprintf('                            this is set to zero.  These x and y shifts\n');
    fprintf('                            are always with respect to the video FRAME.\n');
    fprintf('\n');
    fprintf('''temporal_uncertainty'' t  Specifies the temporal uncertainty\n');
    fprintf('                          (plus or minus, in frames) over which\n');
    fprintf('                          to search.  The processed remains fixed\n');
    fprintf('                          and the original is shifted.  By\n');
    fprintf('                          default, this is set to zero.\n');
    fprintf('\n');
    fprintf('''video_standard'' v    Specifies the video standard of the clips in the\n');
    fprintf('                      clip_dir. Values for v are:\n');
    fprintf('                      ''interlace_lower_field_first''\n');
    fprintf('                      ''interlace_upper_field_first''or .\n');
    fprintf('                      ''progressive'' (default)\n');
    fprintf('                      If the processed clips are interlaced and have a\n');
    fprintf('                      chance of having been ''reframed'', the\n');
    fprintf('                      video_standard option should be used and the\n');
    fprintf('                      spatial_uncertainty y should be set to at least\n');
    fprintf('                      1.  If the clips are interlaced but do NOT have\n');
    fprintf('                      reframing, then the default ''progressive''\n');
    fprintf('                      video_standard can be used. Specifying an\n');
    fprintf('                      interlaced video_standard will cause the program\n');
    fprintf('                      to search reframed conditions in addition to\n');
    fprintf('                      non-reframed conditions.  When reframing occurs,\n');
    fprintf('                      a + or - 0.5 frames is added to the temporal\n');
    fprintf('                      shift (as appropriate).\n');
    fprintf('    Interlace_upper_field_first reframing is explained as follows:\n');
    fprintf('    (field_2 of frame_n+1 = field_1 of frame_n, field_1 of frame_n+1 =\n');
    fprintf('    field_2 of frame_n+1, etc.), where field 1 = lines 2, 4, 6 (late)\n');
    fprintf('                                       field 2 = lines 1, 3, 5 (early)\n');
    fprintf('    Interlace_lower_field_first reframing is explained as follows:\n');
    fprintf('    (field_1 of frame_n+1 = field_2 of frame_n, field_2 of frame_n+1 =\n');
    fprintf('    field_1 of frame_n+1, etc.), where field 1 = lines 2, 4, 6 (early)\n');
    fprintf('                                   field 2 = lines 1, 3, 5 (late)\n');
    fprintf('\n');
    fprintf('''verbose''     Display progress during processing.  An update line is\n');
    fprintf('              printed for spatial and temporal shifts that produce\n');
    fprintf('              an improved PSNR, the average PSNR for the HRC is printed,\n');
    fprintf('              and for progressive video, a 3D graph of PSNR (as a\n');
    fprintf('              function of X and Y shift) is generated for temporal\n');
    fprintf('              shifts that produce a better PSNR (if both x and y\n');
    fprintf('              spatial_uncertainty are non-zero).\n');
    fprintf('\n');
    fprintf('''full_results''     Generate an output file for each clip that has every\n');
    fprintf('                   gain, offset, and PSNR calculations for each spatial\n');
    fprintf('                   and temporal shift that was examined. These files are\n');
    fprintf('                   stored in the same directory that is specified by\n');
    fprintf('                   ''results_file''.\n');
    fprintf('\n');
    fprintf('\n');
    fprintf('RESULTS\n');
    fprintf('Each line in the results_file contains the following information for\n');
    fprintf('each processed clip, in this order, in Comma-separated values (CSV) format:\n');
    fprintf('\n');
    fprintf('test    The test name for the video clip.\n');
    fprintf('scene   The scene name for the video clip.\n');
    fprintf('hrc     The HRC name for the video clip.\n');
    fprintf('yshift  The y shift for maximum PSNR.\n');
    fprintf('xshift  The x shift for maximum PSNR.\n');
    fprintf('tshift  The time shift for maximum PSNR.\n');
    fprintf('gain    The gain*processed+offset for maximum PSNR.\n');
    fprintf('offset\n');
    fprintf('psnr    The maximum PSNR observed over the search range.\n');
    fprintf('\n');
    fprintf('EXAMPLES\n');
    fprintf('These examples illustrate how to call the routine to process the VQEG\n');
    fprintf('MM Phase I test scenes, where test scenes from each subjective\n');
    fprintf('experiment are stored in a unique directory.  These three examples\n');
    fprintf('illustrate the Big YUV format for QCIF, CIF, and VGA.\n');
    fprintf('\n');
    fprintf('psnr_search ''d:\\q01\\'' ''q01'' ''q01_psnr.csv'' ''yuv'' 144 176 ''sroi'' 5 5 140 172 ''spatial_uncertainty'' 1 1 ''temporal_uncertainty'' 8 ''verbose''\n');
    fprintf('psnr_search ''d:\\c01\\'' ''c01'' ''c01_psnr.csv'' ''yuv'' 288 352 ''sroi'' 8 8 281 345 ''spatial_uncertainty'' 1 1 ''temporal_uncertainty'' 8 ''verbose''\n');
    fprintf('psnr_search ''d:\\v01\\'' ''v01'' ''v01_psnr.csv'' ''yuv'' 480 640 ''sroi'' 14 14 467 627 ''spatial_uncertainty'' 1 1 ''temporal_uncertainty'' 8 ''verbose''\n');
    fprintf('\n');
    fprintf('This example illustrates the uncompressed UYVY AVI format for QCIF.\n');
    fprintf('\n');
    fprintf('psnr_search ''d:\\q01\\'' ''q01'' ''q01_psnr.csv'' ''sroi'' 5 5 140 172 ''spatial_uncertainty'' 1 1 ''temporal_uncertainty'' 8 ''verbose''\n');
    fprintf('\n');
    return;
end

% strip off the extra single quotes ''
clip_dir = eval(clip_dir); 
test = eval(test);
results_file = eval(results_file);

% Define the peak signal level
peak = 255.0;

% Add extra \ in clip_dir in case user did not
clip_dir = strcat(clip_dir,'\');

% Validate input arguments and set their defaults
file_type = 'avi';  % default file type, uncompressed UYVY AVI
is_yuv = 0;
is_whole_image = 1;
is_whole_time = 1;
x_uncert = 0;
y_uncert = 0;
t_uncert = 0;
verbose = 0;
fraction_sampled = 1;
video_standard = 'progressive';  % Default
full_results = 0;
dx=1;  % dx, dy, and dt sizes to use for gain and level offset calculations
dy=1;
dt=1;
cnt=1;
while cnt <= length(varargin),
    if strcmpi(eval(char(varargin(cnt))),'yuv') == 1
        rows = str2double(varargin{cnt+1});
        cols = str2double(varargin{cnt+2});
        is_yuv = 1;
        file_type = 'yuv';
        cnt = cnt + 3;
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
        is_whole_time = 0;
        cnt = cnt + 3;
    elseif strcmpi(eval(char(varargin(cnt))), 'fraction_sampled') ==1
        fraction_sampled = str2double(varargin{cnt+1});
        cnt = cnt + 2; 
    elseif strcmpi(eval(char(varargin(cnt))), 'spatial_uncertainty') ==1
        x_uncert = str2double(varargin{cnt+1});
        y_uncert = str2double(varargin{cnt+2});
        cnt = cnt + 3;
    elseif strcmpi(eval(char(varargin(cnt))), 'temporal_uncertainty') ==1
        t_uncert = str2double(varargin{cnt+1});
        cnt = cnt + 2;
    elseif strcmpi(eval(char(varargin(cnt))),'verbose') == 1
        verbose = 1;
        cnt = cnt +1;
    elseif strcmpi(eval(char(varargin(cnt))),'full_results') ==1;
        full_results = 1;
        cnt = cnt + 1;
    elseif strcmpi(eval(char(varargin(cnt))),'video_standard') ==1
        video_standard = eval(varargin{cnt+1});
        cnt = cnt + 2;
    else
        error('Property value passed into psnr_search not recognized');
    end
end

if (fraction_sampled <= 0 || fraction_sampled > 1)
    error('The value of fraction_sampled is invalid.');
end
if (fraction_sampled < .001)
    warning('Fraction_sampled is below the recommended limit.');
end

% If not progressive and user inputs an SROI, they must have an odd top and
% an even bottom.  Otherwise the field ordering will reverse.
if (~strcmpi(video_standard,'progressive') && ~is_whole_image && (~mod(top,2) || mod(bottom,2)))
    error('SROI top must be odd and bottom must be even for interlaced video.');
end

% For interlaced video, multiple the y_uncert by two since it will later be
% divided by two (i.e., field line shifts are searched rather than frame
% lines).
if(~strcmpi(video_standard,'progressive'))
    y_uncert = 2*y_uncert;
end

%  Get a directory listing
files = dir(clip_dir);  % first two files are '.' and '..'
num_files = size(files,1);

% Find the HRCs and their scenes for the specified video test
hrc_list = {};
scene_list = {};
for i=3:num_files
    this_file = files(i).name;
    und = strfind(this_file,'_'); % find underscores and period
    dot = strfind(this_file,'.');
    if(size(und,2)==2) % possible standard naming convention file found
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
    error('No files with standard naming convention found.\n');
end

%Results struct to store results, shifts are how much the original must be
%shifted with respect to the processed
results = struct('test', {}, 'scene', {}, 'hrc', {}, 'yshift', {}, ...
    'xshift', {}, 'tshift', {}, 'gain', {}, 'offset', {}, 'psnr', {});

% Process one HRC at a time to compute average PSNR for that HRC
index = 1;  % index to store results
fid_results = fopen(results_file,'a');  % open results file for appending
fprintf(fid_results,'Test,Scene,HRC,Yshift,Xshift,Tshift,Gain,Offset,PSNR\n');
fclose(fid_results);
for i = 1:num_hrcs
    
    psnr_ave = 0;  % initialize the psnr average summer for this HRC
    this_hrc = hrc_list{i};
    if(strmatch('original',this_hrc,'exact')) % Don't process original
        continue;
    end
    num_scenes = size(scene_list{i},2);  % Number of scenes in this HRC
    
    for j = 1:num_scenes
        
        this_scene = scene_list{i}{j};
        results(index).test = test;
        results(index).scene = this_scene;
        results(index).hrc = this_hrc;
        
        % Read original and processed video files
        if (~is_yuv)  % YUV file parameters not specified, AVI assummed
            % Re-generate the original and processed avi file names
            orig = strcat(clip_dir, test,'_', this_scene, '_', 'original', '.avi');
            proc = strcat(clip_dir, test,'_', this_scene, '_', this_hrc, '.avi');
            [avi_info] = read_avi('Info',orig);
            [avi_info_proc] = read_avi('Info',proc);
            rows = avi_info.Height;
            cols = avi_info.Width;
            % Set/Validate the SROI
            if (is_whole_image) % make SROI whole image less uncertainty
                top = 1+y_uncert;
                left = 1+x_uncert;
                bottom = rows-y_uncert;
                right = cols-x_uncert;
            elseif (top<1 || left<1 || bottom>rows || right>cols)
                error('Requested SROI too large for image size.\n');
            end
            tframes = avi_info.NumFrames;  % total frames in orig file
            tframes_proc = avi_info_proc.NumFrames;
            % Validate that orig and proc have the same number of frames
            if (tframes ~= tframes_proc)
                fprintf('\n%s_%s_%s: orig & proc files have a different number of frames, the longer file will be truncated.\n', test, this_scene, this_hrc);
                tframes = min(tframes,tframes_proc);
            end
            % Set/Validate the time segment to use
            if (is_whole_time) % use whole time segment less uncertainty
                fstart= 1+t_uncert;
                fstop = tframes-t_uncert;
            elseif (fstart<1 || fstop>tframes)
                error('Requested Temporal segment too large for file size.\n');
            end
            % Validate the spatial uncertainty search bounds
            if (left-x_uncert < 1 || right+x_uncert > cols)
                error('Spatial x-uncertainty too large for SROI.\n');
            end
            if (top-y_uncert < 1 || bottom+y_uncert > rows)
                error('Spatial y-uncertainty too large for SROI.\n');
            end
            % Validate the temporal uncertainty search bounds
            if(fstart-t_uncert < 1 || fstop+t_uncert > tframes)
                error('Temporal uncertainty too large for fstart or fstop.\n');
            end
            if(strcmpi(video_standard,'progressive')) %Reads in normal amount of frames
                % Read in video and clear color planes to free up memory
                [y_orig,cb,cr] = read_avi('YCbCr',orig,'frames',fstart-t_uncert,...
                    fstop+t_uncert, 'sroi',top-y_uncert,left-x_uncert,...
                    bottom+y_uncert,right+x_uncert);
                clear cb cr;
                [y_proc,cb,cr] = read_avi('YCbCr',proc,'frames',fstart,fstop,...
                    'sroi',top,left,bottom,right);
                clear cb cr;
            else %Reads in one less processed frame for interlace
                % Read in video and clear color planes to free up memory
                [y_orig,cb,cr] = read_avi('YCbCr',orig,'frames',fstart-t_uncert,...
                    fstop+t_uncert, 'sroi',top-y_uncert,left-x_uncert,...
                    bottom+y_uncert,right+x_uncert);
                clear cb cr;
                [y_proc,cb,cr] = read_avi('YCbCr',proc,'frames',fstart,fstop-1,...
                    'sroi',top,left,bottom,right);
                clear cb cr;
            end

        else  % YUV file
            % Re-generate the original and processed YUV file name
            orig = strcat(clip_dir, test,'_', this_scene, '_', 'original', '.yuv');
            proc = strcat(clip_dir, test,'_', this_scene, '_', this_hrc, '.yuv');
            % Set/Validate the SROI
            if (is_whole_image) % make SROI whole image less uncertainty
                top = 1+y_uncert;
                left = 1+x_uncert;
                bottom = rows-y_uncert;
                right = cols-x_uncert;
            elseif (top<1 || left<1 || bottom>rows || right>cols)
                error('Requested SROI too large for image size.');
            end
            % Find the total frames of the input original file
            [fid, message] = fopen(orig, 'r');
            if fid == -1
                fprintf(message);
                error('Cannot open this clip''s bigyuv file, %s', orig);
            end
            % Find last frame.
            fseek(fid,0, 'eof');
            tframes = ftell(fid) / (2 * rows * cols);
            fclose(fid);
            % Find the total frames of the processed file
            [fid, message] = fopen(proc, 'r');
            if fid == -1
                fprintf(message);
                error('Cannot open this clip''s bigyuv file, %s', proc);
            end
            % Find last frame.
            fseek(fid,0, 'eof');
            tframes_proc = ftell(fid) / (2 * rows * cols);
            fclose(fid);
            % Validate that orig and proc have the same number of frames
            if (tframes ~= tframes_proc)
                fprintf('\n%s_%s_%s: orig & proc files have a different number of frames, the longer file will be truncated.\n', test, this_scene, this_hrc);
                tframes = min(tframes,tframes_proc);
            end
            % Set/Validate the time segment to use
            if (is_whole_time) % use whole time segment less uncertainty
                fstart= 1+t_uncert;
                fstop = tframes-t_uncert;
            elseif (fstart<1 || fstop>tframes)
                error('Requested Temporal segment too large for file size.');
            end
             % Validate the spatial uncertainty search bounds
            if (left-x_uncert < 1 || right+x_uncert > cols)
                error('Spatial x-uncertainty too large for SROI.');
            end
            if (top-y_uncert < 1 || bottom+y_uncert > rows)
                error('Spatial y-uncertainty too large for SROI.');
            end
            % Validate the temporal uncertainty search bounds
            if(fstart-t_uncert < 1 || fstop+t_uncert > tframes)
                error('Temporal uncertainty too large for fstart or fstop.');
            end
            if(strcmpi(video_standard,'progressive')) %Reads in normal amount of frames
                % Read in video and clear color planes to free up memory
                [y_orig,cb,cr] = read_bigyuv(orig,'frames',fstart-t_uncert,...
                    fstop+t_uncert,'size',rows,cols,'sroi',top-y_uncert,...
                    left-x_uncert,bottom+y_uncert,right+x_uncert);
                clear cb cr;
                [y_proc,cb,cr] = read_bigyuv(proc,'frames',fstart,fstop,...
                    'size',rows,cols,'sroi',top,left,bottom,right);
                clear cb cr;
            else %Reads in one less processed frame for interlace
                % Read in video and clear color planes to free up memory
                [y_orig,cb,cr] = read_bigyuv(orig,'frames',fstart-t_uncert,...
                    fstop+t_uncert,'size',rows,cols,'sroi',top-y_uncert,...
                    left-x_uncert,bottom+y_uncert,right+x_uncert);
                clear cb cr;
                [y_proc,cb,cr] = read_bigyuv(proc,'frames',fstart,fstop-1,...
                    'size',rows,cols,'sroi',top,left,bottom,right);
                clear cb cr;
            end
        end
        
        y_orig = double(y_orig);
        y_proc = double(y_proc);
        [nrows, ncols, nsamps] = size(y_proc);
        
        % Compute PSNR for each spatial-temporal shift
        best_psnr = -inf;
        best_xshift = 0;
        best_yshift = 0;
        best_tshift = 0;
        best_gain = 1;
        best_offset = 0;
        if(verbose)
            fprintf('\nTest = %s,   Scene = %s,   HRC = %s\n',test, this_scene, this_hrc);
        end
        
        if(full_results == 1)
            results_file_full = strcat(test,'_',this_scene,'_',this_hrc,'.csv');
            if(size(strfind(results_file,'/'),2) > 0)
                results_dir = strfind(results_file,'/');
                results_dir = results_file(1:results_dir(size(results_dir,2)));
                results_file_full = strcat(results_dir,results_file_full);
            elseif(size(strfind(results_file,'\'),2) > 0)
                results_dir = strfind(results_file,'\');
                results_dir = results_file(1:results_dir(size(results_dir,2)));
                results_file_full = strcat(results_dir,results_file_full);
            end
            fid_results_full = fopen(results_file_full, 'a');
            fprintf(fid_results_full,'Test,Scene,HRC,Yshift,Xshift,Tshift,Gain,Offset,PSNR\n');
            fclose(fid_results_full);
        end
        
        if(fraction_sampled ~= 1)  %Random Sampling
            rand_nums = round(randperm((nrows*ncols*nsamps))); %Randomizes numbers from 1 to nrows*ncols*nsamps
        end
        
        % Interlace lower field first search
        if (strcmpi(video_standard, 'interlace_lower_field_first'))
            
            %split_into_fields -- y_orig
            [row, col, time_orig] = size(y_orig);
            y_orig = reshape(y_orig,2,row/2,col,time_orig);
            late_field_orig_master = squeeze(y_orig(1,:,:,:));
            early_field_orig_master = squeeze(y_orig(2,:,:,:));
            clear y_orig;
            
            % Reshape y_proc for the PSNR calculation
            y_proc = reshape(y_proc,nrows*ncols*nsamps,1);
            
            y_uncert_field = (y_uncert/2);
            for curr_time = -t_uncert:t_uncert
                for curr_h = -x_uncert:x_uncert
                    for curr_v = -y_uncert_field:y_uncert_field

                        if(mod(curr_v,2)) % Odd - reframing
                            late_field_v = ceil(curr_v/2);
                            early_field_v = floor(curr_v/2);

                            % This code assumes lower field first, so first
                            % lower field is discarded & need one extra
                            % field to make up for this at the end.
                            % no change to length in time.
                            
                            lower_field_orig = late_field_orig_master(...
                                (y_uncert_field + 1 + late_field_v):(y_uncert_field + (nrows/2) + late_field_v),...
                                (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                                (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                            
                            upper_field_orig = early_field_orig_master(...
                                (y_uncert_field + 1 + early_field_v):(y_uncert_field + (nrows/2) + early_field_v),...
                                (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                                (t_uncert + 2 + curr_time):(t_uncert + nsamps + curr_time + 1));
                            
                        else % Even - no reframing
                            lower_field_orig = early_field_orig_master((y_uncert_field + 1 + curr_v/2):(y_uncert_field + (nrows/2) + curr_v/2),...
                                (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                                (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                            
                            upper_field_orig = late_field_orig_master(...
                                (y_uncert_field + 1 + curr_v/2):(y_uncert_field + (nrows/2) + curr_v/2),...
                                (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                                (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                            
                        end
                        
                        % Join into frames and reshape into column vector
                        [row1, col1, time1] = size(lower_field_orig);
                        y_orig(1,:,:,:) = upper_field_orig;
                        y_orig(2,:,:,:) = lower_field_orig;
                        clear lower_field_orig upper_field_orig;
                        y_orig = reshape(y_orig, row1*2, col1, time1);  % 3D
                        y_orig = reshape(y_orig, nrows*ncols*nsamps, 1);  % Column Vector
                        
                        % Perform gain and level offset calculation
                        if(fraction_sampled == 1)
                            this_fit = polyfit(y_proc,y_orig,1);
                        else
                            this_fit = polyfit(y_proc(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),...
                                y_orig(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),1);
                        end
                        
                        % Calculate the PSNR
                        this_psnr = 10*(log10(peak*peak)-log10(sum(((this_fit(1)*y_proc + this_fit(2))-...
                                y_orig).^2)/(nrows*ncols*nsamps)));
                        clear y_orig;
                        
                        if(this_psnr > best_psnr)
                            best_psnr = this_psnr;
                            best_yshift = curr_v;
                            best_xshift = curr_h;
                            %  Correct curr_time for reframed video
                            if(mod(best_yshift,2))
                                best_tshift = curr_time + 0.5;
                            else
                                best_tshift = curr_time;
                            end
                            best_gain = this_fit(1);
                            best_offset = this_fit(2);
                            if(verbose)
                                fprintf('dy =%3i, dx =%3i, dt = %1.1f, gain = %5.4f, offset = %5.4f, PSNR = %5.4f\n',...
                                    best_yshift,best_xshift,best_tshift,best_gain,best_offset,best_psnr);
                            end
                        end

                        if(full_results == 1)
                            gain = this_fit(1);
                            offset = this_fit(2);
                            %  Correct curr_time for reframed video
                            if(mod(curr_v,2))
                                this_tshift = curr_time + 0.5;
                            else
                                this_tshift = curr_time;
                            end
                            fid_results_full = fopen(results_file_full, 'a');
                            fprintf(fid_results_full,'%s, %s, %s, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f\n', ...
                                test, this_scene, this_hrc, curr_v, curr_h, this_tshift, gain, offset, this_psnr);
                            fclose(fid_results_full);
                        end
                        
                    end
                end
            end
            
        elseif (strcmpi(video_standard, 'interlace_upper_field_first'))  % Interlace upper field first search
            
            %split_into_fields -- y_orig
            [row, col, time] = size(y_orig);
            y_orig = reshape(y_orig,2,row/2,col,time);
            early_field_orig_master = squeeze(y_orig(1,:,:,:));
            late_field_orig_master = squeeze(y_orig(2,:,:,:));
            clear y_orig;
            
            % Reshape y_proc for the PSNR calculation
            y_proc = reshape(y_proc,nrows*ncols*nsamps,1);
            
            y_uncert_field = (y_uncert/2);
            for curr_time = -t_uncert:t_uncert
                for curr_h = -x_uncert:x_uncert
                    for curr_v = -y_uncert_field:y_uncert_field

                        if(mod(curr_v,2)) % Odd - reframing
                            early_field_v = ceil(curr_v/2);
                            late_field_v = floor(curr_v/2);

                            % This code assumes upper field first, so first
                            % upper field is discarded & need one extra
                            % field to make up for this at the end.
                            % No change to length in time.
                            
                            lower_field_orig = early_field_orig_master(...
                                (y_uncert_field + 1 + early_field_v):(y_uncert_field + (nrows/2) + early_field_v),...
                                (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                                (t_uncert + 2 + curr_time):(t_uncert + nsamps + curr_time + 1));
                            
                            upper_field_orig = late_field_orig_master(...
                                (y_uncert_field + 1 + late_field_v):(y_uncert_field + (nrows/2) + late_field_v),...
                                (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                                (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                            
                        else % Even - no reframing
                            lower_field_orig = late_field_orig_master(...
                                (y_uncert_field + 1 + curr_v/2):(y_uncert_field + (nrows/2) + curr_v/2),...
                                (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                                (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                            
                            upper_field_orig = early_field_orig_master((y_uncert_field + 1 + curr_v/2):(y_uncert_field + (nrows/2) + curr_v/2),...
                                (x_uncert + 1 + curr_h):(x_uncert + ncols + curr_h),...
                                (t_uncert + 1 + curr_time):(t_uncert + nsamps + curr_time));
                            
                        end
                        
                        % Join into frames and reshape into column vector
                        [row1, col1, time1] = size(upper_field_orig);
                        y_orig(1,:,:,:) = upper_field_orig;
                        y_orig(2,:,:,:) = lower_field_orig;
                        clear lower_field_orig upper_field_orig;
                        y_orig = reshape(y_orig, row1*2, col1, time1);  % 3D
                        y_orig = reshape(y_orig, nrows*ncols*nsamps, 1);  % Column Vector
                        
                        % Perform gain and level offset calculation
                        if(fraction_sampled == 1)
                            this_fit = polyfit(y_proc,y_orig,1);
                        else
                            this_fit = polyfit(y_proc(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),...
                                y_orig(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),1);
                        end
                        
                        % Calculate the PSNR
                        this_psnr = 10*(log10(peak*peak)-log10(sum(((this_fit(1)*y_proc + this_fit(2))-...
                                y_orig).^2)/(nrows*ncols*nsamps)));
                        clear y_orig;
                        
                        if(this_psnr > best_psnr)
                            best_psnr = this_psnr;
                            best_yshift = curr_v;
                            best_xshift = curr_h;
                            %  Correct curr_time for reframed video
                            if(mod(best_yshift,2))
                                best_tshift = curr_time + 0.5;
                            else
                                best_tshift = curr_time;
                            end
                            best_gain = this_fit(1);
                            best_offset = this_fit(2);
                            if(verbose)
                                fprintf('dy =%3i, dx =%3i, dt = %1.1f, gain = %5.4f, offset = %5.4f, PSNR = %5.4f\n',...
                                    best_yshift,best_xshift,best_tshift,best_gain,best_offset,best_psnr);
                            end
                        end

                        if(full_results == 1)
                            gain = this_fit(1);
                            offset = this_fit(2);
                            %  Correct curr_time for reframed video
                            if(mod(curr_v,2))
                                this_tshift = curr_time + 0.5;
                            else
                                this_tshift = curr_time;
                            end
                            fid_results_full = fopen(results_file_full, 'a');
                            fprintf(fid_results_full,'%s, %s, %s, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f\n', ...
                                test, this_scene, this_hrc, curr_v, curr_h, this_tshift, gain, offset, this_psnr);
                            fclose(fid_results_full);
                        end

                    end
                end
            end
            
        else % progressive or interlace with no reframing
            cnt_graph = 1;
            current_best_psnr = 0;
            y_proc = reshape(y_proc,nrows*ncols*nsamps,1);
            
            for curr_time = -t_uncert:t_uncert
                for curr_h = -x_uncert:x_uncert
                    for curr_v = -y_uncert:y_uncert
                        
                        if(fraction_sampled == 1)
                            this_fit = polyfit(y_proc,reshape(y_orig(1+curr_v+y_uncert:curr_v+y_uncert+nrows,...
                                1+curr_h+x_uncert:curr_h+x_uncert+ncols,1+curr_time+t_uncert:curr_time+t_uncert+nsamps),...
                                nrows*ncols*nsamps,1),1);
                        else
                            %Reshapes y_orig
                            reshape_y_orig = reshape(y_orig(1+curr_v+y_uncert:curr_v+y_uncert+nrows, 1+curr_h+x_uncert:curr_h+x_uncert+ncols,...
                                1+curr_time+t_uncert:curr_time+t_uncert+nsamps), nrows*ncols*nsamps, 1);
                            
                            % Perform gain and level offset calculation
                            this_fit = polyfit(y_proc(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),...
                                reshape_y_orig(rand_nums(1:round(nrows*ncols*nsamps*fraction_sampled))),1);
                        end
                        
                        % Calculate the PSNR
                        if(fraction_sampled == 1)
                            this_psnr = 10*(log10(peak*peak)-log10(sum(((this_fit(1)*y_proc + this_fit(2))-...
                                reshape(y_orig(1+curr_v+y_uncert:curr_v+y_uncert+nrows,1+curr_h+x_uncert:curr_h+x_uncert+ncols,...
                                1+curr_time+t_uncert:curr_time+t_uncert+nsamps),nrows*ncols*nsamps,1)).^2)/(nrows*ncols*nsamps)));
                        else
                            this_psnr = 10*(log10(peak*peak)-log10(sum(((this_fit(1)*y_proc + this_fit(2))-...
                                reshape_y_orig).^2)/(nrows*ncols*nsamps)));
                            clear reshape_y_orig;
                        end
                        
                        if(this_psnr > best_psnr)
                            best_psnr = this_psnr;
                            best_yshift = curr_v;
                            best_xshift = curr_h;
                            best_tshift = curr_time;
                            best_gain = this_fit(1);
                            best_offset = this_fit(2);
                            if(verbose)
                                fprintf('dy =%3i, dx =%3i, dt =%3i, gain = %5.4f, offset = %5.4f, PSNR = %5.4f\n',...
                                    best_yshift,best_xshift,best_tshift,best_gain,best_offset,best_psnr);
                            end
                        end
                        x(cnt_graph) = curr_h;
                        y(cnt_graph) = curr_v;
                        t(cnt_graph) = curr_time;
                        psnr_graph(cnt_graph) = this_psnr;
                        
                        if(full_results == 1)
                            gain = this_fit(1);
                            offset = this_fit(2);
                            fid_results_full = fopen(results_file_full, 'a');
                            fprintf(fid_results_full,'%s, %s, %s, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f\n', ...
                                test, this_scene, this_hrc, curr_v, curr_h, curr_time, gain, offset, this_psnr);
                            fclose(fid_results_full);
                        end
                        
                        if (cnt_graph==(2*x_uncert+1)*(2*y_uncert+1) && verbose) %If statement for PSNR Graph
                            %  Only call plot_psnr if both x_uncert and
                            %  y_uncert are non-zero, for 3D plot.
                            if ((x_uncert ~= 0) && (y_uncert ~=0))
                                current_best_psnr = plot_psnr(x,y,psnr_graph,(curr_time),current_best_psnr, verbose);
                            end
                            x = [];
                            y = [];
                            psnr_graph = [];
                            cnt_graph = 0;
                        end
                        cnt_graph = cnt_graph + 1;

                    end
                end
            end
            
        end
        
        results(index).yshift = best_yshift;
        results(index).xshift = best_xshift;
        results(index).tshift = best_tshift;
        results(index).gain = best_gain;
        results(index).offset = best_offset;
        results(index).psnr = best_psnr;
        fid_results = fopen(results_file,'a');  % open results file for appending
        fprintf(fid_results,'%s, %s, %s, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f\n', ...
            results(index).test, results(index).scene, results(index).hrc, ...
            results(index).yshift, results(index).xshift, results(index).tshift, ...
            results(index).gain, results(index).offset, results(index).psnr);
        fclose(fid_results);
        
        psnr_ave = psnr_ave+best_psnr;
        index = index+1;
        
    end
    
    % Compute average PSNR for this HRC
        psnr_ave = psnr_ave/(num_scenes);
        if(verbose)
            fprintf('\nHRC = %s, psnr_ave = %5.4f\n\n',this_hrc, psnr_ave);
        end
        close all; %changed
    
end

function current_best_psnr = plot_psnr(x,y,psnr,t, current_best_psnr, verbose) 
if(verbose)
    max_psnr = max(psnr);
    if(max_psnr < current_best_psnr)
        %Leave function
    else
        [X,Y] = meshgrid(linspace(min(x),max(x)), linspace(min(y),max(y)));
        psnr_test = griddata(x,y,psnr,X,Y,'cubic');
        mesh(X,Y,psnr_test);
        hold on;
        xlabel('Spatial X'), ylabel('Spatial Y'), zlabel('PSNR'), title(sprintf('Time Shift (%d)', t));
        contour3(X, Y, psnr_test);
        plot3(x,y,psnr,'.','markersize',10);
        hold off;
        pause(1);
        current_best_psnr = max_psnr;
    end
end



