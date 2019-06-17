function [new_clip_structs] = bvqm_psnr_search_gclips(test_structs, clip_structs, results_file_csv_name, results_dir_path, results_parameter_name, varargin)
% [new_clip_structs] = bvqm_psnr_search_gclips(test_structs, clip_structs, results_file_csv_name, results_dir_path, results_parameter_name, options) 
%
%   Estimate the Y-channel PSNR (Peak Signal to Noise Ratio) of all clips
%   and HRCs (Hypothetical Reference Circuits) in a video test (input
%   argument clip_structs) where the video clips are stored in the
%   specified directory (test_structs).  The video clips must have names
%   that conform to the naming convention test_scene_hrc.yuv, with no extra
%   '_' or '.' in the file names.  "test" is the name of the test, "scene"
%   is the name of the scene, and "hrc" is the name of the HRC.  The name
%   of the original reference clip for the PSNR calculation must be
%   "test_scene_original.yuv".  Please note that 'test_structs' needs to be
%   in the same format as GTests in "Gclips".  Also, the input
%   'clip_structs' needs to be in the same format as "Gclips".  These
%   structures can be generated from pulling from the
%   existing "Gclips" or by using the function "initialize_clip_struct"
%   which will generate both a gclips test and clip structure from clips in
%   a folder.  PSNR results are output to the file name specified by the
%   user (input argument results_file_csv_name).  'results_file_csv_name' can
%   be any name desired as long as it has the file extinsion ".csv".  This
%   will store the results of the PSNR search in a csv format for easy
%   viewing.  The input variable 'results_dir_path' is where all of the
%   results will be stored as 
%   well as some temporary variables which will allow the program to
%   recover from a crash.  These temporary variables will be deleted once
%   the program has fully finished running.  The input argument
%   "results_parameter_name" is the file name of the PSNR model values
%   which will be saved under the 'results_dir_path'.  This file will be a
%   ".mat" format which will allow matlab to easily use the results if
%   desired.  Finally, the results that will be
%   stored under this 'results_dir_path' will be the following.  Assuming that all
%   logging and results flags are called (see the description section
%   below), the following files will be left in the 'results_dir_path':
%   log_file.txt, "results_file_csv_name".csv (where "results_file_csv_name" will be the name
%   that the user entered into that input argument), and "results_parameter_name".mat (where
%   "results_parameter_name" will be the name that the user entered into that input
%   argument).  If the flag "full_results" is used, then there will be a
%   file corresponding to each of the clips that was tested showing the 
%   full PSNR results for every search.  If there appear to be more files in
%   this directory than listed above, this means that the program most
%   likely crashed and needs to be restarted to finish fully running.  Once
%   this is done, those temporary files will be deleted.  For example, if
%   "tshifts.mat" is present after the program is done running, this
%   indicates a crash has orrcured as this file should no longer exist
%   after completion of the program.  Simply restart the program to allow
%   the program to finish running.
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
%   This algorithm now uses "Gclips" (which is the input argument
%   'clip_structs'), unless the option 'uncalibrated' is
%   selected.  This means that while running in 'calibrated' mode (which is
%   the default setting) the search will search spatially and temporally from
%   where "Gclips" thinks it should be aligned.  Meaning, if the
%   "Gclips" for the currect clip states that the spatial shift in the
%   horizontal direction is 1 and a time shift of 4, then when the video
%   clip is read into memory, these shifts will be taken into account and
%   the video clip in memory will be alinged to the shifts that are
%   reported.  Then the program will search plus and minus spatial and
%   temporal uncertanty from that point.  When 'uncalibrated' is selected,
%   the program ignores what is in "Gclips" and reads in the video as
%   if there is no time shifts or spatial shifts.  Please note that PSNR
%   search has been modified to report the results of the program in terms
%   according to what is in "Gclips".  Therefore, x and y shifts, tshifts,
%   gain and offsets have been altered in order to match what is reported
%   in "Gclips".
%
% SYNTAX
%   'new_clip_structs' =  PSNR_SEARCH 'test_structs' 'clip_structs' ...
%                               'results_file_csv_name' 'results_dir_path' 'results_parameter_name' options
%
% DESCRIPTION
%   This function will process all video clips in the user specified
%   clip_structs and test_structs, estimate the Y-channel PSNR of each
%   clip, and then output these PSNR results to the user specified
%   results_file_csv_name and PSNR parameter file.  If the 
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
%   Any or all of the following optional properties may be requested.
%
%
%   'sroi' top left bottom right    Only use the specified spatial region 
%                                   of interest (sroi) for the PSNR
%                                   calculation.  This is the sroi of the
%                                   processed sequence, which remains fixed
%                                   over all spatial shifts.  By default,
%                                   the values in "Gclips" will be used.
%                                   Whether or not the user enters a SROI,
%                                   the values of 'spatial_uncertainty'
%                                   will be taken into account and the SROI
%                                   will be checked so proper searching can
%                                   occur.  Also, if the user enters too
%                                   big or too small of values, the program
%                                   will give an error.  If the flag 
%                                   'calibrated' is used, SROI will be
%                                   loaded from "Gclips" automatically.  
%                                   If 'uncalibrated' is used, SROI will be
%                                   generated from default_sroi and stored
%                                   in the "Gclips" structure.  If
%                                   the 'sroi' flag is used, this over
%                                   writes the SROI that is in "Gclips".
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
%                             (plus or minus, in seconds) over which
%                             to search.  The processed remains fixed
%                             and the original is shifted.  By
%                             default, this is set to zero.
%
%   'verbose'   Display progress during processing.  An update line is
%               printed for spatial and temporal shifts that produce
%               an improved PSNR, the average PSNR for the HRC is printed,
%               and for progressive video, a 3D graph of PSNR (as a
%               function of X and Y shift) is generated for temporal
%               shifts that produce a better PSNR (if both x and y
%               spatial_uncertainty are non-zero).
%
%   'silent'    Does not display any progress during processing.
%               Basically, nothing is printed to the screen during
%               processing.  By default, this setting is active.
%
%   'full_results'   Generate an output file for each clip that has every
%                    gain, offset, and PSNR calculations for each spatial
%                    and temporal shift that was examined. These files are
%                    stored in the directory entered in the input argument
%                    'results_dir_path'.
%
%   'calibrated'    This flag should be used when the data in "GClips" is
%                   vaild.  Meaning that, the spatial and temporal shifts,
%                   spatial scaling, alignment, gain/offset, and valid
%                   region SROI is set correctly inside of "GClips", which
%                   is being passed into this function as the variable
%                   'clip_structs'.  By default, this option is selected.
%
%   'uncalibrated'  This flag should be used when the data in "GClips" is
%                   not valid.  Meaning that, the spatial and temporal
%                   shifts, spatial scaling, alignment, gain/offset, and 
%                   valid region SROI is NOT set correctly inside of
%                   "GClips".  This flag should also be used when the
%                   clipset has not been calibrated before.
%
%   'log_file'      This flag should be used when the user desires an ASCI
%                   results file.  This 'log_file' contains the same
%                   information as the "results_file_csv_name".csv but instead of
%                   being a csv file, it stores the information in an ASCI
%                   format.
%
% 
% RESULTS
%   Each line in the results_file_csv_name contains the following information for
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
%   These examples illustrate how to call the routine.
%
%   new_clip_struct = psnr_search(GTests_E,q01_clipset,'q01_results.csv','E:\temp\','q01_results','yuv',144,176,'sroi',5,5,140,172,'spatial_uncertainty',1,1,'temporal_uncertainty',.25,'verbose','calibrated' 
%   new_clip_struct = psnr_search(GTests_E,q01_clipset,'q01_results.csv','E:\temp\','q01_results','yuv',144,176,'spatial_uncertainty',2,1,'temporal_uncertainty',.3,'verbose','calibrated'
%   new_clip_struct = psnr_search(GTests_E,q01_clipset,'q01_results.csv','E:\temp\','q01_results','yuv',144,176,'spatial_uncertainty',4,4,'temporal_uncertainty',.1,'verbose','uncalibrated' 
%
%   This example illustrates the uncompressed UYVY AVI format for QCIF.
%
%   new_clip_struct = psnr_search(GTests_E,q01_clipset,'q01_results.csv','E:\temp\','q01_results','spatial_uncertainty',1,1,'temporal_uncertainty',.4,'verbose','calibrated'
%

% Define the peak signal level
peak = 255.0;

% Validate input arguments and set their defaults
is_whole_image = 1; 
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
user_set_sroi = 0;
calibrated = 1;
uncalibrated = 0;
completed = 0;
log_file = 0;
while cnt <= length(varargin),
    if strcmpi(char(varargin(cnt)),'sroi') == 1
        top = varargin{cnt+1};
        left = varargin{cnt+2};
        bottom = varargin{cnt+3};
        right = varargin{cnt+4};
        is_whole_image = 0;
        user_set_sroi = 1;
        cnt = cnt + 5;
    elseif strcmpi(char(varargin(cnt)), 'fraction_sampled') ==1
        fraction_sampled = varargin{cnt+1};
        cnt = cnt + 2; 
    elseif strcmpi(char(varargin(cnt)), 'spatial_uncertainty') ==1
        x_uncert = varargin{cnt+1};
        y_uncert = varargin{cnt+2};
        cnt = cnt + 3;
    elseif strcmpi(char(varargin(cnt)), 'temporal_uncertainty') ==1
        t_uncert = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(char(varargin(cnt)),'verbose') == 1
        verbose = 1;
        cnt = cnt +1;
    elseif strcmpi(char(varargin(cnt)),'silent') == 1
        verbose = 0;
        cnt = cnt + 1;
    elseif strcmpi(char(varargin(cnt)),'full_results') ==1;
        full_results = 1;
        cnt = cnt + 1;
    elseif strcmpi(char(varargin(cnt)),'calibrated') == 1
        calibrated = 1;
        cnt = cnt + 1;
    elseif strcmpi(char(varargin(cnt)),'uncalibrated') == 1
        calibrated = 0;
        uncalibrated = 1;
        cnt = cnt + 1;
    elseif strcmpi(char(varargin(cnt)),'log_file') == 1
        log_file = 1;
        cnt = cnt + 1;
    else
        error('Property value passed into psnr_search not recognized');
    end
end

if (fraction_sampled <= 0 || fraction_sampled > 1)
    error('The value of fraction_sampled is invalid.');
end
if (fraction_sampled < .001)
    if(verbose == 1)
        warning('Fraction_sampled is below the recommended limit.');
    end
end

%Assume that the video standard is the same through the whole clipset.
video_standard = clip_structs(1).video_standard;

%Check that this assumption is true
for i = 2:length(clip_structs)
    if(~strcmpi(clip_structs(i).video_standard,video_standard))
        error('Video Standard is not consistent over the clipset.');
    end
end

% For interlaced video, multiple the y_uncert by two since it will later be
% divided by two (i.e., field line shifts are searched rather than frame
% lines).
if(~strcmpi(video_standard,'progressive'))
    y_uncert = 2*y_uncert;
end

%  Get a directory listing and file info
% files = dir(clip_dir);  % first two files are '.' and '..'
num_files = size(clip_structs,2);

%Convert t_uncert from seconds to frames (round using the ceil function)
t_uncert = ceil(clip_structs(1).fps * t_uncert);

%Check that all fps's are the same over the entire clipset
for i = 2:length(clip_structs)
    if(clip_structs(i).fps ~= clip_structs(1).fps)
        error('Frames Per Second is not consistent over the clipset.');
    end
end

%Check that all SROI's are valid for the clipset if interlaced video is
%present. Otherwise the field ordering will reverse.
if(~strcmpi(video_standard,'progressive'))
    for i = 1:length(clip_structs)
        if(~mod(clip_structs(i).cvr.top,2) || mod(clip_structs(i).cvr.bottom,2))
            if(verbose == 1)
                warning('SROI top must be odd and bottom must be even for interlaced video. SROI will be adjusted!');
            end
            temp = adjust_requested_sroi(clip_structs(i),'evenodd');
            clip_structs(i).cvr.top = temp.top;
            clip_structs(i).cvr.left = temp.left;
            clip_structs(i).cvr.bottom = temp.bottom;
            clip_structs(i).cvr.right = temp.right;
        end
    end
end

% Find the HRCs and their scenes for the specified video test
hrc_list = {};
scene_list = {};
hrc_count = 1;
for i = 1:num_files
    %In order to preserve order in hrc_list, unique can not be used.
    %Therefore, a simple algorithm in order to keep duplicates out is
    %below.
    hrc_list_temp(i) = clip_structs(i).hrc;
    hrc_list(hrc_count) = clip_structs(i).hrc;
    temp = strcmp(hrc_list,hrc_list_temp(i));
    if(max(temp(1,1:size(temp,2)-1)))
        %The current HRC is already in the hrc_list.  It needs to be
        %removed.
        hrc_list(hrc_count) = [];
    else
        hrc_count = hrc_count + 1;
    end
    if(i ~= 1)
        if(max(temp(1,1:size(temp,2)-1)))
            where = find(temp(1,1:size(temp,2)-1) == 1);
            scene_list{where} = [scene_list{where} clip_structs(i).scene];
        else
            scene_list{length(scene_list)+1,1} = clip_structs(i).scene;
        end
    else
        scene_list{i,1} = clip_structs(i).scene;
    end
end
% clear hrc_list_temp
hrc_list = hrc_list';
hrc_list_old = hrc_list;
num_hrcs = size(hrc_list,1);
if (num_hrcs == 0)
    error('No files with standard naming convention found.\n');
end

%Need to sort clip_structs by HRC
offset_sorted = sort_clips_by('hrc',clip_structs,test_structs);
counter = 1;
for i = 1:length(offset_sorted)
    for j = 1:length(offset_sorted{i})
        temp_struct(counter) = offset_sorted{i}(j);
        counter = counter + 1;
    end
end
clip_structs = clip_structs(temp_struct);
hrc_list_temp = hrc_list_temp(temp_struct);
%Remove Multiple entries.
for j = 1:length(hrc_list_old)
    temp = find(strcmpi(hrc_list_temp,hrc_list_old(j)));
    if(size(temp) == 1)
        %Only one entery
        continue;
    end
    %Keep the first original and remove the rest
    for i = length(temp):-1:2
        hrc_list_temp(temp(i)) = [];
    end
end
hrc_list = hrc_list_temp';
scene_list_temp = scene_list;
for i = 1:length(hrc_list_old)
    moved_to = find(strcmpi(hrc_list_old(i),hrc_list));
    scene_list_temp(moved_to) = scene_list(i);
end
scene_list = scene_list_temp;
clear temp_struct hrc_list_temp temp moved_to scene_list_temp hrc_lost_old;

crashed = 0;
%Check if program crashed and the user is restarting.  The program should
%pick up where the program crashed.
try
    load([results_dir_path '\' results_parameter_name]);
    crashed = 1;
    try
        %Check that the program has not finished properly and only the
        %results file is left.  If this is true, then if the user wants to
        %re-run the psnr over the clipset then the results file's name
        %needs to be changed or deleted from the "results_dir_path".
        results(1).file_name;
    catch
        %This field "file_name" does not exist in the results file.
        %Therefore, it is assumed that psnr search has completely ran and
        %the program will stop running.
        completed = 1;
    end
catch
    %File does not exist, keep crash set equal to zero
    crashed = 0;
end

if(completed == 1)
    error('PSNR Search has fully ran.  Delete the file "%s" in the results_dir_path "%s" to restart PSNR Search or change the "results_parameter_name" input to a new name.',results_parameter_name,results_dir_path);
else
    clear completed;
end

if(crashed == 1)
    %Find where the programmed crashed.
    for i = 1:length(results)
        if(verbose == 1)
            fprintf('PSNR completed for clip %s. Skipping clip!\n',char(results(i).file_name));
        end
    end
    index = i + 1;
    %Load temp Gclips and start from there
    load([results_dir_path '\' 'temp_gclips.mat']);
    load([results_dir_path '\' 'temp_struct_orig.mat']);
    load([results_dir_path '\' 'tshifts.mat']);
    if(calibrated == 0)
        load([results_dir_path '\' 'temp_gclips_new.mat']);
    end
    wb1=waitbar(0,'PSNR Search is running.  Please wait...','color',[1,1,0.65]);
    ph = findobj(wb1,'type','patch');
    set(ph,'facecolor',[0.3 0.3 1],'edgecolor', [0 0 0]);
else
    %Program is running for the first time
    %Results struct to store results, shifts are how much the original must be
    %shifted with respect to the processed
    results = struct('test', {}, 'scene', {}, 'hrc', {}, 'yshift', {}, ...
        'xshift', {}, 'tshift', {}, 'gain', {}, 'offset', {}, 'psnr', {}, ...
        'file_name', {}, 'hrc_i', {}, 'scene_j', {});
    
    % Process one HRC at a time to compute average PSNR for that HRC
    index = 1;  % index to store results
    
    %Change results_file_csv_name to be in the same path as results_dir_path
    results_file_csv_name = strcat(results_dir_path,'\',results_file_csv_name);
    
    fid_results = fopen(results_file_csv_name,'a');  % open results file for appending
    fprintf(fid_results,'Test,Scene,HRC,Yshift,Xshift,Tshift,Gain,Offset,PSNR\n');
    fclose(fid_results);
    
    if(log_file == 1)
        fid_log_file = fopen(strcat(results_dir_path,'\','log_file.txt'), 'a');
        fprintf(fid_log_file,'Test    Scene    HRC     Vshift    Hshift    Delay     Gain       Offset     PSNR\r\n');
        fclose(fid_log_file);
    end
    
    %Set up waitbar
    wb1=waitbar(0,'PSNR Search is running.  Please wait...','color',[1,1,0.65]);
    ph = findobj(wb1,'type','patch');
    set(ph,'facecolor',[0.3 0.3 1],'edgecolor', [0 0 0])

    count_of_clips = 0;
    where_in_clips = 1;
    %Find how many clips there are
    for z = 1:length(clip_structs)
        if(strcmp(clip_structs(z).hrc,'original'))
            continue;
        else
            count_of_clips = count_of_clips + 1;
        end
    end
    
    %Set up new_clip_structs with the same formatting as gclips
    if(calibrated == 1 && uncalibrated == 0)
        %Test Code
        new_clip_structs = clip_structs;
        for i = 1:length(clip_structs)
            clip_structs(i).luminance_gain = 1;
            clip_structs(i).luminance_offset = 0;
        end
    elseif(calibrated == 0 && uncalibrated == 1)
        %Adjustments need to be made in order to let the program run correctly.
        %These adjustments include resetting some values in the clip structure
        %to their default values.
        new_clip_structs = clip_structs;
        for i = 1:length(clip_structs)
            new_clip_structs(i).spatial.horizontal = 0;
            new_clip_structs(i).spatial.vertical = 0;
            new_clip_structs(i).luminance_gain = 1;
            new_clip_structs(i).luminance_offset = 0;
            new_clip_structs(i).scale.horizontal = 1000;
            new_clip_structs(i).scale.vertical = 1000;
            %This assumes that clip_structs is populated correctly with the
            %values of rows and cols.
            image_size.rows = clip_structs(i).image_size.rows;
            image_size.cols = clip_structs(i).image_size.cols;
            roi = default_sroi(image_size);
            new_clip_structs(i).cvr.top = roi.top;
            new_clip_structs(i).cvr.left = roi.left;
            new_clip_structs(i).cvr.bottom = roi.bottom;
            new_clip_structs(i).cvr.right = roi.right;
        end
        %Call two functions to set the default values for the valid region
        %SROI and alignment.
        new_clip_structs = fix_temporal(test_structs,new_clip_structs,'FirstFrame');
        new_clip_structs_rev = new_clip_structs;
    else
        error('Calibrated status is unknown!');
    end
    
    %Create a temp structure to hold the original struct without the changes
    %that will be made below.
    if(calibrated == 1 && uncalibrated == 0)
        temp_struct_orig = clip_structs;
    elseif(calibrated == 0 && uncalibrated == 1)
        temp_struct_orig = new_clip_structs;
    end
end

%For crash recovery, to get i to the correct value.
move_to_next_i = 0;
correct_value_of_i = 0;

for i = 1:num_hrcs
    if(crashed == 1)
        %Reset i
        i = results(index-1).hrc_i;
        if(i > num_hrcs)
            if(verbose == 1)
                display('PSNR Search has been completed!');
            end
            break;
        end
    end
    if(move_to_next_i == 1)
        %Move i to the correct value to correct from program crash.
        if(i ~= correct_value_of_i)
            continue;
        else
            move_to_next_i = 0;
        end
    end
    
    if(crashed == 0)
        psnr_ave = 0;  % initialize the psnr average summer for this HRC
    end
    this_hrc = hrc_list{i};
    if(strcmp('original',this_hrc)) % Don't process original
        continue;
    end
    num_scenes = size(scene_list{i},2);  % Number of scenes in this HRC
    
    for j = 1:num_scenes
        if(crashed == 1)
            %Reset j
            j_new = results(index-1).scene_j + 1;
            %Reset crashed to zero to continue the test
%             crashed = 0;
            if(j_new > num_scenes)
                %Set i to the correct value.  Since you can not adjust i
                %and let the loop run (the i that is replaced will not
                %increment correctly).
                move_to_next_i = 1;
                correct_value_of_i = i + 1;
                j = j_new;
                crashed = 0;
%                 i = i + 1;
                break;
            else
                %j needs to be adjusted in order to continue running the
                %program properly.
                if(j_new ~= j)
                    continue;
                else
                    %Reset crashed to zero to continue the test
                    crashed = 0;
                    move_to_next_i = 1;
                    correct_value_of_i = i + 1;
                end
            end
        end
        
        this_scene = scene_list{i}{j};
        results(index).test = clip_structs(index).test{1};
        results(index).scene = this_scene;
        results(index).hrc = this_hrc;
        
        % Read original and processed video files
        % Re-generate the original and processed avi file names
        test_num = search_test_list(test_structs,clip_structs(index));
        this_original = find_clip(clip_structs,results(index).test,this_scene,'original');
        
        rows = clip_structs(index).image_size.rows;
        cols = clip_structs(index).image_size.cols;

        if(user_set_sroi == 1)
            %This is not gclips sroi, but it is the user entered sroi.
            %Its called gclips_sroi for the program to run easily.
            gclips_sroi = adjust_requested_sroi(clip_structs(index),'sroi',top,left,bottom,right);
        else
            %Use gclips to determine SROI
            if(calibrated == 1)
                gclips_sroi = adjust_requested_sroi(clip_structs(index),'yxextra',y_uncert,x_uncert);
            else
                gclips_sroi = adjust_requested_sroi(new_clip_structs(index),'yxextra',y_uncert,x_uncert);
            end
            if(gclips_sroi.top ~= 1 || gclips_sroi.left ~= 1 ||...
                    gclips_sroi.bottom ~= rows || gclips_sroi.right ~= cols)
                %Then the SROI has been adjusted and is no longer the
                %whole image.
                is_whole_image = 0;
            end
        end
        
        % Set/Validate the SROI
        if (is_whole_image) % make SROI whole image less uncertainty
            top = 1+y_uncert;
            left = 1+x_uncert;
            bottom = rows-y_uncert;
            right = cols-x_uncert;
            if(calibrated == 1)
                gclips_sroi = adjust_requested_sroi(clip_structs(index),'sroi',top,left,bottom,right);
            else
                gclips_sroi = adjust_requested_sroi(new_clip_structs(index),'sroi',top,left,bottom,right);
            end
        end
        top = gclips_sroi.top;
        left = gclips_sroi.left;
        bottom = gclips_sroi.bottom;
        right = gclips_sroi.right;
        
        
        
        %Since top,left,bottom,right has changed, update gclips
        temp_struct_orig(this_original).cvr.top = top-y_uncert;
        temp_struct_orig(this_original).cvr.left = left-x_uncert;
        temp_struct_orig(this_original).cvr.bottom = bottom+y_uncert;
        temp_struct_orig(this_original).cvr.right = right+x_uncert;
        
        if (top<1 || left<1 || bottom>rows || right>cols)
            error('Requested SROI too large for image size.\n');
        end
        
        %Adjust "GClips"  for the new align start and stop that will allow
        %for the full search.  Also a "fake" "GClips" will be created in
        %order to allow for padding of the original video.
        if(calibrated == 1 && uncalibrated == 0)
            %Use Gclips
            
            %Find the align starts that are correct.
            %orig(align_start - start) <= delta
            %orig(align_stop - stop) >= -delta
            if(temp_struct_orig(this_original).align_start - temp_struct_orig(this_original).loc_start <= t_uncert)
                %Move in align_start by the missing difference to allow for
                %full search.
                how_much = t_uncert - (temp_struct_orig(this_original).align_start - temp_struct_orig(this_original).loc_start);
                clip_structs(this_original).align_start = temp_struct_orig(this_original).align_start + how_much;
                clip_structs(index).align_start = clip_structs(index).align_start + how_much;
                clear how_much;
            end
            if(temp_struct_orig(this_original).align_stop - temp_struct_orig(this_original).loc_stop >= -t_uncert)
                %Move in align_start by the missing difference to allow for
                %full search.
                how_much = t_uncert + (temp_struct_orig(this_original).align_stop - temp_struct_orig(this_original).loc_stop);
                clip_structs(this_original).align_stop = temp_struct_orig(this_original).align_stop - how_much;
                clip_structs(index).align_stop = clip_structs(index).align_stop - how_much;
                clear how_much;
            end
        elseif(calibrated == 0 && uncalibrated == 1)
            %Instead of using GClips, use the newly created clipset.
            
            %orig(align_start - start) <= delta
            %orig(align_stop - stop) >= -delta
            if(temp_struct_orig(this_original).align_start - temp_struct_orig(this_original).loc_start <= t_uncert)
                %Move in align_start by the missing difference to allow for
                %full search.
                how_much = t_uncert - (temp_struct_orig(this_original).align_start - temp_struct_orig(this_original).loc_start);
                new_clip_structs(this_original).align_start = temp_struct_orig(this_original).align_start + how_much;
                new_clip_structs(index).align_start = new_clip_structs(index).align_start + how_much;
                clear how_much;
            end
            if(temp_struct_orig(this_original).align_stop - temp_struct_orig(this_original).loc_stop >= -t_uncert)
                %Move in align_start by the missing difference to allow for
                %full search.
                how_much = t_uncert + (temp_struct_orig(this_original).align_stop - temp_struct_orig(this_original).loc_stop);
                new_clip_structs(this_original).align_stop = temp_struct_orig(this_original).align_stop - how_much;
                new_clip_structs(index).align_stop = new_clip_structs(index).align_stop - how_much;
                clear how_much;
            end
        else
            error('Calibrated status is unknown!');
        end
        
        % Validate the spatial uncertainty search bounds
        if (left-x_uncert < 1 || right+x_uncert > cols)
            error('Spatial x-uncertainty too large for SROI.\n');
        end
        if (top-y_uncert < 1 || bottom+y_uncert > rows)
            error('Spatial y-uncertainty too large for SROI.\n');
        end
        
        %update waitbar
        waitbar(where_in_clips/(count_of_clips+1),wb1);
        pause(0.25);
        where_in_clips = where_in_clips + 1;
        
        if(strcmpi(video_standard,'progressive')) %Reads in normal amount of frames
            % Read in video and clear color planes to free up memory
            if(calibrated == 1 && uncalibrated == 0)
                %Use Gclips
                [y_orig] = read_tslice(test_structs(test_num),temp_struct_orig(this_original),...
                    0,1,...
                    'sroi',top-y_uncert,left-x_uncert,bottom+y_uncert,right+x_uncert,...
                    'align_start',clip_structs(this_original).align_start-t_uncert,...
                    'align_stop', clip_structs(this_original).align_stop+t_uncert,...
                    'all_frames','yxextra',y_uncert,x_uncert);
                [y_proc] = read_tslice(test_structs(test_num),clip_structs(index),...
                    0,1,...
                    'sroi',top,left,bottom,right,'all_frames','aligned');
            elseif(calibrated == 0 && uncalibrated == 1)
                %Instead of using GClips, use the newly created
                %clipset.
                [y_orig] = read_tslice(test_structs(test_num),temp_struct_orig(this_original),...
                    0,1,...
                    'sroi',top-y_uncert,left-x_uncert,bottom+y_uncert,right+x_uncert,...
                    'align_start',new_clip_structs(this_original).align_start-t_uncert,...
                    'align_stop', new_clip_structs(this_original).align_stop+t_uncert,...
                    'all_frames','yxextra',y_uncert,x_uncert);
                [y_proc] = read_tslice(test_structs(test_num),new_clip_structs(index),...
                    0,1,...
                    'sroi',top,left,bottom,right,'all_frames','aligned');
            else
                error('Calibrated status is unknown!');
            end
        else %Reads in one less processed frame for interlace
            % Read in video and clear color planes to free up memory
            if(calibrated == 1 && uncalibrated == 0)
                %Use Gclips
                [y_orig] = read_tslice(test_structs(test_num),temp_struct_orig(this_original),...
                    0,1,...
                    'sroi',top-y_uncert,left-x_uncert,bottom+y_uncert,right+x_uncert,...
                    'align_start',clip_structs(this_original).align_start-t_uncert,...
                    'align_stop', clip_structs(this_original).align_stop+t_uncert,...
                    'all_frames','yxextra',y_uncert,x_uncert);
                [y_proc] = read_tslice(test_structs(test_num),clip_structs(index),...
                    0,1,...
                    'sroi',top,left,bottom,right,...
                    'align_stop',clip_structs(index).align_stop-1,...
                    'all_frames');
                if((size(y_orig,3)-(2*t_uncert)-1) - size(y_proc,3) ~= 0)
                    %Reframing is present, need to re-read y_proc with the
                    %right value of align_stop.  +1 because of the -1
                    %above, this will allow for the correct calculation.
                    how_much = abs((size(y_orig,3)-(2*t_uncert)-1) - size(y_proc,3)) + 1;
                    clear y_proc;
                    [y_proc] = read_tslice(test_structs(test_num),clip_structs(index),...
                        0,1,...
                        'sroi',top,left,bottom,right,...
                        'align_stop',clip_structs(index).align_stop-how_much,...
                        'all_frames');
                    clear how_much;
                end
            elseif(calibrated == 0 && uncalibrated == 1)
                %Instead of using GClips, use the newly created
                %clipset.
                [y_orig] = read_tslice(test_structs(test_num),temp_struct_orig(this_original),...
                    0,1,...
                    'sroi',top-y_uncert,left-x_uncert,bottom+y_uncert,right+x_uncert,...
                    'align_start',new_clip_structs(this_original).align_start-t_uncert,...
                    'align_stop', new_clip_structs(this_original).align_stop+t_uncert,...
                    'all_frames','yxextra',y_uncert,x_uncert);
                [y_proc] = read_tslice(test_structs(test_num),new_clip_structs(index),...
                    0,1,...
                    'sroi',top,left,bottom,right,...
                    'align_stop',new_clip_structs(index).align_stop-1,...
                    'all_frames');
            else
                error('Calibrated status is unknown!');
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
            fprintf('\nTest = %s,   Scene = %s,   HRC = %s\n',clip_structs(index).test{1}, this_scene, this_hrc);
        end
        
        if(full_results == 1)
            results_file_full = strcat(clip_structs(index).test{1},'_',this_scene,'_',this_hrc,'.csv');
            if(size(strfind(results_file_csv_name,'/'),2) > 0)
                results_dir = strfind(results_file_csv_name,'/');
                results_dir = results_file_csv_name(1:results_dir(size(results_dir,2)));
                results_file_full = strcat(results_dir,results_file_full);
            elseif(size(strfind(results_file_csv_name,'\'),2) > 0)
                results_dir = strfind(results_file_csv_name,'\');
                results_dir = results_file_csv_name(1:results_dir(size(results_dir,2)));
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
                                    -1*best_yshift,-1*best_xshift,-1*best_tshift,1/best_gain,(-1*best_offset)/best_gain,best_psnr);
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
                                clip_structs(index).test{1}, this_scene, this_hrc, -1*curr_v, -1*curr_h, -1*this_tshift, 1/gain, (-1*offset)/gain, this_psnr);
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
                                    -1*best_yshift,-1*best_xshift,-1*best_tshift,1/best_gain,(-1*best_offset)/best_gain,best_psnr);
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
                                clip_structs(index).test{1}, this_scene, this_hrc, -1*curr_v, -1*curr_h, -1*this_tshift, 1/gain, (-1*offset)/gain, this_psnr);
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
                                    -1*best_yshift,-1*best_xshift,-1*best_tshift,1/best_gain,(-1*best_offset)/best_gain,best_psnr);
                            end
                        end
                        x(cnt_graph) = -1*curr_h;
                        y(cnt_graph) = -1*curr_v;
                        t(cnt_graph) = -1*curr_time;
                        psnr_graph(cnt_graph) = this_psnr;
                        
                        if(full_results == 1)
                            gain = this_fit(1);
                            offset = this_fit(2);
                            fid_results_full = fopen(results_file_full, 'a');
                            fprintf(fid_results_full,'%s, %s, %s, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f\n', ...
                                clip_structs(index).test{1}, this_scene, this_hrc, -1*curr_v, -1*curr_h, -1*curr_time, 1/gain, (-1*offset)/gain, this_psnr);
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
        
        results(index).yshift = -1*best_yshift;
        results(index).xshift = -1*best_xshift;
        results(index).tshift = -1*best_tshift;
        results(index).gain = (1/best_gain);
        results(index).offset = (-1*best_offset)/best_gain;
        results(index).psnr = best_psnr;
        results(index).file_name = clip_structs(index).file_name;
        results(index).hrc_i = i;
        results(index).scene_j = j;
        fid_results = fopen(results_file_csv_name,'a');  % open results file for appending
        fprintf(fid_results,'%s, %s, %s, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f, %5.4f\n', ...
            results(index).test, results(index).scene, results(index).hrc, ...
            results(index).yshift, results(index).xshift, results(index).tshift, ...
            results(index).gain, results(index).offset, results(index).psnr);
        fclose(fid_results);
        
        if(log_file == 1)
            fid_log_file = fopen(strcat(results_dir_path,'\','log_file.txt'), 'a');
            fprintf(fid_results,'%s    %s     %s    %5.4f    %5.4f    %5.4f    %5.4f    %5.4f    %5.4f\r\n', ...
                results(index).test, results(index).scene, results(index).hrc, ...
                results(index).yshift, results(index).xshift, results(index).tshift, ...
                results(index).gain, results(index).offset, results(index).psnr);
            fclose(fid_log_file);
        end
        
        psnr_ave = psnr_ave+best_psnr;
        
        if(calibrated == 1 && uncalibrated == 0)
            %Get the output "GClips" ready.
            new_clip_structs(index) = clip_structs(index);
        else
            new_clip_structs_rev(index) = new_clip_structs(index);
        end
        
        %Update GClips
        if(calibrated == 1)
            %In order for Read_tslice to read correctly, changes were made to
            %the clip structure so it must be reset to the  original values.
            new_clip_structs(index).align_start = temp_struct_orig(index).align_start;
            new_clip_structs(index).align_stop = temp_struct_orig(index).align_stop;
            new_clip_structs(this_original).align_start = temp_struct_orig(this_original).align_start;
            new_clip_structs(this_original).align_stop = temp_struct_orig(this_original).align_stop;
            if(best_yshift ~= 0)
                %If the shift is not equal to zero, then gclips needs to be
                %updated
                best_yshift = -1*best_yshift; %Recalc for GClips notation.
                new_clip_structs(index).spatial.vertical = new_clip_structs(index).spatial.vertical + best_yshift;
            end
            if(best_xshift ~= 0)
                %If the shift is not equal to zero, then gclips needs to be
                %updated
                best_xshift = -1*best_xshift; %Recalc for GClips notation.
                new_clip_structs(index).spatial.horizontal = new_clip_structs(index).spatial.horizontal + best_xshift;
            end
            %Change into gclips format
            new_clip_structs(index).luminance_gain = (1/best_gain);
            new_clip_structs(index).luminance_offset = (-1*best_offset)/best_gain;
            best_tshift = -1*best_tshift; %Recalc for GClips notation.
            %After recalc for Gclips make sure that "reframing" is proper
            %to store in Gclips by use of the floor function.
            tshifts(index) = floor(best_tshift);
            which_orig(index) = this_original;
        else
            %In order for Read_tslice to read correctly, changes were made to
            %the clip structure so it must be reset to the  original values.
            new_clip_structs_rev(index).align_start = temp_struct_orig(index).align_start;
            new_clip_structs_rev(index).align_stop = temp_struct_orig(index).align_stop;
            new_clip_structs_rev(this_original).align_start = temp_struct_orig(this_original).align_start;
            new_clip_structs_rev(this_original).align_stop = temp_struct_orig(this_original).align_stop;
            if(best_yshift ~= 0)
                %If the shift is not equal to zero, then gclips needs to be
                %updated
                best_yshift = -1*best_yshift; %Recalc for GClips notation.
                new_clip_structs_rev(index).spatial.vertical = new_clip_structs_rev(index).spatial.vertical + best_yshift;
            end
            if(best_xshift ~= 0)
                %If the shift is not equal to zero, then gclips needs to be
                %updated
                best_xshift = -1*best_xshift; %Recalc for GClips notation.
                new_clip_structs_rev(index).spatial.horizontal = new_clip_structs_rev(index).spatial.horizontal + best_xshift;
            end
            %Change into gclips format
            new_clip_structs_rev(index).luminance_gain = (1/best_gain);
            new_clip_structs_rev(index).luminance_offset = (-1*best_offset)/best_gain;
            best_tshift = -1*best_tshift; %Recalc for GClips notation.
            %After recalc for Gclips make sure that "reframing" is proper
            %to store in Gclips by use of the floor function.
            tshifts(index) = floor(best_tshift);
            which_orig(index) = this_original;
        end
        
        
        %Generate a save file for the psnr parameter values
        save(strcat(results_dir_path,'\',results_parameter_name),'results');
        %Save temp Gclips struct just incase of crash
        if(calibrated == 1)
            save(strcat(results_dir_path,'\','temp_gclips.mat'),'new_clip_structs');
        else
            save(strcat(results_dir_path,'\','temp_gclips.mat'),'new_clip_structs');
            save(strcat(results_dir_path,'\','temp_gclips_new.mat'),'new_clip_structs_rev');
        end
        save(strcat(results_dir_path,'\','temp_struct_orig.mat'),'temp_struct_orig');
        save(strcat(results_dir_path,'\','tshifts.mat'),'which_orig','tshifts','psnr_ave','count_of_clips','where_in_clips','results_file_csv_name');
          
        index = index+1;
        
    end
    
    if(j <= num_scenes)    
        % Compute average PSNR for this HRC
            psnr_ave = psnr_ave/(num_scenes);
            if(verbose)
                fprintf('\nHRC = %s, psnr_ave = %5.4f\n\n',this_hrc, psnr_ave);
                %Close figure 1
                close(1);
            end
%             close all; %changed
    end
    if(i == num_hrcs)
        %Have reached the end of the search.  If "break" is not placed in
        %the code the program re calculates an un-needed result.
        if(verbose == 1)
            display('PSNR Search has already been run over all the clips in the clipset.');
        end
        break;
    end
end

%Figure out align starts and stops after all the calculations have been
%completed.
%Find how many originals exist:
how_many_orig = unique(which_orig);
for i = 1:length(how_many_orig)
    %Find which tshifts correspond to each orig
    where_at = find(which_orig == how_many_orig(i));
    temp_tshifts = tshifts(where_at);
    max_pos_tshift = max(temp_tshifts);
    where_max = where_at(find(max_pos_tshift == tshifts(where_at)));
    
    min_pos_tshift = min(temp_tshifts);
    where_min = where_at(find(min_pos_tshift == tshifts(where_at)));
    
    %Min tshift
    if(~isempty(where_min) && min_pos_tshift < 0)
        if(uncalibrated == 1)
            if(length(where_min) > 1)
                for i = 1:length(where_min)
                    new_clip_structs_rev(where_min(i)).align_start = 1;
                    tshifts(where_min(i)) = 0;
                end
                %Change align start of the original
                new_clip_structs_rev(which_orig(where_min(1))).align_start = 1 + abs(min_pos_tshift);
            else
                new_clip_structs_rev(where_min).align_start = 1;
                %Change align start of the original
                new_clip_structs_rev(which_orig(where_min)).align_start = 1 + abs(min_pos_tshift);
                tshifts(where_min) = 0;
            end
            %May need to check other calculated alignments to
            %re-adjust if fix_temporal does not fix it for you.
            new_clip_structs_rev = fix_temporal(test_structs,new_clip_structs_rev,'endpoint');
        else
            if(length(where_min) > 1)
                for i = 1:length(where_min)
                    new_clip_structs(where_min(i)).align_start = new_clip_structs(where_min(i)).align_start + min_pos_tshift;
                    if(new_clip_structs(where_min(i)).align_start < 1 && i == 1)
                        %Cant be below one!  Need to shift original up.
                        %The original only needs to be shifted once.
                        off_by = abs(new_clip_structs(where_min(i)).align_start) + 1;
                        new_clip_structs(where_min).align_start = 1;
                        new_clip_structs(which_orig(where_min)).align_start = new_clip_structs(which_orig(where_min)).align_start + off_by;
                    elseif(new_clip_structs(where_min(i)).align_start < 1 && i > 1)
                        %Since the original has already been shifted it
                        %does not need to be shifted again.  Only the proc
                        %needs to be changed from negative to pos 1.
                        new_clip_structs(where_min).align_start = 1;
                    end
                    tshifts(where_min(i)) = 0;
                end
            else
                new_clip_structs(where_min).align_start = new_clip_structs(where_min).align_start + min_pos_tshift;
                if(new_clip_structs(where_min).align_start < 1)
                    %Cant be below one!  Need to shift original up.
                    off_by = abs(new_clip_structs(where_min).align_start) + 1;
                    new_clip_structs(where_min).align_start = 1;
                    new_clip_structs(which_orig(where_min)).align_start = new_clip_structs(which_orig(where_min)).align_start + off_by;
                end
                tshifts(where_min) = 0;
            end
            %May need to check other calculated alignments to
            %re-adjust if fix_temporal does not fix it for you.
            new_clip_structs = fix_temporal(test_structs,new_clip_structs,'endpoint');
        end
    end
    %Max tshift
    if(~isempty(where_max) && max_pos_tshift > 0)
        if(uncalibrated == 1)
            if(length(where_max) > 1)
                for i = 1:length(where_max)
                    new_clip_structs_rev(where_max(i)).align_start = new_clip_structs_rev(which_orig(where_max(i))).align_start + max_pos_tshift;
                    tshifts(where_max(i)) = 0;
                end
                new_clip_structs_rev = fix_temporal(test_structs,new_clip_structs_rev,'endpoint');
            else
                new_clip_structs_rev(where_max).align_start = new_clip_structs_rev(which_orig(where_max)).align_start + max_pos_tshift;
                new_clip_structs_rev = fix_temporal(test_structs,new_clip_structs_rev,'endpoint');
                tshifts(where_max) = 0;
            end
        else
            if(length(where_max) > 1)
                for i = 1:length(where_max)
                    new_clip_structs(where_max(i)).align_start = new_clip_structs(where_max(i)).align_start + max_pos_tshift;
                    tshifts(where_max(i)) = 0;
                end
                new_clip_structs = fix_temporal(test_structs,new_clip_structs,'endpoint');
            else
                new_clip_structs(where_max).align_start = new_clip_structs(where_max).align_start + max_pos_tshift;
                new_clip_structs = fix_temporal(test_structs,new_clip_structs,'endpoint');
                tshifts(where_max) = 0;
            end
        end
    end
end

for i = 1:length(tshifts)
    if(tshifts(i) == 0)
        continue;
    else
        if(uncalibrated == 1)
            orig_value = new_clip_structs_rev(which_orig(i)).align_start;
            new_clip_structs_rev(i).align_start = orig_value + tshifts(i);
        else
            orig_value = new_clip_structs(i).align_start;
            new_clip_structs(i).align_start = orig_value + tshifts(i);
        end      
    end
end

if(calibrated == 0)
    new_clip_structs = new_clip_structs_rev;
end

%Interlaced video needs to be taken care of by the following call to
%fix_temporal.
new_clip_structs = fix_temporal(test_structs,new_clip_structs,'endpoint');
new_clip_structs = fix_temporal(test_structs,new_clip_structs,'spatial');



%Clear temp Gclips struct from file
delete([results_dir_path '\' 'temp_gclips.mat']);
delete([results_dir_path '\' 'temp_struct_orig.mat']);
delete([results_dir_path '\' 'tshifts.mat']);
if(calibrated == 0)
    delete([results_dir_path '\' 'temp_gclips_new.mat']);
end
results = rmfield(results,{'file_name','hrc_i','scene_j'});
save(strcat(results_dir_path,'\',results_parameter_name),'results');
close(wb1);


function current_best_psnr = plot_psnr(x,y,psnr,t, current_best_psnr, verbose) 
if(verbose)
    max_psnr = max(psnr);
    if(max_psnr < current_best_psnr)
        %Leave function
    else
        [X,Y] = meshgrid(linspace(min(x),max(x)), linspace(min(y),max(y)));
        psnr_test = griddata(x,y,psnr,X,Y,'cubic');
        which_fig = mesh(X,Y,psnr_test);
        hold on;
        xlabel('Spatial X'), ylabel('Spatial Y'), zlabel('PSNR'), title(sprintf('Time Shift (%d)', t));
        contour3(X, Y, psnr_test);
        plot3(x,y,psnr,'.','markersize',10);
        hold off;
        pause(1);
        current_best_psnr = max_psnr;
    end
end



