function status = vfd_feature_loop_mse(test_structs, clip_structs, deg_size, ...
    time_size, feature_base_dir, varargin)
% VFD_FEATURE_LOOP_MSE
%  Loop through a list of clips, compute the variable frame delay (VFD)
%  Mean Squared Error (MSE) between the processed and original, and save
%  the processed and original clip's features in a file.  This is a
%  full reference metric based on the difference between the processed and
%  the original pixels.  To maintain standard naming conventions, the MSE 
%  block features are saved as the Root Mean Square (RMS) of (proc-org).
%  Only the Y image is used for this computation so the features can be
%  used to compute the traditional PSNR based on the Y channel only.
%  
%  VFD feature files for processed clips contain both its features
%  (variable 'data') plus the corresponding VFD-compensated original
%  features (variable 'datao').  VFD feature files for original clips 
%  contain both variables 'data' and 'datao' but these two variables
%  contain identical information.  This is to maintain compatibility with
%  existing parameter and model functions that utilize variable 'data'.
%  Since this is a full reference metric, the 'datao' variables will
%  contain zeros (i.e., orig-orig).
% SYNTAX
%  vfd_feature_loop_mse(test_structs, clip_structs, deg_size, time_size, feature_base_dir)
%  vfd_feature_loop_mse(...,'PropertyName',...);
% DESCRIPTION
%  Compute VFD MSE features
% 
%  This function takes variable 'test_structs' (of the same format as
%  GTests), which describes each video test and the location of the
%  associated video; and variable 'clip_structs' (of the same format as
%  GClips), which lists the clips to be processed (including the
%  originals), 'deg_size' specifies the block size in angular degrees
%  (horizontally & vertically),  'time_size' specifies the length of the
%  block in SECONDS, and 'feature_base_dir' is the path to the directory
%  where the feature's sub-directory should be created.  That
%  sub-directory's name will be this feature's name.
%
%  To calculate the corresponding block sizes (in pixels) for a specified
%  'deg_size', test_structs.viewing_distance (in picture heights) is used.
%  For simplicity, all the image rows are used for the picture height H.
%  For interlaced systems, the block size (in pixels/lines) for a desired
%  'deg_size' is forced to be even so an equal number of lines from each
%  field are included in the block (blocks are calculated from frames, not
%  fields, with this routine).  The nearest even block size is chosen for
%  interlaced systems.
%
%  The following optional arguments may be specified: 
%   'blockmean'     Compute an additional feature using the mean of the sum
%                   of the errors of (proc-orig).  This is like the DC
%                   error of each block, as opposed to the RMS error 
%                   that includes both the DC and AC components.
%   'verbose'   Print progress report to screen.
%   'quiet'     Don't print anything to the screen (default).
%
%  Return 0 if operated correctly; 1 if an error was encountered.

status = 0;  % Normal return status
optional_mean = 0;
verbose = 0;
last_test = '';  % State variable that gives the last test whose vfd information was loaded
aba_t = 8;  % ave_best_aligned threshold for still scene detection

global Y_GLOBAL;  % Large 3D image array to hold video for feature extraction

for cnt = 1:length(varargin),
    if strcmpi(varargin{cnt},'blockmean'),
        optional_mean = 1;
    elseif strcmpi(varargin{cnt},'verbose'),
        verbose = 1;
    elseif strcmpi(varargin{cnt},'quiet'),
        verbose = 0;
    else
        error('optional property not recognized');
    end
end

% Generate the standard names for the features
feature_name = standard_feature_names('Y','degsize',deg_size,'tslice_length_sec',time_size,...
        'block_stat','rms', 'specific','FR_vfd');
feature_name2 = standard_feature_names('Y','degsize',deg_size,'tslice_length_sec',time_size,...
        'block_stat','mean', 'specific','FR_vfd');

name_exists2 = 1;
for loop = 1:max(size(clip_structs)),
    try
        t = clock;
        is_original = 0;  % Boolean that says if this clip is an original clip
        
        %  Pick off the clip's test, scene, and HRC names
        clip_test =  clip_structs(loop).test{1};
        clip_scene = clip_structs(loop).scene{1};
        clip_hrc = clip_structs(loop).hrc{1};
        
        % Find the offset of the test structure for this clip, which gives
        % the path to the files for reading.
        tnum = search_test_list(test_structs, clip_structs(loop));
        
        % Determine if this is an original clip or not
        if (strcmpi(clip_hrc,'original'))
            is_original = 1;
        end
        
        %  Load the appropriate VFD delay information for this test.
        %  This code utilizes the variables 'this_clipset' and
        %  'results' from the VFD mat files.  The align_start and
        %  align_stop of the clip will be checked against the VFD info to
        %  make sure the VFD alignment is synched to the clip_structs alignment.
        if(~strcmpi(last_test,clip_test))  % Only load the VFD results file if it is different from last time
            last_test = clip_test;
            load([feature_base_dir 'vfd_' clip_test '.mat']);
        end
        vfd_ind = find_clip(this_clipset, clip_test, clip_scene, clip_hrc);
        if (clip_structs(loop).align_start ~= this_clipset(vfd_ind).align_start ...
                || clip_structs(loop).align_stop ~= this_clipset(vfd_ind).align_stop)
            fprintf('\tError encountered for clip %s:%s(%s)\n', clip_test, clip_scene, clip_hrc);
            error('Invalid VFD information for this clip: align_start or align_stop')
        end
        
        %  See if the requested feature files already exists and skip 
        %  computation of them if they do.
        [name_exists, file_name] = write_feature_standard(0,feature_base_dir,feature_name,clip_structs(loop),'exists');
        if (optional_mean)
            [name_exists2, file_name2] = write_feature_standard(0,feature_base_dir,feature_name2,clip_structs(loop),'exists');
        end
        if (name_exists && name_exists2)
            if (verbose)
                fprintf('\tfeature already computed for clip %s:%s(%s)\n', clip_test, clip_scene, clip_hrc);
            end
            continue;
        end
        
        %  Output the clip information and time before continuing
        if (verbose)
            fprintf('Clip %d of %d ==> %s:%s(%s) at %d:%d\n', loop, max(size(clip_structs)), clip_test, clip_scene, clip_hrc, t(4), t(5) );
        end
        
        %  Calculate the vsize and hsize for the requested deg_size using
        %  the viewing distance and the image_rows.
        image_rows = clip_structs(loop).image_size.rows;
        image_standard = clip_structs(loop).video_standard;
        clip_vd = test_structs(tnum).viewing_distance;
        vsize = clip_vd*deg_size*image_rows*pi/180;
        if (strcmpi(image_standard,'progressive'))
            vsize = round(vsize);
        else  % Interlaced, set vsize to nearest even number
            vsize = round(vsize/2)*2;
        end
        if (vsize <=1 || vsize >= image_rows)
            error('Invalid deg_size');
        end
        hsize = vsize;
        
        % Find adaptive filter 'extra' size so the feature blocks match the SI & HV feature blocks
        [filter_size, extra] = adaptive_filter (clip_structs(loop).image_size);
        
        % Compute the default adjusted SROI and number of blocks available.
        if (strcmpi(clip_structs(loop).video_standard,'progressive'))
            [sroi,vert,horiz] = adjust_requested_sroi (clip_structs(loop), 'vsize',vsize, 'hsize',hsize, 'extra',extra);
        else
            [sroi,vert,horiz] = adjust_requested_sroi (clip_structs(loop), 'vsize',vsize, 'hsize',hsize, 'extra',extra, 'evenodd');
        end
        nrows = sroi.bottom-sroi.top+1;
        ncols = sroi.right-sroi.left+1;
        
        % Determine the total number of aligned frames that are available
        % (in seconds) and read in the entire clip.
        tslice_length_sec_proc = (clip_structs(loop).align_stop-clip_structs(loop).align_start+1) / clip_structs(loop).fps;
        if (is_reframing_indicated(clip_structs(loop)))
            tslice_length_sec_proc = tslice_length_sec_proc - 1/clip_structs(loop).fps;
        end
        number_tslices_proc = total_tslices(clip_structs(loop),tslice_length_sec_proc);
        if (number_tslices_proc ~= 1)
            error('Error in computing the number of processed aligned frames and their time duration');
        end
        nframes = tslice_conversion(tslice_length_sec_proc, clip_structs(loop).fps);  % Total number of frames in processed clip
        
        % Calculate the temporal size in frames so feature extraction
        % subroutine does not need to be passed the clip fps.
        tsize_frames = time_size*clip_structs(loop).fps;  % could be a fractional number of frames
        
        % Perform VFD correction on the original (unless this clip is an
        % original).
        if (is_original)  % No need for feature extraction so just fill zeros
            
            data = zeros(vert, horiz, floor(nframes/tsize_frames));  % Same formula as feature_mse
            if (optional_mean)
                data2 = data;
            end
            
        else  % make the original look like the processed through VFD correction
            
            %  Check the original's alignment with the VFD alignment
            vfd_ind_orig = find_clip(this_clipset, clip_test, clip_scene, 'original');  % original location in VFD variable 'this_clipset'
            ind_orig = find_clip(clip_structs, clip_test, clip_scene, 'original');  % original location in clip_structs
            if (clip_structs(ind_orig).align_start ~= this_clipset(vfd_ind_orig).align_start ...
                    || clip_structs(ind_orig).align_stop ~= this_clipset(vfd_ind_orig).align_stop)
                fprintf('\tError encountered in VFD correction to original for clip %s:%s(%s)\n', clip_test, clip_scene, clip_hrc);
                error('Invalid VFD information for this original: align_start or align_stop')
            end
            
            % Figure out number of original frames available, unaligned.
            % Will use the same algorithm that was used by the program
            % that calculated VFD since the VFD results are referenced to
            % these original frames numbers.
            tslice_length_sec_orig = (this_clipset(vfd_ind_orig).loc_stop - ...
                this_clipset(vfd_ind_orig).loc_start + 1) / this_clipset(vfd_ind_orig).fps;
            number_tslices_orig = total_tslices(this_clipset(vfd_ind_orig),tslice_length_sec_orig, 'unaligned');
            if number_tslices_orig ~= 1,
                error('Error in computing the number of original unaligned frames and their time duration');
            end
            
            % Set the initial alignment point (in frames) for the VFD
            % correction.  first_align is considered frame or field
            % number 1 in the VFD results file.
            first_align = this_clipset(vfd_ind_orig).align_start - this_clipset(vfd_ind_orig).loc_start + 1;
            
            % Set the interlaced flag for the VFD correction.
            if (strcmpi(this_clipset(vfd_ind_orig).video_standard, 'interlace_lower_field_first') == 1)
                is_interlaced = 1;
                field_first = 1;
                first_align = 2*first_align-1;  % convert to fields
            elseif (strcmpi(this_clipset(vfd_ind_orig).video_standard, 'interlace_upper_field_first') == 1)
                is_interlaced = 1;
                field_first = 2;
                first_align = 2*first_align-1;  % convert to fields
            else
                is_interlaced = 0;  % progressive
            end
            
            % read in entire unaligned original clip:
            y_orig = read_tslice(test_structs(tnum),clip_structs(ind_orig),tslice_length_sec_orig, 1, ...
                'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right, 'unaligned');
            
            % Find the VFD original alignment results for this processed
            % clip.  You must eliminate the originals from this_clipset to
            % find the proper results index.
            this_clipset_proc = this_clipset(find_clip(this_clipset,'*','*','original','not'));
            results_ind = find_clip(this_clipset_proc, clip_test, clip_scene, clip_hrc);
            this_result = results{results_ind};  % Pick off the alignment results for this clip from the retrieved VFD info
            
            %  Determine if the clip fails to meet the ave_best_aligned
            %  threshold for still scene detection.  If so, use first_align
            %  to generate this_result (i.e., use gclips alignment), so
            %  there will be no difference between VFD results and normal
            %  results.
            if (this_result ~= 0)  % Must check for zero results or ave_best_aligned will throw error
                aba = ave_best_aligned(results_fuzzy_mse{results_ind});  % results_fuzzy_mse is loaded from mat file
                if (aba > aba_t)
                    if (is_interlaced)
                        npts = nframes*2;
                    else
                        npts = nframes;
                    end
                    this_result = first_align:first_align+npts-1;  % Overwrite the VFD results with gclips alignment
                    if (verbose)
                        fprintf('%s_%s_%s failed ave_best_aligned threshold test, aba = %f, using gclips alignment\n', clip_test, clip_scene, clip_hrc, aba);
                    end
                end
            end
            
            %  If the VFD alignment algorithm failed, then this_result==0.
            %  In that case, use first_align to generate this_result (i.e.,
            %  use gclips alignment), so there will be no difference
            %  between VFD results and normal results.
            if (this_result == 0)
                if (is_interlaced)
                    npts = nframes*2;
                else
                    npts = nframes;
                end
                this_result = first_align:first_align+npts-1;
                if (verbose)
                    fprintf('%s_%s_%s VFD algorithm failed, using gclips alignment.\n', clip_test, clip_scene, clip_hrc);
                end
            end
            
            % new_orig will be the VFD-corrected original.  Use same rule
            % as read_tslice for image precision
            if (nrows > 650)
                new_orig = zeros(nrows,ncols,nframes,'single');
            else
                new_orig = zeros(nrows,ncols,nframes,'double');
            end
            
            if(is_interlaced)
                
                for j = 1:nframes
                    % Get matching original field for the early processed field
                    orig_frame_num = ceil(this_result(2*j-1)/2);  % The frame number that contains the original field
                    early_field = mod(this_result(2*j-1),2);  % =1 if early field, =0 if late field
                    [yo1 yo2] = split_into_fields(squeeze(y_orig(:,:,orig_frame_num)));
                    if (early_field)
                        switch field_first
                            case(1)
                                this_orig1 = yo1;
                            case(2)
                                this_orig1 = yo2;
                        end
                    else % late field
                        switch field_first
                            case(1)
                                this_orig1 = yo2;
                            case(2)
                                this_orig1 = yo1;
                        end
                    end
                    % Get matching original field for the late processed field
                    orig_frame_num = ceil(this_result(2*j)/2);  % The frame number that contains the original field
                    early_field = mod(this_result(2*j),2);  % =1 if early field, =0 if late field
                    [yo1 yo2] = split_into_fields(squeeze(y_orig(:,:,orig_frame_num)));
                    if (early_field)
                        switch field_first
                            case(1)
                                this_orig2 = yo1;
                            case(2)
                                this_orig2 = yo2;
                        end
                    else % late field
                        switch field_first
                            case(1)
                                this_orig2 = yo2;
                            case(2)
                                this_orig2 = yo1;
                        end
                    end
                    % Joint the two original fields into a frame
                    switch field_first
                        case(1)
                            this_orig = join_into_frames(this_orig1,this_orig2);
                        case(2)
                            this_orig = join_into_frames(this_orig2,this_orig1);
                    end
                    new_orig(:,:,j) = this_orig;
                end
                clear this_orig this_orig1 this_orig2 yo1 yo2;
                
            else  % progressive
                
                for j = 1:nframes
                    new_orig(:,:,j) = y_orig(:,:,this_result(j));
                end
                
            end
            
            clear y_orig;  % The VFD corrected original is now in new_orig
            
            %  Read in the processed clip
            y_proc = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec_proc, 1, ...
                'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
            [nrows,ncols,nframes] = size(y_proc);  % must be a VFD result for each of the fields/frames in nframes
            
            % Compute the difference between the processed and the new VFD corrected original
            Y_GLOBAL = y_proc - new_orig;
            clear y_proc new_orig;
            
            % Compute feature data (rms), data2 (mean, if requested)
            % from global variable Y_GLOBAL.
            if (optional_mean)
                [data data2] = feature_mse(vsize, hsize, tsize_frames);
            else
                [data] = feature_mse(vsize, hsize, tsize_frames);
            end
            
        end
        
        %  For VFD feature file compatibility, create zero arrays for datao
        %  and data2o (if requested).
        [size1,size2,size3] = size(data);
        datao = zeros(size1,size2,size3);
        if (optional_mean)
            data2o = datao;
        end
        
        % Write features all at once to their respective files
        write_feature_standard(data, feature_base_dir, feature_name, clip_structs(loop), 'vfd', datao);
        if (optional_mean)
            write_feature_standard(data2, feature_base_dir, feature_name2, clip_structs(loop), 'vfd', data2o);
        end
        
        % Erase the data variables from memory
        clear data*;
        
    catch exception
        status = 1;
        if (verbose)
            disp(exception.message);
            fprintf('\tSkipping clip.\n');
        end
    end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data data2] = feature_mse(vsize, hsize, tsize_frames)
% FEATURE_MSE
%  Compute MSE block features (rms and mean) given a proc-orig video
%  clip, and the desired block size of vsize by hsize by tsize_frames.  The
%  vsize and hsize must evenly divide the vertical and horizontal size of
%  Y_GLOBAL (a global variable that holds the entire video clip).  Note
%  that tsize_frames may be fractional. 
%
% SYNTAX
%  [data data2] = feature_mse(vsize, hsize, tsize_frames)
%
% DESCRIPTION
%  The function takes a 3D global variable (Y_GLOBAL) that holds the video
%  clip and extracts rms (data) and mean (data2) block features from blocks
%  of size vsize by hsize by tsize_frames.  If the size(Y_GLOBAL) = [nrows,
%  ncols, nframes], then the size of the returned feature arrays data and
%  data2 will be nrows/vsize, ncols/hsize, and floor(nframes/tsize_frames)
%  (i.e., the subroutines returns the maximum number of time slices).  For
%  fractional tsize_frames, some spatial-temporal (ST) blocks will have
%  more pixels than others but a given pixel will only be included in one
%  and only one block. 
%

global Y_GLOBAL;

%  Check the spatial size of Y_GLOBAL for validity
[nrows, ncols, nframes] = size(Y_GLOBAL);
if (rem(nrows,vsize) || rem(ncols,hsize))
    error('Invalid vsize or hsize detected in function feature_mse');
end

% Compute Spatial block sizes
nvert = nrows/vsize;  % Number of blocks in vertical dim
nhoriz = ncols/hsize;  % Number of blocks in horizontal dim

% Compute the number of time slices
nslices = floor(nframes/tsize_frames);  % Maximum number of time slices

%  Dimension output arrays as needed
if (nargout ==2)  % Mean requested
    optional_mean = 1;
    data = zeros(nvert,nhoriz,nslices);
    data2 = zeros(nvert,nhoriz,nslices);
else  % Mean not requested
    optional_mean = 0;
    data = zeros(nvert,nhoriz,nslices);
end

% For each time-slice in this clip, compute the requested block statistic.
for cnt = 1:nslices
    
    beg_frame = 1 + floor((cnt-1)*tsize_frames);
    end_frame = floor(cnt*tsize_frames);
    
    if (~optional_mean)
        data(:,:,cnt) = block_statistic (Y_GLOBAL(:,:,beg_frame:end_frame), vsize, hsize, 'rms');
    else
        [data(:,:,cnt),data2(:,:,cnt)] = block_statistic (Y_GLOBAL(:,:,beg_frame:end_frame), vsize, hsize, 'rms','mean');
    end
    
end

% Has the effect of clearing Y_GLOBAL memory space if the space is 
% not occupied by another variable.
Y_GLOBAL = 0;  
        
end
