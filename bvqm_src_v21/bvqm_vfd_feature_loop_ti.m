function status = vfd_feature_loop_ti(test_structs, clip_structs, deg_size, ...
    time_size, feature_base_dir, varargin)
% VFD_FEATURE_LOOP_TI
%  Loop through a list of clips, compute variable frame delay (VFD) RMS
%  Temporal Information (TI) for Y, Cb, and Cr, and save the processed 
%  and original clip's features in a file.  VFD feature files for
%  processed clips contain both its features (variable 'data') plus the
%  corresponding VFD-compensated original features (variable 'datao').  VFD
%  feature files for original clips contain both variables 'data' and
%  'datao' but these two variables contain identical information.  This is
%  to maintain compatibility with existing parameter and model functions
%  that utilize variable 'data'. 
% SYNTAX
%  vfd_feature_loop_ti(test_structs, clip_structs, deg_size, time_size, feature_base_dir)
%  vfd_feature_loop_ti(...,'PropertyName',...);
% DESCRIPTION
%  Compute VFD Temporal Information (TI) features for Y, Cb, and Cr.
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
%   'blockmean'     Compute feature using block mean, as well as the usual
%                   RMS of each block.
%   'verbose'   Print progress report to screen.
%   'quiet'     Don't print anything to the screen (default).
%
%  Return 0 if operated correctly; 1 if an error was encountered.

status = 0;  % Normal return status
optional_mean = 0;
verbose = 0;
last_test = '';  % State variable that gives the last test whose vfd information was loaded
aba_t = 8;  % ave_best_aligned threshold for still scene detection

global Y_GLOBAL CB_GLOBAL CR_GLOBAL;  % Large 3D image arrays to hold video for feature extraction

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
feature_yrms = standard_feature_names('Y','degsize',deg_size,'tslice_length_sec',time_size,...
        'block_stat','rms', 'specific','vfd_ti');
feature_ymean = standard_feature_names('Y','degsize',deg_size,'tslice_length_sec',time_size,...
        'block_stat','mean', 'specific','vfd_ti');
feature_cbrms = standard_feature_names('Cb','degsize',deg_size,'tslice_length_sec',time_size,...
        'block_stat','rms', 'specific','vfd_ti');
feature_cbmean = standard_feature_names('Cb','degsize',deg_size,'tslice_length_sec',time_size,...
        'block_stat','mean', 'specific','vfd_ti');
feature_crrms = standard_feature_names('Cr','degsize',deg_size,'tslice_length_sec',time_size,...
        'block_stat','rms', 'specific','vfd_ti');
feature_crmean = standard_feature_names('Cr','degsize',deg_size,'tslice_length_sec',time_size,...
        'block_stat','mean', 'specific','vfd_ti');

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
        [yrms_exists, file_name] = write_feature_standard(0,feature_base_dir,feature_yrms,clip_structs(loop),'exists');
        [cbrms_exists, file_name] = write_feature_standard(0,feature_base_dir,feature_cbrms,clip_structs(loop),'exists');
        [crrms_exists, file_name] = write_feature_standard(0,feature_base_dir,feature_crrms,clip_structs(loop),'exists');
        if (optional_mean)
            [ymean_exists, file_name] = write_feature_standard(0,feature_base_dir,feature_ymean,clip_structs(loop),'exists');
            [cbmean_exists, file_name] = write_feature_standard(0,feature_base_dir,feature_cbmean,clip_structs(loop),'exists');
            [crmean_exists, file_name] = write_feature_standard(0,feature_base_dir,feature_crmean,clip_structs(loop),'exists');
        else
            ymean_exists = 1;  % Setting to the exists state so this feature will not be computed
            cbmean_exists = 1;
            crmean_exists = 1;
        end
        if (yrms_exists && ymean_exists && cbrms_exists && cbmean_exists && crrms_exists && crmean_exists)
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
        
        % Determine the total number of aligned frames that are available
        % (in seconds) and read in the entire clip.
        tslice_length_sec = (clip_structs(loop).align_stop-clip_structs(loop).align_start+1) / clip_structs(loop).fps;
        if (is_reframing_indicated(clip_structs(loop)))
            tslice_length_sec = tslice_length_sec - 1/clip_structs(loop).fps;
        end
        number_tslices = total_tslices(clip_structs(loop),tslice_length_sec);
        if (number_tslices ~= 1)
            error('Error in computing the number of processed aligned frames and their time duration');
        end
        [Y_GLOBAL, CB_GLOBAL, CR_GLOBAL] = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec, 1, ...
            'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
        [nrows,ncols,nframes] = size(Y_GLOBAL);  % must be a VFD result for each of the fields/frames in nframes
        
        % Calculate the temporal size in frames so feature extraction
        % subroutine does not need to be passed the clip fps.
        tsize_frames = time_size*clip_structs(loop).fps;  % could be a fractional number of frames
        
        % Compute processed feature data_yrms and data_ymean (if requested)
        % from global variable Y_GLOBAL, and similarly for Cb and Cr.
        if (optional_mean)
            [data_yrms data_cbrms data_crrms data_ymean data_cbmean data_crmean] = feature_ti(vsize, hsize, tsize_frames);
        else
            [data_yrms data_cbrms data_crrms] = feature_ti(vsize, hsize, tsize_frames);
        end
        
        
        % Perform VFD correction on the original (unless this clip is an
        % original) and extract its features.
        
        if (is_original)  % No need for VFD correction of original so just copy
            data_yrmso = data_yrms;
            data_cbrmso = data_cbrms;
            data_crrmso = data_crrms;
            if (optional_mean)
                data_ymeano = data_ymean;
                data_cbmeano = data_cbmean;
                data_crmeano = data_crmean;
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
            tslice_length_sec = (this_clipset(vfd_ind_orig).loc_stop - ...
                this_clipset(vfd_ind_orig).loc_start + 1) / this_clipset(vfd_ind_orig).fps;
            number_tslices = total_tslices(this_clipset(vfd_ind_orig),tslice_length_sec, 'unaligned');
            if number_tslices ~= 1,
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
            [Y_GLOBAL, CB_GLOBAL, CR_GLOBAL] = read_tslice(test_structs(tnum),clip_structs(ind_orig),tslice_length_sec, 1, ...
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
            
            % new_origy will be the VFD-corrected original.  Use same rule
            % as read_tslice for image precision
            if (nrows > 650)
                new_origy = zeros(nrows,ncols,nframes,'single');
                new_origcb = zeros(nrows,ncols,nframes,'single');
                new_origcr = zeros(nrows,ncols,nframes,'single');
            else
                new_origy = zeros(nrows,ncols,nframes,'double');
                new_origcb = zeros(nrows,ncols,nframes,'double');
                new_origcr = zeros(nrows,ncols,nframes,'double');
            end
            
            if(is_interlaced)
                
                for j = 1:nframes
                    % Get matching original field for the early processed field
                    orig_frame_num = ceil(this_result(2*j-1)/2);  % The frame number that contains the original field
                    early_field = mod(this_result(2*j-1),2);  % =1 if early field, =0 if late field
                    [y1 y2] = split_into_fields(squeeze(Y_GLOBAL(:,:,orig_frame_num)));
                    [cb1 cb2] = split_into_fields(squeeze(CB_GLOBAL(:,:,orig_frame_num)));
                    [cr1 cr2] = split_into_fields(squeeze(CR_GLOBAL(:,:,orig_frame_num)));
                    if (early_field)
                        switch field_first
                            case(1)
                                this_origy1 = y1;
                                this_origcb1 = cb1;
                                this_origcr1 = cr1;
                            case(2)
                                this_origy1 = y2;
                                this_origcb1 = cb2;
                                this_origcr1 = cr2;
                        end
                    else % late field
                        switch field_first
                            case(1)
                                this_origy1 = y2;
                                this_origcb1 = cb2;
                                this_origcr1 = cr2;
                            case(2)
                                this_origy1 = y1;
                                this_origcb1 = cb1;
                                this_origcr1 = cr1;
                        end
                    end
                    % Get matching original field for the late processed field
                    orig_frame_num = ceil(this_result(2*j)/2);  % The frame number that contains the original field
                    early_field = mod(this_result(2*j),2);  % =1 if early field, =0 if late field
                    [y1 y2] = split_into_fields(squeeze(Y_GLOBAL(:,:,orig_frame_num)));
                    [cb1 cb2] = split_into_fields(squeeze(CB_GLOBAL(:,:,orig_frame_num)));
                    [cr1 cr2] = split_into_fields(squeeze(CR_GLOBAL(:,:,orig_frame_num)));
                    if (early_field)
                        switch field_first
                            case(1)
                                this_origy2 = y1;
                                this_origcb2 = cb1;
                                this_origcr2 = cr1;
                            case(2)
                                this_origy2 = y2;
                                this_origcb2 = cb2;
                                this_origcr2 = cr2;
                        end
                    else % late field
                        switch field_first
                            case(1)
                                this_origy2 = y2;
                                this_origcb2 = cb2;
                                this_origcr2 = cr2;
                            case(2)
                                this_origy2 = y1;
                                this_origcb2 = cb1;
                                this_origcr2 = cr1;
                        end
                    end
                    % Joint the two original fields into a frame
                    switch field_first
                        case(1)
                            this_origy = join_into_frames(this_origy1,this_origy2);
                            this_origcb = join_into_frames(this_origcb1,this_origcb2);
                            this_origcr = join_into_frames(this_origcr1,this_origcr2);
                        case(2)
                            this_origy = join_into_frames(this_origy2,this_origy1);
                            this_origcb = join_into_frames(this_origcb2,this_origcb1);
                            this_origcr = join_into_frames(this_origcr2,this_origcr1);
                    end
                    new_origy(:,:,j) = this_origy;
                    new_origcb(:,:,j) = this_origcb;
                    new_origcr(:,:,j) = this_origcr;
                end
                clear this_origy this_origy1 this_origy2 y1 y2;
                clear this_origcb this_origcb1 this_origcb2 cb1 cb2;
                clear this_origcr this_origcr1 this_origcr2 cr1 cr2;
                
            else  % progressive
                
                for j = 1:nframes
                    new_origy(:,:,j) = Y_GLOBAL(:,:,this_result(j));
                    new_origcb(:,:,j) = CB_GLOBAL(:,:,this_result(j));
                    new_origcr(:,:,j) = CR_GLOBAL(:,:,this_result(j));
                end
                
            end
            
            Y_GLOBAL = new_origy;  % Y_GLOBAL and new_origy now occupy the same memory space so can't clear new_origy yet
            CB_GLOBAL = new_origcb;
            CR_GLOBAL = new_origcr;
            
            % Compute original feature data_yrmso and data_ymeano (if requested)
            % from global variable Y_GLOBAL, and similarly for Cb and Cr.
            if (optional_mean)
                [data_yrmso data_cbrmso data_crrmso data_ymeano data_cbmeano data_crmeano] = feature_ti(vsize, hsize, tsize_frames);
            else
                [data_yrmso data_cbrmso data_crrmso] = feature_ti(vsize, hsize, tsize_frames);
            end
            
            clear new_origy;  % Has the effect of clearing Y_GLOBAL memory space
            clear new_origcb;
            clear new_origcr;
            
        end
        
        % Write features all at once to their respective files
        write_feature_standard(data_yrms, feature_base_dir, feature_yrms, clip_structs(loop), 'vfd', data_yrmso);
        write_feature_standard(data_cbrms, feature_base_dir, feature_cbrms, clip_structs(loop), 'vfd', data_cbrmso);
        write_feature_standard(data_crrms, feature_base_dir, feature_crrms, clip_structs(loop), 'vfd', data_crrmso);
        if (optional_mean)
            write_feature_standard(data_ymean, feature_base_dir, feature_ymean, clip_structs(loop), 'vfd', data_ymeano);
            write_feature_standard(data_cbmean, feature_base_dir, feature_cbmean, clip_structs(loop), 'vfd', data_cbmeano);
            write_feature_standard(data_crmean, feature_base_dir, feature_crmean, clip_structs(loop), 'vfd', data_crmeano);
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
function [data_yrms data_cbrms data_crrms data_ymean data_cbmean data_crmean] = feature_ti(vsize, hsize, tsize_frames)
% FEATURE_YCBCR
%  Compute Y, Cb, and Cr Temporal Information (TI) block features (rms and
%  mean) given an entire video clip, and the desired block size of vsize by
%  hsize by tsize_frames.  The vsize and hsize must evenly divide the
%  vertical and horizontal size of Y_GLOBAL, CB_GLOBAL, and CR_GLOBAL
%  (global variables that holds the entire video clip).  Note that
%  tsize_frames may be fractional.  
%
% SYNTAX
%  [data_yrms data_cbrms data_crrms data_ymean data_cbmean data_crmean] = feature_ti(vsize, hsize, tsize_frames)
%
% DESCRIPTION
%  The function takes 3D global variables (Y_GLOBAL, CB_GLOBAL, and
%  CR_GLOBAL) that holds the video clip and extracts the Y, Cb, and Cr
%  temporal information RMS (data_yrms, data_cbrms,
%  data_crrms) and mean (data_ymean, data_cbmean, data_crmean) block
%  features from blocks of size vsize by hsize by tsize_frames.  If the
%  size(Y_GLOBAL) = [nrows, ncols, nframes], then the size of the returned
%  feature arrays will be nrows/vsize, ncols/hsize, and
%  floor(nframes/tsize_frames) (i.e., the subroutines returns the maximum 
%  number of time slices).  For fractional tsize_frames, some
%  spatial-temporal (ST) blocks will have more pixels than others but 
%  a given pixel will only be included in one and only one block.  Also,
%  the TI is computed using one extra frame after the block (in time) if
%  available. When the floor(nframes/tsize_frames) = (nframes/tsize_frames), 
%  the last block in the clip will not have a following frame so its TI
%  will be based on one less frame.
%

global Y_GLOBAL CB_GLOBAL CR_GLOBAL;

%  Check the spatial size of Y_GLOBAL for validity
[nrows, ncols, nframes] = size(Y_GLOBAL);
if (rem(nrows,vsize) || rem(ncols,hsize))
    error('Invalid vsize or hsize detected in function feature_ti');
end

% Compute Spatial block sizes
nvert = nrows/vsize;  % Number of blocks in vertical dim
nhoriz = ncols/hsize;  % Number of blocks in horizontal dim

% Compute the number of time slices
nslices = floor(nframes/tsize_frames);  % Maximum number of time slices

%  Dimension output arrays as needed
if (nargout > 3)  % Mean requested
    optional_mean = 1;
    data_yrms = zeros(nvert,nhoriz,nslices);
    data_cbrms = zeros(nvert,nhoriz,nslices);
    data_crrms = zeros(nvert,nhoriz,nslices);
    data_ymean = zeros(nvert,nhoriz,nslices);
    data_cbmean = zeros(nvert,nhoriz,nslices);
    data_crmean = zeros(nvert,nhoriz,nslices);
else  % Mean not requested
    optional_mean = 0;
    data_yrms = zeros(nvert,nhoriz,nslices);
    data_cbrms = zeros(nvert,nhoriz,nslices);
    data_crrms = zeros(nvert,nhoriz,nslices);
end

% For each time-slice in this clip, compute the requested block statistic.
for cnt = 1:nslices
    
    beg_frame = 1 + floor((cnt-1)*tsize_frames);
    end_frame = floor(cnt*tsize_frames);
    
%     %  Add extra frame at the end of block if available (for the diff
%     %  function, since you will loose one frame).
%     if (end_frame < nframes)
%         end_frame = end_frame+1;
%     end
    
    %  Add extra frame at the beginning of block if available (for the diff
    %  function, since you will loose one frame).
    %  This is probably more correct since TI is normally defined as 
    %  motion in current frame w.r.t. previous frame.
    if (beg_frame > 1)
        beg_frame = beg_frame-1;
    end
    
    %  Perform the TI calculation
    ydiff = diff(Y_GLOBAL(:,:,beg_frame:end_frame), 1, 3);  % 1st order difference in 3rd dimension
    cbdiff = diff(CB_GLOBAL(:,:,beg_frame:end_frame), 1, 3);
    crdiff = diff(CR_GLOBAL(:,:,beg_frame:end_frame), 1, 3);
    
    if (~optional_mean)
        data_yrms(:,:,cnt) = block_statistic (ydiff, vsize, hsize, 'rms');
        data_cbrms(:,:,cnt) = block_statistic (cbdiff, vsize, hsize, 'rms');
        data_crrms(:,:,cnt) = block_statistic (crdiff, vsize, hsize, 'rms');
    else
        [data_yrms(:,:,cnt),data_ymean(:,:,cnt)] = block_statistic (ydiff, vsize, hsize, 'rms','mean');
        [data_cbrms(:,:,cnt),data_cbmean(:,:,cnt)] = block_statistic (cbdiff, vsize, hsize, 'rms','mean');
        [data_crrms(:,:,cnt),data_crmean(:,:,cnt)] = block_statistic (crdiff, vsize, hsize, 'rms','mean');
    end
    
end

% Has the effect of clearing Y_GLOBAL memory space if the space is 
% not occupied by another variable (e.g., new_origy) 
Y_GLOBAL = 0;
CB_GLOBAL = 0;
CR_GLOBAL = 0;

end
