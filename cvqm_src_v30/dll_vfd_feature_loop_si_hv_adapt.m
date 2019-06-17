function [datasi, datasio, datahv, datahvo, datahvb, datahvbo, status] = ...
    dll_vfd_feature_loop_si_hv_adapt(clip_structs, deg_size, time_size, viewing_distance)
% VFD_FEATURE_LOOP_SI_HV_ADAPT
%  Take a pair of clips (an original and a processed), compute variable
%  frame delay (VFD) SI, HV, and HVbar, and output the processed clip's
%  features.
% SYNTAX
%  [datasi, datasio, datahv, datahvo, datahvb, datahvbo, status] = ...
%       dll_vfd_feature_loop_si_hv_adapt(clip_structs, deg_size, time_size, viewing_distance)
% DESCRIPTION
%  Compute VFD SI, HV and HVbar features
% 
%  This function takes in a variable 'clip_structs' (of the same format as
%  GClips), which includes information about the pair of clips to be
%  processed, 'deg_size' specifies the block size in angular degrees
%  (horizontally & vertically), 'time_size' specifies the length of the
%  block in SECONDS, and 'viewing_distance' specifies the viewing distance
%  for the data set in PICTURE HEIGHTS.
%
%  For interlaced systems, the block size (in
%  pixels/lines) for a desired 'deg_size' is forced to be even so an equal
%  number of lines from each field are included in the block (blocks are
%  calculated from frames, not fields, with this routine).  The nearest
%  even block size is chosen for interlaced systems.
%
%  Return the feature arrays as well as 0 if operated correctly; 1 if an
%  error was encountered.
% 


status = 0;  % Normal return status


% Y_GLOBAL used to be a global variable

% try
    %  Calculate the vsize and hsize for the requested deg_size using
    %  the assumptions in the header documentation.  The midpoints
    %  between the vertical sizes are used to assign image types (e.g.,
    %  QCIF, CIF, VGA, 720, 1080).
    image_rows = clip_structs(2).image_size.rows;
    image_standard = clip_structs(2).video_standard;
    vsize = viewing_distance*deg_size*image_rows*pi/180;
    if (strcmpi(image_standard,'progressive'))
        vsize = round(vsize);
    else  % Interlaced, set vsize to nearest even number
        vsize = round(vsize/2)*2;
    end
    if (vsize <=1 || vsize >= image_rows)
        error('Invalid deg_size');
    end
    hsize = vsize;

    % find adaptive filter size
    [filter_size, extra] = adaptive_filter (clip_structs(2).image_size);
    
    % Compute the default adjusted SROI and number of blocks available.
    if (strcmpi(clip_structs(2).video_standard,'progressive'))
        [sroi] = adjust_requested_sroi (clip_structs(2), 'vsize',vsize, 'hsize',hsize, 'extra',extra);
    else
        [sroi] = adjust_requested_sroi (clip_structs(2), 'vsize',vsize, 'hsize',hsize, 'extra',extra, 'evenodd');
    end
    
    dll_calib_video('sroi', sroi, extra);
    clip_structs(1).cvr = sroi;
    clip_structs(2).cvr = sroi;

    % Determine the total number of aligned frames that are available
    % (in seconds) and read in the entire clip.
    tslice_length_sec = (clip_structs(2).align_stop-clip_structs(2).align_start+1) / clip_structs(2).fps;
    if (is_reframing_indicated(clip_structs(2)))
        tslice_length_sec = tslice_length_sec - 1/clip_structs(2).fps;
    end
    number_tslices = total_tslices(clip_structs(2),tslice_length_sec);
    if (number_tslices ~= 1)
        error('Error in computing the number of processed aligned frames and their time duration');
    end

    dll_video('set_rewind', 2);

    dll_video('set_tslice', 2, tslice_length_sec);
    Y_GLOBAL = dll_calib_video('tslice', 2);
    dll_video('rewind', 2);


    % Calculate the temporal size in frames so feature extraction
    % subroutine does not need to be passed the clip fps.
    tsize_frames = time_size*clip_structs(2).fps;  % could be a fractional number of frames

    % Compute processed feature data si, hv, and hvb from global
    % variable Y_GLOBAL.
    % Check the spatial size of Y_GLOBAL for validity
    [nrows, ncols, nframes] = size(Y_GLOBAL);
    if (rem(nrows-2*extra,vsize) || rem(ncols-2*extra,hsize))
        error('Invalid vsize or hsize detected in function feature_si_hv_adapt');
    end

    % Compute Spatial block sizes
    nvert = (nrows-2*extra)/vsize;  % Number of blocks in vertical dim
    nhoriz = (ncols-2*extra)/hsize;  % Number of blocks in horizontal dim

    % Compute the number of time slices
    nslices = floor(nframes/tsize_frames);  % Maximum number of time slices

    %  Dimension output arrays as needed
    datasi = zeros(nvert,nhoriz,nslices);
    datahv = zeros(nvert,nhoriz,nslices);
    datahvb = zeros(nvert,nhoriz,nslices);

    % For each time-slice in this clip, compute the requested block statistic.
    for cnt = 1:nslices

        beg_frame = 1 + floor((cnt-1)*tsize_frames);
        end_frame = floor(cnt*tsize_frames);

        % Compute the si_hv adaptive filter on this time slice
        [si,hv,hvb] = filter_si_hv_adapt(Y_GLOBAL(:,:,beg_frame:end_frame), filter_size, extra);

        %  Extract the features
        datasi(:,:,cnt) = block_statistic (si, vsize, hsize, 'std');
        datahv(:,:,cnt) = block_statistic (hv, vsize, hsize, 'mean');
        datahvb(:,:,cnt) = block_statistic (hvb, vsize, hsize, 'mean');

        clear si hv hvb;

    end

    % Has the effect of clearing Y_GLOBAL memory space if the space is 
    % not occupied by another variable (e.g., new_orig) 
    Y_GLOBAL = 0;
    
    tslice_length_sec = (clip_structs(1).align_stop-clip_structs(1).align_start+1) / clip_structs(1).fps;
    if (is_reframing_indicated(clip_structs(1)))
        tslice_length_sec = tslice_length_sec - 1/clip_structs(1).fps;
    end
    number_tslices = total_tslices(clip_structs(1),tslice_length_sec);
    if (number_tslices ~= 1)
        error('Error in computing the number of processed aligned frames and their time duration');
    end
    
    dll_video('set_rewind', 1);
    dll_video('set_tslice', 1, tslice_length_sec);
    Y_GLOBAL = dll_calib_video('tslice', 1);
    dll_video('rewind', 1);
    
    [nrows, ncols, nframes] = size(Y_GLOBAL);
    if (rem(nrows-2*extra,vsize) || rem(ncols-2*extra,hsize))
        error('Invalid vsize or hsize detected in function feature_si_hv_adapt');
    end

    % Compute Spatial block sizes
    nvert = (nrows-2*extra)/vsize;  % Number of blocks in vertical dim
    nhoriz = (ncols-2*extra)/hsize;  % Number of blocks in horizontal dim

    % Compute the number of time slices
    nslices = floor(nframes/tsize_frames);  % Maximum number of time slices

    %  Dimension output arrays as needed
    datasio = zeros(nvert,nhoriz,nslices);
    datahvo = zeros(nvert,nhoriz,nslices);
    datahvbo = zeros(nvert,nhoriz,nslices);

    % For each time-slice in this clip, compute the requested block statistic.
    for cnt = 1:nslices

        beg_frame = 1 + floor((cnt-1)*tsize_frames);
        end_frame = floor(cnt*tsize_frames);

        % Compute the si_hv adaptive filter on this time slice
        [si,hv,hvb] = filter_si_hv_adapt(Y_GLOBAL(:,:,beg_frame:end_frame), filter_size, extra);

        %  Extract the features
        datasio(:,:,cnt) = block_statistic (si, vsize, hsize, 'std');
        datahvo(:,:,cnt) = block_statistic (hv, vsize, hsize, 'mean');
        datahvbo(:,:,cnt) = block_statistic (hvb, vsize, hsize, 'mean');

        clear si hv hvb;

    end
    
% catch
%     status = 1;
%     if (verbose)
%         disp(exception.message);
%         fprintf('\tSkipping clip.\n');
%     end
% end

    
            
        
        
        
        