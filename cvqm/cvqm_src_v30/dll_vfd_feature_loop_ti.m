function [data_yrms, data_yrmso, data_cbrms, data_cbrmso, data_crrms, data_crrmso, ...
    data_ymean, data_ymeano, data_cbmean, data_cbmeano, data_crmean, data_crmeano] = ...
    dll_vfd_feature_loop_ti(clip_structs, deg_size, time_size, viewing_distance)
% VFD_FEATURE_LOOP_TI
%  Take a pair of clips (an original and a processed), compute variable
%  frame delay (VFD) RMS Temporal Information (TI) for Y, Cb, and CR, and
%  output the processed clip's features.
% SYNTAX
%  [data_yrms, data_yrmso, data_cbrms, data_cbrmso, data_crrms, data_crrmso, ...
%      data_ymean, data_ymeano, data_cbmean, data_cbmeano, data_crmean, data_crmeano] = ...
%      vfd_feature_loop_ti(clip_structs, deg_size, time_size, viewing_distance);
% DESCRIPTION
%  Compute VFD Temporal Information (TI) features for Y, Cb, and Cr.
% 
%  This function takes variable 'clip_structs' (of the same format as
%  GClips), which contains information on the pair of clips, 'deg_size'
%  specifies the block size in angular degrees (horizontally & vertically),
%  'time_size' specifies the length of the block in SECONDS, and
%  'viewing_distance' specifies the viewing distance for the data set in
%  PICTURE HEIGHTS.
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


status = 0;
% Y_GLOBAL, CB_GLOBAL, and CR_GLOBAL used to be global variables.

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

dll_calib_video('sroi', sroi, 0);
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
[Y_GLOBAL, CB_GLOBAL, CR_GLOBAL] = dll_calib_video('tslice', 2);
dll_video('rewind', 2);

% Calculate the temporal size in frames so feature extraction
% subroutine does not need to be passed the clip fps.
tsize_frames = time_size*clip_structs(2).fps;  % could be a fractional number of frames

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
optional_mean = 1;
data_yrms = zeros(nvert,nhoriz,nslices);
data_cbrms = zeros(nvert,nhoriz,nslices);
data_crrms = zeros(nvert,nhoriz,nslices);
data_ymean = zeros(nvert,nhoriz,nslices);
data_cbmean = zeros(nvert,nhoriz,nslices);
data_crmean = zeros(nvert,nhoriz,nslices);


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
[Y_GLOBAL, CB_GLOBAL, CR_GLOBAL] = dll_calib_video('tslice', 1);
dll_video('rewind', 1);

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
optional_mean = 1;
data_yrmso = zeros(nvert,nhoriz,nslices);
data_cbrmso = zeros(nvert,nhoriz,nslices);
data_crrmso = zeros(nvert,nhoriz,nslices);
data_ymeano = zeros(nvert,nhoriz,nslices);
data_cbmeano = zeros(nvert,nhoriz,nslices);
data_crmeano = zeros(nvert,nhoriz,nslices);


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
        data_yrmso(:,:,cnt) = block_statistic (ydiff, vsize, hsize, 'rms');
        data_cbrmso(:,:,cnt) = block_statistic (cbdiff, vsize, hsize, 'rms');
        data_crrmso(:,:,cnt) = block_statistic (crdiff, vsize, hsize, 'rms');
    else
        [data_yrmso(:,:,cnt),data_ymeano(:,:,cnt)] = block_statistic (ydiff, vsize, hsize, 'rms','mean');
        [data_cbrmso(:,:,cnt),data_cbmeano(:,:,cnt)] = block_statistic (cbdiff, vsize, hsize, 'rms','mean');
        [data_crrmso(:,:,cnt),data_crmeano(:,:,cnt)] = block_statistic (crdiff, vsize, hsize, 'rms','mean');
    end
    
end

% Has the effect of clearing Y_GLOBAL memory space if the space is 
% not occupied by another variable (e.g., new_origy) 
Y_GLOBAL = 0;
CB_GLOBAL = 0;
CR_GLOBAL = 0;











