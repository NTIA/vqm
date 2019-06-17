function [data, datao, data2, data2o, status] = dll_vfd_feature_loop_mse(...
    clip_structs, deg_size, time_size, viewing_distance)
% VFD_FEATURE_LOOP_MSE
%  Take a pair of clips (an original and a processed), compute the variable
%  frame delay (VFD) Mean Squared Error (MSE) between the processed and the
%  original, and output the processed clip's features.  Only the Y image is
%  used for this computation so the features can be used to compute the
%  traditional PSNR based on the Y channel only.
% SYNTAX
%  [data, datao, data2, data2o, success] = dll_vfd_feature_loop_mse(...
%       clip_structs, deg_size, time_size, viewing_distance);
% DESCRIPTION
%  Compute VFD MSE features
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

% Y_GLOBAL used to be a global variable.

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

nrows = sroi.bottom-sroi.top+1;
ncols = sroi.right-sroi.left+1;

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
y_proc = dll_calib_video('tslice', 2);
dll_video('rewind', 2);

% Calculate the temporal size in frames so feature extraction
% subroutine does not need to be passed the clip fps.
tsize_frames = time_size*clip_structs(2).fps;  % could be a fractional number of frames

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
y_orig = dll_calib_video('tslice', 1);
dll_video('rewind', 1);

Y_GLOBAL = y_proc - y_orig;
clear y_proc y_orig;

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
optional_mean = 1;
data = zeros(nvert,nhoriz,nslices);
data2 = zeros(nvert,nhoriz,nslices);

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

[size1,size2,size3] = size(data);
datao = zeros(size1,size2,size3);
data2o = datao;






