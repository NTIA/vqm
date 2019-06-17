function [data] = model_fastlowbw_features_shift (mode, one, two, three, four, five, six);
% model_fastlowbw_features_shift
%   Compute the processed features for the ITS fast low bandwidth model.
% SYNTAX
%   [data] = model_fastlowbw_features_shift ('clip', otest, one_clip, sroi)
%   model_fastlowbw_features_shift ('memory',  y, cb, cr, fps, filter_size, extra)
%   [data] = model_fastlowbw_features_shift ('eof')
%   model_fastlowbw_features_shift ('clear')
% DESCRIPTION
%   When the first argument is the string 'clip', this function takes a
%   test structure (otest) and ONE clip structure (one_clip) of the 
%   same format as GTests and GClips.  Function will process entire clip
%   (e.g., 1-minute of video).  
%
%   When the first argument is the string 'memory', this function takes one
%   second (EXACTLY) of video already in memory, in the YCbCr colorspace
%   (y, cb, and cr respectively).  All three variables are formatted
%   (row,col,time); and (fps) is the number of frames per second (e.g.,
%   29.97, 30, 25, 15, etc).  'filter_size' and 'extra' are from function
%   'adaptive_filter' when given the current image size.
%   Size of y, cb, & cr spatially should be exactly ROI returned by
%   function model_lowbw_sroi.  See 'eof. 
%
%   When the first argument is the string 'eof', this function completes
%   the feature calculations, empties internal buffers, and returns the
%   feature stream.  See 'memory'.
%
%   When the first argument is the string 'clear', the internal buffers are
%   emptied.  This is done when called with 'eof', also. 
%
%   Warning:  uses 1 set of internal buffers, so don't mix sequences.  
%
%   Return value 'data' is contains for each of the 9 pixel shifts --
%   data(1) .. data(9) -- the following features as structure elements:
%           si_std hv_ratio y_mean cb_mean cr_mean ati_rms
%   See function model_fastlowbw_features for a brief description of each
%   element in this structure.  

    
    

%   

persistent buffer;
persistent buffer_y;

BSIZE = 30;


if strcmpi(mode,'clip') && nargin == 4,
    % error check.
    if length(buffer) ~= 0,
        error('Mixing video clips; must clear buffers with eof call');
    end
    tslice_sec = 1.0;
    [filter_size,extra] = adaptive_filter(two.image_size); 
    for number = 1:total_tslices(two, tslice_sec),
        % read images
        [y, cb, cr] = read_tslice (one, two, tslice_sec, number, 'extra', extra+1, ...
            'hsize', BSIZE, 'vsize', BSIZE, 'sroi', three.top, three.left, three.bottom, three.right );
        model_fastlowbw_features_shift ('memory',  y, cb, cr, two.fps, filter_size, extra);
        clear y cb cr;
    end
    if size(buffer(1).si_std,3) ~= total_tslices(two,tslice_sec),
        error('count somehow off -- fatal length mismatch');
    end
    [data] = model_fastlowbw_features_shift ('eof');
    
elseif strcmpi(mode,'memory') && nargin == 7,
    tis_sec = 0.2;
    [tis_frames] = tslice_conversion(tis_sec, four);
    
    [curr] = ...
        model_fastlowbw_features_memory (one, two, three, four, buffer_y, five, six, BSIZE);
    
    if length(buffer) == 0,
        buffer = curr;
    else
        [row,col,time1] = size(buffer(1).si_std);
        [row,col,time2] = size(buffer(1).ati_rms);
        [row,col,time3] = size(curr(1).ati_rms);
        for loop = 1:9,
            buffer(loop).si_std(:,:,time1+1) = curr(loop).si_std;
            buffer(loop).hv_ratio(:,:,time1+1) = curr(loop).hv_ratio;
            buffer(loop).y_mean(:,:,time1+1) = curr(loop).y_mean;
            buffer(loop).cb_mean(:,:,time1+1) = curr(loop).cb_mean;
            buffer(loop).cr_mean(:,:,time1+1) = curr(loop).cr_mean;

            buffer(loop).ati_rms(:,:,time2+1:time2+time3) = curr(loop).ati_rms;
        end
    end
    [row,col,time] = size(one);
    buffer_y = one(:,:,(time-tis_frames+1):time);
elseif strcmpi(mode,'eof') && nargin == 1,
    
    data = buffer;
    model_fastlowbw_features_shift('clear');
    
elseif strcmpi(mode,'clear'),
    buffer = [];
    buffer_y = [];
else
    error('argument list not recognized');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [curr] = ...
    model_fastlowbw_features_memory (y, cb, cr, fps, buffer_y, filter_size, extra, BSIZE);

HV_THRESHOLD = 4.0;

tis_sec = 0.2;
[tis_frames] = tslice_conversion(tis_sec, fps);

% loop through video
tslice_frames = round(fps);
ati_curr = 1;
[row,col,time] = size(y);
rng1 = (extra+1):(row-extra);
rng2 = (extra+1):(col-extra);
    
% Calculate features.
ym = mean(y,3);
[si, hv, hvb] = filter_si_hv_adapt(ym, filter_size, extra); 

si_std = block_statistic_shift(si, BSIZE, BSIZE, 'std');
hv_mean = block_statistic_shift(hv, BSIZE, BSIZE, 'mean');
hvb_mean = block_statistic_shift(hvb, BSIZE, BSIZE, 'mean');
y_mean = block_statistic_shift(ym(rng1,rng2,:), BSIZE, BSIZE, 'mean');
cb_mean = block_statistic_shift(cb(rng1,rng2,:), BSIZE, BSIZE, 'mean');
cr_mean = block_statistic_shift(cr(rng1,rng2,:), BSIZE, BSIZE, 'mean');

clear si hv hvb;
for loop = 1:9,
    curr(loop).si_std = si_std(loop).std;
    curr(loop).hv_mean = hv_mean(loop).mean;
    curr(loop).hvb_mean = hvb_mean(loop).mean;
    curr(loop).y_mean = y_mean(loop).mean;
    curr(loop).cb_mean = cb_mean(loop).mean;
    curr(loop).cr_mean = cr_mean(loop).mean;
end

% compute Y 0.2s ATI on frames
t1 = 0;
if length(buffer_y) > 0,
    [r,c,t1] = size(buffer_y);
    [r,c,t2] = size(y);
    ati_y(:,:,1:t1) = buffer_y(rng1,rng2,:);
    ati_y(:,:,t1+1:t1+t2) = y(rng1,rng2,:);
else
    ati_y = y;
end

% use 5% of pixels, randomly chosen
ati_y = filter_ati_random(ati_y, tis_frames, 0.05);

% have time-differences; now compute ATI
[rowcol,time] = size(ati_y);
for cnt = 1:time,
    hold = block_statistic( ati_y(:,cnt), rowcol, 1, 'rms');
    for loop = 1:length(curr),
        curr(loop).ati_rms(1,1,ati_curr) = hold;
    end
    ati_curr = ati_curr + 1;
end

for loop = 1:length(curr),
    curr(loop).hv_ratio = max(HV_THRESHOLD, curr(loop).hv_mean) ...
        ./ max(HV_THRESHOLD, curr(loop).hvb_mean);
end

