function [si_std hv_ratio y_mean cb_mean cr_mean ati_rms] = ...
    model_fastlowbw_features (mode, one, two, three, four, five, six);
% MODEL_fastlowbw_features
%   Compute the original features for the ITS fast low bandwidth model.
% SYNTAX
%   [si_std, hv_ratio, y_mean, cb_mean, cr_mean, ati_rms] = ...
%       model_fastlowbw_features ('clip', otest, one_clip, sroi);
% 
%   model_fastlowbw_features ('memory',  y, cb, cr, fps, filter_size, extra);
% 
%   [si_std, hv_ratio, y_mean, cb_mean, cr_mean, ati_rms] = ...
%       model_fastlowbw_features ('eof');
%
%   model_fastlowbw_features ('clear');
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
%   When the first argument is the string 'clear', all internall bufferes
%   are emptied.  This is done when 'eof' is called, also.  
%
%   Warning:  uses 1 set of internal buffers, so don't mix sequences.  
%
%   Return variables when called with 'eof' are the six feature streams:
%   - si_std containing standard deviation of the spatial information (SI)
%   - hv_ratio containing the ratio of HV to HVbar energy
%   - cb_mean containing the average Cb value
%   - cr_mean containing the averge Cr value
%   - ati_rms containing the root mean squared of absolute value of
%     temporal information (TI)

persistent hold_si_std;
persistent hold_hv_ratio;
persistent hold_y_mean;
persistent hold_cb_mean;
persistent hold_cr_mean;
persistent hold_ati_rms;
persistent buffer_y;


if strcmpi(mode,'clip') && nargin == 4,
    % error check.
    if length(hold_si_std) ~= 0,
        error('Mixing video clips; must clear buffers with eof call');
    end
    tslice_sec = 1.0;
    [filter_size,extra] = adaptive_filter(two.image_size); 
    for number = 1:total_tslices(two, tslice_sec),
        % read images
        [y, cb, cr] = read_tslice (one, two, tslice_sec, number, 'extra', extra, ...
            'hsize', 30, 'vsize', 30, 'sroi', three.top, three.left, three.bottom, three.right );
        model_fastlowbw_features ('memory',  y, cb, cr, two.fps, filter_size, extra);
        clear y cb cr;
    end
    if size(hold_si_std,3) ~= total_tslices(two,tslice_sec),
        error('count somehow off -- fatal length mismatch');
    end
    [si_std hv_ratio y_mean cb_mean cr_mean ati_rms] = model_fastlowbw_features ('eof');
    
elseif strcmpi(mode,'memory') && nargin == 7,
    tis_sec = 0.2;
    [tis_frames] = tslice_conversion(tis_sec, four);
    
    [si_std hv_ratio y_mean cb_mean cr_mean ati_rms] = ...
        model_fastlowbw_features_memory (one, two, three, four, buffer_y, five, six);
    
    if length(hold_si_std) == 0,
        hold_si_std = si_std;
        hold_hv_ratio = hv_ratio;
        hold_y_mean = y_mean;
        hold_cb_mean = cb_mean;
        hold_cr_mean = cr_mean;
        hold_ati_rms = ati_rms;
    else
        [row,col,time1] = size(hold_si_std);
        hold_si_std(:,:,time1+1) = si_std;
        hold_hv_ratio(:,:,time1+1) = hv_ratio;
        hold_y_mean(:,:,time1+1) = y_mean;
        hold_cb_mean(:,:,time1+1) = cb_mean;
        hold_cr_mean(:,:,time1+1) = cr_mean;

        [row,col,time1] = size(hold_ati_rms);
        [row,col,time2] = size(ati_rms);
        hold_ati_rms(:,:,time1+1:time1+time2) = ati_rms;
    end
    [row,col,time] = size(one);
    buffer_y = one(:,:,(time-tis_frames+1):time);
elseif strcmpi(mode,'eof') && nargin == 1,
    
    si_std = hold_si_std;
    hv_ratio = hold_hv_ratio;
    y_mean = hold_y_mean;
    cb_mean = hold_cb_mean;
    cr_mean = hold_cr_mean;
    ati_rms = hold_ati_rms;
    model_fastlowbw_features('clear');
    
elseif strcmpi(mode,'clear'),
    hold_si_std = [];
    hold_hv_ratio = [];
    hold_y_mean = [];
    hold_cb_mean = [];
    hold_cr_mean = [];
    hold_ati_rms = [];
    buffer_y = [];
else
    error('argument list not recognized');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [si_std hv_ratio y_mean cb_mean cr_mean ati_rms] = ...
    model_fastlowbw_features_memory (y, cb, cr, fps, buffer_y, filter_size, extra);

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
[si_std] = block_statistic(si, 30, 30, 'std');
[hv_mean] = block_statistic(hv, 30, 30, 'mean');
[hvb_mean] = block_statistic(hvb, 30, 30, 'mean');
[r,c,t] = size(y);
[y_mean] = block_statistic(ym(rng1, rng2, :), 30, 30, 'mean');
clear si hv hvb;

% Compute Cb, Cr
[cb_mean] = block_statistic(cb(rng1,rng2,:), 30, 30, 'mean');
[cr_mean] = block_statistic(cr(rng1,rng2,:), 30, 30, 'mean');

% compute YCbCr 0.2s ATI on frames
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
    ati_rms(1,1,ati_curr) = block_statistic( ati_y(:,cnt), rowcol, 1, 'rms');
    ati_curr = ati_curr + 1;
end

hv_ratio = max(HV_THRESHOLD, hv_mean) ./ max(HV_THRESHOLD, hvb_mean);

