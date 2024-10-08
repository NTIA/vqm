function [vqm, hv_loss_par, hv_gain_par,si_loss_par, si_gain_par, ...
    color_comb_par, noise_par, error_par] = model_lowbw_parameters (...
    si_proc, hv_feat_proc, y_proc, cb_proc, cr_proc, ati_proc, ...
    si_orig, hv_feat_orig, y_orig, cb_orig, cr_orig, ati_orig, ...
    fps, part_si_min, part_si_max, part_hv_min, part_hv_max, ...
    part_c_min, part_c_max, part_c, part_ati_min, part_ati_max, ...
    part_ati, code_ati, do_test_print, TIME_DELTA);
% MODEL_LOWBW_PARAMETERS
% [vqm, hv_loss_par, hv_gain_par,si_loss_par, si_gain_par, ...
%     color_comb_par, noise_par, error_par] = model_lowbw_parameters (...
%     si_proc, hv_feat_proc, y_proc, cb_proc, cr_proc, ati_proc, ...
%     si_orig, hv_feat_orig, y_orig, cb_orig, cr_orig, ati_orig, ...
%     fps, part_si_min, part_si_max, part_hv_min, part_hv_max, ...
%     part_c_min, part_c_max, part_c, part_ati_min, part_ati_max, ...
%     part_ati, code_ati, do_test_print);
% [vqm, hv_loss_par, hv_gain_par,si_loss_par, si_gain_par, ...
%     color_comb_par, noise_par, error_par] = model_lowbw_parameters (...
%     si_proc, hv_feat_proc, y_proc, cb_proc, cr_proc, ati_proc, ...
%     si_orig, hv_feat_orig, y_orig, cb_orig, cr_orig, ati_orig, ...
%     fps, part_si_min, part_si_max, part_hv_min, part_hv_max, ...
%     part_c_min, part_c_max, part_c, part_ati_min, part_ati_max, ...
%     part_ati, code_ati, do_test_print, TIME_DELTA);
% If 'TIME_DELTA' not defined, then 10 seconds used (10).
%
% The updates from the previous model include *** in the comments field.
% SW 12/28/05

if ~exist('TIME_DELTA'),
    TIME_DELTA = 10;
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate parameters %
%%%%%%%%%%%%%%%%%%%%%%%%
[hv_loss_par, hv_gain_par,si_loss_par, si_gain_par, ...
    color_comb_par, noise_par, error_par] = calc_parameters (...
    si_proc, hv_feat_proc, y_proc, cb_proc, cr_proc, ati_proc, ...
    si_orig, hv_feat_orig, y_orig, cb_orig, cr_orig, ati_orig, ...
    fps, part_si_min, part_si_max, part_hv_min, part_hv_max, ...
    part_c_min, part_c_max, part_c, part_ati_min, part_ati_max, ...
    part_ati, code_ati, TIME_DELTA);

%% clear extra variables
clear ati_orig ati_proc cb_orig cb_proc col cr_orig cr_proc hv_feat_orig;
clear hv_feat_proc part_ati part_ati_max part_ati_min part_c part_c_max;
clear part_c_min part_hv_max part_hv_min part_si_max part_si_min;
clear range row si_orig si_proc y_orig y_orig2 y_proc y_proc2;

%%%%%%%%%%%%%%%%%%%%%
% Fix results   %
%%%%%%%%%%%%%%%%%%%%%
hv_loss_par = [hv_loss_par(1) hv_loss_par'];
hv_gain_par = [hv_gain_par(1) hv_gain_par'];
si_loss_par = [si_loss_par(1) si_loss_par'];
color_comb_par = [color_comb_par(1) color_comb_par'];

% si_gain, noise_par, & error_par don't need a fix -- no macro blocks!

%%%%%%%%%%%%%%%%%%%%%
% calculate VQM     %
%%%%%%%%%%%%%%%%%%%%%
if do_test_print,
    % for testing purposes only!
    fprintf('print out parameters & VQM for testing\n');
    save raw_pars.mat hv_loss_par hv_gain_par si_loss_par si_gain_par color_comb_par noise_par error_par;
end

% ***Weights updated
DC_WEIGHT = 0.0;
HV_LOSS_WEIGHT = 0.68265798470360;
HV_GAIN_WEIGHT = 0.38028815869787;
SI_LOSS_WEIGHT = 0.35526265483436;
SI_GAIN_WEIGHT = 0.69490465292299;
% COLOR_WEIGHT = 1.75166641317702;
COLOR_WEIGHT = 1.00000000000000;
NOISE_WEIGHT = 0.06657446367248;
ERROR_WEIGHT = 0.02751797957319;

% upsample (1 sample per sec --> 2 samples per sec) and  weight each parameter.
hv_loss_par = my_interp(2, hv_loss_par) * HV_LOSS_WEIGHT;
hv_gain_par = my_interp(2, hv_gain_par) * HV_GAIN_WEIGHT;
si_loss_par = my_interp(2, si_loss_par) * SI_LOSS_WEIGHT;
si_gain_par = my_interp(2, si_gain_par) * SI_GAIN_WEIGHT;
color_comb_par = my_interp(2, color_comb_par) * COLOR_WEIGHT;
noise_par = noise_par * NOISE_WEIGHT;
error_par = error_par * ERROR_WEIGHT;

% clear a few variables
clear *par1 *par2 *orig *proc;
clear col video_dir ans code_ati oclip otest part* range row time;

% even out lengths, if needed.  This should never happen if the video
% sequences are really 1-minute long each.
hold = [length(hv_loss_par) length(hv_gain_par) length(si_loss_par) length(si_gain_par) ...
    length(color_comb_par) length(noise_par) length(error_par)];
if min(hold) ~= max(hold),
%    fprintf('WARNING:  time-length of features being rectified\n');
    hold = min(hold);
    hv_loss_par = hv_loss_par(1:hold);
    hv_gain_par = hv_gain_par(1:hold);
    si_loss_par = si_loss_par(1:hold);
    si_gain_par = si_gain_par(1:hold);
    color_comb_par = color_comb_par(1:hold);
    noise_par = noise_par(1:hold);
    error_par = error_par(1:hold);
end

%  Final VQM has clipping at low end and crushing at high end
vqm = DC_WEIGHT + hv_loss_par + hv_gain_par + si_loss_par + si_gain_par + color_comb_par + noise_par + error_par;

% Clip VQM for values less than zero
vqm(find(vqm < 0)) = 0.0;  % No quality improvements allowed

% Compress VQM values greater than 1 using standard crushing function
c = 0.5;
crush = find(vqm > 1);
vqm(crush) = (1 + c)*vqm(crush) ./ (c + vqm(crush));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data2] = my_interp(times, data1);
% perform upsampling (1 sample per second to 2 samples per second)
% put extra sample at the beginning (a repeat of the first sample) rather
% than at the end.

if times ~= 2,
    error('Only implemented for 2');
end

hold = length(data1);
data2(1) = data1(1);
data2(2:2:hold*2,1) = data1;
data2(3:2:hold*2,1) = ( data1(2:hold) + data1(1:hold-1) ) / 2;
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hv_loss_par, hv_gain_par,si_loss_par, si_gain_par, ...
    color_comb_par, noise_par, error_par] = calc_parameters (...
    si_proc, hv_feat_proc, y_proc, cb_proc, cr_proc, ati_proc, ...
    si_orig, hv_feat_orig, y_orig, cb_orig, cr_orig, ati_orig, ...
    fps, part_si_min, part_si_max, part_hv_min, part_hv_max, ...
    part_c_min, part_c_max, part_c, part_ati_min, part_ati_max, ...
    part_ati, code_ati, TIME_DELTA);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model calculation.  Change from collapse all frames, to running-collapse
% of TIME_DELTA seconds at once.  Remove multiple-clip functionality.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ATI_SEC = 0.2; % width of ATI time-difference
ATI_WIDTH  = tslice_conversion (ATI_SEC, fps);

%%%%%%%%%%%%%%%%%%%%%
% calculate HV pars %
%%%%%%%%%%%%%%%%%%%%%

% Zero parameter values whose hv_feat_orig < part_hv(1).
hvgain_low = find(hv_feat_orig < part_hv_min);
hvloss_low = find(hv_feat_orig < 0.435);  % 0.435 > part_hv(1), better threshold for hv_loss
% Zero parameter values whose hv_feat_orig > part_hv(code_hv_size-1).
hvloss_high = find(hv_feat_orig > part_hv_max);
hvgain_high = find(hv_feat_orig > 1.90);  % 1.90 < hv_part(code_hv_size-1), better threshold for hv_gain

% Implement the hv_loss block parameter
hv_loss = (hv_feat_proc-hv_feat_orig)./hv_feat_orig;  % ratio loss calculation
hv_loss = min(hv_loss,0);  % negative part calculation
hv_loss(hvloss_low) = 0.0;
hv_loss(hvloss_high) = 0.0;

% Implement the hv_gain block parameter
hv_gain = log10(hv_feat_proc./hv_feat_orig);  % log10 gain calculation
hv_gain = max(hv_gain,0);  % positive part calculation
hv_gain(hvgain_low) = 0.0;
hv_gain(hvgain_high) = 0.0;

% If low spatial detail, reduce weight:
% Linear weight reduction, weight goes from weight_low to 1 as spatial goes
% from a to b
[frows,fcols,ftime] = size(si_orig);
si_weight = ones(frows, fcols, ftime);
weight_low = 0;
a = 5;
b = 25;
weight_low_slope = (1-weight_low)/(b-a);
si_weight(find(si_orig < a)) = weight_low;
low = find(si_orig >= a & si_orig < b);
si_weight(low) = weight_low_slope*(si_orig(low)-a) + weight_low;
hv_loss = hv_loss.*si_weight;
clear si_weight;
    
    % Pick off the valid region for this clip, using same rules as
    % read_feature.

    hv_loss_par_clip = squeeze(hv_loss(:,:,:));
    hv_gain_par_clip = squeeze(hv_gain(:,:,:));
     
    % y weighting function by macroblocks for less bandwidth
    weight_high = 0;
    c = 175;
    d = 255;
    weight_high_slope = (weight_high-1)/(d-c);
    y_orig_clip = squeeze(y_orig(:,:,:));
    [pr, pc, pt] = size(y_orig_clip);  
    y_weight = ones(pr, pc, pt);

    % Quantize y_orig average macroblock values which will be RR info
    high = find(y_orig_clip > c & y_orig_clip <= d);
    y_weight(high) = weight_high_slope*(y_orig_clip(high)-c) + 1;
    y_weight(find(y_orig_clip > d)) = weight_high;
    
    % MB(3,3,2)below1%_STminkowski(1,1.5) for hv_loss
    hv_loss_par_clip = hv_loss_par_clip.*y_weight;
    
    hv_loss_par_clip = st_collapse('below1%', hv_loss_par_clip, 'OverlapMacroBlock', 3, 3, 2);
    hv_loss_par_clip = running_collapse ('minkowski(1,1.5)', hv_loss_par_clip, TIME_DELTA-1, '3D');
    hv_loss_par = hv_loss_par_clip;
    
    % MB(3,3,2)above99%_STminkowski(1.5,3) for hv_gain
    hv_gain_par_clip = hv_gain_par_clip.*y_weight;
    
    % *** Added:  Clip the low end and subtract this clipping level
    hv_gain_par_clip = max(hv_gain_par_clip, 0.06) - 0.06;
    
    hv_gain_par_clip = st_collapse('above99%', hv_gain_par_clip, 'OverlapMacroBlock', 3, 3, 2);
    hv_gain_par_clip = running_collapse ('minkowski(1.5,3)', hv_gain_par_clip, TIME_DELTA-1, '3D');
    hv_gain_par = hv_gain_par_clip;
    
clear hv_loss hv_gain y_orig;

% *** Changed clipping level from 0.06 to 0.08
% Clip the low end and subtract this clipping level
hv_loss_clip = 0.08;
hv_loss_par = max(hv_loss_par, hv_loss_clip) - hv_loss_clip;

% *** Deleted, clipping now done prior to S-T collapse:
% hv_gain_clip = 0.06;
% hv_gain_par = max(hv_gain_par, hv_gain_clip) - hv_gain_clip;

% % Compress hv_gain values greater than training using crushing function.
% *** High clipping level changed from 0.74 to 0.75, which is approximately
% the new maximum from the training data.
hv_gain_max = 0.75;  % from training data
crush_hv = find(hv_gain_par > hv_gain_max);  % these values will be crushed
hv_gain_par(crush_hv) = (hv_gain_max + 0.25)*hv_gain_par(crush_hv) ./ (0.25 + hv_gain_par(crush_hv));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SI loss and gain parameters
%  Y_si13_32x32_1s_std_8_ratio_loss_MB(3,3,2)minkowski(1,2)_STminkowski(1.5,2.5)
%  Y_si13_32x32_1s_std_8_log_gain_above95%tail_minkowski(2,3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
si_t = part_si_min;  % si clipping threshold

% Processed features, Original features were already loaded for hv_loss
% si_proc = read_feature(feat_base_dir, si_name, sclips, 'processed', frows, fcols, ftime);

% Parameters outside of the si_orig quantizer range will be zeroed.
si_high = find(si_orig > part_si_max);

%  Clip si at low end
si_proc = max(si_proc,si_t);
si_orig = max(si_orig,si_t);

% Implement the si_loss parameter
si_loss = (si_proc-si_orig)./si_orig;  % ratio loss calculation
si_loss = min(si_loss,0);  % negative part calculation
si_loss(si_high) = 0.0;

% Implement the si_gain parameter
si_gain = log10(si_proc./si_orig);  % log gain calculation
si_gain = max(si_gain,0);  % positive part calculation
si_gain(si_high) = 0.0;

clear si_orig si_proc si_high;

    % Pick off the valid region for this clip, using same rules as
    % read_feature.
    si_loss_par_clip = squeeze(si_loss(:,:,:));
    si_gain_par_clip = squeeze(si_gain(:,:,:));
    
    si_loss_par_clip = si_loss_par_clip.*y_weight;
    
    si_loss_par_clip = st_collapse('minkowski(1,2)', si_loss_par_clip, 'OverlapMacroBlock', 3, 3, 2);
    si_loss_par_clip = running_collapse ('minkowski(1.5,2.5)', si_loss_par_clip, TIME_DELTA-1, '3D');
    si_loss_par = si_loss_par_clip;
    
    % above95%tail_minlowski(1.5,2) for si_gain
    [pr, pc, pt] = size(si_gain_par_clip);  % Find size after picking off valid
    
    % *** Added:  Clip the low end and subtract this clipping level
    si_gain_par_clip = max(si_gain_par_clip, 0.1) - 0.1;
    
    si_gain_par_clip = st_collapse('above95%tail', reshape(si_gain_par_clip, pr*pc,pt));
    si_gain_par_clip = running_collapse ('minkowski(1.5,2)', si_gain_par_clip, TIME_DELTA, '3D');
    si_gain_par = si_gain_par_clip;
    
clear si_gain si_loss;

% *** Changed clipping level from 0.1 to 0.12
% Clip the low end and subtract this clipping level
si_loss_clip = 0.12;
si_loss_par = max(si_loss_par, si_loss_clip) - si_loss_clip;

% *** Deleted: clipping now done prior to S-T collapse
% si_gain_clip = 0.03;
% si_gain_par = max(si_gain_par, si_gain_clip) - si_gain_clip;

% % Compress si_gain values greater than training using crushing function.
% *** High clipping level changed from 0.38 to 0.48, which is approximately
% the new maximum from the training data.
si_gain_max = 0.48;  % from training data
crush_si = find(si_gain_par > si_gain_max);  % these values will be crushed
si_gain_par(crush_si) = (si_gain_max + 0.25)*si_gain_par(crush_si) ./ (0.25 + si_gain_par(crush_si));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cb & Cr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters outside of the quantizer ranges will be zeroed, including the 
% bin that quantizes to zero.
cb_zero_high = find(cb_orig <= part_c_min | cb_orig >= part_c_max | cb_orig == 0);
cb_orig(cb_zero_high) = 0.0;
cb_proc(cb_zero_high) = 0.0;
clear cb_zero__high;
cr_zero_high = find(cr_orig <= part_c_min | cr_orig >= part_c_max | cr_orig == 0);
cr_orig(cr_zero_high) = 0.0;
cr_proc(cr_zero_high) = 0.0;
clear cb_zero_high cr_zero_high;

% Color extreme parameter with modifications to Euclidean distance
% Better Euclidean distance for new color extreme parameter
p = 1.0;  % Normal Euclidean uses p = 2.0, this absolute diff works better
r = 0.5;
w = 1.5;
color_euclid = (abs(cb_proc-cb_orig).^p + (w*abs(cr_proc-cr_orig)).^p).^r;
clear cb_proc cr_proc cb_orig cr_orig;

    % Pick off the valid region for this clip, using same rules as
    % read_feature.
    color_extreme_par_clip = squeeze(color_euclid(:,:,:));

    color_spread_par_clip = color_extreme_par_clip;
    
    % mb(3,3,2)above99%_STminkowski(0.5,1) for color_extreme
    color_extreme_par_clip = st_collapse('above99%', color_extreme_par_clip, 'OverlapMacroBlock',3,3,2);
    color_extreme_par_clip = running_collapse ('minkowski(0.5,1)', color_extreme_par_clip, TIME_DELTA-1, '3D');
    color_extreme_par = color_extreme_par_clip;
    
    % MB(3,3,2)minkowski(2,4)_ST90% for color_spread
    color_spread_par_clip = st_collapse('minkowski(2,4)', color_spread_par_clip, 'OverlapMacroBlock',3,3,2);
    color_spread_par_clip = running_collapse ('90%', color_spread_par_clip, TIME_DELTA-1, '3D');
    color_spread_par = color_spread_par_clip;
    
clear color_euclid;

% *** Old values multiplied by 1.75166641317702*1.00262368141757
% color_comb_par = -0.35186*color_spread_par + 0.39384*color_extreme_par;
color_comb_par = -0.617958*color_spread_par + 0.691686*color_extreme_par;

%  Clip the color comb parameter
%  *** Old value of 0.065 multiplied by 1.75166641317702*1.00262368141757
color_comb_clip = 0.114;  %  This will be clipped to prevent quality improvements
color_comb_par = max(color_comb_par, color_comb_clip) - color_comb_clip;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ati_nt = part_ati(24);  % ati clipping threshold for noise, approximately 5
ati_et = part_ati(57);  % ati clipping threshold for error, approximately 12

% Clip processed features if they exceed the maximum quantizer value
code_ati_size = size(code_ati,2);
ati_proc(find(ati_proc > part_ati(code_ati_size-1))) = code_ati(code_ati_size);

% Do the ATI parameters at each temporal alignment from plus to minus 0.4 seconds.
try_num = 0;
ati_search = floor(ceil(fps) * 0.4);
whole_length = min(length(ati_proc), length(ati_orig));
ati_proc = ati_proc((ati_search+1):(whole_length-ati_search));
whole_ati_orig = ati_orig;

for try_time = -ati_search:ati_search,
    try_num = try_num + 1;
    clear ati_nproc ati_norig;
    ati_orig = whole_ati_orig((ati_search+1+try_time):(whole_length-ati_search+try_time));
    
    %  Clip ati at low end for noise calculation
    ati_nproc = max(ati_proc,ati_nt);
    ati_norig = max(ati_orig,ati_nt);
    
    % Implement the ati gain parameter that will be used for noise
    ati_ngain = (ati_nproc-ati_norig)./ati_norig;  % ratio gain calculation
    ati_ngain = max(ati_ngain,0);  % positive part calculation

    % between25%50% for noise_par
    noise_par_clip = running_collapse ('between25%50%', ati_ngain, TIME_DELTA * ceil(fps) - ATI_WIDTH, '3D');
    
    % pick off samples, one every half second
    pattern = ceil( length(noise_par_clip):-ceil(fps)/2:1 );
    % reverse order & extend each end by one.
    pattern = [ pattern(length(pattern)) pattern(length(pattern):-1:1) pattern(1)];
    % error check the length of pattern.  This code will only be triggered
    % if some really weird frame rate is chosen.  Then, instead of adding
    % up to exactly 1sec, the discards might round to a different number.
    if ftime * 2 ~= length(pattern),
        while length(pattern) < ftime * 2,
            pattern = [ pattern(1) pattern ];
        end
        while length(pattern) > ftime * 2,
            pattern = pattern(2:length(pattern));
        end
    end
    
    noise_par(:,try_num) = noise_par_clip(pattern);
    
    % Pick off the valid region for this clip, using same rules as
    % read_feature.
    this_proc = ati_proc;
    this_orig = ati_orig;
    
    % 7-point Max filter
    this_proc = max_filterw(squeeze(this_proc),7);
    this_orig = max_filterw(squeeze(this_orig),7);
    
    %  Clip ati at low end for error calculation
    this_proc = max(this_proc,ati_et);
    this_orig = max(this_orig,ati_et);

    % Implement the ati gain parameter that will be used for errors
    ati_egain = (this_proc-this_orig)./this_orig;  % ratio gain calculation
    ati_egain = max(ati_egain,0);  % positive part calculation
  
    % above90% for error_par
    error_par_clip = running_collapse ('above90%', ati_egain, TIME_DELTA * ceil(fps) - ATI_WIDTH, '3D');
    error_par(:,try_num) = error_par_clip(pattern);
    
end

% take minimum at each point in time, of the temporal shifts (+- 2F)
noise_par = min(noise_par')';
error_par = min(error_par')';

clear ati_nproc ati_norig;
clear ati_ngain ati_proc ati_orig;
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataout] = macro_interp(datain);
%  Takes an input 3d array dimensioned as (r,c,t) and adds interpolated 
%  time samples halfway between.  So the output array will have dimension 
%  (r,c,2t-1).

[r,c,t] = size(datain);

dataout = zeros(r,c,2*t-1);

for i = 1:2*t-1
    if (floor(i/2) == i/2) % interpolate
        dataout(:,:,i) = (datain(:,:,i/2)+datain(:,:,1+i/2))/2;
    else % don't interpolate
        dataout(:,:,i) = datain(:,:,(i+1)/2);
    end
end


