function [one, two, three] = dll_model(control, varargin)
% DLL_MODEL
%   Complete VQM model calculations. 
%
% To Initialize:
%   [model_tslice_sec, model_planes] = dll_model('initialize', model_name, durration, fn);
%         'model_name' is the name of the model to be run: 
%           'Low'        Low-Bandwidth Model
%           'General'    NTIA General Model
%           'Developers' Developer's Model
%           'Fast'       Fast Low-Bandwidth Model
%         'fn' is 1 for original and 2 for processed -- either is okay -- where
%         function dll_video has been initialized on this computer for (fn).
%         'fn' presumed for following 'tslice' and 'get' calls, until next 'initialize'.
%         'durration' is the 
%
% To Calculate features for next time-slice 
%       where y (if 'model_planes' == 'y') 
%       or y, cb, cr, & fps (if 'model_planes' == 'ycbcr'):
%   [ready_for_vqm] = dll_model('tslice', y);
%   [ready_for_vqm] = dll_model('tslice', y, cb, cr, fps);
%
% To Get features
%   [features] = dll_model('get');
%
% To Complete VQM model calculations.
%   [vqm, pars, par_names] = dll_model('vqm', source_features, proc_features);
%       'source_features' is the 'features' return value from dll_features(fn=1)
%           for general & developer's models.  For Lowbw & Fast models,
%           'source_features' is the file name containing compressed
%           features.
%       'proc_features' is the 'features' return value from dll_features(fn=2)
%       Function 'dll_features' must already have been run & retreived with
%       dll_model('get') for fn=1 (source_features) and fn=2 (processed
%       featues).

  
persistent data;


one = NaN;
two = NaN;
three = NaN;

if strcmp(control,'initialize'),
    if strcmp(varargin{1},'Developers'),
        [data,one,two] = model_run_Callback_initialize_developers(varargin{2}, varargin{3});
    elseif strcmp(varargin{1}, 'General'),
        [data,one,two] = model_run_Callback_initialize_general(varargin{2}, varargin{3});
    elseif strcmp(varargin{1}, 'Low'),
        [data,one,two] = model_run_Callback_initialize_lowbw(varargin{2}, varargin{3});
    elseif strcmp(varargin{1}, 'Fast'),
        [data,one,two] = model_run_Callback_initialize_fastlowbw(varargin{2}, varargin{3});
    else
        error('model name not recognized');
    end
    data.model = varargin{1};
    
elseif strcmp(control,'tslice'),
    if strcmp(data.model, 'Developers'),
        [data, one] = model_run_Callback_feature_developers(data, varargin{1});
    elseif strcmp(data.model,'General'),
        [data, one] = model_run_Callback_feature_general(data, varargin{1}, varargin{2}, varargin{3});
    elseif strcmp(data.model,'Low'),
        [data, one] = model_run_Callback_feature_lowbw(data, varargin{1}, varargin{2}, varargin{3}, varargin{4});
    elseif strcmp(data.model,'Fast'),
        [data, one] = model_run_Callback_feature_fastlowbw(data, varargin{1}, varargin{2}, varargin{3}, varargin{4});
    end
    
elseif strcmp(control,'vqm'),
    if (strcmp(data.model,'Low') || strcmp(data.model,'Fast') ),
        file_name = varargin{1};
        [orig_features.si_orig, orig_features.part_si_min, orig_features.part_si_max, ...
             orig_features.hv_feat_orig, orig_features.part_hv_min, orig_features.part_hv_max, orig_features.y_orig, ...
             orig_features.cb_orig, orig_features.cr_orig, orig_features.part_c_min, orig_features.part_c_max, orig_features.part_c, ...
             orig_features.ati_orig, orig_features.part_ati_min, orig_features.part_ati_max, orig_features.part_ati, orig_features.code_ati ] ...
             = model_lowbw_compression('uncompress', file_name);
    else
        orig_features = varargin{1};
    end

    [fps] = dll_video('fps');  

    if strcmp(data.model, 'Developers'),
        [one, two, three] = model_run_Callback_vqm_developers(varargin{2}, orig_features);
    elseif strcmp(data.model,'General'),
        [one, two, three] = model_run_Callback_vqm_general(varargin{2}, orig_features);
    elseif strcmp(data.model,'Low'),
        [one, two, three] = model_run_Callback_vqm_lowbw(varargin{2}, orig_features, fps);
    elseif strcmp(data.model,'Fast'),
        [one, two, three] = model_run_Callback_vqm_fastlowbw(varargin{2}, orig_features, fps);
    end
elseif strcmp(control,'get'),
    one = data;
end


%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, model_tslice_sec, model_planes] = model_run_Callback_initialize_developers(durration, fn);

model_tslice_sec = 0.6;
model_planes = 'y';

data.tslice_total = floor(durration / model_tslice_sec);
data.tslices = 0;
data.avg_prev_tslice = 0;

% get default SROI
[temp.image_size.rows,temp.image_size.cols] = dll_video('size',fn);
[temp.cvr] = dll_calib_video('pvr');
hv_size = 8;
extra = 6;
[sroi,vert,horiz] = adjust_requested_sroi (temp, ...
    'vsize',hv_size, 'hsize',hv_size, 'extra',extra);
dll_calib_video('sroi', sroi, extra);

%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,model_tslice_sec, model_planes] = model_run_Callback_initialize_general(durration, fn);

model_tslice_sec = 0.2;
model_planes = 'ycbcr';

data.tslice_total = floor(durration / model_tslice_sec);
data.tslices = 0;

% get default SROI
[temp.image_size.rows,temp.image_size.cols] = dll_video('size',fn);
[temp.cvr] = dll_calib_video('pvr');
hv_size = 8;
extra = 6;
[sroi,vert,horiz] = adjust_requested_sroi (temp, ...
    'vsize',hv_size, 'hsize',hv_size, 'extra',extra);
dll_calib_video('sroi', sroi, extra);

%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,model_tslice_sec, model_planes] = model_run_Callback_initialize_lowbw(durration, fn);

% figure out side & control option.

data.destination = (fn == 2);

model_tslice_sec = 1.0;
model_planes = 'ycbcr';

data.tslice_total = floor(durration / model_tslice_sec);
data.tslices = 0;
if data.destination,
    model_lowbw_features_shift('clear');
else
    model_lowbw_features('clear');
end

% set lowbw SROI
[image_size.rows,image_size.cols] = dll_video('size',fn);
[filter_size, extra] = adaptive_filter(image_size);
[pvr] = dll_calib_video('pvr');
[valid, cvr, sroi] = model_lowbw_sroi(extra, pvr.top, pvr.left, pvr.bottom, pvr.right);
dll_calib_video('sroi', sroi, extra+1);


%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,model_tslice_sec, model_planes] = model_run_Callback_initialize_fastlowbw(durration, fn);

% figure out side & control option.

data.destination = (fn == 2);

model_tslice_sec = 1.0;
model_planes = 'ycbcr';

data.tslice_total = floor(durration / model_tslice_sec);
data.tslices = 0;
if data.destination,
    model_fastlowbw_features_shift('clear');
else
    model_fastlowbw_features('clear');
end

% set lowbw SROI
[image_size.rows,image_size.cols] = dll_video('size',fn);
[filter_size, extra] = adaptive_filter(image_size);
[pvr] = dll_calib_video('pvr');
[valid, cvr, sroi] = model_lowbw_sroi(extra, pvr.top, pvr.left, pvr.bottom, pvr.right);
dll_calib_video('sroi', sroi, extra+1);


%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, ready_for_vqm] = model_run_Callback_feature_developers(data, tslice);
% process next Time-slice

if data.tslices == data.tslice_total,
    data.tslices = 0;
    data.si_feature = [];
    data.hv_feature = [];
    data.hvb_feature = [];
    data.ati_feature = [];
end

% si&hv&ati, preaveraged.
avg_tslice = mean(tslice,3);
[si,hv,hvb] = filter_si_hv(avg_tslice);

% Compute ATI.  discard extra border.
if data.tslices > 0,
	[row,col] = size(avg_tslice);
	ati = abs(data.avg_tslice_prior - avg_tslice);
	ati=ati(7:row-6,7:col-6);
end
data.avg_tslice_prior = avg_tslice;

% shrink rows/cols if not power of 8.
[row,col] = size(si);
if floor(row/8) ~= ceil(row/8) | floor(col/8) ~= ceil(col/8),
    row = 8 * floor(row/8);
    col = 8 * floor(col/8);
    si=si(1:row,1:col);
    hv=hv(1:row,1:col);
    hvb=hvb(1:row,1:col);
    if data.tslices > 0,
        ati=ati(1:row,1:col);
    end
end

% take mean or standard deviation of block.
data.si_feature(:,:,data.tslices+1) = block_statistic(si,8,8,'std');
data.hv_feature(:,:,data.tslices+1) = block_statistic(hv,8,8,'mean');
data.hvb_feature(:,:,data.tslices+1) = block_statistic(hvb,8,8,'mean');
if data.tslices > 0,
    data.ati_feature(:,:,data.tslices) = block_statistic(ati,8,8,'std');
end

% update number of tslices destination.
data.tslices = data.tslices + 1;

if data.tslices == data.tslice_total,
    ready_for_vqm = 1;
else
    ready_for_vqm = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, ready_for_vqm] = model_run_Callback_feature_general(data, y, cb, cr);
% process next Time-slice

[rows,cols,time] = size(y);
if data.tslices == data.tslice_total,
    data.tslices = 0;
    data.si_feature = [];
    data.hv_feature = [];
    data.hvb_feature = [];
    data.ati_feature = [];
    data.cont_feature = [];
    data.cb_feature = [];
    data.cr_feature = [];
end

% si&hv perceptual filter.
[si,hv,hvb] = filter_si_hv(y);
pause(0.02);

% ATI perceptual filter.
if data.tslices == 0;
    ati = filter_ati(y);
else
    yplus = cat(3, data.prior_frame, y);
    [ati] = filter_ati(yplus);
end
data.prior_frame = y(:,:,time);
pause(0.02);

% shrink ATI & Y & Cb & Cr same amount as as SI/HV/HVbar
ati=ati(7:rows-6,7:cols-6,:);
y=y(7:rows-6,7:cols-6,:);
cb=cb(7:rows-6,7:cols-6,:);
cr=cr(7:rows-6,7:cols-6,:);

% shrink rows/cols if not power of 8.
[rows,cols,time] = size(y);
if floor(rows/8) ~= ceil(rows/8) | floor(cols/8) ~= ceil(cols/8),
    rows = 8 * floor(rows/8);
    cols = 8 * floor(cols/8);
    si=si(1:rows,1:cols,:);
    hv=hv(1:rows,1:cols,:);
    hvb=hvb(1:rows,1:cols,:);
    ati=ati(1:rows,1:cols,:);
    y=y(1:rows,1:cols,:);
    cb=cb(1:rows,1:cols,:);
    cr=cr(1:rows,1:cols,:);
    pause(0.02);
end

% take mean or standard deviation of blocks.
data.si_feature(:,:,data.tslices+1) = block_statistic(si,8,8,'std');
data.hv_feature(:,:,data.tslices+1) = block_statistic(hv,8,8,'mean');
data.hvb_feature(:,:,data.tslices+1) = block_statistic(hvb,8,8,'mean');
data.ati_feature(:,:,data.tslices+1) = block_statistic(ati,4,4,'std');
data.cont_feature(:,:,data.tslices+1) = block_statistic(y,4,4,'std');
for cnt=1:time,
    data.cb_feature(:,:,data.tslices*time+cnt) = block_statistic(cb(:,:,cnt),8,8,'mean');
    data.cr_feature(:,:,data.tslices*time+cnt) = block_statistic(cr(:,:,cnt),8,8,'mean');
end

% update number of tslices destination.
data.tslices = data.tslices + 1;

if data.tslices == data.tslice_total,
    ready_for_vqm = 1;
else
    ready_for_vqm = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, ready_for_vqm] = model_run_Callback_feature_lowbw(data, y, cb, cr, fps);
% process next Time-slice

if data.tslices == data.tslice_total,
    data.tslices = 0;
end


% cut out SROI +/- 6 pixels
[rows,cols,time] = size(y);
image_size.rows = rows;
image_size.cols = cols;
[filter_size,extra] = adaptive_filter(image_size);
[valid, pvr, sroi] = model_lowbw_sroi(extra, 1, 1, rows, cols);
if ~valid,
    report_Callback('add', 'Valid Region too small.  Low Bandwidth Model cannot execute.');
    stop_Callback('button');
    return;
end

[image_size.rows,image_size.cols,junk] = size(y);
y = y(pvr.top:pvr.bottom, pvr.left:pvr.right,:);
cb = cb(pvr.top:pvr.bottom, pvr.left:pvr.right,:);
cr = cr(pvr.top:pvr.bottom, pvr.left:pvr.right,:);

% compute features. 
[filter_size, extra] = adaptive_filter(image_size);
if data.destination,
    model_lowbw_features_shift ('memory',  y, cb, cr, fps, filter_size, extra);
else
    % discard one pixel on all sides, then compute features
    [row,col,time] = size(y);
    y = y(2:row-1, 2:col-1, :);
    cb = cb(2:row-1, 2:col-1, :);
    cr = cr(2:row-1, 2:col-1, :);
    model_lowbw_features ('memory',  y, cb, cr, fps, filter_size, extra);
end


% update number of tslices destination.
data.tslices = data.tslices + 1;

if data.tslices == data.tslice_total,
    if data.destination,
        [data.data] = model_lowbw_features_shift ('eof');
    else
        [data.si_std data.hv_ratio data.y_mean data.cb_mean data.cr_mean data.ati_rms] ...
            = model_lowbw_features ('eof');
    end
    ready_for_vqm = 1;
else
    ready_for_vqm = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, ready_for_vqm] = model_run_Callback_feature_fastlowbw(data, y, cb, cr, fps);
% process next Time-slice

if data.tslices == data.tslice_total,
    data.tslices = 0;
end


% cut out SROI +/- 6 pixels
[rows,cols,time] = size(y);
image_size.rows = rows;
image_size.cols = cols;
[filter_size,extra] = adaptive_filter(image_size);
[valid, pvr, sroi] = model_lowbw_sroi(extra, 1, 1, rows, cols);
if ~valid,
    report_Callback('add', 'Valid Region too small.  Low Bandwidth Model cannot execute.');
    stop_Callback('button');
    return;
end

[image_size.rows,image_size.cols,junk] = size(y);
y = y(pvr.top:pvr.bottom, pvr.left:pvr.right,:);
cb = cb(pvr.top:pvr.bottom, pvr.left:pvr.right,:);
cr = cr(pvr.top:pvr.bottom, pvr.left:pvr.right,:);

% compute features. 
[filter_size, extra] = adaptive_filter(image_size);
if data.destination,
    model_fastlowbw_features_shift ('memory',  y, cb, cr, fps, filter_size, extra);
else
    % discard one pixel on all sides, then compute features
    [row,col,time] = size(y);
    y = y(2:row-1, 2:col-1, :);
    cb = cb(2:row-1, 2:col-1, :);
    cr = cr(2:row-1, 2:col-1, :);
    model_fastlowbw_features ('memory',  y, cb, cr, fps, filter_size, extra);
end


% update number of tslices destination.
data.tslices = data.tslices + 1;

if data.tslices == data.tslice_total,
    if data.destination,
        [data.data] = model_fastlowbw_features_shift ('eof');
    else
        [data.si_std data.hv_ratio data.y_mean data.cb_mean data.cr_mean data.ati_rms] ...
            = model_fastlowbw_features ('eof');
    end
    ready_for_vqm = 1;
else
    ready_for_vqm = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%
function [vqm_value, pars, par_names] = model_run_Callback_vqm_developers(proc, src);


si_loss = compare_feature(src.si_feature, proc.si_feature, 'ratio_loss', 'MinThreshold', 6);
[r,c,t] = size(si_loss);
si_loss = reshape(si_loss,r*c,t);
si_loss = 0.03 + min( -0.03, st_collapse( 'mean', st_collapse('below5%', si_loss)));
si_loss = -0.6289 * si_loss;

hv_loss = compare_dual_feature(src.hv_feature, src.hvb_feature, proc.hv_feature, proc.hvb_feature, 'divide', 'MinThreshold', 3, 'compare', 'ratio_loss');
[r,c,t] = size(hv_loss);
hv_loss = reshape(hv_loss,r*c,t);
hv_loss = -0.06 + max( 0.06, ( st_collapse( '10%', st_collapse( 'below5%', hv_loss))).^2);
hv_loss = 0.2305 * hv_loss;

hv_gain = compare_dual_feature(src.hv_feature, src.hvb_feature, proc.hv_feature, proc.hvb_feature, 'divide', 'MinThreshold', 3, 'compare', 'log_gain');
[r,c,t] = size(hv_gain);
hv_gain = reshape(hv_gain,r*c,t);
hv_gain = st_collapse( 'mean', st_collapse( 'above95%', hv_gain));
hv_gain = 0.1551 * hv_gain;

ati_gain = compare_feature(src.ati_feature, proc.ati_feature, 'log_gain', 'MinThreshold', 1);
[r,c,t] = size(ati_gain);
ati_gain = reshape(ati_gain,r*c,t);
ati_gain = st_collapse( '10%', st_collapse( 'mean', ati_gain));
ati_gain = 1.0587 * ati_gain;

ati_loss = compare_feature(src.ati_feature,proc.ati_feature, 'ratio_loss', 'MinThreshold', 3);
[r,c,t] = size(ati_loss);
ati_loss = reshape(ati_loss,r*c,t);
ati_loss = st_collapse( '10%', st_collapse( 'below5%', ati_loss));
ati_loss = -0.1444 * ati_loss;

% compute model
pars = [si_loss hv_loss hv_gain ati_gain ati_loss];
vqm_value = sum(pars);

% Clip VQM for values less than zero
if vqm_value < 0,
	vqm_value = 0.0;  % No quality improvements allowed
end
% Compress VQM values greater than 1 using standard crushing function
c = 0.5;
if vqm_value > 1,
    vqm_value = (1 + c)*vqm_value ./ (c + vqm_value);
end

par_names = { 'si_loss' 'hv_loss' 'hv_gain' 'ati_gain' 'ati_loss' };

%%%%%%%%%%%%%%%%%%%%%%%%%
function [vqm_value, pars, par_names] = model_run_Callback_vqm_general(proc, src);

% SI loss
si_loss = compare_feature(src.si_feature, proc.si_feature, 'ratio_loss', 'MinThreshold', 12);
[r,c,t] = size(si_loss);
si_loss = reshape(si_loss,r*c,t);
si_loss = st_collapse( '10%', st_collapse('below5%', si_loss));
si_loss = -0.2097 * si_loss;

% HV loss
hv_loss = compare_dual_feature(src.hv_feature, src.hvb_feature, proc.hv_feature, proc.hvb_feature, ...
    'divide', 'MinThreshold', 3, 'compare', 'ratio_loss');
[r,c,t] = size(hv_loss);
hv_loss = reshape(hv_loss,r*c,t);
hv_loss = -0.06 + max( 0.06, ( st_collapse( 'mean', st_collapse( 'below5%', hv_loss))).^2);
hv_loss = 0.5969 * hv_loss;

% HV gain
hv_gain = compare_dual_feature(src.hv_feature, src.hvb_feature, proc.hv_feature, proc.hvb_feature, ...
    'divide', 'MinThreshold', 3, 'compare', 'log_gain');
[r,c,t] = size(hv_gain);
hv_gain = reshape(hv_gain,r*c,t);
hv_gain = st_collapse( 'mean', st_collapse( 'above95%', hv_gain));
hv_gain = 0.2483 * hv_gain;

% first color metric
color_compare = compare_dual_feature(src.cb_feature, src.cr_feature, proc.cb_feature, proc.cr_feature, ...
    'euclid', 'Weight', 1.5);
[r,c,t] = size(color_compare);
color_compare = reshape(color_compare,r*c,t);
color1 = -0.6 + max(0.6, st_collapse( '10%', st_collapse( 'std', color_compare)));
color1 = 0.0192 * color1;

% SI gain
si_gain = compare_feature(src.si_feature, proc.si_feature, 'log_gain', 'MinThreshold', 8);
[r,c,t] = size(si_gain);
si_gain = reshape(si_gain,r*c,t);
si_gain = min(0.14, -0.004 + max(0.004, st_collapse( 'mean', st_collapse('mean', si_gain))));
si_gain = -2.3416 * si_gain;

% contrast * ati
contati = compare_dual_feature(src.cont_feature, src.ati_feature, proc.cont_feature, proc.ati_feature, ...
    'multiply', 'MinThreshold', 3, 'compare', 'ratio_gain', 'LogBefore');
[r,c,t] = size(contati);
contati = reshape(contati,r*c,t);
contati = st_collapse( '10%', st_collapse('mean', contati));
contati = 0.0431 * contati;

% second color metric
color2 = st_collapse( 'std', st_collapse( 'above99%tail', color_compare));
color2 = 0.0076 * color2;

%
pars = [si_loss hv_loss hv_gain color1 si_gain contati color2];
vqm_value = sum(pars);

% Clip VQM for values less than zero
if vqm_value < 0,
	vqm_value = 0.0;  % No quality improvements allowed
end
% Compress VQM values greater than 1 using standard crushing function
c = 0.5;
if vqm_value > 1,
    vqm_value = (1 + c)*vqm_value ./ (c + vqm_value);
end

par_names = { 'si_loss' 'hv_loss' 'hv_gain' 'color1' 'si_gain' 'contati' 'color2' };


%%%%%%%%%%%%%%%%%%%%%%%%%
function [vqm_value, pars, par_names] = model_run_Callback_vqm_lowbw(proc, src, fps);

    [row,col,TIME_DELTA] = size(src.si_orig);
    for loop = 1:9,
        % calculate model
        do_test_print = 0;
        [data(loop).vqm, data(loop).hv_loss_par, data(loop).hv_gain_par,data(loop).si_loss_par, data(loop).si_gain_par, ...
            data(loop).color_comb_par, data(loop).noise_par, data(loop).error_par] = ...
            model_lowbw_parameters (...
                proc.data(loop).si_std, proc.data(loop).hv_ratio, proc.data(loop).y_mean, proc.data(loop).cb_mean, proc.data(loop).cr_mean, proc.data(loop).ati_rms, ...
                src.si_orig, src.hv_feat_orig, src.y_orig, src.cb_orig, src.cr_orig, src.ati_orig, ...
                fps, src.part_si_min, src.part_si_max, src.part_hv_min, src.part_hv_max, ...
                src.part_c_min, src.part_c_max, src.part_c, src.part_ati_min, src.part_ati_max, ...
                src.part_ati, src.code_ati, 0, TIME_DELTA);
    end

    % select smallest average VQM score shift
    for shift=1:9,
        vqm_mean(shift) = mean(data(shift).vqm);
    end
    [junk shift] = min(vqm_mean);
    row = floor((shift-1)/3)-1;
    col = mod(shift-1,3)-1;
    
    % keep only the last sample; above function produces an entire
    % time-history and we only want the ending number.
    num = length(data(shift).vqm);
    
    % select & copy data to return.
    pars = [data(shift).hv_loss_par(num), data(shift).hv_gain_par(num), ...
                data(shift).si_loss_par(num), data(shift).si_gain_par(num), ...
                data(shift).color_comb_par(num), data(shift).noise_par(num), data(shift).error_par(num), ...
                row, col];

    par_names = {'hv_loss' 'hv_gain' 'si_loss' 'si_gain' 'color_comb' 'noise' 'error' 'vshift' 'hshift'};

    vqm_value = data(shift).vqm(num);

    
%%%%%%%%%%%%%%%%%%%%%%%%%
function [vqm_value, pars, par_names] = model_run_Callback_vqm_fastlowbw(proc, src, fps);

    [row,col,TIME_DELTA] = size(src.si_orig);
    for loop = 1:9,
        % calculate model
        do_test_print = 0;
        [data(loop).vqm, data(loop).hv_loss_par, data(loop).hv_gain_par,data(loop).si_loss_par, data(loop).si_gain_par, ...
            data(loop).color_comb_par, data(loop).noise_par, data(loop).error_par] = ...
            model_fastlowbw_parameters (...
                proc.data(loop).si_std, proc.data(loop).hv_ratio, proc.data(loop).y_mean, proc.data(loop).cb_mean, proc.data(loop).cr_mean, proc.data(loop).ati_rms, ...
                src.si_orig, src.hv_feat_orig, src.y_orig, src.cb_orig, src.cr_orig, src.ati_orig, ...
                fps, src.part_si_min, src.part_si_max, src.part_hv_min, src.part_hv_max, ...
                src.part_c_min, src.part_c_max, src.part_c, src.part_ati_min, src.part_ati_max, ...
                src.part_ati, src.code_ati, 0, TIME_DELTA);
    end

    % select smallest average VQM score shift
    for shift=1:9,
        vqm_mean(shift) = mean(data(shift).vqm);
    end
    [junk shift] = min(vqm_mean);
    row = floor((shift-1)/3)-1;
    col = mod(shift-1,3)-1;
    
    % keep only the last sample; above function produces an entire
    % time-history and we only want the ending number.
    num = length(data(shift).vqm);
    
    % select & copy data to return.
    pars = [data(shift).hv_loss_par(num), data(shift).hv_gain_par(num), ...
                data(shift).si_loss_par(num), data(shift).si_gain_par(num), ...
                data(shift).color_comb_par(num), data(shift).noise_par(num), data(shift).error_par(num), ...
                row, col];

    par_names = {'hv_loss' 'hv_gain' 'si_loss' 'si_gain' 'color_comb' 'noise' 'error' 'vshift' 'hshift'};

    vqm_value = data(shift).vqm(num);

