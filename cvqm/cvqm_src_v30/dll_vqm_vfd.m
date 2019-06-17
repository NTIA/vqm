function nn_model = dll_vqm_vfd(spatial_shift, luminance, scale, viewing_distance, delay)
% DLL_VQM_VFD
%  This function computes the Video Quality Metric with Variable Frame
%  Delay (VQM_VFD) for an original and processed clip in CVQM.
%
%  The function performs the following:
%  1. Creates mock GClips structures.
%  2. Perform VFD feature extraction.
%  3. Perform VFD parameter calculation.
%  4. Perform VFD VQM calculation (which uses the Neural Network Toolbox).

% Populate mock Gclips structures.
[rows, cols, fps] = dll_video('size', 1);
clipset(1).image_size.rows = rows;
clipset(1).image_size.cols = cols;
clipset(1).fps = fps;
clipset(1).loc_start = 1;
clipset(1).align_start = 1;
clipset(1).loc_stop = dll_video('total_frames', 1);
clipset(1).align_stop = dll_video('total_frames', 1);
clipset(1).video_standard = dll_video('get_video_standard',1);
clipset(1).cvr = dll_calib_video('pvr');
clipset(1).spatial.vertical = 0;
clipset(1).spatial.horizontal = 0;
clipset(1).luminance_gain = 1;
clipset(1).luminance_offset = 0;
clipset(1).scale.horizontal = 1000;
clipset(1).scale.vertical = 1000;
clipset(1).hrc = 'original';

[rows, cols, fps] = dll_video('size', 2);
clipset(2).image_size.rows = rows;
clipset(2).image_size.cols = cols;
clipset(2).fps = fps;
clipset(2).loc_start = 1;
clipset(2).align_start = 1;
clipset(2).loc_stop = dll_video('total_frames', 2);
clipset(2).align_stop = dll_video('total_frames', 2);
clipset(2).video_standard = dll_video('get_video_standard',1);
clipset(2).cvr = dll_calib_video('pvr');
clipset(2).spatial.vertical = spatial_shift.vertical;
clipset(2).spatial.horizontal = spatial_shift.horizontal;
clipset(2).luminance_gain = luminance.gain;
clipset(2).luminance_offset = luminance.offset;
clipset(2).scale.horizontal = scale.horizontal;
clipset(2).scale.vertical = scale.vertical;
clipset(2).hrc = 'processed';

if delay < 0  % added 8/5/11
    clipset(2).align_stop = clipset(2).align_stop + delay;
end


%  Define the spatial degrees and time extent of the ST blocks
deg_size = 0.4;  % in angular degrees, a square block spatially
time_size = 0.2;  % in seconds

% Extract SI and HV features
[datasi, datasio, datahv, datahvo, datahvb, datahvbo] = ...
    dll_vfd_feature_loop_si_hv_adapt(clipset, deg_size, time_size, viewing_distance);

% Extract Y (CONT) features
[data_cont, data_conto, data_contm, data_contmo] = dll_vfd_feature_loop_cont(clipset, deg_size, time_size, viewing_distance);

% Extract TI features
[data_yrms, data_yrmso, data_crbrms, data_cbrmso, data_crrms, data_crrmso, ...
    data_ymean, data_ymeano, data_cbmean, data_cbmeano, data_crmean, data_crmeano] = ...
    dll_vfd_feature_loop_ti(clipset, deg_size, time_size, viewing_distance);

% Extract Mean Square Error (MSE) features
[data_mse, data_mseo, data_msem, data_msemo] = dll_vfd_feature_loop_mse(clipset, deg_size, time_size, viewing_distance);

% HV loss calculation
hv_loss = dll_vfd_clippar_loop_exp_hv_loss(datahv, datahvo, datahvb, datahvbo, data_yrms, data_contm);

% HV gain calculation
hv_gain = dll_vfd_parameter_dual_loop(datahv, datahvo, datahvb, datahvbo, ...
    {'rms'}, {'rms'}, 'MinThreshold', 3, 'divide', 'compare', 'log_gain');

% SI loss calculation
si_loss = dll_vfd_parameter_loop(datasi, datasio, 'ratio_loss', {'mean'}, {'above90%'}, 'MinThreshold', 12);

% SI gain calculation
si_gain = dll_vfd_parameter_loop(datasi, datasio, 'log_gain', {'above98%tail'}, {'rms'}, 'MinThreshold', 8);

% TI gain calculation
ti_gain = dll_vfd_parameter_loop(data_yrms, data_yrmso, 'log_gain', {'above95%tail'}, {'above95%tail'}, 'MinThreshold', 3, '3D');  % Use 3D collapsing

% RMSE calculation
vfd_rmse = dll_vfd_parameter_loop(data_mse, data_mseo, 'minus_gain', {'mean'}, {'mean'}, '3D');  % Use 3D collapsing

% vfd_par1 calculation
vfd_par1 = dll_vfd_clippar_loop_par1(clipset);

% vfd_xpar calculation (vfd_par1*psnr_vfd)
psnr_vfd = dll_vfd_clippar_loop_psnr(data_mse);
vfd_xpar = vfd_par1.*psnr_vfd;

nn_obj = [hv_loss;hv_gain;si_loss;si_gain;ti_gain;vfd_rmse;vfd_par1;vfd_xpar];

load nn_8par;
nn_model = nn_8par(nn_obj);


