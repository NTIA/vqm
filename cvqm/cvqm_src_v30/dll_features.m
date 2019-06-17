function [features] = dll_features (model_name, fn, durration, compressed_file);
% DLL_FEATURES
%   Calculate features for a model.
% SYNTAX
%   [features] = dll_features(model_name, fn, durration);
%   [features] = dll_features(model_name, fn, durration, compressed_features_file);
% DESCRIPTION
%   Calculate original or processed features needed to calculate one model.
%   Function 'dll_video' must be initialized for (fn).  'model_name' is the
%   name of the model to be run: 
%       'Low'        Low-Bandwidth Model
%       'General'    NTIA General Model
%       'Developers' Developer's Model
%   'fn' is 1 for original and 2 for processed.  
%   'durration' is the durration of the video sequence for which the
%   features are to be calculated, in seconds (from 5 to 30 seconds)
%   'compressed_features_file' is the name of a file where the compressed
%   features should be written.  Currently, this option is only available
%   for fn=1 and model_name='Low' (low-bandwidth model, original features).
%
%   Return variable 'features' is a structure holding the uncompressed
%   feature data. 
%

warning off MATLAB:max:mixedSingleDoubleInputs

% initialize model.
[model_tslice_sec, model_planes] = dll_model('initialize', model_name, durration, fn);
dll_video('set_tslice', fn, model_tslice_sec);

% run tslices through features
if strcmp(model_planes, 'y'),
    ready_for_vqm = 0;
    while ~ready_for_vqm,
        [y] = dll_calib_video('tslice', fn);
        [ready_for_vqm] = dll_model('tslice', y);
    end
elseif strcmp(model_planes, 'ycbcr'),
    [fps] = dll_video('fps', fn);  
    ready_for_vqm = 0;
    while ~ready_for_vqm,
        [y, cb, cr] = dll_calib_video('tslice', fn);
        [ready_for_vqm] = dll_model('tslice', y, cb, cr, fps);
    end
end

% retrieve the features.
[features] = dll_model('get');

if (strcmp(model_name,'Low') | strcmp(model_name,'Fast')) & fn==1 & exist('compressed_file','var'),
    model_lowbw_compression ('compress', compressed_file, features.si_std, features.hv_ratio, ...
        features.y_mean, features.cb_mean, features.cr_mean, features.ati_rms );
end
