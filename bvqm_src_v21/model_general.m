function [pars, status] = model_general(test_structs, clip_structs, feature_base_dir, varargin);
% MODEL_GENERAL
%   Compute General model on a clip structure.
% SYNTAX
%     [pars] = model_general(test_structs, clip_structs, feature_base_dir);
% DESCRIPTION
%   Calculate the General model on the specified clips.   
%       'test_structs' describes the video tests, formatted like GTests
%       'clip_structs' describges the video clips, formatted like GClips
%       'feature_base_dir' is used to write / store features.
%   Returns the following:
%       'pars', a parameter structure with the general model & parameters
%       'status' is 0 if operated correctly; 1 if an error was encountered
%
%  The following optional parameters are also available.  
%   'verbose'   Print progress report to screen.
%   'quiet'     Don't print anything to the screen (default).
%       
%  Note: Code used to compute from features to general model is idetical to
%  that used by IVQM.  This code has been validated using VQM_pc.

status = 0;
verbose = 0;
for cnt = 1:length(varargin),
    if strcmp(lower(varargin{cnt}),'verbose'),
        verbose = 1;
    elseif strcmp(lower(varargin{cnt}),'quiet'),
        verbose = 0;
    else
        error('optional property not recognized');
    end
end

try

    if verbose,
        verbose_string = 'verbose';
    else
        verbose_string = 'quiet';
    end

    % calculate features
    feature_loop_si_hv(test_structs,clip_structs, 8,8,0.2,feature_base_dir,verbose_string);
    feature_loop_coherent_color(test_structs,clip_structs, 8,8,1/30,feature_base_dir,verbose_string);
    feature_loop_cont(test_structs,clip_structs, 4,4,0.2,feature_base_dir,verbose_string);
    feature_loop_ati(test_structs,clip_structs, 4,4,0.2,feature_base_dir,verbose_string);

    feature_base_dir = [feature_base_dir '/'];

    % figure out the list of clip names.  Fill in the names when computing
    % metrics.
    %pars.clip_name = cell(num_clips);
    ccnt = 1;

    % Loop through all clips, sorted by scene
    offsets = sort_clips_by('scene',clip_structs, test_structs);
    for loop = 1:length(offsets),
        % pick off the original for this scene.  Skip if no original defined.
        curr_offsets = offsets{loop};
        if ~strcmp(clip_structs(curr_offsets(1)).hrc{1},'original'),
            if verbose,
                fprintf('Warning:  Original missing for %s:%s.  Skipping these clips\n', ...
                    clip_structs(curr_offsets(1)).test{1}, clip_structs(curr_offsets(1)).scene{1});
            end
            continue;
        end

        % load original features. 
        clip_name = sprintf('%s_%s_%s.mat', ...
            clip_structs(curr_offsets(1)).test{1}, ...
            clip_structs(curr_offsets(1)).scene{1}, ...
            clip_structs(curr_offsets(1)).hrc{1});

        load( [feature_base_dir 'feature_Y_si13_8x8_0.2s_std/' clip_name ] );
        src.si_feature = data;
        load( [feature_base_dir 'feature_Y_hv13_8x8_0.2s_mean/' clip_name] );
        src.hv_feature = data;
        load( [feature_base_dir 'feature_Y_hvbar13_8x8_0.2s_mean/' clip_name] );
        src.hvb_feature = data;
        load( [feature_base_dir 'feature_Y_cont_4x4_0.2s_std/' clip_name] );
        src.cont_feature = data;
        load( [feature_base_dir 'feature_Y_ati_4x4_0.2s_std/' clip_name] );
        src.ati_feature = data;
        load( [feature_base_dir 'feature_color_coher_cb_8x8_0.03333s_mean/' clip_name] );
        src.cb_feature = data;
        load( [feature_base_dir 'feature_color_coher_cr_8x8_0.03333s_mean/' clip_name] );
        src.cr_feature = data;

        % loop through all processed versions of this clip.
        for cnt = 2:length(curr_offsets),
            % load processed features. 
            clip_name = sprintf('%s_%s_%s.mat', ...
                clip_structs(curr_offsets(cnt)).test{1}, ...
                clip_structs(curr_offsets(cnt)).scene{1}, ...
                clip_structs(curr_offsets(cnt)).hrc{1});
            load( [feature_base_dir 'feature_Y_si13_8x8_0.2s_std/' clip_name ] );
            proc.si_feature = data;
            load( [feature_base_dir 'feature_Y_hv13_8x8_0.2s_mean/' clip_name] );
            proc.hv_feature = data;
            load( [feature_base_dir 'feature_Y_hvbar13_8x8_0.2s_mean/' clip_name] );
            proc.hvb_feature = data;
            load( [feature_base_dir 'feature_Y_cont_4x4_0.2s_std/' clip_name] );
            proc.cont_feature = data;
            load( [feature_base_dir 'feature_Y_ati_4x4_0.2s_std/' clip_name] );
            proc.ati_feature = data;
            load( [feature_base_dir 'feature_color_coher_cb_8x8_0.03333s_mean/' clip_name] );
            proc.cb_feature = data;
            load( [feature_base_dir 'feature_color_coher_cr_8x8_0.03333s_mean/' clip_name] );
            proc.cr_feature = data;

            % calculate model
            do_test_print = 0;
            [vqm, clip_pars, par_names] = model_run_Callback_vqm_general(proc, src);

            % fill in the name of this clip
            pars.clip_name{ccnt} = sprintf('%s_%s_%s', clip_structs(curr_offsets(cnt)).test{1}, ...
                    clip_structs(curr_offsets(cnt)).scene{1}, clip_structs(curr_offsets(cnt)).hrc{1});
            pars.inlsa_mos(ccnt) = clip_structs(curr_offsets(cnt)).inlsa_mos;
            pars.mos(ccnt) = clip_structs(curr_offsets(cnt)).mos;

            num = length(vqm);
            pars.data(:,ccnt) = [vqm clip_pars];
            ccnt = ccnt + 1;
        end
    end

    % record name of each parameter
    pars.par_name{1} = 'NTIA_General_Model';
    if ~exist('par_names')
        par_names = { 'si_loss' 'hv_loss' 'hv_gain' 'color1' 'si_gain' 'contati' 'color2' };
    end

    for cnt = 2:length(par_names)+1,
        pars.par_name{cnt} = par_names{cnt-1};
    end
catch
    status = 1;
    pars = [];
    if verbose,
        fprintf('\n%s\n', lasterr);
        fprintf('\tAborting.\n'); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vqm_value, pars, par_names] = model_run_Callback_vqm_general(proc, src);
% proc.si_feature
% proc.hv_feature
% proc.hvb_feature
% proc.ati_feature
% proc.cont_feature
% src.cb_feature

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

