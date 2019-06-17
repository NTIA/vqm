function par = vfd_clippar_loop_exp_hv_loss(test_structs, clip_structs, feature_base_dir)
% VFD_CLIPPAR_LOOP_EXP_HV_LOSS
%  Computes the vfd_exp_hv_loss_below5%_mean_square_clip_0.06 parameter
%  (i.e., par) of the VQM_VFD model. This function takes variable (1)
%  'test_structs' (of the same format as GTests), which describes each
%  video test and the location of the associated video, and variable (2) 
%  'clip_structs' (of the same format as GClips), which specifies the set
%  of video clips, and (3) the feature directory 'feature_base_dir', which
%  specifies the directory location of the features. The returned par will
%  be an ITS parameter structure that contains the parameters for all the
%  processed clips in clip_structs. The sort_clips_by function is used to
%  sort the clips before parameter calculation, and the original clips are
%  skipped if they are encountered.
% SYNTAX
%  par = vfd_clippar_loop_exp_hv_loss(test_structs, clip_structs, feature_base_dir)
% DESCRIPTION
%  This script uses the features generated by vfd_feature_loop_si_hv_adapt,
%  vfd_feature_loop_ti, and vfd_feature_loop_cont. In particular, the
%  following features must reside in the feature_base_dir:
%  1. feature_Y_vfd_hvA_0.4deg_0.2s_mean
%  2. feature_Y_vfd_hvbarA_0.4deg_0.2s_mean
%  3. feature_Y_vfd_ti_0.4deg_0.2s_rms
%  4. feature_Y_vfd_0.4deg_0.2s_mean
%  The hv_loss ST BLOCK parameters are weighted by quadratic motion and
%  intensity functions of the processed clip before spatial and temporal
%  collasping. The reason being that low/high motion and low/high intensity
%  blocks should have less weight than mid-level motion and intensity. The
%  general model's minimum thresholds, ST collapsing functions, parameter
%  squaring, and final parameter clipping thresholds are all used.
%

%  Define the features used for the block parameters.
%  This exp_hv_loss parameter uses 0.4deg_0.2s ST blocks.
f_hv = 'feature_Y_vfd_hvA_0.4deg_0.2s_mean';
f_hvb = 'feature_Y_vfd_hvbarA_0.4deg_0.2s_mean';
f_mot = 'feature_Y_vfd_ti_0.4deg_0.2s_rms';  % The average motion of a block (Y channel)
f_int = 'feature_Y_vfd_0.4deg_0.2s_mean';  % The average intensity of a block (Y channel)

%  Define the minimum threshold on HV and HVB features before dividing,
%  same as general model.
minthreshold = 3;

%  Define the spatial and temporal collasping functions, same as general
%  model.
spatial = 'below5%';
temporal = 'mean';

%  Define the final parameter clipping threshold after squaring.
clipthreshold = 0.06;

% Define the par_name for the returned parameter structure.
par.par_name = {'vfd_exp_hv_loss_below5%_mean_square_clip_0.06'};

% Parameters for quadratic intensity weighting (x, y): (0, ci), (di, 1),
% (2di, ci):  y = ai*x.^2 + bi*x + ci
% The amount of intensity reduction is limited to li.
di = 100.0;  % Peak location of intensity function : (di, 1.0)
ci = 0.64;  % Intensity reduction at 0 and 2di (1.0 at di)
li = 0.40;  % Limit on intensity reduction
ai = (ci-1)/di^2;
bi = 2*(1-ci)/di;

% Parameters for quadractic motion weighting
dm = 23.0;  % Peak location of motion function : (dm, 1.0)
cm = 0.75;  % Motion reduction at 0 and 2dm (1.0 at dm)
lm = 0.3;  % Limit on motion reduction
am = (cm-1)/dm^2;
bm = 2*(1-cm)/dm;

% Loop through all clips, sorted alphabetically
offsets = sort_clips_by('none', clip_structs, test_structs);
clip_structs = clip_structs(offsets);
ccnt = 1;  % processed clip counter for returned par array
for loop = 1:length(offsets)
    
    % Skip if an original clip was included.
    if strcmpi(clip_structs(loop).hrc{1},'original'),
        continue;
    end
    
    % Load intensity features
    data = zeros(1,0);
    datao = zeros(1,0);
    name = sprintf('%s/%s/%s_%s_%s.mat', feature_base_dir, f_int, ...
        clip_structs(loop).test{1}, ...
        clip_structs(loop).scene{1}, ...
        clip_structs(loop).hrc{1});
    if ~exist(name, 'file'),
        fprintf('Skipping clip %s:%s(%s), intensity feature filename does not exist.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    load( name );
    % Also skip this clip if no data or datao exists.
    [proc_r,proc_c,proc_t] = size(data);
    if proc_r*proc_c*proc_t == 0,
        fprintf('Skipping clip %s:%s(%s), no processed intensity features.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    [orig_r,orig_c,orig_t] = size(datao);
    if orig_r*orig_c*orig_t == 0,
        fprintf('Skipping clip %s:%s(%s), no original intensity features.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    intensity = data;
    
    %  Load motion features
    data = zeros(1,0);
    datao = zeros(1,0);
    name = sprintf('%s/%s/%s_%s_%s.mat', feature_base_dir, f_mot, ...
        clip_structs(loop).test{1}, ...
        clip_structs(loop).scene{1}, ...
        clip_structs(loop).hrc{1});
    if ~exist(name, 'file'),
        fprintf('Skipping clip %s:%s(%s), motion feature filename does not exist.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    load( name );
    % Also skip this clip if no data or datao exists.
    [proc_r,proc_c,proc_t] = size(data);
    if proc_r*proc_c*proc_t == 0,
        fprintf('Skipping clip %s:%s(%s), no processed motion features.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    [orig_r,orig_c,orig_t] = size(datao);
    if orig_r*orig_c*orig_t == 0,
        fprintf('Skipping clip %s:%s(%s), no original motion features.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    motion = data;
    
    % Load HV features
    data = zeros(1,0);
    datao = zeros(1,0);
    name = sprintf('%s/%s/%s_%s_%s.mat', feature_base_dir, f_hv, ...
        clip_structs(loop).test{1}, ...
        clip_structs(loop).scene{1}, ...
        clip_structs(loop).hrc{1});
    if ~exist(name, 'file'),
        fprintf('Skipping clip %s:%s(%s), HV feature filename does not exist.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    load( name );
    % Also skip this clip if no data or datao exists.
    [proc_r,proc_c,proc_t] = size(data);
    if proc_r*proc_c*proc_t == 0,
        fprintf('Skipping clip %s:%s(%s), no processed HV features.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    [orig_r,orig_c,orig_t] = size(datao);
    if orig_r*orig_c*orig_t == 0,
        fprintf('Skipping clip %s:%s(%s), no original HV features.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    hv = data;
    hvo = datao;
    
    % Load HVB features
    data = zeros(1,0);
    datao = zeros(1,0);
    name = sprintf('%s/%s/%s_%s_%s.mat', feature_base_dir, f_hvb, ...
        clip_structs(loop).test{1}, ...
        clip_structs(loop).scene{1}, ...
        clip_structs(loop).hrc{1});
    if ~exist(name, 'file'),
        fprintf('Skipping clip %s:%s(%s), HVB feature filename does not exist.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    load( name );
    % Also skip this clip if no data or datao exists.
    [proc_r,proc_c,proc_t] = size(data);
    if proc_r*proc_c*proc_t == 0,
        fprintf('Skipping clip %s:%s(%s), no processed HVB features.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    [orig_r,orig_c,orig_t] = size(datao);
    if orig_r*orig_c*orig_t == 0,
        fprintf('Skipping clip %s:%s(%s), no original HVB features.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    hvb = data;
    hvbo = datao;
    clear data datao;
    
    %  Calculate the intensity/motion block weighting for the pars.
    wint = ai*intensity.^2 + bi*intensity + ci;
    wint = max(wint,li);  % limit intensity reduction to li
    wmot = am*motion.^2 + bm*motion + cm;
    wmot = max(wmot,lm);  % limit motion reduction to lm
    weight = wint.*wmot;
    clear intensity motion wint wmot;
    
    %  Calculate the hv_loss BLOCK parameter, same as general model.
    this_par = compare_dual_feature(hvo, hvbo, hv, hvb, 'minthreshold', minthreshold, 'divide', 'compare','ratio_loss');
    this_par = this_par.*weight;
    clear weight;
    [trows, tcols, tsamps] = size(this_par);
    this_par = reshape(this_par, trows*tcols, tsamps);
    this_par = st_collapse(spatial, this_par);
    this_par = st_collapse(temporal, this_par);
    this_par = max(this_par.^2, clipthreshold) - clipthreshold;
    
    % Fill the clip_structs information into the returned parameter.
    par.clip_name{ccnt} = sprintf('%s_%s_%s', clip_structs(loop).test{1}, ...
        clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
    par.inlsa_mos(ccnt) = clip_structs(loop).inlsa_mos;
    par.mos(ccnt) = clip_structs(loop).mos;
    par.data(1, ccnt) = this_par;
    
    % Update clip counter.
    ccnt = ccnt + 1;
    
end
