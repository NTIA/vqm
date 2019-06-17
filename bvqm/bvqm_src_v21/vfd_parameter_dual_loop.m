function [pars] = vfd_parameter_dual_loop(test_structs, clip_structs, feature_base_dir,...
    feature_name1, feature_name2, feature_name3, spatial, temporal, varargin)
% VFD_PARAMETER_DUAL_LOOP
%   Given two Variable Frame Delay (VFD) features computed for a list of
%   clips, compute and return variant parameters.  Available functions
%   include Euclidean distance, and multiplying/dividing feature streams to
%   form a single feature stream.  The returned pars array is sorted
%   alphabetically by clip.
% SYNTAX
%  [pars] = vfd_parameter_dual_loop(test_structs, clip_structs, feature_base_dir,
%               feature_name1, feature_name2, feature_name3, spatial, temporal)
%  [pars] = vfd_parameter_dual_loop(...,'PropertyName', PropertyValue, ...);
%  [pars] = vfd_parameter_dual_loop(...,'PropertyName', ...);
% DESCRIPTION
%  Given two VFD features, compute parameters with different spatial-temporal
%  collapsing functions.  Ideally, the feature should have been computed
%  for all clips listed in 'clip_structs'.
%
%  This function takes variable 'test_structs' (of the same format as
%  GTests), which describes each video test and the location of the
%  associated video; and variable 'clip_structs' (of the same format as
%  GClips).  clip_structs may contain original clips but these are skipped
%  since the processed VFD feature files contain the original clip's
%  VFD features (i.e., original after VFD correction is made to look like
%  the processed clip).
%  'feature_base_dir' is the path to the directory where the feature's 
%  sub-directory were created.  The two features' names ('feature_name1')
%  are sub-directories in that directory, containing one file for each
%  clip.  'feature_name3' is the fictitious feature name to be used when
%  constructing the parameter's name (usually a combination of
%  feature_name1 and feature_name2), with 'feature_' prefix auto-removed.
%  'Spatial' and 'temporal' are cell arrays, listing the spatial and
%  temporal collapsing functions (respectively) to be applied to this
%  feature.  All combinations of the above are computed.
%
%  See function 'compare_dual_feature()' for the list of mandatory and
%  optional properties.
%
%  The following optional Properties may also be specified.
%   'MacroBlock', row, col, time,
%              Use macro-blocks of size (row,col,time) instead of standard
%              spatail then temporal collapsing.  After macroblocks, use 3D
%              collapsing of all remaining values.
%   '3D'       Collapse all spatial and temporal values simultaneously,
%              using the 3D collapsing functions listed in 'temporal'.
%              Cannot be used with 'MacroBlock'.  Cell array 'spatial' will
%              be ignored. 
%   'verbose'
%
%  The list of parameters is returned in struct 'pars', which contains the
%  following elements:
%   'par_name'  -- cell array listing the name of each parameter.
%   'clip_name' -- cell array listing the name of each clip.
%   'inlsa_mos' -- double array listing each clip's inlsa MOS.
%   'mos'       -- double array listing each clip's MOS.
%   'data'      -- matrix, indexed "(par, clip)", containing parameters.
%
% Example call for HV/HVbar, loss:
%
%  spatial = {'below5%' 'mean'};
%  temporal = {'10%' 'mean'};
% vfd_hv_loss = vfd_parameter_dual_loop(GTests_E, sd_Clips, 'e:/features/', ...
%         'feature_Y_vfd_hvA_0.4deg_0.2s_mean', 'feature_Y_vfd_hvbarA_0.4deg_0.2s_mean', ...
%         'feature_Y_vfd_hvA_0.4deg_0.2s_mean', spatial, temporal, ...
%         'MinThreshold', 3, 'divide', 'compare', 'ratio_loss');
%
% Example call for coherent color:
%
%  vfd_coher_color = vfd_parameter_dual_loop(GTests_E, sd_Clips,'e:/features/',...
%    'feature_Cb_vfd_0.4deg_0.2s_mean', 'feature_Cr_vfd_0.4deg_0.2s_mean', ...
%    'feature_coher_color_vfd_0.4deg_0.2s_mean', {'std'}, {'10%'}, 'weight', 1.5, 'euclid');
%

is_threshold = 0;
is_euclid = 0;
is_multiply = 0;
is_divide = 0;
is_logbefore = 0;
is_weight = 0;
compare = NaN;
is_verbose = 0;

is_macroblock = 0;
mb_row = 1;
mb_col = 1;
mb_time = 1;
is_3D = 0;

% parse varargin for naming convention only.  List of properties here MUST
% exactly match function compare_dual_feature().  Local properties won't be
% passed on to compare_dual_feature.

cnt = 1;
while cnt <= nargin-8,
    if strcmpi(varargin{cnt},'euclid')
        is_euclid = 1;
        varargin_pass{cnt} = varargin{cnt};
    elseif strcmpi(varargin{cnt},'multiply')
        is_multiply = 1;
        varargin_pass{cnt} = varargin{cnt};
    elseif strcmpi(varargin{cnt},'verbose')
        is_verbose = 1;
    elseif strcmpi(varargin{cnt},'logbefore')
        is_logbefore = 1;
        varargin_pass{cnt} = varargin{cnt};
    elseif strcmpi(varargin{cnt},'divide')
        is_divide = 1;
        varargin_pass{cnt} = varargin{cnt};
    elseif strcmpi(varargin{cnt},'compare')
        compare = varargin{cnt+1};
        varargin_pass{cnt} = varargin{cnt};
        cnt = cnt + 1;
        varargin_pass{cnt} = varargin{cnt};
    elseif strcmpi(varargin{cnt},'MinThreshold')
        is_threshold = 1;
        threshold = varargin{cnt+1};
        varargin_pass{cnt} = varargin{cnt};
        cnt = cnt + 1;
        varargin_pass{cnt} = varargin{cnt};
    elseif strcmpi(varargin{cnt},'Weight')
        is_weight = 1;
        weight = varargin{cnt+1};
        varargin_pass{cnt} = varargin{cnt};
        cnt = cnt + 1;
        varargin_pass{cnt} = varargin{cnt};
    elseif strcmpi(varargin{cnt},'MacroBlock'),
        is_macroblock = 1;
        is_3D = 0;
        mb_row = varargin{cnt+1};
        mb_col = varargin{cnt+2};
        mb_time = varargin{cnt+3};
        cnt = cnt + 3;
    elseif strcmpi(varargin{cnt},'3D'),
        is_macroblock = 0;
        is_3D = 1;
        if length(spatial) > 1,
            warning('Input argument ''spatial'' is undefined when property ''3D'' selected!  ''spatial'' will be ignored.');
        end
        clear spatial;
        spatial{1} = 'mean';
    else
        error('Property ''%s'' not recognized', varargin{cnt});
    end
    cnt = cnt + 1;
end

varargin = varargin_pass;

% make base parameter name.
if is_euclid,
    if is_threshold,
        par_name_piece = sprintf('%d_euclid', threshold);
    else
        par_name_piece = 'euclid';
    end
else
    if is_multiply,
        par_name_piece = 'multiply';
    elseif is_divide
        par_name_piece = 'divide';
    end
    par_name_piece = sprintf('%d_%s_%s', threshold, par_name_piece, compare);
end

% Pick 'feature_' off of feature_name3 if it is present, i.e., don't
% include in the construction of the parameter names!
if strncmp(feature_name3,'feature_',8),
    feature_name3 = feature_name3(9:length(feature_name3));
end

% figure out the list of parameters
num_parameters = length(spatial) * length(temporal);

% figure out the list of parameter names.
pars.par_name = cell(1,num_parameters);
pcnt = 1;
for space = 1:length(spatial),
    for time = 1:length(temporal),
        if is_macroblock,
            pars.par_name{pcnt} = sprintf('%s_%s_MB(%g,%g,%g)%s_ST%s', ...
                feature_name3, par_name_piece, mb_row, mb_col, mb_time, spatial{space}, temporal{time});
        elseif is_3D,
            pars.par_name{pcnt} = sprintf('%s_%s_ST%s', ...
                feature_name3, par_name_piece, temporal{time});
        else
            pars.par_name{pcnt} = sprintf('%s_%s_%s_%s', ...
                feature_name3, par_name_piece, spatial{space}, temporal{time});
        end
        pcnt = pcnt + 1;
    end
end

ccnt = 1;  % processed clip counter for returned pars array

% Loop through all clips, sorted alphabetically
offsets = sort_clips_by('none',clip_structs, test_structs);
clip_structs = clip_structs(offsets);
for loop = 1:length(offsets),
    
    % Skip if an original clip was included.
    if strcmpi(clip_structs(loop).hrc{1},'original'),
        continue;
    end
    
    % Load clip data, feature1;
    data = zeros(1,0);
    datao = zeros(1,0);
    name = sprintf('%s/%s/%s_%s_%s.mat',feature_base_dir,feature_name1, ...
        clip_structs(loop).test{1}, ...
        clip_structs(loop).scene{1}, ...
        clip_structs(loop).hrc{1});
    if ~exist(name, 'file'),
        fprintf('Skipping clip %s:%s(%s), feature filename does not exist.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    load( name );
    
    % Skip this clip if no data or datao exists.
    [proc_r,proc_c,proc_t] = size(data);
    if proc_r*proc_c*proc_t == 0,
        fprintf('Skipping clip %s:%s(%s), no processed features1.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    [orig_r,orig_c,orig_t] = size(datao);
    if orig_r*orig_c*orig_t == 0,
        fprintf('Skipping clip %s:%s(%s), no original features1.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    
    % Check to make sure original and processed dimensions match
    if orig_r ~= proc_r || orig_c ~= proc_c || orig_t ~= proc_t,
        fprintf('Skipping clip %s:%s(%s) processed feature1 size (%dx%dx%d) doesn''t match\noriginal feature1 size (%dx%dx%d)\n', ...
            clip_structs(loop).test{1}, clip_structs(loop).scene{1}, ...
            clip_structs(loop).hrc{1}, proc_r, proc_c, proc_t, orig_r, orig_c, orig_t);
        continue;
    end
    data1 = data;
    data1_orig = datao;
    
    % load clip data, feature2;
    data = zeros(1,0);
    datao = zeros(1,0);
    load( sprintf('%s/%s/%s_%s_%s.mat',feature_base_dir,feature_name2, ...
        clip_structs(loop).test{1}, ...
        clip_structs(loop).scene{1}, ...
        clip_structs(loop).hrc{1}));
    
    % Skip this clip if no data or datao exists.
    [proc_r,proc_c,proc_t] = size(data);
    if proc_r*proc_c*proc_t == 0,
        fprintf('Skipping clip %s:%s(%s), no processed features2.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    [orig_r,orig_c,orig_t] = size(datao);
    if orig_r*orig_c*orig_t == 0,
         fprintf('Skipping clip %s:%s(%s), no original features2.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    if orig_r ~= proc_r || orig_c ~= proc_c || orig_t ~= proc_t,
        fprintf('Skipping clip %s:%s(%s) processed feature2 size (%dx%dx%d) doesn''t match\noriginal feature2 size (%dx%dx%d)\n', ...
            clip_structs(loop).test{1}, clip_structs(loop).scene{1}, ...
            clip_structs(loop).hrc{1}, proc_r, proc_c, proc_t, orig_r, orig_c, orig_t);
        continue;
    end
    data2 = data;
    data2_orig = datao;
    
    % Fill the clip_structs information into the returned parameter
    pars.clip_name{ccnt} = sprintf('%s_%s_%s', clip_structs(loop).test{1}, ...
        clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
    pars.inlsa_mos(ccnt) = clip_structs(loop).inlsa_mos;
    pars.mos(ccnt) = clip_structs(loop).mos;
    
    if is_verbose,
        fprintf('%s\n', pars.clip_name{ccnt});
    end
    
    % Compare orig & processed.
    data = compare_dual_feature(data1_orig, data2_orig, data1, data2, varargin);
    
    if ~is_macroblock,
        % Reshape from 3D into 2D
        data = reshape(data,proc_r*proc_c,proc_t);
    end
    
    % Loop through all metrics & compute them.
    pcnt = 1;
    for space = 1:length(spatial),
        % Spatially collapse.
        if is_macroblock,
            collapse1 = st_collapse(spatial{space}, data, 'MacroBlock', mb_row, mb_col, mb_time);
        elseif is_3D,
            % Do nothing
        else
            collapse1 = st_collapse(spatial{space}, data);
        end
        % Loop through temporal collapsing functions.
        for time = 1:length(temporal),
            % Temporally collapse
            if is_macroblock,
                collapse2 = st_collapse(temporal{time}, collapse1, '3D');
            elseif is_3D,
                collapse2 = st_collapse(temporal{time}, data, '3D');
            else
                collapse2 = st_collapse(temporal{time}, collapse1);
            end
            
            % Record this parameter value for this clip.
            pars.data(pcnt, ccnt) = collapse2;
            pcnt = pcnt + 1;
        end
    end
    
    % Update clip counter.
    ccnt = ccnt + 1;
    
    % Clear arrays
    clear data*;
    
    pause(0.01);

end
