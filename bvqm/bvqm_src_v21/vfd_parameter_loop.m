function [pars] = vfd_parameter_loop(test_structs, clip_structs, feature_base_dir,...
    feature_name, compare, spatial, temporal, varargin)
% VFD_PARAMETER_LOOP
%  Given one Variable Frame Delay (VFD) feature computed for a list of
%  clips, compute and return variant parameters.  User specifies one
%  feature clipping threshold, one comparison function (see
%  compare_feature), and multiple spatial-temporal collapsing functions.
% SYNTAX
%  [pars] = vfd_parameter_loop(test_structs, clip_structs, feature_base_dir,
%               feature_name, compare, spatial, temporal);
%  [pars] = vfd_parameter_loop(...,'PropertyName', PropertyValue, ...);
% DESCRIPTION
%  Given one VFD feature, compute parameters with different spatial-temporal
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
%  sub-directory were created.  The feature's name ('feature_name') is a
%  sub-directory's name in that directory, containing one file for each
%  clip.
%  'compare' is the comparison function to be used, as specified by
%  function compare_feature().  Only one comparison function may be
%  specified.
%  'Spatial' and 'temporal' are cell arrays, listing the spatial and
%  temporal collapsing functions (respectively) to be applied to this
%  feature.  All combinations of the above are computed.
%
%  The following optional Properties may be specified.
%
%    'MinThreshold', #    Minimum threshold to be applied to features.
%                         This property must be specified with comparison
%                         functions log loss, log gain, ratio loss, and
%                         ratio gain.
%   'MacroBlock', row, col, time,
%              Use macro-blocks of size (row,col,time) instead of standard
%              spatail then temporal collapsing.  After macroblocks, use 3D
%              collapsing of all remaining values.
%   '3D'       Collapse all spatial and temporal values simultaneously,
%              using a single 3D collapsing functions listed in 'temporal'.
%              Cannot be used with 'MacroBlock'.  Cell array 'spatial' will
%              be ignored. 
%   'verbose'       Print out current clip name
%
%  The list of parameters is returned in struct 'pars', which contains the
%  following elements:
%   'par_name'  -- cell array listing the name of each parameter.
%   'clip_name' -- cell array listing the name of each clip.
%   'inlsa_mos' -- double array listing each clip's inlsa MOS.
%   'mos'       -- double array listing each clip's MOS.
%   'data'      -- matrix, indexed "(par, clip)", containing parameters.
%
% Example call for SI loss:
%
%  spatial = {'below5%' 'mean'};
%  temporal = {'10%' 'mean'};
% vfd_si_loss = vfd_parameter_loop(GTests_E, sd_Clips, 'e:/features/', ...
%         'feature_Y_vfd_siA_0.4deg_0.2s_std', 'ratio_loss', spatial, temporal, 'MinThreshold', 12);
%

is_threshold = 0;

is_macroblock = 0;
mb_row = 1;
mb_col = 1;
mb_time = 1;
is_3D = 0;
is_verbose = 0;

cnt = 1;
while cnt <= nargin - 7,
    if strcmpi(varargin{cnt},'MinThreshold')
        is_threshold = 1;
        min_threshold = varargin{cnt+1};
        cnt = cnt + 1;
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
    elseif strcmpi(varargin{cnt},'verbose'),
        is_verbose = 1;
    else
        if ischar(varargin{cnt}),
            error('Property ''%s'' not recognized', varargin{cnt});
        else
            error('Fatal function argument list error.  Optional property not recognized.');
        end    
    end
    cnt = cnt + 1;
end

% error checking.
if strcmpi(compare,'ratio_gain') || strcmpi(compare,'ratio_loss') || strcmpi(compare,'log_gain') || strcmpi(compare,'log_loss'),
    if ~is_threshold,
        error('A minimum threshold must be specified with comparison functions log loss, log gain, ratio loss, and ratio gain.');
    end
end

% figure out the list of parameters
num_parameters = length(spatial) * length(temporal);

% Pick 'feature_' off of the feature name if it is present, don't include
% in the construction of parameter names!
if strncmp(feature_name,'feature_',8),
    hold_feature_name = feature_name(9:length(feature_name));
else
    hold_feature_name = feature_name;
end

% figure out the list of parameter names.
pars.par_name = cell(1,num_parameters);
pcnt = 1;
for space = 1:length(spatial),
    for time = 1:length(temporal),
        holdn = hold_feature_name;
        if is_threshold,
            holdn = sprintf('%s_%d', holdn, min_threshold);
        end
        holdn = sprintf('%s_%s', holdn, compare);
        if is_macroblock,
            holdn = sprintf('%s_MB(%g,%g,%g)%s_ST%s', holdn, mb_row, mb_col, mb_time, spatial{space}, temporal{time});
        elseif is_3D,
            holdn = sprintf('%s_ST%s', holdn, temporal{time});
        else
            holdn = sprintf('%s_%s_%s', holdn, spatial{space}, temporal{time});
        end
        pars.par_name{pcnt} = holdn;
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
    
    % Load clip data;
    data = zeros(1,0);
    datao = zeros(1,0);
    name = sprintf('%s/%s/%s_%s_%s.mat',feature_base_dir,feature_name, ...
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
        fprintf('Skipping clip %s:%s(%s), no processed features.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    [orig_r,orig_c,orig_t] = size(datao);
    if orig_r*orig_c*orig_t == 0,
        fprintf('Skipping clip %s:%s(%s), no original features.\n', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
        continue;
    end
    
    % Check to make sure original and processed dimensions match
    if orig_r ~= proc_r || orig_c ~= proc_c || orig_t ~= proc_t,
        fprintf('Skipping clip %s:%s(%s) processed feature size (%dx%dx%d) doesn''t match\noriginal feature size (%dx%dx%d)\n', ...
            clip_structs(loop).test{1}, clip_structs(loop).scene{1}, ...
            clip_structs(loop).hrc{1}, proc_r, proc_c, proc_t, orig_r, orig_c, orig_t);
        continue;
    end
    
    % Fill the clip_structs information into the returned parameter
    pars.clip_name{ccnt} = sprintf('%s_%s_%s', clip_structs(loop).test{1}, ...
        clip_structs(loop).scene{1}, clip_structs(loop).hrc{1});
    pars.inlsa_mos(ccnt) = clip_structs(loop).inlsa_mos;
    pars.mos(ccnt) = clip_structs(loop).mos;
    
    if is_verbose,
        fprintf('%s\n', pars.clip_name{ccnt});
    end
    
    % Compare orig & processed.
    if is_threshold,
        data = compare_feature(datao, data, compare, 'MinThreshold', min_threshold);
    else
        data = compare_feature(datao, data, compare);
    end
    
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
            % Do nothing.
        else
            collapse1 = st_collapse(spatial{space}, data);
        end
        % Loop through temporal collapsing functions.
        for time = 1:length(temporal),
            % temporally collapse
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
