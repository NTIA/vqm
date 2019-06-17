function par = dll_vfd_parameter_loop(data, datao, compare, spatial, temporal, varargin)
% VFD_PARAMETER_LOOP
%  Given a pair of clips (an original and a processed) with a Variable
%  Frame Delay feature already computed for the pair, compute and return
%  variant parameters.  User specifies a data feature previously
%  extracted from the clips, one feature clipping threshold, one comparison
%  function (see compare_feature), and multiple spatial-temporal collapsing
%  functions.  
%
%  This function works in conjunction with other dll functions,
%  so many pieces of information come from various dll functions.
% SYNTAX
%  [pars] = dll_vfd_parameter_loop(data, datao, feature_name, compare, spatial, temporal);
%  [pars] = dll_vfd_parameter_loop(...,'PropertyName', PropertyValue, ...);
% DESCRIPTION
%  Given one VFD feature (user specified), compute parameters with
%  different spatial-temporal collapsing functions.
%
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
%  The requested parameter will be returned.
%
%
is_threshold = 0;

is_macroblock = 0;
mb_row = 1;
mb_col = 1;
mb_time = 1;
is_3D = 0;

cnt = 1;
while cnt <= nargin - 5,
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

[orig_r,orig_c,orig_t] = size(datao);
[proc_r,proc_c,proc_t] = size(data);

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
        par = collapse2;
    end
end





