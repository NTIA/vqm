function pars = dll_vfd_parameter_dual_loop(data1, data1_orig, data2, data2_orig, ...
    spatial, temporal, varargin)
% VFD_PARAMETER_DUAL_LOOP
%   Given two Variable Frame Delay (VFD) features computed for a pair of
%   clips, compute and return variant parameters.  Available functions
%   include Euclidean distance, and multiplying/dividing feature streams to
%   form a single feature stream.  The returned pars array is sorted
%   alphabetically by clip.  This function works in conjunction with other
%   dll functions for aditional information.
% SYNTAX
%  [pars] = dll_vfd_parameter_dual_loop(data1, data1_orig, data2, ...
%       data2_orig, spatial, temporal);
%  [pars] = dll_vfd_parameter_dual_loop(...,'PropertyName', PropertyValue, ...);
%  [pars] = dll_vfd_parameter_dual_loop(...,'PropertyName', ...);
% DESCRIPTION
%  Given two VFD features, compute parameters with different spatial-temporal
%  collapsing functions.
%
%  The first four inputs of this function are the arrays created by the
%  various dll feature extraction functions.
%
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
%  The returned value is the parameter the user specifies.
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

cnt = 1;
while cnt <= nargin-6,
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

% Compare orig & processed.
data = compare_dual_feature(data1_orig, data2_orig, data1, data2, varargin);

[proc_r,proc_c,proc_t] = size(data2);
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
        pars(pcnt) = collapse2;
        pcnt = pcnt + 1;
    end
end





