function [par] = compare_dual_feature (orig1, orig2, proc1, proc2, varargin)
% COMPARE_DUAL_FEATURE
%  Perform comparison function on two original & two processed features and
%  return the block parameters.
%
% SYNTAX
%  [par] = compare_dual_feature (orig1, orig2, proc1, proc2);
%  [par] = compare_dual_feature (..., 'PropertyName', PropertyValue, ...);
%  [par] = compare_dual_feature (..., 'PropertyName', ...);
%
% DEFINITION
%  Given two original ('orig1', 'orig2') and two processed ('proc1',
%  'proc2') spatial-temporal features, and a set of comparision function
%  options, return a set of block parameters ('par') by performing the
%  point-by-point comparisions.  The orig and proc features can either be:
%   (1) all are 3D arrays (rows, cols, time) of features associated with
%       one processed video sequence and the associated original.
%   (2) all are 4D arrays (rows, cols, time, clips) of features associated
%       with two or more processed video sequences and thir originals.
%   (3) all are cell arrays of length equal to the number of clips 
%  These formats are output by the function read_feature.
%
%  Required Arguments:  The comparision function argument list must specify
%  one of the following PropertyNames: 'euclid', 'multiply' or 'divide'.
%
%  'euclid'    Compare original and processed features using the Euclidean
%              distance between the vectors formed by (feature1,feature2).
%  'divide'    As per the hv13 feature, divide feature1 by feature2,
%              prior to comparing the original and processed.  The
%              'MinThreshold' and 'compare' options must also be present
%              (see below). 
%  'multiply'  As per the contrast-ati feature, multiply feature1 by
%              feature2, prior to comparing the original and processed.
%              The 'MinThreshold' and 'compare' options must also be
%              present (see below).
%
%  The following optional Properties may be specified.
%
%  'Weight', value        Weight (multiply) the second feature by the
%                         specified value before performing any further
%                         calculations.  For example, this option is useful
%                         for weighting the red chrominance component with
%                         respect to the blue before performing the
%                         'euclid' operation.
%  'MinThreshold', value  Specifies a minimum threshold value to be 
%                         applied to all the features before performing 
%                         the required computations are performed.
%  'LogBefore'            Do log10 before the 'multiply' or 'divide'.
%                         Useful when performing log_gain or log_loss for
%                         the contrast-ati parameter (see footnote on page
%                         54 of NTIA Report 02-392).  When the LogBefore
%                         option is specified, you must specify either the
%                         log_gain or log_loss 'compare' func.  This
%                         routine will then convert the log_gain or
%                         log_loss 'compare' func to minus_gain 
%                         or minus_loss, respectively, before calling
%                         compare_feature.
%  'compare', 'func'      After the two feature streams have been combined
%                         into one feature stream, this option specifies
%                         the comparison function 'func' to be used to
%                         generate the block parameter.  See function 
%                         compare_feature for the list of available
%                         functions.
%
% Example properties for the HV/HVbar loss parameter:
%         ..., 'MinThreshold', 3, 'divide', 'compare', 'ratio_loss');
%
% Example properties for contrast_ati log gain parameter:
%       ..., 'MinThreshold', 3, 'LogBefore', 'multiply', 'compare', 'log_gain');
%
% Example properties for the coher_color parameter where Cb is feature 1
% and Cr is feature 2: 
%       ..., 'Weight', 1.5, 'euclid');
%

%  Determine if features are 4D or cell arrays
if (iscell(orig1) && iscell(orig2) && iscell(proc1) && iscell(proc2)) % true if all orig and proc are cell arrays
    cell_input = 1;
    nclips = length(proc1);
    nclipso = length(orig1);
    if (nclips ~= nclipso)
        error('proc and orig cell arrays must be have same number of clips');
    end
elseif (~iscell(orig1) && ~iscell(orig2) && ~iscell(proc1) && ~iscell(proc2)) % true if all orig and proc are 4D arrays
    cell_input = 0;
    [nrows,ncols,ntimes,nclips] = size(proc1);
    [nrowso,ncolso,ntimeso,nclipso] = size(orig1);
    if (nrows ~= nrowso || ncols ~= ncolso || ntimes ~= ntimeso || nclips ~= nclipso)
        error('proc and orig 3D or 4D arrays must have same dimensions');
    end
else
    error ('Both orig and proc must be either 3D arrays, 4D arrays or cell arrays');
end

is_threshold = 0;
is_euclid = 0;
is_multiply = 0;
is_divide = 0;
is_logbefore = 0;
is_weight = 0;
is_compare = 0;

% If this function was called from another function, let that func. pass in
% its entire varargin as one variable.
if nargin == 5 && iscell(varargin{1}),
    varargin = varargin{1};
    addto = length(varargin) - 1;
else
    addto = 0;
end

cnt = 1;
while cnt <= nargin+addto-4,
    if strcmpi(varargin{cnt},'euclid')
        is_euclid = 1;
    elseif strcmpi(varargin{cnt},'multiply')
        is_multiply = 1;
    elseif strcmpi(varargin{cnt},'logbefore')
        is_logbefore = 1;
    elseif strcmpi(varargin{cnt},'divide')
        is_divide = 1;
    elseif strcmpi(varargin{cnt},'compare')
        is_compare = 1;
        compare = varargin{cnt+1};
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'MinThreshold')
        is_threshold = 1;
        threshold = varargin{cnt+1};
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'Weight')
        is_weight = 1;
        weight = varargin{cnt+1};
        cnt = cnt + 1;
    else
        error('compare_dual_feature property ''%s'' not recognized', varargin{cnt});
    end
    cnt = cnt + 1;
end

% Check arguments.
if is_euclid + is_multiply + is_divide ~= 1,
    error('Must specify only one of the following properties: ''euclid'', ''multiply'' or ''divide'' ');
elseif ~is_euclid && ~is_compare,
    error('Must specify ''compare'' function');
elseif (is_multiply || is_divide) && ~is_threshold,
    error('Must specify a minimum threshold when using ''multiply'' or ''divide''');
end

if (~cell_input) % Code for 3D and 4D arrays: no looping required
    
    % weight original feature2, if requested.
    if is_weight,
        orig2 = orig2 * weight;
        proc2 = proc2 * weight;
    end
    
    % Compare original & processed.
    
    % Compare and combine at once, Euclidean distance.
    if is_euclid,
        if is_threshold,
            orig1 = max(orig1,threshold);
            orig2 = max(orig2,threshold);
            proc1 = max(proc1,threshold);
            proc2 = max(proc2,threshold);
        end
        par = sqrt((orig1-proc1).^2 + (orig2-proc2).^2);
        
    % If doing log_gain or log_loss comparison func and logbefore is specified,
    % do the log10 first and then the multiply or divide. Change the
    % comparision function to minus_gain or minus_loss before calling
    % compare_feature.
    elseif (is_multiply || is_divide) && strncmpi(compare,'log',3) && is_logbefore,
        if is_multiply,
            orig = log10(max(orig1,threshold)) .* log10(max(orig2,threshold));
            proc = log10(max(proc1,threshold)) .* log10(max(proc2,threshold));
        elseif is_divide,
            orig = log10(max(orig1,threshold)) ./ log10(max(orig2,threshold));
            proc = log10(max(proc1,threshold)) ./ log10(max(proc2,threshold));
        else
            error('compare_dual_feature ''LogBefore'' option requires either ''multiply'' or ''divide''');
        end
        if strcmpi(compare,'log_gain'),
            par = compare_feature(orig, proc, 'minus_gain');
        elseif strcmpi(compare,'log_loss'),
            par = compare_feature(orig, proc, 'minus_loss');
        else
            error('compare_dual_feature ''LogBefore'' option requires ''compare'' func log_gain or log_loss');
        end
        
    % All other cases.  Doing multiply or divide of feature streams and
    % then the comparison function.
    else
        if is_multiply,
            orig = max(orig1,threshold) .* max(orig2,threshold);
            proc = max(proc1,threshold) .* max(proc2,threshold);
        elseif is_divide,
            orig = max(orig1,threshold) ./ max(orig2,threshold);
            proc = max(proc1,threshold) ./ max(proc2,threshold);
        else
            error('Unrecognized option in compare_dual_feature');
        end
        par = compare_feature(orig, proc, compare);
    end

else   % Cell input, must loop over clips since feature dimension can change
    
    for i = 1:nclips
        
        this_orig1 = orig1{i};
        this_orig2 = orig2{i};
        this_proc1 = proc1{i};
        this_proc2 = proc2{i};
        
        % weight original feature2, if requested.
        if is_weight,
            this_orig2 = this_orig2 * weight;
            this_proc2 = this_proc2 * weight;
        end
        
        % Compare original & processed.
        
        % Compare and combine at once, Euclidean distance.
        if is_euclid,
            if is_threshold,
                this_orig1 = max(this_orig1,threshold);
                this_orig2 = max(this_orig2,threshold);
                this_proc1 = max(this_proc1,threshold);
                this_proc2 = max(this_proc2,threshold);
            end
            par{i} = sqrt((this_orig1-this_proc1).^2 + (this_orig2-this_proc2).^2);
            
        % If doing log_gain or log_loss comparison func and logbefore is specified,
        % do the log10 first and then the multiply or divide. Change the
        % comparision function to minus_gain or minus_loss before calling
        % compare_feature.
        elseif (is_multiply || is_divide) && strncmpi(compare,'log',3) && is_logbefore,
            if is_multiply,
                this_orig = log10(max(this_orig1,threshold)) .* log10(max(this_orig2,threshold));
                this_proc = log10(max(this_proc1,threshold)) .* log10(max(this_proc2,threshold));
            elseif is_divide,
                this_orig = log10(max(this_orig1,threshold)) ./ log10(max(this_orig2,threshold));
                this_proc = log10(max(this_proc1,threshold)) ./ log10(max(this_proc2,threshold));
            else
                error('compare_dual_feature ''LogBefore'' option requires either ''multiply'' or ''divide''');
            end
            if strcmpi(compare,'log_gain'),
                par{i} = compare_feature(this_orig, this_proc, 'minus_gain');
            elseif strcmpi(compare,'log_loss'),
                par{i} = compare_feature(this_orig, this_proc, 'minus_loss');
            else
                error('compare_dual_feature ''LogBefore'' option requires ''compare'' func log_gain or log_loss');
            end
            
        % All other cases.  Doing multiply or divide of feature streams and
        % then the comparison function.
        else
            if is_multiply,
                this_orig = max(this_orig1,threshold) .* max(this_orig2,threshold);
                this_proc = max(this_proc1,threshold) .* max(this_proc2,threshold);
            elseif is_divide,
                this_orig = max(this_orig1,threshold) ./ max(this_orig2,threshold);
                this_proc = max(this_proc1,threshold) ./ max(this_proc2,threshold);
            else
                error('Unrecognized option in compare_dual_feature');
            end
            
            par{i} = compare_feature(this_orig, this_proc, compare);
            
        end
        
    end
    
end

end
