function [par] = compare_feature (orig, proc, compare, varargin)
% COMPARE_FEATURE
%  Perform comparison function on original & processed features and return
%  the block parameters.
%
% SYNTAX
%  [par] = compare_feature (orig, proc, compare);
%  [par] = compare_feature (..., 'PropertyName', PropertyValue, ...);
%
% DEFINITION
%  Given a set of original ('orig') and processed ('proc') spatial-temporal
%  features, and a comparision function ('compare'), return a set of block
%  parameters ('par') by performing a point-by-point comparision.  The
%  'orig' and 'proc' features can either be:
%   (1) 3D arrays (rows, cols, time) of features associated with
%       one processed video sequence and the associated original,
%   (2) 4D arrays (rows, cols, time, clips) of features associated
%       with two or more processed video sequences and thir originals, or
%   (3) cell arrays of length equal to the number of clips 
%  These formats are output by the function read_feature.
%
%  The implemented comparision functions 'compare' are as follows:
%    'ratio_gain'    positive part of ((proc - orig)/orig)
%    'ratio_loss'    negative part of ((proc - orig)/orig)
%    'log_gain'      positive part of log10(proc/orig)
%    'log_loss'      negative part of log10(proc/orig)
%    'minus_gain'    positive part of (proc - orig)
%    'minus_loss'    negative part of (proc - orig)
%
%  The following optional property name / property flag pairs may be
%  specified:
%  'MinThreshold', value    Clip 'orig' and 'proc' at the given minimum
%                           threshold value before computing the requested
%                           comparison function.
%  'MaxThreshold', value    Clip 'orig' and 'proc' at the given maximum
%                           threshold value before computing the requested
%                           comparison function.
%  

%  Determine if features are 3D/4D or cell arrays
if (iscell(orig) && iscell(proc)) % true if both orig and proc are cell arrays
    cell_input = 1;
    nclips = length(proc);
    nclipso = length(orig);
    if (nclips ~= nclipso)
        error('proc and orig cell arrays must be have same number of clips');
    end
elseif (~iscell(orig) && ~iscell(proc)) % true if both orig and proc are 4D arrays
    cell_input = 0;
    [nrows,ncols,ntimes,nclips] = size(proc);
    [nrowso,ncolso,ntimeso,nclipso] = size(orig);
    if (nrows ~= nrowso || ncols ~= ncolso || ntimes ~= ntimeso || nclips ~= nclipso)
        error('proc and orig 4D arrays must have same dimensions');
    end
else
    error ('Both orig and proc must be either 4D arrays or cell arrays');
end

minthres = 0;
maxthres = 0;
cnt = 1;
while cnt <= nargin - 3,
    if strcmpi(varargin(cnt),'minthreshold'),
        minthres = 1;
        minvalue = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt),'maxthreshold'),
        maxthres = 1;
        maxvalue = varargin{cnt+1};
        cnt = cnt + 2;
    else
        error('compare_feature optional property name "%s" is not recognized', varargin{cnt});
    end
end

if (~cell_input)  % Code for #D/4D arrays, no looping required

    if(minthres)  %  Apply min threshold if requested
        orig = max(orig, minvalue);
        proc = max(proc, minvalue);
    end
    
    if(maxthres)  %  Apply max threshold if requested
        orig = min(orig, maxvalue);
        proc = min(proc, maxvalue);
    end
    
    % Perform the comparision function
    if strcmpi(compare,'log_gain'),
        par = max( log10(proc ./ orig), 0);
    elseif strcmpi(compare, 'log_loss'),
        par = min( log10(proc ./ orig), 0);
    elseif strcmpi(compare,'ratio_gain'),
        par = max( (proc - orig) ./ orig, 0);
    elseif strcmpi(compare, 'ratio_loss'),
        par = min( (proc - orig) ./ orig, 0);
    elseif strcmpi(compare,'minus_gain'),
        par = max( (proc - orig), 0);
    elseif strcmpi(compare, 'minus_loss'),
        par = min( (proc - orig), 0);
    else
        error('compare_feature comparison function "%s" is not recognized',compare);
    end

else  % Cell input, must loop over clips since feature dimension can change
    
    for i = 1:nclips
        
        this_orig = orig{i};
        this_proc = proc{i};
        
        if(minthres)  %  Apply min threshold if requested
            this_orig = max(this_orig, minvalue);
            this_proc = max(this_proc, minvalue);
        end
        
        if(maxthres)  %  Apply max threshold if requested
            this_orig = min(this_orig, maxvalue);
            this_proc = min(this_proc, maxvalue);
        end
        
        % Perform the comparision function
        if strcmpi(compare,'log_gain'),
            this_par = max( log10(this_proc ./ this_orig), 0);
        elseif strcmpi(compare, 'log_loss'),
            this_par = min( log10(this_proc ./ this_orig), 0);
        elseif strcmpi(compare,'ratio_gain'),
            this_par = max( (this_proc - this_orig) ./ this_orig, 0);
        elseif strcmpi(compare, 'ratio_loss'),
            this_par = min( (this_proc - this_orig)./ this_orig, 0);
        elseif strcmpi(compare,'minus_gain'),
            this_par = max( (this_proc - this_orig), 0);
        elseif strcmpi(compare, 'minus_loss'),
            this_par = min( (this_proc - this_orig), 0);
        else
            error('compare_feature comparison function "%s" is not recognized',compare);
        end
        
        par{i} = this_par;
    
    end
    
end

end

