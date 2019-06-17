function [data] = st_collapse(request, raw_data, varargin)
% ST_COLLAPSE
%  Compute the requested spatial or temporal collapsing function.
% SYNTAX
%  [data] = st_collapse(request, raw_data)
%  [data] = st_collapse(...,'PropertyName',...);
% DESCRIPTION
%  Compute the requested spatial or temporal collapsing function on the
%  FIRST dimension of the array or matrix, 'raw_data', and return the
%  results in 'data'.  The available precentile functions are:
%   'mean', 'std', 'rms', 'min', 'max', '10%', '25%', '50%', '90%', 
%   'above90%', 'above95%', 'above99%', 'above90%tail', 'above95%tail', 'above98%tail', 'above99%tail', 
%   'below5%', 'below10%', 'below1%', 'below1%tail, 'below2%tail', 'below5%tail', 'below10%tail'
%   'below25%', 'above75%', 'below2%', 'above98%', 'below50%tail'
%   'between25%50%'
%           [ The meanings of the above are as defined in "Video Quality
%             Measurement Techniques" NTIA Technical Report 02-392. ]
%   and 'minkowski(P,R)'.
%           [ minkowski = mean(abs(raw_data).^P).^(1/R) ]
%           Where 'P' and 'R' are replaced with the actual values to be
%           used.  For example, 'minkowski(1.8,2.8)' or 'minkowski(6,7)'
%
% The following optional parameters are also available.  All of these
% options are mutually exclusive.
%
%   'MacroBlock', row, col, time, 
%           Apply the requested function to macro blocks, of size (row,col,time).
%           If the region does not evenly divide, center spatially, and
%           abut with the end of the raw data temporally.
%   'OverlapMacroBlock', row, col, time, 
%           Apply the requested function to macro blocks, of size (row,col,time).
%           Unlike option 'MacroBlock', blocks overlap rather than
%           abutting.  Thus, returned data will be smaller only by
%           (row-1) rows,(col-1) columns and (time-1) in time.
%   'SlideMacroBlock', row, col, time, 
%           Apply the requested function to macro blocks, of size (row,col,time).
%           Blocks will overlap in time and abut spatially.  
%   '3D',   Apply the requested function simultaneously to all
%           dimensions.  Thus, convert all of the data into a 1D array and
%           apply the colllapsing function to that 1D array.
%   '1D',   Apply the requested function to only the first
%           dimensions.  For example, the variable 'raw_data' should be 2D
%           (spatial,temporal) or 1D (temporal).  This is the default
%           behavior.
%
%
% WARNING: Unless 'MacroBlock' or '3D' arguments are selected, the
% requested function will be applied to the first dimension only! 
%

collapse_3d = 0;
collapse_macroblock = 0;
collapse_overlap_macroblock = 0;
collapse_slide_macroblock = 0;
mb_row = 1;
mb_col = 1;
mb_time = 1;

cnt = 1;
while cnt <= nargin-2,
    switch lower(varargin{cnt}),
        case { 'macroblock' },
            collapse_macroblock = 1;
            mb_row = varargin{cnt+1};
            mb_col = varargin{cnt+2};
            mb_time = varargin{cnt+3};
            cnt = cnt + 4;
        case { 'overlapmacroblock' },
            collapse_overlap_macroblock = 1;
            mb_row = varargin{cnt+1};
            mb_col = varargin{cnt+2};
            mb_time = varargin{cnt+3};
            cnt = cnt + 4;
        case { 'slidemacroblock' },
            collapse_slide_macroblock = 1;
            mb_row = varargin{cnt+1};
            mb_col = varargin{cnt+2};
            mb_time = varargin{cnt+3};
            cnt = cnt + 4;
        case { '3d' },
            collapse_3d = 1;
            cnt = cnt + 1;
        case { '1d' },
            collapse_3d = 0;
            cnt = cnt + 1;
        otherwise
            error('Function st_collapse, optional argument "%s" not recognized', varargin{cnt});
    end
    
    if cnt <= nargin-2,
        error('Optional arguments are mutually exclusive.  Choose one only.');
    end
end

% Overlapping Macroblock Case.  Recurse to calculate.
if collapse_overlap_macroblock,
    [row_raw,col_raw,time_raw] = size(raw_data);
    data = zeros(row_raw-mb_row+1, col_raw-mb_col+1, time_raw-mb_time+1);
    [row,col,time] = size(data);
    
    for cnt1 = 1:mb_row,
        temp = floor((row_raw - cnt1 + 1) / mb_row) * mb_row;
        rng1 = cnt1:(cnt1-1+temp);
        if length(rng1) < mb_row,
            continue;
        end
        for cnt2 = 1:mb_col,
            temp = floor((col_raw - cnt2 + 1) / mb_col) * mb_col;
            rng2 = cnt2:(cnt2-1+temp);
            if length(rng2) < mb_col,
                continue;
            end
            for cnt3 = 1:mb_time,
                temp = floor((time_raw - cnt3 + 1) / mb_time) * mb_time;
                rng3 = cnt3:(cnt3-1+temp);
                if length(rng3) < mb_time,
                    continue;
                end
                data(cnt1:mb_row:row, cnt2:mb_col:col, cnt3:mb_time:time) = ...
                    st_collapse(request, raw_data(rng1,rng2,rng3), 'macroblock', ...
                    mb_row,mb_col,mb_time);
            end
        end
    end
    return;
end

% Slide Macroblock Case.  Recurse to calculate.
if collapse_slide_macroblock,
    [row_raw,col_raw,time_raw] = size(raw_data);
    row = floor(row_raw / mb_row);
    col = floor(col_raw / mb_col);
    time = time_raw - mb_time + 1;
    data = zeros(row,col,time);
    
    for cnt = 1:mb_time,
        temp = floor((time_raw - cnt + 1) / mb_time) * mb_time;
        rng = cnt:(temp+cnt-1);
        if length(rng) < mb_time,
            continue;
        end
        data(:,:, cnt:mb_time:time) = ...
            st_collapse(request, raw_data(:,:,rng), 'macroblock', ...
            mb_row,mb_col,mb_time);
    end
    return;
end

% if wanting to collapse a 3D structure all at once, reshape into an array.
if collapse_3d,
    [r,c,t] = size(raw_data);
    raw_data = reshape(raw_data, r*c*t,1);
end

% reshape for macroblocks.
if collapse_macroblock,
    [r,c,t] = size(raw_data);

    % error check.
    if floor(r / mb_row) == 0 | floor(c / mb_col) == 0 | floor(t / mb_time) == 0,
        error('st_collapse, macroblock too large.');
    elseif mb_row == 1 & mb_col == 1 & mb_time == 1,
        error('st_collapse, macro-block size must be greater than one.');
    end
    
    % figure out how to center the macroblocks
    [r,c,t] = size(raw_data);
    r_start = floor( mod(r,mb_row) / 2 ) + 1;
    r_stop = r_start + floor(r / mb_row)*mb_row - 1;
    
    c_start = floor( mod(c,mb_col) / 2 ) + 1;
    c_stop = c_start + floor(c / mb_col)*mb_col - 1;

    t_start = mod(t,mb_time) + 1;
    t_stop = t;
    
    % copy over macro-block data.
    raw_data = raw_data(r_start:r_stop,c_start:c_stop,t_start:t_stop);
    
    % reorganize first dimension
    [r,c,t] = size(raw_data);
    raw_data = reshape(raw_data, mb_row, r / mb_row, c, t);
    raw_data = permute(raw_data, [1 3 4 2]);
    raw_data = reshape(raw_data, mb_row * mb_col, c / mb_col, t, r / mb_row);
    raw_data = permute(raw_data, [1 3 4 2]);
    raw_data = reshape(raw_data, mb_row * mb_col * mb_time, t / mb_time, r / mb_row, c / mb_col);
    raw_data = permute(raw_data, [1 3 4 2]);    
end

% error checking.
[r,c,t] = size(raw_data);
if t > 1 & ~collapse_macroblock,
    error('Function st_collapse requires 2D or 1D arrays, only.  See help information warning.');
end

% Apply requested function.
above = 0;
below = 0;
tail = 0;

if r == 1,
    % Special case.  ST-collapse over a singleton dimension.  This is
    % always either the same as the input or zero.  
    if strcmp(request,'std') | length(findstr(request,'tail')) > 0,
        data = 0 * raw_data;
    else
        data = raw_data;
    end
    
% handle the usual cases.
elseif strcmp(request,'mean'),
    data = mean(raw_data);
elseif strcmp(request,'std'),
    data = std(raw_data);
elseif strcmp(request,'rms'),
    data = sqrt(mean(raw_data.^2));
elseif strcmp(request,'min'),
    data = min(raw_data);
elseif strcmp(request,'max'),
    data = max(raw_data);
elseif strncmp(request, 'minkowski(', 10),
    [mink n] = sscanf(request(11:length(request)), '%f,%f');
    if n ~= 2,
        error('Cannot parse minkowski P and R values in string "%s"', request(10:length(request)));
    end
    data = mean(abs(raw_data).^mink(1)).^(1.0/mink(2));
elseif strcmp(request,'between25%50%'),
    percentile1 = 0.25;
    percentile2 = 0.50;

	% if 1D but wrong direction vector, transpose it.
	if ndims(raw_data) == 2 & r == 1,
        raw_data = raw_data';
    end
	
	% compute percentile functions
	[rows,cols] = size(raw_data);
	
	want1 = 1 + round((rows-1) * percentile1);
	want2 = 1 + round((rows-1) * percentile2);
	
	temp = sort(raw_data, 1);
    data = mean(temp(want1:want2,:,:,:),1);
else
    if strcmp(request,'10%'),
        percentile = 0.10;
	elseif strcmp(request,'25%'),
        percentile = 0.25;
	elseif strcmp(request,'50%'),
        percentile = 0.50;
	elseif strcmp(request,'75%'),
        percentile = 0.75;
	elseif strcmp(request,'90%'),
        percentile = 0.90;
	elseif strcmp(request,'above95%'),
        percentile = 0.95;
        above = 1;
	elseif strcmp(request,'below5%'),
        percentile = 0.05;
        below = 1;
	elseif strcmp(request,'above99%'),
        percentile = 0.99;
        above = 1;
	elseif strcmp(request,'below1%'),
        percentile = 0.01;
        below = 1;
	elseif strcmp(request,'above98%'),
        percentile = 0.98;
        above = 1;
	elseif strcmp(request,'below2%'),
        percentile = 0.02;
        below = 1;
	elseif strcmp(request,'above90%'),
        percentile = 0.90;
        above = 1;
	elseif strcmp(request,'below10%'),
        percentile = 0.10;
        below = 1;
	elseif strcmp(request,'above75%'),
        percentile = 0.75;
        above = 1;
	elseif strcmp(request,'below25%'),
        percentile = 0.25;
        below = 1;
	elseif strcmp(request,'above95%tail'),
        percentile = 0.95;
        above = 1;
        tail = 1;
	elseif strcmp(request,'below5%tail'),
        percentile = 0.05;
        below = 1;
        tail = 1;
	elseif strcmp(request,'below50%tail'),
        percentile = 0.50;
        below = 1;
        tail = 1;
	elseif strcmp(request,'below2%tail'),
        percentile = 0.02;
        below = 1;
        tail = 1;
	elseif strcmp(request,'above98%tail'),
        percentile = 0.98;
        above = 1;
        tail = 1;
	elseif strcmp(request,'above99%tail'),
        percentile = 0.99;
        above = 1;
        tail = 1;
	elseif strcmp(request,'below1%tail'),
        percentile = 0.01;
        below = 1;
        tail = 1;
	elseif strcmp(request,'above90%tail'),
        percentile = 0.90;
        above = 1;
        tail = 1;
	elseif strcmp(request,'below10%tail'),
        percentile = 0.10;
        below = 1;
        tail = 1;
	else
        error('ERROR: percentile function "%s" not recognized by function compute_percentile', request);
    end
	
	% if 1D but wrong direction vector, transpose it.
	if ndims(raw_data) == 2 & r == 1,
        raw_data = raw_data';
    end
	
	% compute percentile functions
	[rows,cols] = size(raw_data);
	
	want = 1 + round((rows-1) * percentile);
	%fprintf('r=%d, c=%d, percentile %f, want=%d\n', rows, cols, percentile, want);
	
	temp = sort(raw_data, 1);
	if ~below & ~above & ~tail,
        data = temp(want,:,:,:);
	elseif above & ~tail,
        data = mean(temp(want:rows,:,:,:),1);
	elseif below & ~tail,
        data = mean(temp(1:want,:,:,:),1);
	elseif above & tail,
        if want == rows,
            % special case, can't do tail.
            data = temp(want,:,:,:) * 0;
        else
            data = mean(temp(want:rows,:,:,:),1) - temp(want,:,:,:);
        end
	elseif below & tail,
        if want == 1,
            % special case, can't do tail.
            data = temp(want,:,:,:) * 0;
        else
            data = temp(want,:,:,:) - mean(temp(1:want,:,:,:),1);
        end
	end;
end	
	
% get rid of extra dimension.
[a,b,c,d] = size(data);
data = reshape(data,b,c,d);
