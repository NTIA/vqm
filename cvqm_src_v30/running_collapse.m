function [data] = running_collapse (request, raw_data, delta, property);
% RUNNING_COLLAPSE
%   Perform a "running" spatial-temporal collapse over a long sequence.
% SYNTAX
%   [data] = running_collapse (request, raw_data, delta, property);
% DESCRIPTION
%   Take a request & raw_data, as defined by function st_collapse.
%   'property' is one of st_collapse's optional properties.  Instead of
%   calling st_collapse directly, divide it into shorter arrays, where the
%   last non-unit dimension is of length 'delta'.  Return the "running"
%   collapse in 'data'.  The length of data will be identical to the last
%   dimension of 'raw_data'.  The leading values (1 to time-1) will use
%   shorter time lengths.
%
%   'raw_data' must be either 1D or 3D.  See also function st_collapse.m

[row,col,time] = size(raw_data);
if time == 1,
    % 1D
    if col ~= 1,
        error('raw_data cannot be 2D');
    end
    time = row;
    data = zeros(time,1);
    for cnt = 1:time,
        start = max(1, cnt - delta + 1);
        stop = cnt;
        data(cnt) = st_collapse(request, raw_data(start:stop), property); 
    end
else
    % 3D
    data = zeros(time,1);
    for cnt = 1:time,
        start = max(1, cnt - delta + 1);
        stop = cnt;
        data(cnt) = st_collapse(request, raw_data(:,:,start:stop), property);
    end
end
