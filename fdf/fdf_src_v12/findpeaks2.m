function [mags, locs] = findpeaks2(x, varargin);
% FINDPEAKS2
%   Find the peaks in vector x.
% SYNTAX
%   [mags, locs] = findpeaks2(x, 'threshold', thres);
% DESCRIPTION
%   This function will find all peaks in vector x and return their
%   magnitudes (mags) and locations (locs).  Note that the first and last
%   element in the vector will never be declared a peak.
%
%   Any or all of the following optional properties may be requested.
%
%   'threshold',thres    Specifies the threshold for the peak finder.  The
%                        current peak must be at least this much above the
%                        neighboring values to be declared a peak.  The
%                        default value for thres is zero.
%

% Validate input arguments and set their defaults
thres = 0;
cnt = 1;
while cnt <= length(varargin),
    if ~isstr(varargin{cnt}),
        error('Property value passed into findpeaks2 is not recognized');
    end
    if strcmpi(varargin(cnt),'threshold') == 1
        thres = varargin{cnt+1};
        cnt = cnt + 2;
    else
        error('Property value passed into findpeaks2 not recognized');
    end
end
n = length(x);
d = diff(x);
locs = 1+find(d(2:n-1)<-1.0*thres & d(1:n-2)>thres);
if (length(locs) ~= 0)
    mags = x(locs);
else
    locs = [];
    mags = [];
end


