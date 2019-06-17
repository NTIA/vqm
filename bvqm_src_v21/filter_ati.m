function ati = filter_ati(y, frames_wide)
% FILTER_ATI
%  Compute absolute value, temporal information (ATI) filter
% SYNTAX
%  ati = filter_ati(y);
%  ati = filter_ati(y, frames_wide);
% DEFINITION
%  Compute absolute value, temporal information (ATI) filter on a time-slice of
%  images.  Time-slice must include at least two images temporally!
%  Variable 'y' can either be a three dimensional time-slice (row,col,time); 
%  or a cell array of length {time}, where each cell contains one frame (row,col).
%  The second option will result in slightly faster run-times.
%
%  If 'frames_wide' is present, do a 'frames_wide' wide ATI filter rather
%  than using sequential images.

if nargin == 1,
    frames_wide = 1;
end
if frames_wide < 1,
    error('''frames_wide'' must be a positive integer');
end

if isa(y,'double') | isa(y,'single'),
    [r,c,t] = size(y);
    if t < frames_wide + 1;,
        error('ERROR:  Too few frames in time to compute ATI');
    end
    
    hold1 = y(:,:,1:t-frames_wide);
    hold2 = y(:,:,frames_wide+1:t);
    ati = abs(hold1 - hold2);

elseif isa(y,'cell'),
    t = length(y);
    if t < frames_wide + 1;,
        error('ERROR:  Too few frames in time to compute ATI');
    end
    [r,c] = size(y{1});

    ati = zeros(r,c,t-frames_wide);
    for loop = 1:t-frames_wide,
        ati(:,:,loop) = abs (y{loop} - y{loop+frames_wide});
    end
else
    error('data type of y is unknown');
end