function ati = filter_ati_random(y, frames_wide, factor)
% FILTER_ATI_RANDOM
%  Compute absolute value, temporal information (ATI) filter.  Use a random
%  sub-sampling of pixels.
% SYNTAX
%  ati = filter_ati_random(y, frames_wide, factor);
% DEFINITION
%  Compute absolute value, temporal information (ATI) filter on a time-slice of
%  images.  Time-slice must include at least two images temporally!
%  Variable 'y'  must be a three dimensional time-slice (row,col,time); 
%
%  Do a 'frames_wide' wide ATI filter (e.g., 'frames_wide' = 2 to use
%  sequential images).  'factor' indicates the fraction of pixels to be used
%  where 0 < factor < 1.0
%
%  Returned variable 'ati' has dimensions (pixels, 1, time) and is
%  appropriate for use in image collapsing but not blocks (e.g.,
%  st_collapse.)


if frames_wide < 1,
    error('''frames_wide'' must be a positive integer');
end
if factor <= 0 | factor >= 1,
    error('''factor'' must be between zero and one');
end

[rows,cols,time] = size(y);
if time < frames_wide + 1;,
    error('ERROR:  Too few frames in time to compute ATI');
end

need = round(rows * cols * factor);
list_row = round(rand(1,need) * (rows) + 0.5);
list_col = round(rand(1,need) * (cols) + 0.5);
list = (list_col-1)*rows + list_row;

y = reshape(y, rows * cols, 1, time);

hold1 = y(list,1:time-frames_wide);
hold2 = y(list,frames_wide+1:time);
ati = abs(hold1 - hold2);

