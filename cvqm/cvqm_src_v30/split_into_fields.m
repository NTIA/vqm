function [one, two] = split_into_fields(y);
% SPLIT_INTO_FIELDS
%  Splits one frame into two fields.
% SYNTAX
%  [one two] = split_into_fields(y);
% DESCRIPTION
%  Split frame 'y' into field one ('one') and field two ('two').  Field two
%  contains the top line of the image; field one contains the second line
%  of the image.  For NTSC, field one occurs before field two in time.  For
%  PAL, the reverse is the case.  Y can be a time-slice of frames.
%
%  If image 'y' contains an odd number of rows, then the last (bottom) row
%  will be eliminated.  See also function 'join_into_frames'

[row, col, time] = size(y);
if mod(row,2),
    y = y(1:row-1,1:col,1:time);
    [row, col, time] = size(y);
end
y = reshape(y,2,row/2,col,time);
two = squeeze(y(1,:,:,:));
one = squeeze(y(2,:,:,:));
