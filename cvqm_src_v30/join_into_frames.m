function [y] = join_into_frames(one, two)
% JOIN_INTO_FRAMES
%  Join two fields into one frame.
% SYNTAX
%  [y] = join_into_frames(one, two);
% DESCRIPTION
%  Create frame 'y' from field one ('one') and field two ('two').  Field two
%  contains the top line of the image; field one contains the second line
%  of the image.  For NTSC, field one occurs before field two in time.  For
%  PAL, the reverse is the case.  Y can be a time-slice of frames.
%
%  If image 'one' and 'two' are not identically sized, then an error will occur.
%  See also function 'split_into_fields'

[row, col, time] = size(one);
[row2, col2, time2] = size(two);

if row ~= row2 || col ~= col2 || time ~= time2,
    error('two fields must be identically sized to be joined into a frame');
end

y(1,:,:,:) = two;
y(2,:,:,:) = one;
y = reshape(y, row*2, col, time);

