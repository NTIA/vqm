function [valid, cvr, sroi] = model_lowbw_sroi(extra, top, left, bottom, right);
% MODEL_LOWBW_SROI
%   Check image-size validity for lowbw or fastlowbw model; and return spatial region
%   of interest (SROI)
% SYNTAX
%   [valid, cvr, sroi] = model_lowbw_sroi(extra, top, left, bottom, right);
% DESCRIPTION
%   'Extra' is the number of pixels needed on all sides for filtering.
%   ONE EXTRA pixel will be needed for shifting.
%   (top, left, bottom, right) are the  valid region coordinates.
%   Within above coordinates, find where the low bandwidth model should be run. 
%
%   Return whether the model can validly be used (valid == 1) or not 
%   (valid == 0); the common valid region (cvr) including just
%   the extra +1 pixel border; and the spatial region of interest without
%   the extra pixels & lines (sroi).   Return variables 'cvr' and 'sroi' 
%   are structures, whose elements are top, left, bottom, and right.
%
%   Presume a block-size of 30x30, macro-block size of 3x3.
%   One extra pixel required in all directions, for spatial shift search. 

% extra pixel for shifting +- 1 in all directions.
extra = extra + 1;

if nargin ~= 5,
    error('number of input arguments to model_lowbw_sroi invalid.');
end
    
% if left or top are even, add one
top = top + (1 - mod(top,2));
left = left + (1 - mod(left,2));

% if bottom or right are odd, subtract one
right = right - mod(right,2);
bottom = bottom - mod(bottom,2);

% compute number of rows available, after discard border
num_rows = 30 * floor( ((bottom - top + 1) - extra*2)/30 );
extra_rows = floor(((bottom - top + 1) - num_rows)/2);

num_cols = 30 * floor( ((right - left + 1) - extra*2)/30 );
extra_cols = floor(((right - left + 1) - num_cols)/2);

% figure if valid
if num_rows / 30 < 3 | num_cols / 30 < 3,
    sroi = [];
    cvr = [];
    valid = 0;
    return;
end

% figure out SROI & CVR
sroi.top = top + extra_rows;
sroi.left = left + extra_cols;
sroi.bottom = sroi.top + num_rows - 1;
sroi.right = sroi.left + num_cols - 1;

valid = 1;

cvr.top = sroi.top-extra;
cvr.left = sroi.left-extra;
cvr.bottom = sroi.bottom+extra;
cvr.right = sroi.right+extra;

