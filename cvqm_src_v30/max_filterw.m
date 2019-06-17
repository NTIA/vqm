function [result] = max_filterw(feat,w)
% MAX_FILTERW
%   Perform a max filter of width 'w' on one or more 1-D arrays.
% SYNTAX
%   [result] = max_filterw(feat,w)
% DESCRIPTION
%  'feat' is either a 1-dimension array, where the first dimension contains 
%         the data (i.e., a column vector), or a 2-dimension array, in which 
%         case the function operates independently on each column.  Put
%         another way, 'feat' is either a column vector or a matrix of column
%         vectors.
%  'w' is the filter width, which must be odd. 
% 
%  Each point is replaced by the maximum of itself and the adjacent points
%  (in that column) +/- ((w-1)/2) neighboring points on either side.  This is
%  aplied to each column vector independently.  At the end points, zero
%  buffering is used (e.g., the beginning point presumes that zeros occured
%  previously). 

if mod(w,2) == 0 | w <= 1,
    error('Argument ''w'' must be an odd number greater than one');
end

[nr,nc] = size(feat);

% mat0 holds the zero padded, centered array
mat0 = zeros(nr+(w-1),nc);
mat0(1+floor(w/2):nr+floor(w/2),:) = feat;

% do negative shifts
result = mat0;
for i = -1:-1:-1*floor(w/2)
    result = max(result,circshift(mat0,i));
end

% do positive shifts
for i = 1:floor(w/2)
    result = max(result,circshift(mat0,i));
end

result = result(1+floor(w/2):nr+floor(w/2),:);
    
