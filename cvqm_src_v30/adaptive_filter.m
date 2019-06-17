function [filter_size, extra] = adaptive_filter (image_size);
% ADAPTIVE_FILTER
%   Given an image size, return adaptive filter size.
%   This function adapts the filter size of the spatial filters SI and HV
%   described in SPIE 1999 paper
% SYNTAX
%   [filter_size, extra] = adaptive_filter (image_size);
% DESCRIPTION
%   Given an image size (image_size.rows, image_size.cols), return the
%   optionimal SI & HV filter length ('filter_size') and the number of extra
%   pixels needed on all sides of the image ('extra').
%
%   Filter size adjusts automatically for the image size as follows:
%       QCIF, QSIF      SI5
%       CIF, SIF        SI9
%       VGA, 601, HDTV  SI13

% find adaptive filter size
if image_size.rows <= 216,
    filter_size = 5;
    extra = 2;
elseif image_size.rows <= 384,
    filter_size = 9;
    extra = 4;
else
    filter_size = 13;
    extra = 6;
end

        