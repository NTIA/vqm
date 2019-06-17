function [roi] = default_sroi (image_size)
% DEFAULT_SROI
%  Return the default spatial region of interest (SROI) for a given image
%  size.
% SYNTAX
%  [roi] = default_sroi (image_size)
% DESCRIPTION
%  [roi] = default_sroi (image_size); takes an image size structure with
%  two elements, 'image_size.rows' and 'image_size.cols', and returns the
%  default SROI for that image size, 'roi'.  The returned variable is also
%  a structure, with four elements, 'roi.top', 'roi.left', 'roi.bottom',
%  and 'roi.right'.
% REMARKS
%  If an image size is requested that does not NTSC / 525-line or PAL / 625-line,
%  then the default SROI encompasses the entire image.
%  

if (image_size.rows == 486 | image_size.rows == 480) & image_size.cols == 720,
    % NTSC / 525-line
    roi.top =    21;
    roi.left =   25;
    roi.bottom = 20+448;
    roi.right =  24+672;
elseif image_size.rows == 576 & image_size.cols == 720,
    % PAL / 625-line
    roi.top =    17;
    roi.left =   25;
    roi.bottom = 16+544;
    roi.right =  24+672;
elseif image_size.rows == 720 & image_size.cols == 1280, 
	% initialize maximum valid region.
	roi.top = 7;
	roi.left = 17;
	roi.bottom = image_size.rows - 6;
	roi.right = image_size.cols - 16;
elseif image_size.rows == 1080 & image_size.cols == 1920, 
	% initialize maximum valid region.
	roi.top = 7;
	roi.left = 17;
	roi.bottom = image_size.rows - 6;
	roi.right = image_size.cols - 16;
else
    roi.top = 1;
    roi.left = 1;
    roi.bottom = image_size.rows;
    roi.right = image_size.cols;
end
