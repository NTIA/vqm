function [roi] = dll_default_vr (fn)
% DLL_DEFAULT_VR
%  Return the default valid region for a given image size.
% SYNTAX
%  [roi] = dll_default_vr (fn)
% DESCRIPTION
%  This function takes fn, initialized in dll_video, and returns the default valid 
%  region for that image size, 'roi'.  The returned variable is also
%  a structure, with four elements, 'roi.top', 'roi.left', 'roi.bottom',
%  and 'roi.right'.

[image_size.rows,image_size.cols] = dll_video('size',fn);

if (image_size.rows == 486 | image_size.rows == 480) & image_size.cols == 720,
    % NTSC / 525-line
    roi.top =    19;
    roi.left =   23;
    roi.bottom = image_size.rows - 18;
    roi.right =  image_size.cols - 22;
elseif image_size.rows == 576 & image_size.cols == 720,
    % PAL / 625-line
    roi.top =    15;
    roi.left =   23;
    roi.bottom = image_size.rows - 14;
    roi.right =  image_size.cols - 22;
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
