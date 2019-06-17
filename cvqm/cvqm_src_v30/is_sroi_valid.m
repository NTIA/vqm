function [status, message] =  is_sroi_valid (roi, clip_struct, varargin)
% IS_SROI_VALID
%  Check whether a spatial region of interest is valid for this 
%  specific clip.  Optional check for even / odd coordinates.
% SYNTAX
%  [status, message] = is_sroi_valid(roi, clip_struct);
%  [...] = is_sroi_valid(..., 'Flag', ...);
% DESCRIPTION
%  [...] = is_sroi_valid(roi, clip_struct); checks the validity of 
%  region of interest, 'roi', given the spatial registration, processed 
%  valid region, and image size specified in 'clip_struct' (of the same 
%  format as GClips).
%  'roi' is a region of insterest, specified as a structure with four
%  elements: top, left, bottom, and right.  
%  'status' returned is 1 if ROI is valid, 0 if ROI is invalid.
%  'message' returns a string describing the problem, if any.
% Optional flags available are:
%  'verbose' print out to screen cause of ROI's invalidity.
%  'quiet'   Don't print anything to the screen.
%  'evenodd' check whether the top-left coordinate values are odd
%           and the bottom-right coordinate values are even.
% EXAMPLE
%  roi.top = 21; roi.left = 21; roi.bottom = 468; roi.right = 700;
%  if is_sroi_valid(roi, GClips), ...
% REMARKS
%  Fully checked for defects.

verbose = 0;
evenodd = 0;
for cnt = 3:nargin,
    if strcmp(varargin{cnt-2},'verbose') == 1,
        verbose = 1;
    elseif strcmp(varargin{cnt-2},'quiet') == 1,
        verbose = 0;
    elseif strcmp(varargin{cnt-2},'evenodd') == 1,
        evenodd = 1;
    else
        Error('Flag not recognized by is_sroi_valid');
    end
end

status = 1;
message = [];

% make sure the top-left coordinate is odd & the bottom-right coordinate
% even.
if evenodd,
	if mod(roi.top, 2) == 0 | mod(roi.left,2) == 0,
        status = 0;
        message = [ message sprintf('Error:  top and left coordinates must be odd\n')];
	end
	if mod(roi.bottom,2) == 1 | mod(roi.right,2) == 1,
        status = 0;
        message = [ message sprintf('Error:  bottom and right coordinates must be even\n')];
	end
end

% check that numbers are valid overall. Top-left most coordinate of the 
% image is (1,1).
if roi.top >= roi.bottom | roi.left >= roi.right,
    status = 0;
    message = [ message ...
        sprintf('Error:  top must be larger than bottom, right must be larger than left\n')];
end
if roi.top < 1 | roi.left < 1,
    status = 0;
    message = [ message sprintf('Error:  top and left must each be at least 1\n')];
end

% check that roi is within the image size
if roi.bottom > clip_struct.image_size.rows,
    status = 0;
    message = [ message sprintf('Error:  region of interest contains more lines than the image.\n')];
end
if roi.right > clip_struct.image_size.cols,
    status = 0;
    message = [ message sprintf('Error:  region of interest contains more columns than the image\n')];
end
if clip_struct.spatial.horizontal >= 0,
    % processed image was shifted to the right.  Must be shifted back to
    % the left.  This will leave an invalid portion on the right hand
    % side of the image.
    if clip_struct.image_size.cols - roi.right < clip_struct.spatial.horizontal,
        status = 0;
        message = [ message ...
            sprintf('Error:  region of interest includes right side of the image.\n')];
        message = [ message ...
            sprintf('Correcting the processed video''s spatial shift will invalidate this area.\n')];
    end
else
    % Processed image was shifted to the left.
    if roi.left - 1 < -clip_struct.spatial.horizontal,
        status = 0;
        message = [ message ...
            sprintf('Error:  region of interest includes left side of the image.\n')];
        message = [ message ...
            sprintf('Correcting the processed video''s spatial shift will invalidate this area.\n')];
    end
end
if clip_struct.spatial.vertical >= 0,
    % processed image was shifted to the down.  Must be shifted back
    % up.  This will leave an invalid portion on the bottom
    % of the image.
    if clip_struct.image_size.rows - roi.bottom < ...
            clip_struct.spatial.vertical + is_reframing_indicated(clip_struct),
        status = 0;
        message = [ message ...
            sprintf('Error:  region of interest includes bottom of the image.\n')];
        message = [ message ...
            sprintf('Correcting the processed video''s spatial shift will invalidate this area.\n')];
    end
else
    % Processed image was shifted up.
    if roi.top - 1 < -clip_struct.spatial.vertical + is_reframing_indicated(clip_struct),
        status = 0;
        message = [ message ...
            sprintf('Error:  region of interest includes top of the image.\n')];
        message = [ message ...
            sprintf('Correcting the processed video''s spatial shift will invalidate this area.\n')];
    end
end

% check if ROI is within CVR
[status2, message2]= is_within_cvr(roi,clip_struct,'verbose');

if status2 == 0,
    message = [ message message2 ];
    status = 0;
end

if verbose,
    fprintf('%s', message);
end

