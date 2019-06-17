function [sroi,vert,horiz] = adjust_requested_sroi (struct, varargin)
% ADJUST_REQUESTED_SROI
%  Adjust the requested Spatial Region of Interest (SROI) as specified.
% SYNTAX
%  [sroi] = adjust_requested_sroi (struct)
%  [sroi] = adjust_requested_sroi (...,'PropertyName',PropertyValue,...);
%  [sroi,vert,horiz] = adjust_requested_sroi (...);
% DESCRIPTION
%  Given clip 'struct' (in the same format as GClips or Gsscqe), return the 
%  adjusted spatial region of interest ('sroi').  Return variable is a
%  structure, with four elements, 'roi.top', 'roi.left', 'roi.bottom',
%  and 'roi.right'.  NOTE: This sub-routine uses elements image_size, cvr,
%  and video_standard from the input variable 'struct'.
%
%  'sroi',...       Requested spatial region of interest (SROI), overriding
%                   default values of SROI.  Must be followed by 4 values, 
%                   specifying the region of interest, in the order: top, 
%                   bottom, left, right.  SROI will be adjusted.  Default
%                   SROI given by function 'default_sroi'.
%  'hsize', value,  Horizontal size of S-T blocks.  SROI must evenly divide
%                   by this value, horizontally.  Default is 1.
%  'vsize', value,  Vertically size of S-T blocks.  SROI must evenly divide
%                   by this value, vertizontally.  Default is 1.
%  'extra', value,  This many valid pixels are required on all sides of the
%                   SROI, extra pixels for filtering.  Returned SROI will
%                   NOT include those extra pixels!
%  'yxextra', y, x, Number 'y' indicates the number of valid pixels required 
%                   on top and bottom sides of the SROI; number 'x' indicates
%                   the number of valid pixels required on the left and right
%                   of the SROI.  Returned SROI will NOT include those extra pixels!
%  'evenodd',       Force the top-left coordinate to be odd, and force the
%                   bottom-right coordinate to be even.  Note that the
%                   default hsize and vsize will effectively be 2 instead
%                   of 1.  Default is to not have this optional restriction.
%
%  Optional return arguments 'vert' and 'horiz' will, if present, be filled
%  with the number of abutting blocks that fit vertically and
%  horizontally within the SROI.
% REMARKS
%  The top-left coordinate will be odd, and the bottom-right coordinate
%  even.  Thus, an equal number of pixels will be used from both fields.
%
%  Functionality tested.

% read values from struct that can be over written by variable argument
% list.
roi = default_sroi(struct.image_size);
xextra = 0;
yextra = 0;
hsize = 1;
vsize = 1;
evenodd = 0;

% parse variable argument list (property values)
cnt = 1;
while cnt <= nargin - 1,
    if strcmpi(varargin(cnt),'sroi') == 1,
        roi.top = varargin{cnt+1};
        roi.left = varargin{cnt+2};
        roi.bottom = varargin{cnt+3};
        roi.right = varargin{cnt+4};
        cnt = cnt + 5;
    elseif strcmpi(varargin(cnt),'hsize') == 1,
        hsize = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt),'vsize') == 1,
        vsize = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt),'extra') == 1,
        xextra = varargin{cnt+1};
        yextra = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmpi(varargin(cnt),'yxextra') == 1,
        yextra = varargin{cnt+1};
        xextra = varargin{cnt+2};
        cnt = cnt + 3;
    elseif strcmpi(varargin(cnt),'evenodd') == 1,
        evenodd = 1;
        cnt = cnt + 1;
    else
        error('Property value passed into adjust_requestetd_sroi not recognized');
    end
end

% Check minimum argument values.
if (hsize <= 0 || vsize <= 0)
    error('hsize and vsize must be greater than 0');
end
if xextra < 0 || yextra < 0,
    error('Number of extra pixels for filtering must be zero or positive');
end

%  Check to make sure block size is even for non-progressive (interlaced) video
if (isfield(struct,'video_standard'))
    if (~strcmpi(struct.video_standard,'progressive') && ((hsize > 2 && mod(hsize,2)) || (vsize > 2 && mod(vsize,2))) )
        error('hsize and vsize must be even for interlaced 4:2:2 video');
    end
end

% make sure are within CVR and have the extra pixels for filtering.
if roi.top < struct.cvr.top + yextra,
    roi.top = struct.cvr.top + yextra;
end
if roi.left < struct.cvr.left + xextra,
    roi.left = struct.cvr.left + xextra;
end
if roi.bottom > struct.cvr.bottom - yextra,
    roi.bottom = struct.cvr.bottom - yextra;
end
if roi.right > struct.cvr.right - xextra,
    roi.right = struct.cvr.right - xextra;
end


if evenodd,
    % We agreed to remove this restriction on Dec 13, 2005.
    % make sure top-left coordinates are odd, and bottom-right even
    if mod(roi.top,2) == 0,
        roi.top = roi.top + 1;
    end
    if mod(roi.left,2) == 0,
        roi.left = roi.left + 1;
    end
    if mod(roi.bottom,2) ~= 0,
        roi.bottom = roi.bottom - 1;
    end
    if mod(roi.right,2) ~= 0,
        roi.right = roi.right - 1;
    end
    
    if vsize == 1 && hsize== 1,
        % Called with the default block size is 1-pixels.
        % Pretend block size of 2 instead of 1, so that the even/odd
        % adjustment below works properly (i.e., adjust region to be
        % divisible by vsize and hsize)
        vsize = 2;
        hsize = 2;
    else
        % Called with an actual block size.
        % This block size must be divisible by 2 (both hsize and vsize)
        if mod(vsize,2),
            error('When ''evenodd'' flag is selected in ''adjust_requested_sroi'', vsize must be even');
        end
        % need hsize divisible by 2
        if mod(hsize,2),
            error('When ''evenodd'' flag is selected in ''adjust_requested_sroi'', hsize must be even');
        end
    end
end


% make sure region evenly divides by vsize & vsize.
while mod((roi.bottom - roi.top + 1),vsize),
    if roi.top < struct.image_size.rows - roi.bottom,
        if evenodd,
            roi.top = roi.top + 2;
        else
            roi.top = roi.top + 1;
        end
    else
        if evenodd
            roi.bottom = roi.bottom - 2;
        else
            roi.bottom = roi.bottom - 1;
        end
    end
end
while mod((roi.right - roi.left + 1),hsize),
    if roi.left < struct.image_size.cols - roi.right,
        if evenodd,
            roi.left = roi.left + 2;
        else
            roi.left = roi.left + 1;
        end
    else
        if evenodd,
            roi.right = roi.right - 2;
        else
            roi.right = roi.right - 1;
        end
    end
end

sroi = roi;
vert = (sroi.bottom-sroi.top+1)/vsize;
horiz = (sroi.right-sroi.left+1)/hsize;

