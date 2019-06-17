function [status, message] = is_within_cvr(roi, clip_struct,varargin)
% IS_WITHIN_CVR
%  checks whether the given region of interest lies within the common
%  valid region.
% SYNTAX
%  [status, message] = is_within_cvr(roi, clip_struct);
%  [...] = is_within_cvr(...,'Flag',...);
% DESCRIPTION
%  is_within_cvr(clip_struct,roi) checks whether the region of interest,
%  roi, lies entirely within the common valid region specified in
%  clip_struct (of the same format as GClips).
%  'roi' is a region of insterest, specified as a structure with four
%  elements: top, left, bottom, and right.  
%  'status' returned is 1 if ROI is within the CVR, 0 otherwise.
%  'message' returns a string describing the issue, if any.
% Optional flags available are:
%  'verbose' print out to screen cause of ROI's invalidity.
% REMARKS
%  Functionality fully tested.

status = 1;
message = [];
verbose = 0;
for cnt = 3:nargin,
    if strcmp(varargin{cnt-2},'verbose') == 1,
        verbose = 1;
    else
        error('Flag not recognized by is_within_cvr');
    end
end


% check against cvr.
if clip_struct.cvr.top > roi.top,
    status = 0;
    message = [ message sprintf(...
        'ERROR:  top line of region of interest must be within the common valid region\n')];
end
if clip_struct.cvr.left > roi.left,
    status = 0;
    message = [ message sprintf(...
        'ERROR:  left pixel of region of interest must be within the common valid region\n')];
end
if clip_struct.cvr.bottom < roi.bottom,
    status = 0;
    message = [ message sprintf(...
        'ERROR:  bottom line of region of interest must be within the common valid region\n')];
end
if clip_struct.cvr.right < roi.right,
    status = 0;
    message = [ message sprintf(...
        'ERROR:  right pixel of region of interest must be within the common valid region\n')];
end

if verbose,
    fprintf('%s', message);
end
