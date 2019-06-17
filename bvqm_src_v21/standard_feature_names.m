function [subdir] = standard_feature_names(varargin)
% STANDARD_FEATURE_NAMES
%  Compute the standard name (i.e., subdirectory name) for a feature. 
% SYNTAX
%  [subdir] = standard_feature_names( );
%  [subdir] = standard_feature_names(..., PropertyName', PropertyValue, ...);
%  [subdir] = standard_feature_names(..., PropertyName', ...);
% DESCRIPTION
%  Feature are named using a standard naming convention.  This standard
%  name will be used to name a subdirectory, within which each clip's data
%  is likewise given a standard name.  That functionality, however, is
%  outside of the scope of this function.
%
%  The following properties further specify the feature.
%  Some minimum set of properties must be specified to get a unique
%  name.  Some properties have associated values, while others do not.
%
%  'preaverage'     Images were pre-averaged.
%  'color'          Feature used color (Cb & Cr) information.
%  'YCbCr'          Feature used luminance and color information.
%  'Cr'             Feature used Cb color information.
%  'Cb'             Feature used Cr color information.
%  'Y'              Feature used luminance information, only.
%  'specific',value Variable 'value' contains a text string describes the
%                   calculations that make this feature unique.  (e.g.,
%                   hv13 or si13)
%  'sliding'        Present when S-T blocks slide (e.g., overlap in time)
%                   by one frame at a time.  By default, S-T blocks abut.
%  'fullimage'      Present when the S-T block size contains the entire
%                   valid region of each frame separately.
%  'fullfield'      Present when the S-T block size contains the entire
%                   valid region of each field separately (frames if
%                   progressive).
%  'vsize',value    Vertical size of the S-T block, in frame lines.
%  'hsize',value    Horizontal size of the S-T block, in pixels.
%  'degsize',value  A square block size specified in angular degrees.  This
%                   option is not compatible with the vsize and hsize
%                   options (it will override these options).
%  'tslice_length_sec',value    Temporal extent of the S-T block was 'value', given in
%                   seconds.
%  'block_stat',svalue   Statistical function used to extract feature from each
%                       S-T region, passed as a string (i.e., 'rms', 'std', 'mean')
%  'sscqe'               SSCQE feature.
%
%  The name of the sub-directory, which name identifies this feature, is
%  returned in 'subdir'. 

preaverage = 0;
color = 0;
specific = 0;
sliding = 0;
full_image = 0;
hsize = 1;
vsize = 1;
degsize = 0;  % =0 if normal vsize x hsize, =1 if degsize option is input
tslice_length_sec = 1;
block_stat = 0;
sscqe = 0;

cnt = 1;
while cnt <= nargin,
    if strcmpi(varargin{cnt},'preaverage'),
        preaverage = 1;
    elseif strcmpi(varargin{cnt},'y'),
        color = 0;
    elseif strcmpi(varargin{cnt},'color'),
        color = 1;
    elseif strcmpi(varargin{cnt},'ycbcr'),
        color = 2;
    elseif strcmpi(varargin{cnt},'cb'),
        color = 3;
    elseif strcmpi(varargin{cnt},'cr'),
        color = 4;
    elseif strcmpi(varargin{cnt},'specific'),
        specific = 1;
        specific_name = varargin{cnt+1};
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'sliding'),
        sliding = 1;
    elseif strcmpi(varargin{cnt},'fullimage'),
        full_image = 1;
    elseif strcmpi(varargin{cnt},'fullfield'),
        full_image = 2;
    elseif strcmpi(varargin{cnt},'vsize'),
        vsize = varargin{cnt+1};
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'hsize'),
        hsize = varargin{cnt+1};
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'degsize'),
        dsize = varargin{cnt+1};
        degsize = 1;
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'tslice_length_sec'),
        tslice_length_sec = varargin{cnt+1};
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'block_stat'),
        block_stat = 1;
        block_stat_name = varargin{cnt+1};
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'sscqe'),
        sscqe = 1;
        cnt = cnt + 1;
    else
        error('write_feature_standard:  property "%s" not recognized', varargin{cnt});
    end
    
    cnt = cnt + 1;
end


subdir = 'feature';
if sscqe,
    subdir = [ subdir '_sscqe' ];
end

if preaverage,
    subdir = [ subdir '_avg' num2str(tslice_length_sec,4) 's' ];
end
if color == 1,
    subdir = [subdir '_color'];
elseif color == 2,
    subdir = [subdir '_YCbCr'];
elseif color == 3,
    subdir = [subdir '_Cb'];
elseif color == 4,
    subdir = [subdir '_Cr'];
else
    subdir = [subdir '_Y'];
end
if specific,
    subdir = [subdir '_' specific_name];
end
if sliding,
    subdir = [subdir '_sliding' ];
end
if full_image == 1,
    subdir = [subdir '_image' ];
elseif full_image == 2,
    subdir = [subdir '_field' ];
elseif preaverage,
    if (~degsize)
        subdir = [subdir '_' int2str(vsize) 'x' int2str(hsize) ];
    else
        subdir = [subdir '_' num2str(dsize,4) 'deg' ];
    end
else
    if (~degsize)
        subdir = [subdir '_' int2str(vsize) 'x' int2str(hsize) '_' num2str(tslice_length_sec,4) 's' ];
    else
        subdir = [subdir '_' num2str(dsize,4) 'deg' '_' num2str(tslice_length_sec,4) 's' ];
    end
end

if block_stat,
    subdir = [subdir '_' block_stat_name ];
end

