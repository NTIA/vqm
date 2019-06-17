function [name_exists, name] = write_feature_standard(data, is_path, feature_name, clip_struct, varargin)
% WRITE_FEATURE_STANDARD
%  Write out features to a file with a standard name.
% SYNTAX
%  [name_exists name] = write_feature_standard(data, is_path, feature_name, clip_struct);
%  [name_exists name] = write_feature_standard(..., PropertyName', PropertyValue, ...);
%  [name_exists name] = write_feature_standard(..., PropertyName', ...);
% DESCRIPTION
%  Write out a feature stream ('data') to file into a sub-directory
%  'feature_name', with the is_path & drive specified by 'is_path'.  Variable
%  'feature_name' should be created using function 'standard_feature_names'
%  The feature stream ('data') should encompas all or part of the features
%  for one clip, clip_struct, of the same structure as GClips.  That 
%  information will be used to produce the file name.  Feature variable 
%  written will have the name 'data', to simplify the task of later 
%  routines automatically reading many such files.  The following 
%  properties are available:
%
%  'vfd', datao     Include an original variable frame delay (VFD) feature
%                   array 'datao' in the output feature file.
%  'overwrite'      Overwrite any existing feature files without asking.
%  'number',value   Number of this segment.  First segment in series must be
%                   #1.  If this property is not specified, then the entire
%                   time history of features is presumed to all be in 'data'.
%                   (Used to split the time history into segments that fit
%                   in memory without swapping to disk.)
%  'exists'         Just check on the existence of a feature file.  Don't 
%                   any write data.
%
%  By default, the full name of the feature file, 'name', will  be returned.
%  Return variable 'name_exists' is set to true (2) if that file existed
%  before this function call, and false (0) if it did not already exist.

write = 1;
overwrite = 0;
number = 0;
vfd = 0;

cnt = 1;
while cnt <= nargin-4,
    if strcmpi(varargin{cnt},'overwrite'),
        overwrite = 1;
    elseif strcmpi(varargin{cnt},'exists'),
        write = 0;
    elseif strcmpi(varargin{cnt},'number'),
        number = varargin{cnt+1};
        if number ~= floor(number) || number < 1,
            error('segment number must be an integer >= 1');
        end
        cnt = cnt + 1;
    elseif strcmpi(varargin{cnt},'vfd'),
        datao = varargin{cnt+1};
        vfd = 1;
        cnt = cnt +1;
    else
        error('write_feature_standard:  property not recognized');
    end
    
    cnt = cnt + 1;
end

% Remove end slashes from feature_name and is_path
if strcmp(feature_name(length(feature_name)),'/') || strcmp(feature_name(length(feature_name)),'\')
    feature_name = feature_name(1:(length(feature_name)-1));
end
if strcmp(is_path(length(is_path)),'/') || strcmp(is_path(length(is_path)),'\')
    is_path = is_path(1:(length(is_path)-1));
end

% Generate feature file name
name = [ is_path '/' feature_name '/' clip_struct.test{1} '_' clip_struct.scene{1} ];
name = [name '_' clip_struct.hrc{1} ];
if number > 0,
    name = [name '_' int2str(number) ];
end

name_mat = [name '.mat'];
name_exists = exist(name_mat, 'file');

% Create directory if needed
if (write && ~exist([is_path '/' feature_name],'dir'))
    mkdir(is_path, feature_name);
end

% Write feature file
if write,
    if (name_exists && ~overwrite),
        fprintf('Feature file "%s" exists.  File not overwritten.\n', name_mat);
        return;
    end
    if (vfd)
        save (name_mat, 'data', 'datao');
    else
        save (name_mat, 'data');
    end
end


