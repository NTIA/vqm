function [status, message] = check_clips(clip_struct, test_struct, varargin)
% CHECK_CLIPS
%  Check the validity of a clip structure:  all required fields present and
%  values are internally consistent.
% SYNTAX
%  [status, message] = check_clips(clip_struct, test_struct);
%  [...] = check_clips(...,'PropertyName',...);
% DESCRIPTION
%  Check internal consistency of clip_struct (of the same format as GClips)
%  and test_struct (of the same format as GTests).
%  Returned 'status' is false (0) when all data checks out as okay.
%  1 indicates a structural error / missing field.  2 indicates data
%  inconsistent with expected values.  3 indicates extreme values that are
%  likely erroneous.
%  Returned 'message' contains a string with a text description of any
%  problem.  Warnings may occur when status == 0.
%
%  Optional properties:
%   'verbose'   Print inconsistency issues to screen
%   'quiet'     Print nothing to the screen.
%

verbose = 0;
for cnt = 1:length(varargin),
    if strcmp(lower(varargin{cnt}),'verbose'),
        verbose = 1;
    elseif strcmp(lower(varargin{cnt}),'quiet'),
        verbose = 0;
    else
        error('optional property not recognized');
    end
end

if verbose,
    verbose_string = 'verbose';
else
    verbose_string = 'quiet';
end


fields(1) = {'test'};
fields(2) = {'scene'};
fields(3) = {'hrc'};
fields(4) = {'inlsa_mos'};
fields(5) = {'mos'};
fields(6) = {'stdev'};
fields(7) = {'viewers'};
fields(8) = {'image_size'};
fields(9) = {'fps'};
fields(10) = {'video_standard'};
fields(11) = {'subj_system'};
fields(12) = {'cvr'};
fields(13) = {'spatial'};
fields(14) = {'luminance_gain'};
fields(15) = {'luminance_offset'};
fields(16) = {'file_name'};
fields(17) = {'align_start'};
fields(18) = {'align_stop'};
fields(19) = {'loc_start'};
fields(20) = {'loc_stop'};
fields(21) = {'hrc_definition'};
fields(22) = {'scene_definition'};

status = 0;
message = [];
for cnt = 1:max(size(fields)),
    if ~isfield( clip_struct(1),fields{cnt}),
        message = [ message sprintf('Error:  Structure missing field ''%s''\n', fields{cnt})];
        status = 1;
    end
end
if status,
    bool = 0;
    if verbose,
        fprintf('%s\n', message);
    end
    return;
end
if verbose,
    message = sprintf('All mandatory top-level fields are present.\n\n');
end

% Go on to check validity of data elements that have limitations.
for loop = 1:max(size(clip_struct)),
    clip_name = sprintf('Clip [%d] %s:%s(%s) File: %s', loop, clip_struct(loop).test{1}, ...
        clip_struct(loop).scene{1}, clip_struct(loop).hrc{1}, ...
        clip_struct(loop).file_name{1});
    
    if ~isfield(clip_struct(loop).image_size,'rows') || ...
            ~isfield(clip_struct(loop).image_size,'cols'),
        new_issue = sprintf('Error:  image_size must consist of fields ''rows'' and ''cols''\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 1;
        break;
    end
    if ~isfield(clip_struct(loop).scale,'horizontal') || ...
            ~isfield(clip_struct(loop).scale,'vertical'),
        new_issue = sprintf('Error:  scale must consist of fields ''horizontal'' and ''vertical''\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 1;
        break;
    end
    if ~isfield(clip_struct(loop).spatial,'horizontal') || ...
            ~isfield(clip_struct(loop).spatial,'vertical'),
        new_issue = sprintf('Error:  spatial must consist of fields ''horizontal'' and ''vertical''\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 1;
        break;
    end

    if ~isfield(clip_struct(loop).cvr,'top') || ...
            ~isfield(clip_struct(loop).cvr,'left') || ...
            ~isfield(clip_struct(loop).cvr,'bottom') || ...
            ~isfield(clip_struct(loop).cvr,'right'),
        new_issue = sprintf('Error:  cvr must consist of fields ''top'', ''left'', ''bottom'' and ''right''\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 1;
        break;
    end
    if ~strcmp(clip_struct(loop).video_standard,'progressive') && ...
           ~strcmp(clip_struct(loop).video_standard,'interlace_lower_field_first') && ...
           ~strcmp(clip_struct(loop).video_standard,'interlace_upper_field_first'),
        new_issue = sprintf('Error:  value for ''video_standard'' is invalid\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if clip_struct(loop).image_size.rows < 2 || clip_struct(loop).image_size.cols < 2,
        new_issue = sprintf('Error:  image_size must be at least 2x2\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if clip_struct(loop).stdev < 0,
        new_issue = sprintf('Error:  std is negative\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    
    if ~(clip_struct(loop).spatial.horizontal >= -100 && clip_struct(loop).spatial.horizontal <= 100),
        new_issue = sprintf('Error:  Horizontal spatial shift must be between -100 and 100\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if ~(clip_struct(loop).spatial.vertical >= -100 && clip_struct(loop).spatial.vertical <= 100),
        new_issue = sprintf('Error:  Vertical spatial shift must be between -100 and 100\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if ~(clip_struct(loop).scale.horizontal >= 500 && clip_struct(loop).scale.horizontal <= 1500),
        new_issue = sprintf('Error:  Horizontal spatial scaling must be between 500 and 1500\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if ~(clip_struct(loop).scale.vertical >= 500 && clip_struct(loop).scale.vertical <= 1500),
        new_issue = sprintf('Error:  Vertical spatial scaling must be between 500 and 1500\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if ~(clip_struct(loop).luminance_gain >= 0.1 || clip_struct(loop).luminance_gain <= 10),
        new_issue = sprintf('Error:  Luminance gain invalid.\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if ~(clip_struct(loop).luminance_offset >= -150 || clip_struct(loop).luminance_offset <= 150),
        new_issue = sprintf('Error:  Luminance offset invalid.\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    [status2, message2] = is_sroi_valid(clip_struct(loop).cvr, clip_struct(loop),verbose_string);
    if ~status2,
        new_issue = sprintf('Error:  Common Valid Region is invalid\n%s', message2);
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if clip_struct(loop).loc_start < 1,
        new_issue = sprintf('Error:  Location of the first valid frame (loc_start) is invalid.  Must be 1+\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if clip_struct(loop).align_start < clip_struct(loop).loc_start,
        new_issue = sprintf('Error:  align_start must be at or after first valid frame (loc_start)\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if clip_struct(loop).loc_stop <= clip_struct(loop).loc_start,
        new_issue = sprintf('Error:  last valid frame (loc_stop) must be at or after first valid frame (loc_start)\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end  
    if clip_struct(loop).align_stop > clip_struct(loop).loc_stop,
        new_issue = sprintf('Error:  align_stop must be at or before last valid frame (loc_stop)\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if ~isnan(clip_struct(loop).align_start) || ~isnan(clip_struct(loop).align_stop),
        if strcmp(clip_struct(loop).hrc{1},'original'),
            if isnan(clip_struct(loop).align_start) || isnan(clip_struct(loop).align_stop) || ...
                    mod(clip_struct(loop).align_start,1) || mod(clip_struct(loop).align_stop,1),
                new_issue = sprintf('Error: align start and align stop should be integers\n');
                message = cc_append_message(message, clip_name, new_issue);
                status = 2;
            end
        elseif isnan(clip_struct(loop).align_start) || isnan(clip_struct(loop).align_stop) || ...
                mod(clip_struct(loop).align_start*2,1) || mod(clip_struct(loop).align_stop*2,1),
            new_issue = sprintf('Error: align start and align stop should be integers\n');
            message = cc_append_message(message, clip_name, new_issue);
            status = 2;
        end
    end
    if floor(clip_struct(loop).loc_start)~=clip_struct(loop).loc_start || ...
            floor(clip_struct(loop).loc_stop)~=clip_struct(loop).loc_stop || ...
            isnan(clip_struct(loop).loc_start) || isnan(clip_struct(loop).loc_stop),
        new_issue = sprintf('Error: first and last valid frame must both be integers\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end
    if (clip_struct(loop).cvr.right - clip_struct(loop).cvr.left + 1) < (32*2+12) || ...
            (clip_struct(loop).cvr.bottom - clip_struct(loop).cvr.top + 1) < (32*2+12),
        new_issue = sprintf('Error: Common Valid Region (CVR) too small for models.\n');
        message = cc_append_message(message, clip_name, new_issue);
        status = 2;
    end

    % look for low priority problems.
    if clip_struct(loop).fps ~= 30 && clip_struct(loop).fps ~= 29.97 && ...
            clip_struct(loop).fps ~= 60 && clip_struct(loop).fps ~= 59.94 && ...
            clip_struct(loop).fps ~= 50 && ...
            clip_struct(loop).fps ~= 24 && ......
            clip_struct(loop).fps ~= 25 && clip_struct(loop).fps ~= 15 && ...
            clip_struct(loop).fps ~= 10,
        
        if (clip_struct(loop).fps <= 30 && clip_struct(loop).fps >= 29.95) || ...
                (clip_struct(loop).fps <= 60 && clip_struct(loop).fps >= 59.93) || ...
                (clip_struct(loop).fps <= 24 && clip_struct(loop).fps >= 23.97),
            % do ignore if close to 30, close to 24, or close to 60
        else
            new_issue = sprintf('Warning:  Unusual Frame Rate may cause problems.\nfps expected to be 60, 59.94, 50, 30, 29.97, 25, 24, 23.98, 15, or 10.\n');
            message = cc_append_message(message, clip_name, new_issue);
            if status == 0, status = 3; end;
        end
            
    end
    if clip_struct(loop).spatial.horizontal < -30 || ...
            clip_struct(loop).spatial.horizontal > 30 || ...
            clip_struct(loop).spatial.vertical < -30 || ...
            clip_struct(loop).spatial.vertical > 30,
        new_issue = sprintf('Warning:  Spatial registration values exteremely large.\nExpecting values between -30 and 30.  Check validity.\n');
        message = cc_append_message(message, clip_name, new_issue);
        if status == 0, status = 3; end;
    end
    if clip_struct(loop).luminance_gain < 0.6 || clip_struct(loop).luminance_gain > 1.6,
        new_issue = sprintf('Warning:  Luminance gain suspicous.\nExpecting values between 0.6 and 1.6.  Check validity.\n');
        message = cc_append_message(message, clip_name, new_issue);
        if status == 0, status = 3; end;
    end
    if clip_struct(loop).luminance_offset < -80 || clip_struct(loop).luminance_offset > 80,
        new_issue = sprintf('Warning:  Luminance offset suspicous.\nExpecting values between -80 and 80.  Check validity.\n');
        message = cc_append_message(message, clip_name, new_issue);
        if status == 0, status = 3; end;
    end
    if clip_struct(loop).scale.horizontal < 900 || clip_struct(loop).scale.horizontal > 1100,
        new_issue = sprintf('Warning:  Horizontal spatial scaling suspicious.\nExpecting values between 900 and 1100.  Check validity.\n');
        message = cc_append_message(message, clip_name, new_issue);
        if status == 0, status = 3; end;
    end
    if clip_struct(loop).scale.vertical < 900 && clip_struct(loop).scale.vertical > 1100,
        new_issue = sprintf('Warning:  Vertical spatial scaling suspicious.\nExpecting values between 900 and 1100.  Check validity.\n');
        message = cc_append_message(message, clip_name, new_issue);
        if status == 0, status = 3; end;
    end
    if ~isnan(clip_struct(loop).align_start) && ~isnan(clip_struct(loop).align_stop),
        if mod(clip_struct(loop).align_start*2,2) || mod(clip_struct(loop).align_stop*2,2),
            new_issue = sprintf('Warning: 0.5 fractions for align start and align stop only valid for reduced reference calibration after initial temporal registration and before spatial registration\n');
            message = cc_append_message(message, clip_name, new_issue);
            if status == 0, status = 3; end;
        end
    end
    
    % stop looping if found problems.
    if status == 1 || status == 2,
        break;
    end
end

if status == 2,
    message = [message  sprintf('\nOther clips not checked\n')];
end

if status,
    bool = 0;
    if verbose,
        fprintf('%s', message);
    end
    if status == 1 || status == 2,
        return;
    end
end

if verbose,
    message = [message sprintf('Clip values appear reasonable.\n\n')];
end

% look for larger structural problems
first_time = 1;
all_offsets = sort_clips_by('scene', clip_struct, test_struct);
for loop = 1:length(all_offsets),
    % all of these clips have the same source scene.  Find the original.
    offsets = all_offsets{loop};

    clip_name = sprintf('%s:%s(original)', clip_struct(offsets(1)).test{1}, ...
        clip_struct(offsets(1)).scene{1});

    % make sure original exists for clip.
    if ~strcmp(clip_struct(offsets(1)).hrc{1},'original'),
        status = 2;
        message = [message sprintf('Error:  %s\nOriginal clip missing for this scene.\n', clip_name)];
        break;
    elseif length(offsets) > 1 && strcmp(clip_struct(offsets(2)).hrc{1},'original'),
        status = 2;
        message = [message sprintf('Error:  Clips %d & %d, %s\nTwo or more original clips defined for one scene\n', ...
            offsets(1), offsets(2), clip_name)];
        break;
    elseif length(offsets) == 1,
        if status == 0, status = 3; end;
        message = [message sprintf('Warning:  %s\nProcessed clips missing for this scene.\n', clip_name)];
    end
    
    % find length of the original aligned segment.  If not defined, go to
    % the next clip
    first_length = clip_struct(offsets(1)).align_stop - ...
                clip_struct(offsets(1)).align_start;
    if isnan(clip_struct(offsets(1)).align_start) || ...
            isnan(clip_struct(offsets(1)).align_stop),
        continue;
    end

    % Now, check that other clips are of the same length.
    for cnt = 2:length(offsets),
        clip_name = sprintf('Clip %d: %s:%s(%s) File: %s', ...
            offsets(cnt), clip_struct(offsets(cnt)).test{1}, ...
            clip_struct(offsets(cnt)).scene{1}, clip_struct(offsets(cnt)).hrc{1}, ...
            clip_struct(offsets(cnt)).file_name{1});
        
        % If aligned segment not defined, skip this one.
        if isnan(clip_struct(offsets(cnt)).align_start) || ...
                isnan(clip_struct(offsets(cnt)).align_stop),
            continue;
        end
        
        % find the length of this clip's aligned segment.  Account for
        % reframing, if this clip reframes.
        this_length = clip_struct(offsets(cnt)).align_stop - ...
            clip_struct(offsets(cnt)).align_start;
        if mod(clip_struct(offsets(cnt)).align_stop*2,2), % 0.5 for reframing
            ; % reframing, but length should already be equal -- 0.5 on start & stop
        elseif is_reframing_indicated(clip_struct(offsets(cnt))),
            this_length = this_length - 1;
        end
        
        % make sure length matches previously seen.
        if this_length ~= first_length,
            if first_time,
                message = [message sprintf('Aligned segment of each processed clip must match length of aligned original segment.\n')];
                first_time = 0;
            end
            
            message = [message sprintf('\nError:  %s\nProcessed clip is %g frames long, original (#%d) is %g frames long\n', ...
                clip_name, this_length, offsets(1), first_length)];
            if is_reframing_indicated(clip_struct(offsets(cnt))),
                message = [message sprintf('Re-framing reduced the processed clip''s aligned length by one (1) frame.\n')];
            end
            status = 2;
        end
    end
end

if status,
    bool = 0;
    if verbose,
        fprintf('%s', message);
    end
    if status == 1 || status == 2,
        return;
    end
end

if verbose && status == 0,
    message = [message sprintf('Clips in each test appears internally consistent.\n')];
end

first_time = 1;
all_offsets = sort_clips_by('scene', clip_struct, test_struct);
for loop = 1:length(all_offsets),
    % all of these clips have the same source scene.  Find the original.
    offsets = all_offsets{loop};

    orig_clip_name = sprintf('Clip %d: %s:%s(%s) File: %s', ...
        offsets(1), clip_struct(offsets(1)).test{1}, ...
        clip_struct(offsets(1)).scene{1}, clip_struct(offsets(1)).hrc{1}, ...
        clip_struct(offsets(1)).file_name{1});
        

    % Now, check that other clips are of the same length.
    for cnt = 2:length(offsets),
        clip_name = sprintf('Clip %d: %s:%s(%s) File: %s', ...
            offsets(cnt), clip_struct(offsets(cnt)).test{1}, ...
            clip_struct(offsets(cnt)).scene{1}, clip_struct(offsets(cnt)).hrc{1}, ...
            clip_struct(offsets(cnt)).file_name{1});
        
        % check common valid region size
        if clip_struct(offsets(cnt)).cvr.top ~= clip_struct(offsets(1)).cvr.top || ...
                clip_struct(offsets(cnt)).cvr.left ~= clip_struct(offsets(1)).cvr.left || ...
                clip_struct(offsets(cnt)).cvr.bottom ~= clip_struct(offsets(1)).cvr.bottom || ...
                clip_struct(offsets(cnt)).cvr.right ~= clip_struct(offsets(1)).cvr.right,
            status = 2;
            message = [message sprintf('Original %s\nProcessed ', orig_clip_name, clip_name)];
            message = [message sprintf('Error:  Common Valid Regions differ.  Valid Region must be identical')];
        end
    end
end

if verbose,
    message = [message sprintf('Clips common valid regions appear internally consistent.\n')];
end

if verbose,
    fprintf('%s', message);
end


function message = cc_append_message(message, clip_name, new_issue);

message = [message sprintf('%s\n', clip_name) sprintf('%s\n', new_issue) ];
