function [new_clip_structs, status, unfiltered_clip_structs] = luminance_gain_offset(test_structs, clip_structs, varargin)
% LUMINANCE_GAIN_OFFSET
%   Calculate luminance gain & offset for each video clip.
%   2004 ITU standard, General Model FRTV Calibration
% SYNTAX
%  [new_clip_structs, status, unfiltered_clip_structs] = luminance_gain_offset(test_structs, clip_structs);
%  [...] = luminance_gain_offset(...,'PropertyName',PropertyValue,...)
%  [new_clip_structs, status] = luminance_gain_offset(...,'NoHRCFilter',...);
% DEFINITION
%   Calculate luminance gain & offset.  Spatial registration & CVR must already
%   be correctly set for each processed clip.  Return un-filtered results
%   in 'unfiltered_clip_structs'
%
%   Optional properties:
%   'Frequency',value,       Examine frames at a frequency of 'value',
%                            specified in seconds.  By default, this
%                            will be one second.
%   'Uncertainty',value,     Assume temporal registration uncertainty of
%                            'value' seconds.  By default, 1/2 second.
%   'unaligned',             Used temporally unaligned video sequences.
%                            This is the default behavior.
%   'aligned',               Used temporally ALIGNED video sequences.
%                            Skip unaligned clips.  Useful when
%                            'first frames align' assumption is very wrong.
%   'verbose',               Print results
%   'quiet',                 Print nothing to the screen.  Default behavior.
%   'NoHRCFilter'            Don't median filter luminance gain & offset by
%                            HRC.  This option results in LESS accurate
%                            estimates; the estimate for each clip depends
%                            upon that clip's video, only.
%
% return value 'status' contains the number of video clips for which
% gain/offset results could not be computed.  Those will default to gain=1,
% offset=0;  If no errors are encountered, status == 0;

status = 0;

% handle optional properties.
frequency = 1.0;
uncert_sec = 0.5;
verbose = 0;
hrc_filter = 1;
flag = 'unaligned';
cnt = 1;
while cnt <= nargin-2,
    if strcmp(lower(varargin{cnt}),'frequency'),
        frequency = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmp(lower(varargin{cnt}),'uncertainty'),
        uncert_sec = varargin{cnt+1};
        cnt = cnt + 2;
    elseif strcmp(lower(varargin{cnt}),'verbose'),
        verbose = 1;
        cnt = cnt + 1;
    elseif strcmp(lower(varargin{cnt}),'quiet'),
        verbose = 0;
        cnt = cnt + 1;
    elseif strcmp(lower(varargin{cnt}),'unaligned'),
        flag = 'unaligned';
        cnt = cnt + 1;
    elseif strcmp(lower(varargin{cnt}),'aligned'),
        flag = 'aligned';
        cnt = cnt + 1;
    elseif strcmp(lower(varargin{cnt}),'nohrcfilter'),
        hrc_filter = 0;
        cnt = cnt + 1;
    else
        error('Property Name not recognized.  Aborting.');
    end
end

% sort clips by scene.
offsets = sort_clips_by('scene', clip_structs, test_structs);

% set subsampling size.
hsize = 16;
vsize = 16;

% clear any luminance gain/offset settings.
for i=1:length(clip_structs),
    clip_structs(i).luminance_gain = 1;
    clip_structs(i).luminance_offset = 0;
end

% replicate clip structure.
new_clip_structs = clip_structs;

% Loop through each scene.
for cnt = 1:length(offsets),
    curr_offsets = offsets{cnt};
    clip = curr_offsets(1);
    
    t = clock;
    if verbose,
        fprintf('Scene %d of %d ==> %s:%s at %d:%d\n', cnt, length(offsets), ...
            clip_structs(clip).test{1}, clip_structs(clip).scene{1}, t(4), t(5) );
    end

    % Find the offset of the test structure for this clip.
    tnum = search_test_list(test_structs, clip_structs(clip));
    
    % Compute the default adjusted SROI and number of blocks available
    [sroi,vert,horiz] = adjust_requested_sroi (clip_structs(clip), ...
                        'vsize',vsize,'hsize',hsize, 'evenodd');
    tslice_length_sec = 1.0 / clip_structs(clip).fps;
    number_tslices = total_tslices(clip_structs(clip),tslice_length_sec, flag);
    
    % Check if this is an interlace or progressive sequence
    % Allocate memory to hold all features for this clip.
    if strcmp(clip_structs(clip).video_standard,'progressive'),
        is_progressive = 1;
        src = zeros(vert*horiz, number_tslices);
    elseif strcmp(clip_structs(clip).video_standard,'interlace_lower_field_first') | strcmp(clip_structs(clip).video_standard,'interlace_upper_field_first'),
        is_progressive = 0;
        src = zeros(vert*horiz, number_tslices*2);
    else
        error('video standard not recognized');
    end
    
    % read in each original field (or frame), and sub-sample it.
    num = 0;
    for frame=1:number_tslices,
        y = read_tslice( test_structs(tnum), clip_structs(clip), tslice_length_sec, frame,flag, ...
            'sroi', sroi.top, sroi.left, sroi.bottom, sroi.right);
        if is_progressive,
            num = num + 1;
            src(:,num) = reshape( block_statistic(y,hsize,vsize,'mean'), vert*horiz, 1);
        else
            [one,two] = split_into_fields(y);
            
            num = num + 1;
            src(:,num) = reshape( block_statistic(one,hsize/2,vsize,'mean'), vert*horiz, 1);
            
            num = num + 1;
            src(:,num) = reshape( block_statistic(two,hsize/2,vsize,'mean'), vert*horiz, 1);
        end
    end
    
    % display_xyt(src,'zoom4');
    % pause;
    
    % loop through each processed version of this scene.
    for loop = 2:length(curr_offsets),
        clip = curr_offsets(loop);
        
        % read in processed fields (or frames), at the requested interval, and sub-sample them.
        all_gain = [];
        all_offset = [];
        add_frames = round(clip_structs(clip).fps * frequency);
        uncert = round(clip_structs(clip).fps * uncert_sec);
        number_tslices = total_tslices(clip_structs(clip),tslice_length_sec, flag);
        for frame=(1+add_frames):add_frames:(number_tslices-add_frames),
            y = read_tslice( test_structs(tnum), clip_structs(clip), tslice_length_sec, frame,flag, ...
                'sroi', sroi.top, sroi.left, sroi.bottom, sroi.right);
            if is_progressive,
                one = reshape( block_statistic(y,hsize,vsize,'mean'), vert*horiz, 1);
                num = frame;
            else
                [one,two] = split_into_fields(y);
                one = reshape( block_statistic(one,hsize/2,vsize,'mean'), vert*horiz, 1);
                two = reshape( block_statistic(two,hsize/2,vsize,'mean'), vert*horiz, 1);
                num = frame*2-1;
            end
            
            [valid, y_gain, y_offset] = luminance_gain_offset_search(src, one, num, uncert, verbose);
            if valid,
                all_gain = [all_gain y_gain];
                all_offset = [all_offset y_offset];
            end
            
            if ~is_progressive,
                [valid, y_gain, y_offset] = luminance_gain_offset_search(src, two, num+1, uncert, verbose);
                if valid,
                    all_gain = [all_gain y_gain];
                    all_offset = [all_offset y_offset];
                end
            end
        end
        % median filter & record result for this clip.
        if length(all_gain) > 0,
            new_clip_structs(clip).luminance_gain = median(all_gain);
            new_clip_structs(clip).luminance_offset = median(all_offset);
        else
            new_clip_structs(clip).luminance_gain = 1.0;
            new_clip_structs(clip).luminance_offset = 0.0;
            status = status + 1;
            if verbose,
                fprintf('\tGain&offset calculation failed for clip %s:%s(%s)n', new_clip_structs(clip).test{1}, new_clip_structs(clip).scene{1}, ...
                    new_clip_structs(clip).hrc{1});
            end
        end
        if verbose,
            fprintf('\t%s:%s(%s)\tgain %f, offset %f\n', new_clip_structs(clip).test{1}, new_clip_structs(clip).scene{1}, ...
                new_clip_structs(clip).hrc{1}, new_clip_structs(clip).luminance_gain , new_clip_structs(clip).luminance_offset);
        end
    end
end

% median filter over HRC
if hrc_filter,
    unfiltered_clip_structs = new_clip_structs;
    % sort clips by hrc.
    offsets = sort_clips_by('hrc', clip_structs, test_structs);

    % Loop through each hrc.
    for cnt = 1:length(offsets),
        curr_offsets = offsets{cnt};

        % find median gain&offset
        if strcmp(clip_structs(curr_offsets(1)).hrc,'original'),
            gain = 1.0;
            offset = 0.0;
        else
            gain = median([new_clip_structs(curr_offsets).luminance_gain]);
            offset = median([new_clip_structs(curr_offsets).luminance_offset]);
        end

        % and record it.
        for i=1:length(curr_offsets),
            new_clip_structs(curr_offsets(i)).luminance_gain = gain;
            new_clip_structs(curr_offsets(i)).luminance_offset = offset;
        end
        if verbose,
            fprintf('HRC %s\tgain %f, offset %f\n', new_clip_structs(curr_offsets(1)).hrc{1}, gain, offset);
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [valid, y_gain, y_offset] = luminance_gain_offset_search(src, deg, num, uncert, verbose);
%

[a, src_length] = size(src);
std_of_diff = zeros(1, src_length);
std_of_diff(:) = NaN;

done = 0;
while ~done,
	for i=max(1,num-uncert) : min(src_length,num+uncert),
        % if haven't filled in this delay's temporal registration, do it
        % now.
        if isnan(std_of_diff(i)),
            std_of_diff(i) = std(src(:,i) - deg);
        end
    end
    
    [a,where] = min(std_of_diff);
    
    % if at end point, return failure.
    if where == 1 | where == src_length 
        valid = 0;
        y_gain = 1;
        y_offset = 0;
        return;
    elseif where == num-uncert | where == num+uncert,
        %fprintf('\tat edge of uncertainty.  extending\n');
        uncert = uncert + 15;
    else
        done = 1;
    end
        
end

%  Check added 4/26/10 to make sure that there is some variance in both the
%  src and deg images to prevent failure of the linear fit
if (std(deg)==0 || std(src(:,where))==0)
    valid = 0;
    y_gain = 1;
    y_offset = 0;
    return;
end

% compute initial gain via linear regression
y = deg;
x = [ones(length(y),1) src(:,where)];
b = x\y;
r = y - x*b;

done = 0;
prev_b = b;
while ~done,
	epsilon = 0.1;
	cost = 1.0 ./ (abs(r) + epsilon); % cost vector, reciprocal of errors
	cost = cost ./ sqrt(sum(cost)); % normalize for unity norm
	cost = (cost.^2);
    temp = eye(length(y),length(y));
    for i=1:length(y),
        temp(i,i) = temp(i,i) * cost(i);
    end
    cost = temp;
	
	b = inv(x'*cost*x)*x'*cost*y;
    r = y - x*b;
    
    if abs(prev_b(2) - b(2)) < 0.0001,
        done = 1;
    else
        prev_b = b;
    end
end

valid = 1;
y_gain = b(2);
y_offset = b(1);

if verbose,
    fprintf('\t\talign to %d, gain %f & offset %f \n', where, y_gain, y_offset);
end


