function [new_clip_structs, status] = temporal_registration(test_structs, clip_structs,varargin);
% TEMPORAL_REGISTRATION 
%   Automatic computation of temporal registration.
%   2004 ITU standard, General Model FRTV Calibration
% SYNTAX
%  [new_clip_structs, status] = temporal_registration(test_structs, clip_structs);
%  [...] = temporal_registration(...,'PropertyName',...)
% DESCRIPTION
%  Perform temporal registration on each processed clip in 'clip_structs'.
%  Frame based algorithm.
%
%  Optional properties:
%   'Uncertainty',value,     Assume temporal registration uncertainty of
%                                   'value' seconds.  By default, 1/2 second.
%   'unaligned',                Start with temporally unaligned video sequences.
%                                   This is the default behavior.
%   'aligned',                  Start with temporally ALIGNED video sequences.
%                                   Useful when 'first frames align' assumption is very 
%                                   wrong.
%   'verbose',                  Print status information to screen
%   'quiet',                    Print no information to screen.  Default behavior.
%   'log',file1,file2,          Create log files.  Echo to 'file1' the
%                                   default screen summary.  Write
%                                   histograms ('verbose') to 'file2'.
%
%  Return the new (temorally aligned) clip structure.  Also return status,
%  which has the following fields defined:
%   'status.error'              0 if okay; 1 if a fatal error was encountered.
%                               2 if original is missing for a scene.  
%                               3 if multiple original sequences are
%                               defined for one scene.  If status.error is
%                               not 0, then the other structure fields will 
%                               be unavailable.  See also function "lasterr"
%   'status.uncertainty'        Number of video clips that should be re-run
%                               with a larger uncertainty.
%   'status.ambiguous'          Number of video clips with ambiguous alignments
%   'status.still'              Number of video clips that came from still scenes

try
    % initialize status.
    status.error = 1;
    status.uncertainty = 0;
    status.ambiguous = 0;
    status.still = 0;
    
    % handle optional properties.
    frequency = 1.0;
    uncert_sec = 0.5;
    verbose = 0;
    flag = 'unaligned';
    log = 0;

    cnt = 1;
    while cnt <= nargin-2,
        if strcmp(lower(varargin{cnt}),'uncertainty'),
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
        elseif strcmp(lower(varargin{cnt}),'log'),
            if nargin-1 < cnt+2,
                error('Specify log file names');
            end
            log = 1;
            logfile1 = fopen(varargin{cnt+1},'w');
            if logfile1 == -1,
                error('Could not open log file "%s"', varargin{cnt+1});
            end
            logfile2 = fopen(varargin{cnt+2},'w');
            if logfile2 == -2,
                error('Could not open log file "%s"', varargin{cnt+2});
            end
            cnt = cnt + 3;
        else
            error('Property Name not recognized.  Aborting.');
        end
    end

    % sort clips by scene.
    offsets = sort_clips_by('scene', clip_structs, test_structs);


    % set subsampling size & other constants.
    hsize = 16;
    vsize = 16;
    STILL_THRESHOLD = 0.002;
    HFW = 3;
    BELOW_WARN = 0.9;
    DELTA = 4;


    % replicate clip structure.
    new_clip_structs = clip_structs;

    % keep track of statistics
    uncert_too_small = 0;
    align_ambiguous = 0;
    still_scene = 0;

    % Loop through each scene.
    for cnt = 1:length(offsets),
        curr_offsets = offsets{cnt};
        clip = curr_offsets(1);

        t = clock;
        if verbose,
            fprintf('Scene %d of %d ==> %s:%s at %d:%d\n', cnt, length(offsets), ...
                clip_structs(clip).test{1}, clip_structs(clip).scene{1}, t(4), t(5) );
        end
        if log,
            fprintf(logfile1,'Scene %d of %d ==> %s:%s\n', cnt, length(offsets), ...
                clip_structs(clip).test{1}, clip_structs(clip).scene{1});
        end

        % Find the offset of the test structure for this clip.
        tnum = search_test_list(test_structs, clip_structs(clip));

        % Compute the default adjusted SROI and number of blocks available
        [sroi,vert,horiz] = adjust_requested_sroi (clip_structs(clip), ...
            'vsize',vsize,'hsize',hsize, 'evenodd');
        tslice_length_sec = 1.0 / clip_structs(clip).fps;
        src_number_tslices = total_tslices(clip_structs(clip),tslice_length_sec, flag);

        % Check if this is an interlace or progressive sequence
        % Allocate memory to hold all features for this clip.
        if strcmp(clip_structs(clip).video_standard,'progressive'),
            is_progressive = 1;
            src = zeros(vert*horiz, src_number_tslices);
        elseif strcmp(clip_structs(clip).video_standard,'interlace_lower_field_first') | strcmp(clip_structs(clip).video_standard,'interlace_upper_field_first'),
            is_progressive = 0;
            src = zeros(vert*horiz, src_number_tslices*2);
        else
            error('video standard not recognized');
        end
        
        % make sure an original sequence is defined.
        if strcmp(clip_structs(clip).hrc{1},'original') ~= 1,
            status.error = 2;
            error(sprintf('Original video sequence missing for scene "%s"', clip_structs(clip).scene{1}), 'error');
        end
        % make sure an original sequence is unique.
        if length(curr_offsets) == 1,
            % no HRC defined.  Silly but okay.
            if verbose,
                warning(sprintf('No processed video sequences defined for scene "%s"', clip_structs(clip).scene{1}), 'warning');
            end
            continue;
        end
        % make sure an original sequence is unique.
        if strcmp(clip_structs(curr_offsets(2)).hrc{1},'original') == 1,
            status.error = 3;
            error(sprintf('Multiple original video sequence defined for scene "%s"', clip_structs(clip).scene{1}), 'error');
        end

        % read in each original field (or frame), and sub-sample it.
        num = 0;
        for frame=1:src_number_tslices,
            y = read_tslice( test_structs(tnum), clip_structs(clip), tslice_length_sec, frame,flag, ...
                'sroi', sroi.top, sroi.left, sroi.bottom, sroi.right);
            if is_progressive,
                src_one(:,frame) = reshape( block_statistic(y,hsize,vsize,'mean'), vert*horiz, 1);
            else
                [one,two] = split_into_fields(y);

                src_one(:,frame) = reshape( block_statistic(one,hsize/2,vsize,'mean'), vert*horiz, 1);
                src_two(:,frame) = reshape( block_statistic(two,hsize/2,vsize,'mean'), vert*horiz, 1);
            end
        end

        % normalize original images.  Don't amplify noise.
        for frame=1:src_number_tslices,
            src_one(:,frame) = src_one(:,frame) ./ max( 1.0, std(src_one(:,frame)));
            if ~is_progressive,
                src_two(:,frame) = src_two(:,frame) ./ max( 1.0, std(src_two(:,frame)));
            end
        end

        % loop through each processed version of this scene.
        for loop = 2:length(curr_offsets),
            clip = curr_offsets(loop);

            if verbose,
                fprintf('\tClip %d of %d ==> %s:%s(%s)', loop-1, length(curr_offsets)-1, ...
                    clip_structs(clip).test{1}, clip_structs(clip).scene{1}, clip_structs(clip).hrc{1} );
            end
            if log,
                fprintf(logfile1, '\tClip %d of %d ==> %s:%s(%s)', loop-1, length(curr_offsets)-1, ...
                    clip_structs(clip).test{1}, clip_structs(clip).scene{1}, clip_structs(clip).hrc{1} );
                fprintf(logfile2, '\n\n------------------------------------------------------------------------------\n\n%s:%s(%s)', ...
                    clip_structs(clip).test{1}, clip_structs(clip).scene{1}, clip_structs(clip).hrc{1} );
            end

            % read in processed fields (or frames), at the requested interval, and sub-sample them.
            number_tslices = total_tslices(clip_structs(clip),tslice_length_sec, flag);
            all_align = [];

            uncert = round(clip_structs(clip).fps * uncert_sec);
            all_std_of_diff = zeros(1,uncert*2+1);

            for frame=1+uncert:number_tslices-uncert,
                % stop if run out of source images.
                if frame + uncert > src_number_tslices,
                    break;
                end

                % read in next frame
                y = read_tslice( test_structs(tnum), clip_structs(clip), tslice_length_sec, frame,flag, ...
                    'sroi', sroi.top, sroi.left, sroi.bottom, sroi.right);

                % split into fields if needed.  Record number /offset relative
                % to source array of images.
                if is_progressive,
                    one = reshape( block_statistic(y,hsize,vsize,'mean'), vert*horiz, 1);
                else
                    [one,two] = split_into_fields(y);
                    one = reshape( block_statistic(one,hsize/2,vsize,'mean'), vert*horiz, 1);
                    two = reshape( block_statistic(two,hsize/2,vsize,'mean'), vert*horiz, 1);
                end


                % normalize & find temporal registration for this image.
                one =one ./ max(1.0, std(one));
                [still, src_num, std_of_diff] = temporal_registration_search(src_one, one, frame, uncert, STILL_THRESHOLD);
                if ~still,
                    all_align = [ all_align (frame-src_num) ];
                    all_std_of_diff = all_std_of_diff + std_of_diff;
                end

                if ~is_progressive,
                    % normalize & find temporal registration for this image.
                    two = two ./ max(1.0, std(two));
                    [still, src_num, std_of_diff] = temporal_registration_search(src_two, two, frame, uncert, STILL_THRESHOLD);
                    if ~still,
                        all_align = [ all_align (frame-src_num) ];
                        all_std_of_diff = all_std_of_diff + std_of_diff;
                    end
                    % fprintf('\tframe %d field1 offset %d; field2 offset %d\n', frame, all_align(length(all_align)-1), frame-src_num);
                end
            end

            % still check.  If still scene, record default alignment & go to
            % next processed clip.  Take mean of all_std_of_diff AFTER the
            % min/max, so that you only need to divide once (not 2*uncert+1 times).
            if length(all_align) == 0 | (max(all_std_of_diff) - min(all_std_of_diff)) / length(all_align) < STILL_THRESHOLD,
                if verbose,
                    fprintf('\tStill Scene\n');
                end
                if log,
                    fprintf(logfile1,'\tStill Scene\n');
                end
                clip_delay(loop) = 0;
                still_scene = still_scene + 1;
                continue;
            end

            % histogram, 
            hist_unsmooth = histc(all_align,-uncert:uncert);

            % smooth
            k=0:(2*HFW);
            fk = 0.5 + 0.5 * cos(pi * (k - HFW) / (1+HFW));
            fk = fk ./ sum(fk);
            hist_smooth = conv(fk,hist_unsmooth);

            % take off extra length, so left only with smoothed portion,
            % length is now (uncert-HFW)*2+1
            hist_smooth = hist_smooth((HFW*2+1):(length(hist_smooth) - HFW*2));

            %
            max_H_value = max(hist_unsmooth);
            [max_SH_value,max_SH_offset] = max(hist_smooth);
            clip_delay(loop) = uncert + 1 - (max_SH_offset + HFW);

            % check if uncert was large enough.
            bins_missing = [1:HFW (length(hist_unsmooth) - HFW):length(hist_unsmooth)];
            if max(hist_unsmooth(bins_missing)) > max_H_value * BELOW_WARN,
                if verbose,
                    fprintf('\tUncertanty too small.');
                end
                if log,
                    fprintf(logfile1, '\tUncertanty too small.');
                    fprintf(logfile2, '\tUncertanty too small.');
                end
                uncert_too_small = uncert_too_small + 1;
            end
            bins_away = [1:(max_SH_offset-DELTA-1) (max_SH_offset + DELTA+1):length(hist_smooth)];
            if hist_smooth(bins_away) > max_SH_value * BELOW_WARN,
                if verbose,
                    fprintf('\ttemporal registration ambiguous');
                end
                if log,
                    fprintf(logfile1,'\ttemporal registration ambiguous');
                    fprintf(logfile2,'\ttemporal registration ambiguous');
                end
                align_ambiguous = align_ambiguous + 1;
            end

            if log,
                fprintf(logfile2,'\tNew delay %d\n\n', max_SH_offset+HFW-uncert-1);
                fprintf(logfile2,'Histogram\ndelay, votes, smoothed votes\n');
                temp = zeros(2*uncert+1,1);
                temp(:) = NaN;
                temp((HFW+1):(length(temp)-HFW)) = hist_smooth;
                for i=1:2*uncert+1,
                    if i == max_SH_offset+HFW,
                        fprintf(logfile2,'*');
                    end
                    fprintf(logfile2,'\t%d\t%d\t%f\n', i-uncert-1, hist_unsmooth(i), temp(i));
                end
                fprintf(logfile1, '\tNew delay %d\n', max_SH_offset+HFW-uncert-1);
            end
            if verbose,
                fprintf('\tNew delay %d\n', max_SH_offset+HFW-uncert-1);
            end
        end

        % 
        clip_delay(1) = 0;
        orig_extra_start = max(clip_delay(2:length(curr_offsets)));
        clip = curr_offsets(1);
        if strcmp(flag,'aligned'),
            new_clip_structs(clip).align_start = new_clip_structs(clip).align_start + orig_extra_start;
        else
            new_clip_structs(clip).align_start = new_clip_structs(clip).loc_start + orig_extra_start;
        end
        clip_delay = orig_extra_start - clip_delay;
        
        %%% problem = clip_delay is supposed to be all >= 0, and original
        %%% here is -1; could have been less? check other clips to be sure!
        while min(clip_delay) < 0,
            clip_delay = clip_delay + 1;
        end

        % record overlapping range for all clips of this scene.
        % Cut off the beginning of the clip, to cause clips to align
        for loop = 1:length(curr_offsets),
            clip = curr_offsets(loop);
            if strcmp(flag,'aligned'),
                new_clip_structs(clip).align_start = clip_structs(clip).align_start + clip_delay(loop);
                new_clip_structs(clip).align_stop = clip_structs(clip).align_stop + clip_delay(loop);
            else % unaligned
                new_clip_structs(clip).align_start = clip_structs(clip).loc_start + clip_delay(loop);
                new_clip_structs(clip).align_stop = clip_structs(clip).loc_stop + clip_delay(loop);
            end
            is_length(loop) = new_clip_structs(clip).loc_stop - new_clip_structs(clip).align_start;
        end

        % shorten clip as needed by file length
        use_length = min(is_length);
        for loop = 1:length(curr_offsets),
            clip = curr_offsets(loop);
            new_clip_structs(clip).align_stop = new_clip_structs(clip).align_start + use_length;
        end

        % clear variables
        clear src_one src_two extra_start extra_stop;
        clear one two clip curr_offsets loop too_few_start too_few_stop;
        clear clip_delay number_tslices src_number_tslices std_of_diff all_std_of_diff;

    end

    status.error = 0;
    status.uncertainty = uncert_too_small;
    status.ambiguous = align_ambiguous;
    status.still = still_scene;

    if verbose,
        fprintf('\n\n%d clips should be re-run with a larger uncertainty\n', uncert_too_small);
        fprintf('%d clips had amibugous alignments\n', align_ambiguous);
        fprintf('%d clips came from still scenes\n', still_scene);
    end
    if log,
        fprintf(logfile1, '\n\n%d clips should be re-run with a larger uncertainty\n', uncert_too_small);
        fprintf(logfile1, '%d clips had amibugous alignments\n', align_ambiguous);
        fprintf(logfile1, '%d clips came from still scenes\n', still_scene);
    end

    if strcmp(flag,'aligned'),
        if verbose,
            fprintf('\n\nWarning:  segments selected will imperfectly match those originally specified.\n');
            fprintf('A manual check is recommended, if the intent was to select an exact portion\n');
            fprintf('of the available video.\n');
        end
        if log,
            fprintf(logfile1, '\n\nWarning:  segments selected will imperfectly match those originally specified.\n');
            fprintf(logfile1, 'A manual check is recommended, if the intent was to select an exact portion\n');
            fprintf(logfile1, 'of the available video.\n');
        end
    end
    
    new_clip_structs = fix_temporal(test_structs, new_clip_structs, 'EndPoint');


    if log,
        fclose(logfile1);
        fclose(logfile2);
    end
catch
    if status.error == 0,
        status.error = 1;
    end
    new_clip_structs = [];
    if verbose,
        fprintf('\n%s\n', lasterr);
        fprintf('\tAborting.\n'); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [still, src_num, std_of_diff] = temporal_registration_search(src, deg, num, uncert, still_threshold);
%

% compare processed & original images.
std_of_diff = std(src(:,(num-uncert):(num+uncert)) - repmat(deg,1,2*uncert+1));

% find minimum stdev of difference.
[a,where] = min(std_of_diff);

% detect whether still
if max(std_of_diff) - min(std_of_diff) < still_threshold,
    still = 1;
else
    still = 0;
end

src_num = (where - uncert - 1) + num;