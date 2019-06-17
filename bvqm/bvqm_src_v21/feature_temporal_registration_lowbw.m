function status = feature_temporal_registration_sequence(test_structs, clip_structs,...
    feature_base_dir, varargin)
% FEATURE_TEMPORAL_REGISTRATION_SEQUENCE
%  Loop through a list of clips, compute field-based/frame-based metrics used by the
%  sequence-based temporal registration algorithm, and save each clip's
%  features in a separate file. 
% SYNTAX
%  feature_temporal_registration_sequence(test_structs, clip_structs,feature_base_dir)
%  feature_temporal_registration_sequence(...,'PropertyName',...);
% DESCRIPTION
%  Compute sequence-based temporal registration feature
% 
%  This function takes variable 'test_structs' (of the same format as
%  GTests), which describes each video test and the location of the
%  associated video; and variable 'clip_structs' (of the same format as
%  GClips), which lists the clips to be processed. 'feature_base_dir' is
%  the path to the directory where the feature's sub-directory should be
%  created.  
%
%  The following optional parameters are also available.  
%   'verbose'   Print progress report to screen.
%   'quiet'     Don't print anything to the screen (default).
%
%  Return 0 if operated correctly; 1 if an error was encountered.

status = 0;

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

feature_name_y_mean = standard_feature_names('Y', 'fullfield',...
        'block_stat','mean', 'specific','cont');
feature_name_ti2 = standard_feature_names('Y', 'fullfield',...
        'block_stat','rms', 'specific','cont', 'specific','ti2');
feature_name_ti10 = standard_feature_names('Y', 'fullfield',...
        'block_stat','rms', 'specific','cont', 'specific','ti10');

        
for loop = 1:max(size(clip_structs)),
    try,
        t = clock;
        [name_exists_1] = write_feature_standard(0,feature_base_dir,feature_name_y_mean,clip_structs(loop),'exists');
        [name_exists_2] = write_feature_standard(0,feature_base_dir,feature_name_ti2,clip_structs(loop),'exists');
        [name_exists_4] = write_feature_standard(0,feature_base_dir,feature_name_ti10,clip_structs(loop),'exists');
        if name_exists_1 && name_exists_2 && name_exists_4,
            if verbose,
                fprintf('\tfeature already computed for clip %s:%s(%s)\n',...
                    clip_structs(loop).test{1}, clip_structs(loop).scene{1}, ...
                    clip_structs(loop).hrc{1});
            end
            continue;
        end
        if verbose,
            fprintf('Clip %d of %d ==> %s:%s(%s) at %d:%d\n', loop, max(size(clip_structs)), ...
                clip_structs(loop).test{1}, clip_structs(loop).scene{1}, ...
                clip_structs(loop).hrc{1}, t(4), t(5) );
        end
	
        % Find the offset of the test structure for this clip.
        tnum = search_test_list(test_structs, clip_structs(loop));
        tis_frames = 1;
        tslice_length_sec = 1 / clip_structs(loop).fps;
        
        % Compute the default adjusted SROI and number of blocks available, so
        % can allocate memory to hold all features for this clip.
        [sroi,vert,horiz] = adjust_requested_sroi (clip_structs(loop));
        number_tslices = total_tslices(clip_structs(loop),tslice_length_sec);
        if ~strcmp(clip_structs(loop).video_standard,'progressive'),
            time = number_tslices * 2;
            data_y = zeros(1, 1, time);
            data_ti2 = zeros(1, 1, time-2);
            data_ti10 = zeros(1, 1, time-10);
            
            % For each time-slice in this clip, read the time-slice; perceptually filter
            % (none); and compute block statistic.
            buffer = [];
            curr = 1;
            for cnt = 1:number_tslices,
                [y] = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec,cnt, ...
                    'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
                if strcmp(clip_structs(loop).video_standard,'interlace_lower_field_first'),
                    [one,two] = split_into_fields(y);
                else
                    [two,one] = split_into_fields(y);
                end
                [rows,cols] = size(one);

                [data_y_mean(:,:,cnt*2-1)]  = block_statistic (one, rows, cols, 'mean');
                [data_y_mean(:,:,cnt*2)]  = block_statistic (two, rows, cols, 'mean');

                if cnt >= 2,
                    use = curr - 2;
                    if use < 1,
                        use = use + 10;
                    end
                    [data_ti2(:,:,cnt*2-3)]  = block_statistic (one - buffer(:,:,use), rows, cols, 'rms');
                    [data_ti2(:,:,cnt*2-2)]  = block_statistic (two - buffer(:,:,use+1), rows, cols, 'rms');
                end

                if cnt >= 6,
                    use = curr;
                    [data_ti10(:,:,cnt*2-11)]  = block_statistic (one - buffer(:,:,use), rows, cols, 'rms');
                    [data_ti10(:,:,cnt*2-10)]  = block_statistic (two - buffer(:,:,use+1), rows, cols, 'rms');
                end

                buffer(:,:,curr) = one;
                buffer(:,:,curr+1) = two;
                curr = curr + 2;
                if curr == 11,
                    curr = 1;
                end
            end
        else
            time = number_tslices;
            data_y = zeros(1, 1, time);
            data_ti2 = zeros(1, 1, time-1);
            data_ti10 = zeros(1, 1, time-5);
            % For each time-slice in this clip, read the time-slice; perceptually filter
            % (none); and compute block statistic.
            buffer = [];
            curr = 1;
            for cnt = 1:number_tslices,
                [y] = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec,cnt, ...
                    'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
                [rows,cols] = size(y);

                [data_y_mean(:,:,cnt)]  = block_statistic (y, rows, cols, 'mean');

                if cnt >= 2,
                    use = curr - 1;
                    if use < 1,
                        use = use + 5;
                    end
                    [data_ti2(:,:,cnt)]  = block_statistic (y - buffer(:,:,use), rows, cols, 'rms');
                end

                if cnt >= 6,
                    use = curr;
                    [data_ti10(:,:,cnt)]  = block_statistic (y - buffer(:,:,use), rows, cols, 'rms');
                end

                buffer(:,:,curr) = y;
                curr = curr + 1;
                if curr == 6,
                    curr = 1;
                end
            end
        end
        
        
        % Write features all at once to a single file.
        write_feature_standard(data_y_mean, feature_base_dir,feature_name_y_mean, clip_structs(loop));
        write_feature_standard(data_ti2, feature_base_dir,feature_name_ti2, clip_structs(loop));
        write_feature_standard(data_ti10, feature_base_dir,feature_name_ti10, clip_structs(loop));
        
        % Erase the large variable 'data' from memory to prevent problems.
        clear data_y_mean data_ti2 data_t10;
    catch
        status = 1;
        if verbose,
            fprintf('\n%s\n', lasterr);
            fprintf('\tSkipping clip.\n'); 
        end
    end
end
