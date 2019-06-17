function status = feature_frame_clipping(test_structs, clip_structs,...
    feature_base_dir, varargin)
% FEATURE_LOOP_CLIPPING
%  Loop through a list of clips, compute percent of clipped pixels.
% SYNTAX
%  feature_frame_clipping(test_structs, clip_structs,feature_base_dir)
%  feature_frame_clipping(...,'PropertyName',...);
% DESCRIPTION
%  Compute feature indicating percent of clipped pixels.
% 
%  This function takes variable 'test_structs' (of the same format as
%  GTests), which describes each video test and the location of the
%  associated video; and variable 'clip_structs' (of the same format as
%  GClips), which lists the clips to be processed. 'feature_base_dir' is
%  the path to the directory where the feature's sub-directory should be
%  created.  That sub-directory's name will be this feature's name.
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

feature_name_white = standard_feature_names('Y', 'fullimage',...
        'block_stat','fraction', 'specific','clipwhite');
feature_name_black = standard_feature_names('Y', 'fullimage',...
        'block_stat','fraction', 'specific','clipblack');

        
for loop = 1:max(size(clip_structs)),
    try,
        t = clock;
        [name_exists_w] = write_feature_standard(0,feature_base_dir,feature_name_white,clip_structs(loop),'exists');
        [name_exists_b] = write_feature_standard(0,feature_base_dir,feature_name_black,clip_structs(loop),'exists');
        if name_exists_w && name_exists_b,
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
        
        % Compute the default adjusted SROI and number of blocks available, so
        % can allocate memory to hold all features for this clip.
        tslice_length_sec = 1 / clip_structs(loop).fps; % one image at a time, exactly
        [sroi,vert,horiz] = adjust_requested_sroi (clip_structs(loop));
        number_tslices = total_tslices(clip_structs(loop),tslice_length_sec);
        data_white = zeros(1, 1, number_tslices);
        data_black = zeros(1, 1, number_tslices);
        
        % get rid of gain & offset, if any
        clip_structs(loop).luminance_gain = 1.0;
        clip_structs(loop).luminance_offset = 0.0;
        
        % For each time-slice in this clip, read the time-slice; perceptually filter
        % (none); and compute block statistic.
        [tis_frames] = tslice_conversion(tslice_length_sec, clip_structs(loop).fps);
        for cnt = 1:number_tslices,
            [y] = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec,cnt, ...
                'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
            [rows,cols,time] = size(y);

            y = reshape(y, rows*cols*time,1);

            data_white(:,:,cnt) = length(find(y >= 254)) / (rows*cols*time);
            data_black(:,:,cnt) = length(find(y <= 1)) / (rows*cols*time);
            
            clear y;
        end
        
        % Write features all at once to a single file.
        write_feature_standard(data_white, feature_base_dir,feature_name_white, clip_structs(loop));
        write_feature_standard(data_black, feature_base_dir,feature_name_black, clip_structs(loop));
        
        % Erase the large variable 'data' from memory to prevent problems.
        clear data_black data_white;
    catch
        status = 1;
        if verbose,
            fprintf('\n%s\n', lasterr);
            fprintf('\tSkipping clip.\n'); 
        end
    end
end
