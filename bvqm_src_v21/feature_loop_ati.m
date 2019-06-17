function status = feature_loop_ati(test_structs, clip_structs,vert_size, horiz_size,...
    tslice_length_sec,feature_base_dir, varargin)
% FEATURE_LOOP_ATI
%  Loop through a list of clips, compute ATI features, and save each clip's
%  features in a separate file.
% SYNTAX
%  feature_loop_ati(test_structs, clip_structs,vert_size, horiz_size,tslice_length_sec,feature_base_dir)
%  feature_loop_ati(...,'PropertyName',...);
% DESCRIPTION
%  Comute ATI feature
% 
%  This function takes variable 'test_structs' (of the same format as
%  GTests), which describes each video test and the location of the
%  associated video; and variable 'clip_structs' (of the same format as
%  GClips), which lists the clips to be processed. 'vert_size' and 'horiz_size'
%  specify the block size, horizontally & vertically.  'tslice_length_sec'
%  specifies the length of the block in SECONDS, and 'feature_base_dir' is
%  the path to the directory where the feature's sub-directory should be
%  created.  That sub-directory's name will be this feature's name.
%
%  The following optional parameters are also available.  
%   'verbose'   Print progress report to screen.
%   'quiet'     Don't print anything to the screen (default).
%
%  Return 0 if operated correctly; 1 if an error was encountered.
% NOTES
%  This feature loop function holds one time-slice of images in memory all at
%  once; and also holds in memory the 3D feature matrix associated with the 
%  current video clip.

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

hsize = horiz_size;
vsize = vert_size;
feature_name = standard_feature_names('Y',...
        'vsize',vsize,'hsize',hsize,'tslice_length_sec',tslice_length_sec,...
        'block_stat','std', 'specific','ati');

for loop = 1:max(size(clip_structs)),
    try,
        t = clock;
        [name_exists, file_name] = write_feature_standard(0,feature_base_dir,feature_name,clip_structs(loop),'exists');
        if name_exists,
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
        [sroi,vert,horiz] = adjust_requested_sroi (clip_structs(loop), ...
                            'vsize',vsize,'hsize',hsize);
        number_tslices = total_tslices(clip_structs(loop),tslice_length_sec);
        data = zeros(vert, horiz, number_tslices);
        
        % error check.
        [tslice_length_frames] = tslice_conversion(tslice_length_sec, clip_structs(loop).fps);
        if tslice_length_frames == 1,
            error('Time slice length must encompars two or more frames for all clips.\nCondition failed for clip %d -- %s:%s(%s)', ...
                loop, clip_structs(loop).test{1}, clip_structs(loop).scene{1}, ...
                clip_structs(loop).hrc{1});
        end
         
        % For each time-slice in this clip, read the time-slice; perceptually filter
        % (none); and compute block statistic.
        for cnt = 1:number_tslices,
            y = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec,cnt, ...
                'vsize',vsize,'hsize',hsize, 'TI', ...
                'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
            ati = filter_ati(y);
            data(:,:,cnt) = block_statistic (ati, vsize, hsize, 'std');
            
            clear ati y;
        end
        
        % Write features all at once to a single file.
        write_feature_standard(data, feature_base_dir,feature_name, clip_structs(loop));
        
        % Erase the large variable 'data' from memory to prevent problems.
        clear data;
    catch
        status = 1;
        if verbose,
            fprintf('\n%s\n', lasterr);
            fprintf('\tSkipping clip.\n');
        end
    end
end
