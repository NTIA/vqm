function status = feature_loop_coherent_color(test_structs, clip_structs,vert_size,horiz_size,...
    tslice_length_sec,feature_base_dir,varargin)
% FEATURE_LOOP_COHERENT_COLOR
%  Loop through a list of clips, compute coherent color features, and save each clip's
%  features in a separate file.
% SYNTAX
%  feature_loop_coherent_color(test_structs, clip_structs,vert_size,horiz_size,tslice_length_sec,feature_base_dir)
%  feature_loop_coherent_color(...,'PropertyName',...)
% DESCRIPTION
%  Comute coherent color features
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
%  The following optional arguments may be specified: 
%   'BlockStd'      Compute feature using block standard deviation, as well
%                   as the usualmean of each block.
%   'verbose'   Print progress report to screen.
%   'quiet'     Don't print anything to the screen (default).
%
%  Return 0 if operated correctly; 1 if an error was encountered.
% NOTES
%  This feature loop function holds one time-slice of images in memory all at
%  once; and also holds in memory the 3D feature matrix associated with the 
%  current video clip.


status = 0;

optional_std = 0;
verbose = 0;

for cnt = 1:length(varargin),
    if strcmp(lower(varargin{cnt}),'blockstd'),
        optional_std = 1;
    elseif strcmp(lower(varargin{cnt}),'verbose'),
        verbose = 1;
    elseif strcmp(lower(varargin{cnt}),'quiet'),
        verbose = 0;
    else
        error('optional property not recognized');
    end
end

hsize = horiz_size;
vsize = vert_size;
feature_name1 = standard_feature_names('color',...
        'vsize',vsize,'hsize',hsize,'tslice_length_sec',tslice_length_sec,...
        'block_stat','mean', 'specific','coher_cb');
feature_name2 = standard_feature_names('color',...
        'vsize',vsize,'hsize',hsize,'tslice_length_sec',tslice_length_sec,...
        'block_stat','mean', 'specific','coher_cr');

feature_name1s = standard_feature_names('color',...
        'vsize',vsize,'hsize',hsize,'tslice_length_sec',tslice_length_sec,...
        'block_stat','std', 'specific','coher_cb');
feature_name2s = standard_feature_names('color',...
        'vsize',vsize,'hsize',hsize,'tslice_length_sec',tslice_length_sec,...
        'block_stat','std', 'specific','coher_cr');

name_exists1s = 1;
name_exists2s = 1;
for loop = 1:max(size(clip_structs)),
    try,
        % check if features already computed for this clip.
        t = clock;
        [name_exists1, file_name] = write_feature_standard(0,feature_base_dir,feature_name1,clip_structs(loop),'exists');
        [name_exists2, file_name] = write_feature_standard(0,feature_base_dir,feature_name2,clip_structs(loop),'exists');
        if optional_std,
            [name_exists1s, file_name] = write_feature_standard(0,feature_base_dir,feature_name1s,clip_structs(loop),'exists');
            [name_exists2s, file_name] = write_feature_standard(0,feature_base_dir,feature_name2s,clip_structs(loop),'exists');
        end
        if name_exists1 & name_exists2 & name_exists1s & name_exists2s,
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
        data1 = zeros(vert, horiz, number_tslices);
        data2 = data1;
        if optional_std,
            data1s = data1;
            data2s = data1;
        end
        
        % For each time-slice in this clip, read the time-slice; perceptually filter
        % (none); and compute block statistic.
        for cnt = 1:number_tslices,
            [y cb cr] = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec,cnt, ...
                'vsize',vsize,'hsize',hsize, ...
                'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
            clear y;
            
            if ~optional_std,
                data1(:,:,cnt) = block_statistic (cb, vsize, hsize, 'mean');
                data2(:,:,cnt) = block_statistic (cr, vsize, hsize, 'mean');
            else
                [data1(:,:,cnt),data1s(:,:,cnt)] = block_statistic (cb, vsize, hsize, 'mean','std');
                [data2(:,:,cnt),data2s(:,:,cnt)] = block_statistic (cr, vsize, hsize, 'mean','std');
            end

            clear cb cr y;
        end
        
        % Write each feature all at once to a single file.
        write_feature_standard(data1, feature_base_dir,feature_name1, clip_structs(loop));
        write_feature_standard(data2, feature_base_dir,feature_name2, clip_structs(loop));
        if optional_std,
            write_feature_standard(data1s, feature_base_dir,feature_name1s, clip_structs(loop));
            write_feature_standard(data2s, feature_base_dir,feature_name2s, clip_structs(loop));
        end
        
        % Erase the large variables 'data' from memory to prevent problems.
        clear data*;
    catch
        status = 1;
        if verbose,
            fprintf('\n%s\n', lasterr);
            fprintf('\tSkipping clip.\n');
        end
    end
end
