function status = feature_loop_color_si_hv(test_structs, clip_structs,...
    vert_size,horiz_size, tslice_length_sec,feature_base_dir, varargin)
% FEATURE_LOOP_COLOR_SI_HV
%  Loop through a list of clips; compute SI, HV and HVbar; and saves each 
%  clip's features in a separate file.
% SYNTAX
%  feature_loop_color_si_hv(test_structs, clip_structs,vert_size,horiz_size,tslice_length_sec,feature_base_dir)
%  feature_loop_color_si_hv(...,'PropertyName',...);
% DESCRIPTION
%  This function computes std(SI) and mean(HV) and mean(HVbar)
%  on combined (Y, Cb and Cr) color planes, combining results
%  This function takes variable 'test_structs' (of the same format as
%  GTests), which describes each video test and the location of the
%  associated video; and variable 'clip_structs' (of the same format as
%  GClips), which lists the clips to be processed. 'tslice_length_sec'
%  specifies the length of the time-slice in SECONDS, and 'feature_base_dir' is
%  the path to the directory where the feature's sub-directory should be
%  created.  That sub-directory's name will be this feature's name.  Also
%  takes the size of the block (vert_size x horiz_size) in pixels, horizontally 
%  and vertically.
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
feature_si_name = standard_feature_names('YCbCr',...
        'vsize',vsize,'hsize',hsize,'tslice_length_sec',tslice_length_sec,...
        'block_stat','std', 'specific','si13');
feature_hv_name = standard_feature_names('YCbCr',...
        'vsize',vsize,'hsize',hsize,'tslice_length_sec',tslice_length_sec,...
        'block_stat','mean', 'specific','hv13');
feature_hvb_name = standard_feature_names('YCbCr',...
        'vsize',vsize,'hsize',hsize,'tslice_length_sec',tslice_length_sec,...
        'block_stat','mean', 'specific','hvbar13');

for loop = 1:max(size(clip_structs)),
    try
        t = clock;
        [name_si_exists, file_name] = write_feature_standard(0,feature_base_dir,...
            feature_si_name,clip_structs(loop),'exists');
        [name_hv_exists, file_name] = write_feature_standard(0,feature_base_dir,...
            feature_hv_name,clip_structs(loop),'exists');
        [name_hvb_exists, file_name] = write_feature_standard(0,feature_base_dir,...
            feature_hvb_name,clip_structs(loop),'exists');
        if name_si_exists & name_hv_exists & name_hvb_exists,
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
                            'vsize',vsize,'hsize',hsize, 'extra',6);
        number_tslices = total_tslices(clip_structs(loop),tslice_length_sec);

        data_si = zeros(vert, horiz, number_tslices);
        data_hv = zeros(vert, horiz, number_tslices);
        data_hvb = zeros(vert, horiz, number_tslices);

        % For each time-slice in this clip, read the time-slice; perceptually filter
        % (si13 filter here); and compute block statistic.
        for cnt = 1:number_tslices,
            [y,cb,cr] = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec,cnt, ...
                'vsize',vsize,'hsize',hsize, 'extra',6, ...
                'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
            [si,hv,hvb] = filter_si_hv(y);

            [sic,hvc,hvbc] = filter_si_hv(cb);
            si = si + sic;
            hv = hv + hvc;
            hvb = hvb + hvbc;
            clear sic hvc hvbc;

            [sic,hvc,hvbc] = filter_si_hv(cr);
            si = si + sic;
            hv = hv + hvc;
            hvb = hvb + hvbc;
            clear sic hvc hvbc;

            data_si(:,:,cnt) = block_statistic (si, vsize, hsize, 'std');
            data_hv(:,:,cnt) = block_statistic (hv, vsize, hsize, 'mean');
            data_hvb(:,:,cnt) = block_statistic (hvb, vsize, hsize, 'mean');

            clear si hv hvb y cb cr sic hvc hvbc;
        end

        % Write each features all at once to a single file.
        write_feature_standard(data_si, feature_base_dir,feature_si_name, clip_structs(loop));
        write_feature_standard(data_hv, feature_base_dir,feature_hv_name, clip_structs(loop));
        write_feature_standard(data_hvb, feature_base_dir,feature_hvb_name, clip_structs(loop));

        % Erase the large variable 'data' from memory to prevent problems.
        clear data_si data_hv data_hvb;
    catch
        status = 1;
        if verbose,
            fprintf('\n%s\n', lasterr);
            fprintf('\tSkipping clip.\n');
        end
    end
end
