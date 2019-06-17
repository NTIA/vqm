function [pars, status] = ...
    model_fastlowbw(test_structs, clip_structs, feature_base_dir, tmp_nam, varargin);
% MODEL_fastlowbw
%   Compute fast lowbw/sec model on a clip structure.
% SYNTAX
%     [pars] = model_fastlowbw(test_structs, clip_structs, feature_base_dir, tmp_nam);
% DESCRIPTION
%   Calculate the lowbw model on the specified clips.  
%       'test_structs' describes the video tests, formatted like GTests
%       'clip_structs' describges the video clips, formatted like GClips
%       'feature_base_dir' is used to write / store features, calculated
%           using the lowbw model's SROI
%       'tmp_nam' is a name of a temporary file that can be used & erased
%   Returns the following:
%       'pars', a parameter structure with the lowbw model & parameters
%       'status' is 0 if operated correctly; 1 if an error was encountered
%
%  The following optional parameters are also available.  
%   'verbose'   Print progress report to screen.
%   'quiet'     Don't print anything to the screen (default).

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

try
    if verbose,
        verbose_string = 'verbose';
    else
        verbose_string = 'quiet';
    end

    % calculate features
    status = model_fastlowbw_feature_loop(test_structs, clip_structs, feature_base_dir, verbose_string);

    feature_base_dir = [feature_base_dir '/'];

    % figure out the list of clip names.  Fill in the names when computing
    % metrics.
    %pars.clip_name = cell(num_clips);
    ccnt = 1;

    % Loop through all clips, sorted by scene
    offsets = sort_clips_by('scene',clip_structs, test_structs);
    for loop = 1:length(offsets),
        % pick off the original for this scene.  Skip if no original defined.
        curr_offsets = offsets{loop};
        if ~strcmp(clip_structs(curr_offsets(1)).hrc{1},'original'),
            if verbose,
                fprintf('Warning:  Original missing for %s:%s.  Skipping these clips\n', ...
                    clip_structs(curr_offsets(1)).test{1}, clip_structs(curr_offsets(1)).scene{1});
            end
            continue;
        end

        % load original features. 
        clip_name = sprintf('%s_%s_%s.mat', ...
            clip_structs(curr_offsets(1)).test{1}, ...
            clip_structs(curr_offsets(1)).scene{1}, ...
            clip_structs(curr_offsets(1)).hrc{1});
        load( [feature_base_dir 'feature_fastlowbw_model/' clip_name ] );
        orig_data = data;

        % compress original features
        model_lowbw_compression ('compress', tmp_nam, orig_data.si_std, orig_data.hv_ratio, orig_data.y_mean, orig_data.cb_mean, orig_data.cr_mean, orig_data.ati_rms );
        [sio_data, part_si_min, part_si_max, ...
             hvo_ratio, part_hv_min, part_hv_max, lo_data, ...
             cbo_data, cro_data, part_c_min, part_c_max, part_c, ...
             atio_data, part_ati_min, part_ati_max, part_ati, code_ati ] = model_lowbw_compression('uncompress', tmp_nam);
        delete(tmp_nam);

        % loop through all processed versions of this clip.
        for cnt = 2:length(curr_offsets),
            % load processed features. 
            clip_name = sprintf('%s_%s_%s.mat', ...
                clip_structs(curr_offsets(cnt)).test{1}, ...
                clip_structs(curr_offsets(cnt)).scene{1}, ...
                clip_structs(curr_offsets(cnt)).hrc{1});
            load( [feature_base_dir 'feature_fastlowbw_model/' clip_name] );
            proc_data = data;

            % find out length of clip in seconds.
            [a,b,TIME_DELTA] = size(proc_data(5).si_std);

            for loop = 1:9,
                % calculate model
                do_test_print = 0;
                [data(loop).vqm, data(loop).hv_loss_par, data(loop).hv_gain_par,data(loop).si_loss_par, data(loop).si_gain_par, ...
                    data(loop).color_comb_par, data(loop).noise_par, data(loop).error_par] = ...
                    model_fastlowbw_parameters (...
                        data(loop).si_std, data(loop).hv_ratio, data(loop).y_mean, data(loop).cb_mean, data(loop).cr_mean, data(loop).ati_rms, ...
                        sio_data, hvo_ratio, lo_data, cbo_data, cro_data, atio_data, ...
                        clip_structs(curr_offsets(cnt)).fps, part_si_min, part_si_max, part_hv_min, part_hv_max, ...
                        part_c_min, part_c_max, part_c, part_ati_min, part_ati_max, ...
                        part_ati, code_ati, do_test_print, TIME_DELTA);
            end
            
            % select smallest average VQM score shift
            for shift=1:9,
                vqm_mean(shift) = mean(data(shift).vqm);
            end
            [junk shift] = min(vqm_mean);
            col = floor((shift-1)/3)-1;
            row = mod(shift-1,3)-1;            

            % fill in the name of this clip
            pars.clip_name{ccnt} = sprintf('%s_%s_%s', clip_structs(curr_offsets(cnt)).test{1}, ...
                    clip_structs(curr_offsets(cnt)).scene{1}, clip_structs(curr_offsets(cnt)).hrc{1});
            pars.inlsa_mos(ccnt) = clip_structs(curr_offsets(cnt)).inlsa_mos;
            pars.mos(ccnt) = clip_structs(curr_offsets(cnt)).mos;

            num = length(data(loop).vqm);
            pars.data(:,ccnt) = [data(shift).vqm(num), data(shift).hv_loss_par(num), data(shift).hv_gain_par(num), ...
                data(shift).si_loss_par(num), data(shift).si_gain_par(num), ...
                data(shift).color_comb_par(num), data(shift).noise_par(num), data(shift).error_par(num), ...
                row, col];
            ccnt = ccnt + 1;
        end
    end
    
    if ~exist('pars','var'),
        error('Due to previous errors, No parameter data available');
    end

    % record name of each parameter
    pars.par_name = { 'VQM_fastlowbw', 'hv_loss', 'hv_gain', 'si_loss', 'si_gain', 'color_comb', 'noise', 'error', 'shiftH', 'shiftV'};
catch
    status = 1;
    pars = [];
    if verbose,
        fprintf('\n%s\n', lasterr);
        fprintf('\tAborting.\n'); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
function [si, hv, hvb, cb, cr, lum] = m10k_fix_sroi(si, hv, hvb, cb, cr, lum);
% Make sure all 30x30_1s features have same SROI -- evenly divisible by 96!

[r1,c1,t1] = size(si);
[r2,c2,t2] = size(hv);
[r3,c3,t3] = size(hvb);
[r4,c4,t4] = size(cb);
[r5,c5,t5] = size(cr);

if nargin ~= 6,
    r = min([r1 r2 r3 r4 r5]);
    c = min([c1 c2 c3 c4 c5]);
    t = min([t1 t2 t3 t4 t5]);
else
    [r6,c6,t6] = size(lum);
    r = min([r1 r2 r3 r4 r5 r6]);
    c = min([c1 c2 c3 c4 c5 c6]);
    t = min([t1 t2 t3 t4 t5 t6]);
end


while mod(r,3) ~= 0,
    r = r - 1;
end
while mod(c,3) ~= 0,
    c = c - 1;
end
while mod(t,2) ~= 0,
    t = t - 1;
end

si = si(1:r,1:c,1:t);
hv = hv(1:r,1:c,1:t);
hvb = hvb(1:r,1:c,1:t);
cb = cb(1:r,1:c,1:t);
cr = cr(1:r,1:c,1:t);
if nargin == 6,
    lum = lum(1:r, 1:c, 1:t);
end
