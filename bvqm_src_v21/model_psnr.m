function [psnr,status] = model_psnr(test_structs, clip_structs, varargin)
% MODEL_PSNR
%   Compute standard PSNR on a clip structure (i.e., constant delay).
% SYNTAX
%   [pars,status] = model_psnr(test_structs, clip_structs);
%   [...] = model_psnr(...,'OptionalFlag',...);
% DESCRIPTION
%   Calculate the PSNR model on the specified clips.   
%       'test_structs' describes the video tests, formatted like GTests
%       'clip_structs' describges the video clips, formatted like GClips
%   Returns the following:
%       'pars', a parameter structure with the PSNR
%       'status' is 0 if operated correctly; 1 if an error was encountered
%
%  The following optional parameters are also available.  
%   'verbose'   Print progress report to screen.
%   'quiet'     Don't print anything to the screen (default).
%   'waitbar'   Show progress on a wait-bar.
%   'NoWaitbar' Do not put up a wait-bar.
%   '255'       Use 255 (computer white) for peak signal.
%   '235'       Use 235 (ITU-R Rec. BT-601 white) for peak signal.  This is
%               the default.
%   'byframe'   Process file pairs frame-by-frame.  Use for HDTV and
%               low-memory systems, where an attempt to load two scenes
%               into memory could cause problems. 
%   'byscene'   Load entirity of both file-pairs into memory and process.  
%               Faster but requires much memory.  Switch to By-Frame
%               processing if any problems occur (e.g., out of memory).
%               This is the default behavior. 


status = 0;
verbose = 0;
byscene = 1;
use_waitbar = 1;
peak_signal = 235;
for cnt = 1:length(varargin),
    if strcmpi(varargin{cnt},'verbose'),
        verbose = 1;
    elseif strcmpi(varargin{cnt},'quiet'),
        verbose = 0;
    elseif strcmpi(varargin{cnt},'235'),
        peak_signal = 235;
    elseif strcmpi(varargin{cnt},'255'),
        peak_signal = 255;
    elseif strcmpi(varargin{cnt},'waitbar'),
        use_waitbar = 1;
    elseif strcmpi(varargin{cnt},'nowaitbar'),
        use_waitbar = 0;
    elseif strcmpi(varargin{cnt},'byscene'),
        byscene = 1;
    elseif strcmpi(varargin{cnt},'byframe'),
        byscene = 0;
    else
        error('optional property not recognized');
    end
end

% get rid of scenes with no aligned segment defined.
new_loop = 1;
for loop = 1:length(clip_structs),
    if ~isnan(clip_structs(loop).align_start),
        new_clip_structs(new_loop) = clip_structs(loop);
        new_loop = new_loop + 1;
    end
end
clip_structs = new_clip_structs;

% initialize waitbar.

if use_waitbar,
    waitbar_handle = waitbar(0.0, 'Initializing', 'Name', 'PSNR Model');
    wait_add = 1/length(clip_structs);
end

psnr.clip_name = [];
psnr.inlsa_mos = [];
psnr.mos = [];
psnr.data = [];
psnr.par_name{1} = sprintf('psnr_%d', peak_signal);

curr = 1;
for loop = 1:max(size(clip_structs)),
    try
        if strcmp(clip_structs(loop).hrc{1}, 'original'),
            continue;
        end

        if verbose,
            fprintf('Clip %d of %d ==> %s:%s(%s)\n', loop, length(clip_structs), ...
                clip_structs(loop).test{1}, clip_structs(loop).scene{1}, ...
                clip_structs(loop).hrc{1});
        end
        if use_waitbar,
            waitbar(wait_add * (loop-1), waitbar_handle, ...
                sprintf('Clip ==> %s:%s(%s)', clip_structs(loop).test{1}, ...
                  clip_structs(loop).scene{1}, clip_structs(loop).hrc{1}), ...
                  'Name', 'PSNR Model');
            pause(0.2);
        end

        % Find the offset of the test structure for this clip.
        tnum = search_test_list(test_structs, clip_structs(loop));

        % Compute the default adjusted SROI and number of blocks available, so
        % can allocate memory to hold all features for this clip.
        [sroi,vert,horiz] = adjust_requested_sroi (clip_structs(loop));

        % find original for this PVS
        orig_offset = find_original(clip_structs, loop);

        if byscene,
            try
                tslice_length_sec = (clip_structs(loop).align_stop - ...
                    clip_structs(loop).align_start + 1) / clip_structs(loop).fps;
                if is_reframing_indicated(clip_structs(loop)),
                    tslice_length_sec = tslice_length_sec - 1/clip_structs(loop).fps;
                end
                number_tslices = total_tslices(clip_structs(loop),tslice_length_sec);

                if number_tslices ~= 1,
                    byscene = 0;
                    error('computation of number tslices wrong');
                end

                % For each time-slice in this clip, 
                sse = 0;
                yp = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec, 1, ...
                    'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
                yo = read_tslice(test_structs(tnum),clip_structs(orig_offset),tslice_length_sec, 1, ...
                    'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
                
                yp = (yp - yo).^2;
                
                sse = mean(mean(mean(yp)));
                
                clear yo yp;
            catch
                byscene = 0;
            end
        end
        
        if ~byscene,
            % figure out number of frames available.
            tslice_length_sec = 1.0 / clip_structs(loop).fps;
            number_tslices = total_tslices(clip_structs(loop),tslice_length_sec);

            % For each time-slice in this clip, 
            sse = 0;
            for cnt = 1:number_tslices,
                yp = read_tslice(test_structs(tnum),clip_structs(loop),tslice_length_sec, cnt, ...
                    'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);
                yo = read_tslice(test_structs(tnum),clip_structs(orig_offset),tslice_length_sec, cnt, ...
                    'sroi',sroi.top,sroi.left,sroi.bottom,sroi.right);

                yp = (yp - yo).^2;
                sse = sse + mean(mean(yp));

                clear yo yp;
            end
        end

        % finish computation of psnr
        mse = sse / number_tslices;
        psnr.data(curr) = 10.0 * log10( (peak_signal).^2 ./ mse);
        
        % 8-bit images, so max value 48
        if psnr.data(curr) > 48,
            psnr.data(curr) = 48;
        end

        % record results for this PVS
        psnr.clip_name{curr} = sprintf('%s_%s_%s', clip_structs(loop).test{1}, clip_structs(loop).scene{1}, ...
                clip_structs(loop).hrc{1});
        psnr.inlsa_mos(curr) = clip_structs(loop).inlsa_mos;
        psnr.mos(curr) = clip_structs(loop).mos;
        curr = curr + 1;

    catch
        status = 1;
        if verbose,
            fprintf('\n%s\n', lasterr);
            fprintf('\tSkipping clip.\n');
        end
    end
end

if use_waitbar,
    waitbar(1.0, waitbar_handle, 'Done', 'Name', 'PSNR Model');
    pause(0.2);
    close(waitbar_handle);
end



