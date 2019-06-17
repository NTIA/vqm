function dll_vfd_print(results, file_name, delay)
% DLL_VFD_PRINT
%  This function uses the results from dll_vfd as well as any temporal delay
%  previously calculated and creates a file with the corresponding original
%  index for every processed index (corrected for temporal delay).  The
%  three inputs are the 'results' of the vfd correction, computed by a
%  function like dll_vfd; 'file_name', the name of the file where the user
%  wants the information stored; and 'delay', the amount of temporal delay
%  calculated by a calibration option.

scan_type = dll_video('get_video_standard', 1);
reframe_flag = dll_calib_video('get_reframe');
% Translate results to use the orig and proc file indexing
if (strcmpi(scan_type,'progressive'))
    if delay >= 0
        proc_indices = delay + (1:length(results));
        orig_indices = results;
    else
        proc_indices = 1:length(results);
        orig_indices = abs(delay) + results;
    end
else % interlaced
    if delay >= 0
        proc_indices = 2*delay + (1:length(results));
        orig_indices = results;
    else
        proc_indices = 1:length(results);
        orig_indices = 2*abs(delay) + results;
    end
    if reframe_flag
        proc_indices = proc_indices + 1;
    end
end


% Save results
fid_results = fopen(file_name,'w');  % open results file for appending

if (strcmpi(scan_type,'progressive'))
    fprintf(fid_results,'File Name, Matching Frame Indices\n');
else % interlaced
    fprintf(fid_results,'File Name, Matching Field Indices\n');
end

npts = length(proc_indices);

fprintf(fid_results,'processed, ');
for i = 1:npts-1
    fprintf(fid_results,'%f, ',proc_indices(i));
end
fprintf(fid_results,'%f\n',proc_indices(i+1));

fprintf(fid_results,'original, ');
for i = 1:npts-1
    fprintf(fid_results,'%f, ',orig_indices(i));
end
fprintf(fid_results,'%f\n',orig_indices(i+1));
fclose(fid_results);

close all;

end