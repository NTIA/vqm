function [jclips,jtests, status] = bvqm_pc_parse(jclips,jtests)
% parse or truncate clips to have maximum length of 15sec
% if anything goes wrong, loop back to the beginning!

max_len = max([jclips.loc_stop] - [jclips.loc_start] + 1);
min_len = min([jclips.loc_stop] - [jclips.loc_start] + 1);
if max_len / jclips(1).fps <= 15,
    status.is_parsed = 0;
    return;
end


button = questdlg('Some files contain more than 15 seconds of video.  Do you want to truncate all video files to 15 second length, or parse all video files into shorter sequences of the same length?', ...
    '15sec Maximum', 'Truncate','Parse','Truncate');
while strcmp(button,''),
    button = questdlg('Please select "truncate" or "parse".  This question cannot be skipped or exited.  Some files contain more than 15 seconds of video.  Do you want to truncate all video files to 15 second length, or parse all video files into shorter sequences of the same length?', ...
        '15sec Maximum', 'Truncate','Parse','Truncate');
end

% handle truncation.  Easy.
if isequal(button,'Truncate'),
    status.is_parsed = 1;
    for jc=1:length(jclips),
        jclips(jc).loc_stop = min(jclips(jc).loc_stop,15*jclips(jc).fps);
    end
else
    % parse instead.  More complicated.
    is_max = min(floor(max_len / jclips(1).fps), 15); 
    prompt = { sprintf('Enter parse length in seconds (5 to %d), an integer', is_max), ...
        'Enter parse shift in seconds (1 minimum), an integer'}; 
    the_answer = inputdlg( prompt, 'Parse', 1, {'10','5'});
    
    % if 'cancel' or 'x' pressed, go to previous message.
    if length(the_answer) == 0,
        [jclips,jtests, status] = bvqm_pc_parse(jclips,jtests);
        return;
    end
    % check if return values are sensible
    if length(the_answer) ~= 2 || length(str2num(the_answer{1})) == 0 || length(str2num(the_answer{1})) == 0,
        uiwait(msgbox('Error reading parse request', 'error', 'error'));
        [jclips,jtests, status] = bvqm_pc_parse(jclips,jtests);
        return;
    end
    
    parse_length = str2num(the_answer{1});
    parse_overlap = str2num(the_answer{2});
    if round(parse_length) ~= parse_length || round(parse_overlap) ~= parse_overlap,
        uiwait(msgbox('Invalid parse request.  Parse length and shift must be integers.', 'error', 'error'));
        [jclips,jtests, status] = bvqm_pc_parse(jclips,jtests);
        return;
    end
    if parse_length < 1 || parse_length > 15 || parse_overlap < 1,
        uiwait(msgbox('Invalid parse request.  Parse length must be between 1sec and 15 sec, and the shift must be at least 1sec', 'error', 'error'));
        [jclips,jtests, status] = bvqm_pc_parse(jclips,jtests);
        return;
    end
    
    if min_len < parse_length,
        uiwait(msgbox('Some clips are shorter than parse length request.  Those clips will be discarded','parse','warn'));
    end
    
    % go from seconds to frames.
    status.is_parsed = 2;
    status.parse_length = parse_length;
    status.parse_overlap = parse_overlap;

    parse_length = round(parse_length * jclips(1).fps);
    parse_overlap = round(parse_overlap * jclips(1).fps);

    % get new list of clips
    jcc = 1;
    for jc=1:length(jclips),
        start_next = jclips(jc).loc_start;
        stop_next = start_next + parse_length - 1;
        while stop_next <= jclips(jc).loc_stop,
            jclips_new(jcc) = jclips(jc);
            jclips_new(jcc).loc_start = start_next;
            jclips_new(jcc).loc_stop = stop_next;
            jclips_new(jcc).scene = { sprintf('%s%d', jclips(jc).scene{1}, round(stop_next/jclips(1).fps ) ) };
            
            jcc = jcc + 1;
            start_next = start_next + parse_overlap;
            stop_next = stop_next + parse_overlap;
            continue;
        end
    end
    
    % update Tests.
    jtests_new = jtests;
    jtests_new.hrc = unique( [jclips_new.hrc] );
    jtests_new.scenes = unique( [jclips_new.scene] );
    
    % check if this is valid!
    [clip_err, clip_msg]=check_clips(jclips_new, jtests_new, 'quiet');
    if clip_err ~= 0,
        uiwait( msgbox( sprintf('Parsing failed. %s', clip_msg) ) );
        [jclips,jtests, status] = bvqm_pc_parse(jclips,jtests);
        return;
    end
    
    % return new jclips, jtests parsed.
    offsets = sort_clips_by('none',jclips_new, jtests_new);
    jclips = jclips_new(offsets);
    jtests = jtests_new;

end
