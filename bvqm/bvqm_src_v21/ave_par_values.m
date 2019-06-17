function [pars_out] = ave_par_values(pars_in, ave_query)
% AVE_PAR_VALUES
%  Average the mos and parameter values in a parameter structure (of the same format
%  as GPars) and return a new parameter structure with these average parameter values.
% SYNTAX
%  [pars_out] = ave_par_values(pars_in, ave_query);
% DESCRIPTION
%  "par_in" and "pars_out" are parameter structures, of the same format as GPars.
%  "ave_query" controls how the averaging is performed and is one of the following:
%       'scene'  average across scenes to produce a single value for each HRC
%       'hrc'  average across hrcs to produce a single value for each scene
%
%  The portion of the clip name in pars_out that contains the quantity being averaged 
%  over will be replaced by "AVE".  For instance, if the input clip name is 
%  cable_basket_c13, the output clip name will be cable_AVE_c13 if ave_query = 'scene'
%  and cable_basket_AVE if ave_query = 'hrc'.
%
%  Averaging is never performed over multiple tests, only within an individual test.
%

[num_pars, num_clips] = size(pars_in.data); 

%  Parse the clip names into three arrays for test, scene, and hrc
[tests, scenes, hrcs] = pars_parse(pars_in);

%  Set up the common output arguments for averaging across scenes or hrcs
ave_num_clips = 1;
ave_clip_count = 1;
ave_data = pars_in.data(:,1);
ave_mos = pars_in.mos(1);
ave_inlsa_mos = pars_in.inlsa_mos(1);

%  Loop through all the clips, starting from the second clip, and sum data
%  and mos values into their appropriate summers.

switch lower(ave_query)
    case {'scene'} % average across scenes
		ave_names = [tests(1); {'AVE'}; hrcs(1)];
    case {'hrc'}  % average across hrcs
		ave_names = [tests(1); scenes(1); {'AVE'}];
    otherwise
        disp('Unknown ave_query')
        pars_out = 0;
        return;
end

for i_clip = 2:num_clips
    switch lower(ave_query)
        case {'scene'}
			this_clip = [tests(i_clip); {'AVE'}; hrcs(i_clip)];
        case {'hrc'}
			this_clip = [tests(i_clip); scenes(i_clip); {'AVE'}];
	end
	new_clip = 1;  %  Boolean variable to see if this is a new HRC or scene, set to true
	for i_ave_clip = 1:ave_num_clips
		this_ave_clip = ave_names(:,i_ave_clip);
		if (strcmp(this_clip,this_ave_clip))
            new_clip = 0;  %  Not a new clip
            ave_clip_count(i_ave_clip) = ave_clip_count(i_ave_clip)+1;
            ave_data(:,i_ave_clip) = ave_data(:,i_ave_clip) + pars_in.data(:,i_clip);
            ave_mos(i_ave_clip) = ave_mos(i_ave_clip) + pars_in.mos(i_clip);
            ave_inlsa_mos(i_ave_clip) = ave_inlsa_mos(i_ave_clip) + pars_in.inlsa_mos(i_clip);
		end
	end
   if (new_clip == 1)
      ave_num_clips = ave_num_clips + 1;
      ave_names = [ave_names this_clip];
      ave_clip_count = [ave_clip_count 1];
      ave_data = [ave_data pars_in.data(:,i_clip)];
      ave_mos = [ave_mos pars_in.mos(i_clip)];
      ave_inlsa_mos = [ave_inlsa_mos pars_in.inlsa_mos(i_clip)];
   end
end

%  Compute the averages
for i_ave_clip = 1:ave_num_clips
    ave_data(:,i_ave_clip) = ave_data(:,i_ave_clip)/ave_clip_count(i_ave_clip);
    ave_mos(i_ave_clip) = ave_mos(i_ave_clip)/ave_clip_count(i_ave_clip);
    ave_inlsa_mos(i_ave_clip) = ave_inlsa_mos(i_ave_clip)/ave_clip_count(i_ave_clip);
end

%  Organize output information
pars_out.par_name = pars_in.par_name;  % parameter names stay the same
new_clip_name = cell(1,ave_num_clips);
for i = 1:ave_num_clips
    new_clip_name{i} = strcat(ave_names{1,i},'_',ave_names{2,i},'_',ave_names{3,i});
end
pars_out.clip_name = new_clip_name;
pars_out.inlsa_mos = ave_inlsa_mos;
pars_out.mos = ave_mos;
pars_out.data = ave_data;
