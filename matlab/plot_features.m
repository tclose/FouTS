function plot_features(features_dir, varargin)
%  
% PURPOSE: Plots histograms of prior distributions collated by extract_features.cpp
%  
% ARGUMENTS: 
%  
%           features_directory     The directory containing the features text files.
%  
% OPTIONS (name, description, type, default):
%  
%          -all           Plots all the features. Overrides all other flags
%                              bool
%                              0
%  
%          -v0_norm       Plot the norm of the v0 vector
%                              bool
%                              0
%  
%          -v1_norm       Plot the norm of the v1 vector
%                              bool
%                              0
%  
%          -v2_norm       Plot the norm of the v2 vector
%                              bool
%                              0
%  
%          -dot_v1_v2     Plot the dot product of the v1 and v2 vectors
%                              bool
%                              0
%  
%          -angle_v1_v2   Plot the angle between the v1 and v2 vectors
%                              bool
%                              0
%  
%          -dot_v0_v1     Plot the dot product of the v0 and v1 vectors
%                              bool
%                              0
                             
                             
                             
  if ~isdir(features_dir)
    error(['Argument 1 (''' features_dir ''') is not a directory, it should be a directory containing the features text files.']);
  end

  description = 'Plots histograms of prior distributions collated by extract_features.cpp';
  
  arguments = {'features_directory', 'The directory containing the features text files.'};

  options = {...
            'plot_all     ', 0, 'bool', 'Plots all the features. Overrides all other flags';...
            'v0_norm      ', 0, 'bool', 'Plot the norm of the v0 vector';...
            'v1_norm      ', 0, 'bool', 'Plot the norm of the v1 vector';...
            'v2_norm      ', 0, 'bool', 'Plot the norm of the v2 vector';...
            'dot_v1_v2    ', 0, 'bool', 'Plot the dot product of the v1 and v2 vectors';...
            'angle_v1_v2  ', 0, 'bool', 'Plot the angle between the v1 and v2 vectors';...
            'dot_v0_v1    ', 0, 'bool', 'Plot the dot product of the v0 and v1 vectors';...
            'd2_v0_v1_dot ', 0, 'bool', 'Plot the dot product of the v0 and v1 vectors against norm of v0 on a d2 histogram';...                                    
            'd2_v1_v2_dot ', 0, 'bool', 'Plot the dot product of the v1 and v2 vectors against norm of v2 on a d2 histogram';...
            'd2_v0_v1_norm', 0, 'bool', 'Plot the v0 and v1 norms against each other in 2d histogram';...
            'd2_v1_v2_norm', 0, 'bool', 'Plot the v1 and v2 norms against each other in 2d histogram';...                        
            'd2_v0_v2_norm', 0, 'bool', 'Plot the v0 and v2 norms against each other in 2d histogram'...                                    
            };
                                    
  supplied_options = parse_options(options, varargin);

  for option_i = 1:size(supplied_options,1)
    if strcmp(supplied_options{option_i,3}, 'string')
      eval([supplied_options{option_i,1} ' = ''' supplied_options{option_i,2} ''';']);
    else
      
      if ~isstr(supplied_options{option_i,2})
        supplied_options{option_i,2} = mat2str(supplied_options{option_i,2});
      end
      
      eval([supplied_options{option_i,1} ' = ' supplied_options{option_i,2} ';']);
    end
  end


  if help_display
    display_help_message(description, arguments, options);
    return;
  end

num_options_selected = 0;

for option_i = 1:size(options,1)
  
  if eval(options{option_i,1})
    num_options_selected = num_options_selected + 1;
  end
  
  if plot_all
    eval([options{option_i} '= 1']);
  end
  
end

if ~num_options_selected
  error('No options selected, use ''-help'' to display possible options.');
end
 
% if ~plot_all && ~v0_norm && ~v1_norm && ~v2_norm && ~dot_v1_v2 && ~dot_v0_v1 && ~angle_v1_v2 && ~d2_v0_v1_norm && ~d2_v1_v2_norm && ~d2_v0_v2_norm && ~d2_v1_v2_dot && ~d2_v0_v1_dot
%   error('No options selected, use ''-help'' to display possible options.');
% end
  
% if plot_all
%   v0_norm = 1;
%   v1_norm = 1;
%   v2_norm = 1;
%   dot_v1_v2 = 1;
%   angle_v1_v2 = 1;
%   dot_v0_v1 = 1;  
%   d2_v0_v1_dot = 1;  
%   d2_v1_v2_dot = 1;  
%   d2_v0_v1_norm = 1;  
%   d2_v1_v2_norm = 1;      
%   d2_v0_v2_norm = 1;    
% end

if  v0_norm || d2_v0_v1_dot || d2_v0_v2_norm || d2_v0_v1_norm
  v0_norm_data = load([ features_dir '/v0_norm.txt']);
end

if  v1_norm || d2_v0_v1_norm || d2_v1_v2_norm
  v1_norm_data = load([ features_dir '/v1_norm.txt']);
end

if  v2_norm || d2_v1_v2_dot || d2_v0_v2_norm || d2_v1_v2_norm
  v2_norm_data = load([ features_dir '/v2_norm.txt']);
end

if  dot_v1_v2 || d2_v1_v2_dot
  dot_v1_v2_data = load([ features_dir '/dot_v1_v2.txt']);
end

if  angle_v1_v2
  angle_v1_v2_data = load([ features_dir '/angle_v1_v2.txt']);
end

if  dot_v0_v1 || d2_v0_v1_dot
  dot_v0_v1_data = load([ features_dir '/dot_v0_v1.txt']);
end
  




if  v0_norm
  f = figure();
  set(f,'Units','normalized');  
  set(f, 'Position', [0.3 0.25 0.4 0.5]);
  set(f, 'Name', strrep(features_dir,'_',' '));  

  hist(v0_norm_data,100);
  title('v0 norm');    

end

if  v1_norm
  f = figure();
  set(f,'Units','normalized');  
  set(f, 'Name', strrep(features_dir,'_',' '));  
  hist(v1_norm_data,100);  
  set(f, 'Position', [0.3 0.25 0.4 0.5]);  
  title('v1 norm');  

end

if  v2_norm
  f = figure();
  set(f,'Units','normalized');  
  set(f, 'Name', strrep(features_dir,'_',' '));  
  hist(v2_norm_data,100);  
  set(f, 'Position', [0.3 0.25 0.4 0.5]);  
  title('v2 norm');  

end

if  dot_v1_v2
  f = figure();
  set(f,'Units','normalized');  
  set(f, 'Name', strrep(features_dir,'_',' '));  
  hist(dot_v1_v2_data,100);
  set(f, 'Position', [0.3 0.25 0.4 0.5]);  
  title('dot v1 v2');  

end

if  angle_v1_v2
  f = figure();
  set(f,'Units','normalized');  
  hist(angle_v1_v2_data,100);  
  set(f, 'Position', [0.3 0.25 0.4 0.5]);  
  title('angle v1 v2');  

end

if  dot_v0_v1
  f = figure();
  set(f,'Units','normalized');  
  set(f, 'Name', strrep(features_dir,'_',' '));  
  hist(dot_v0_v1_data,100);
  set(f, 'Position', [0.3 0.25 0.4 0.5]);  
  title('dot v0 v1');  

end

if d2_v0_v1_dot
  f = figure();
  set(f,'Units','normalized');  
  set(f, 'Name', strrep(features_dir,'_',' '));  
  title('dot v0 v1 - v0 norm');  
  set(f, 'Position', [0.3 0.25 0.4 0.5]);
  hist2d(dot_v0_v1_data, v0_norm_data,100,100);
  
%   cameratoolbar('Show');
%   cameratoolbar('SetMode','orbit');    

end

if d2_v1_v2_dot
  f = figure();
  set(f,'Units','normalized');  
  
  title('dot v1 v2 - v2 norm');  
  set(f, 'Position', [0.3 0.25 0.4 0.5]);  
  set(f, 'Name', strrep(features_dir,'_',' '));  
  hist2d(dot_v1_v2_data, v2_norm_data,100,100);
  
%   cameratoolbar('Show');
%   cameratoolbar('SetMode','orbit');    

end

if d2_v0_v2_norm
  f = figure();
  set(f,'Units','normalized');  
  set(f, 'Name', strrep(features_dir,'_',' '));  
  title('v0 norm - v2 norm');  
  set(f, 'Position', [0.3 0.25 0.4 0.5]);
  hist2d(v0_norm_data, v2_norm_data,100,100);
  
%   cameratoolbar('Show');
%   cameratoolbar('SetMode','orbit');    

end

if d2_v1_v2_norm
  f = figure();
  set(f,'Units','normalized');  
  
  title('v1 norm - v2 norm');  
  set(f, 'Position', [0.3 0.25 0.4 0.5]);  
  set(f, 'Name', strrep(features_dir,'_',' '));  
  hist2d(v1_norm_data, v2_norm_data,100,100);
  
%   cameratoolbar('Show');
%   cameratoolbar('SetMode','orbit');    

end

if d2_v0_v1_norm
  f = figure();
  set(f,'Units','normalized');  
  
  title('v0 norm - v2 norm');  
  set(f, 'Position', [0.3 0.25 0.4 0.5]);  
  set(f, 'Name', strrep(features_dir,'_',' '));  
  hist2d(v0_norm_data, v1_norm_data,100,100);
  
%   cameratoolbar('Show');
%   cameratoolbar('SetMode','orbit');    

end

