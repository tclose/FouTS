function values = parse_properties_row(line)

  line = deblank(line);

  if isempty(line)
    error('No keys found in extended properties file.');
    
  elseif isempty(regexp(line,'\t', 'once'))
    values{1} = line;
    
  else  
    values = regexp(line,'\t','split');
  end

end