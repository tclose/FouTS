function ext = file_extension(filename)
  dots_i = strfind(filename,'.');
  
  if dots_i
    ext = filename(dots_i(end)+1:end);
  else
    ext = '';
  end
end