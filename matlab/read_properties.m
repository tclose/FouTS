function props = read_properties (filename)
% function props = read_properties (filename)
% Reads properties from a file and returns them in a struct.

f = fopen (filename, 'r');
if (f<1) 
  disp (['could not open file ' filename ]);
  return
end

while 1
  L = fgetl(f);
  if ~ischar(L), break, end;
  L = strtrim(L);
  if strcmp(L, 'END'), break, end;
  d = strfind (L,':');
  if isempty(d)
    disp (['invalid line in header: ''' L ''' - ignored']);
  else
    key = strrep(lower(strtrim(L(1:d(1)-1))), ' ', '_');
    value = strtrim(L(d(1)+1:end));
    eval(['props.' key ' = value;']);
  end
end
fclose(f);




function S = split_strings (V, delim)
  S = {};
  while size(V,2) > 0
    [R, V] = strtok(V,delim);
    S{end+1} = R;
  end

