function read_extracted_features(dir_name)

fid = fopen([ dir_name '/v1_norm.dat']);

v1_norm = load([ dir_name '/v1_norm.dat']);

v2_norm = fopen([ dir_name '/v2_norm.dat']);

v0_norm = fopen([ dir_name '/v0_norm.dat']);

dot_v1_v2 = fopen([ dir_name '/dot_v1_v2.dat']);



% 
% fid = fopen([ dir_name '/v1_norm.dat']);
% 
% [w,v,machine_format] = fopen(fid);
% 
% v1_norm = fread(fid, inf, 'float32=>double',machine_format);
% fclose(fid);
% 
% 
% fid = fopen([ dir_name '/v2_norm.dat']);
% v2_norm = fread(fid, inf, 'float32=>double',machine_format);
% fclose(fid);
% 
% fid = fopen([ dir_name '/v0_norm.dat']);
% v0_norm = fread(fid, inf, 'float32=>double',machine_format);
% fclose(fid);
% 
% fid = fopen([ dir_name '/dot_v1_v2.dat']);
% dot_v1_v2 = fread(fid, inf, 'float32=>double',machine_format);
% fclose(fid);

hist(v1_norm)