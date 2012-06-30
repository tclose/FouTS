% partial_track.m
% Old file: smooth_tracks.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads the TCK file generated with MRTRIX and generates a file with smoothed
% versions of the tracks.
% The smoothing method used is MOVING and the SPAN window must be chosen. 
% Modified to read a TCK and generate a new TCK file with part of the tracks
% contained (FC-31/10/08)
% Select: input and output files, Ntracks. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%%%%%%%%%%%%%%%%%
% Parameters to select
%%%%%%%%%%%%%%%%%
% Ntracks = 10000;      % Number of trakcs to be selected from the the TCK input 
Ntracks = input('Number of tracks?= ');
%%%%%%%%%%%%%%%%%

% Select TCK file to be smoothed
[fname,pname] = uigetfile('*.tck','Select TCK file');  
fidtck = fopen([pname fname],'r');
% fidtck = fopen('/home/fercala/DWI150/150lisa/curvature/20K/wholebrain3.tck','r');

% Select OUTPUT  file name for the smoothed TCK file 
[fname2,pname2] = uiputfile('*.tck','Select OUTPUT file');  
fidout = fopen([pname2 fname2],'w');
% fidout = fopen('/home/fercala/DWI150/150lisa/curvature/20K/wholebrain3sm.tck','w');

% Find file offset
temp ='abcd';
i=1;
while 1
    temp =  fgetl(fidtck);
    if strcmp(temp(1:4),'file'); break; end;
    i=i+1;
end;
[temp1 offset] = strtok(temp);
[temp1 offset] = strtok(offset);
offset = str2double(offset);

% Copy header to output file
status1 = fseek(fidtck,0,-1);    % rewind file
headerinfo = fread(fidtck,offset);  % header info
fwrite(fidout,headerinfo);

% Find coordinates of tracks (one at a time)
status = fseek(fidtck,offset,-1);

ptemp = zeros(3,1);
% k = 1;

% while (~isinf(sum(ptemp)))      %repeat until Inf is found (i.e. EOF)
for k=1:Ntracks;
    clear p
    i=1;
     while 1;
        ptemp = fread(fidtck,3,'float32','ieee-le');
        if isinf(sum(ptemp)); 
            fwrite(fidout,ptemp,'float32','ieee-le');   % end file with Inf
            fclose(fidout);
            break;
        elseif ~isfinite(sum(ptemp)); break; 
        end;
        p(i,:) = ptemp;
        i = i+1;
    end;
    if isinf(sum(ptemp)); break; end;
%     figure(10), plot3(squeeze(p(:,1)),squeeze(p(:,2)),squeeze(p(:,3)),'-k'), grid on, hold on
    fwrite(fidout,p','float32','ieee-le');
    fwrite(fidout,ptemp,'float32','ieee-le');   % add NaN at end of track

%     k=k+1;
    disp(k);
    
end;    % end of TRACKS loop

fclose('all')




