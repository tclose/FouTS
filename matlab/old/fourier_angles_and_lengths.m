filename = '../chunked_fourier_d4/090403192509-9-8-11.tck'

strands = load_mrtrix_strands(filename);

for (strand_i = 1:size(strands,1))
	strand = strands{strand_i,1};
	
	f0 = strand(1,:)';
	f1 = strand(2,:)';
	f2 = strand(3,:)';
%	f3 = strand(4,:)';

%	length_f0 = norm(f0)
	length_f1 = norm(f1);
	length_f2 = norm(f2);
	
	ang = 180/pi * acos(dot(f1,f2) /  (length_f1 * length_f2)) - 90;
	
	disp(['f1 length = ' num2str(length_f1) ',     f2_length = ' num2str(length_f2) ',     f2/f1_length = ' num2str(length_f2/length_f1) ',     f1/f2 angle = ' num2str(ang) ]);
	
end