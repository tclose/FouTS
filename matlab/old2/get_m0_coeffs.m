function m0 = get_m0_coeffs(SH)

reshape(SH, length(SH),1);


total_index = 0;
i = 0;
indices = [];

while (total_index < length(SH))   
    
    
    block_size = 2 * i;
    
    new_index = total_index + block_size + 1;
           
    
    indices = [indices; new_index];
    
    total_index = new_index + block_size;
    i = i + 1;
end

m0 = SH(indices);