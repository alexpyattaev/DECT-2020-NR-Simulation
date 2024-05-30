
% Follows section 5.3 in the spec
function [N_TB_bits, C] = Transport_block_size(tm_mode, mcs, N_PDC_REs, Z)

    % copy all relevant variables
    N_bps = mcs.N_bps;  
    
    % Total number of PDC bits carried by data symbols
    N_PDC_bits = floor(tm_mode.N_SS * N_PDC_REs * N_bps * mcs.R_t / mcs.R_b);
    
    % CRC length
    L = 24;
    
    if N_PDC_bits <= 512
        M = 8;
    elseif N_PDC_bits <= 1024
        M = 16;
    elseif N_PDC_bits <= 2048
        M = 32;
    else
        M=64;
    end


    %Calculate largest multiple of M that is less than N_PDC_bits
    N_M = floor(N_PDC_bits/M) * M;
    
    if N_M <= Z
        N_TB_bits = N_M - L;
        C = 1;
    else
        C = ceil((N_M-L)/Z);
        N_TB_bits = N_M - (C+1)*L;
    end
end

