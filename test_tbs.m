
validation_file = fopen('validation_data_tbs_mu_1_beta_1_nss_1.txt');

line_format = "Z: %d, mu: %d, b: %d, n_ss: %d, packet_len: %d, mcs: %d, chb: %d, tbs: %d, C_plus: %d, C_minus: %d, K_plus: %d, K_minus: %d, cr: %f\n";
validation_data = fscanf(validation_file, line_format, [13, Inf])';

validation_pkt_len = validation_data(:, 5);
validation_mcs = validation_data(:, 6);
validation_tbs = validation_data(:, 8);

redundancy_version=0

% Current configuration is for 1.7MHz channel
% these variables need to be set before creating tx and rx
mac_meta_tx.u = 1;                                  % mu = 1, 2, 4 or 8
% Choice of u sets the subcarrier spacing
mac_meta_tx.b = 1;                                  % beta = 1, 2, 4, 8, 12 or 16
numerology = lib_4_physical_layer_principles.numerologies(mac_meta_tx.u,mac_meta_tx.b);


val_idx = 1;
for packet_length_subslots = [1:16, 18:2:32]
    for mcs_index = 0:11
        if packet_length_subslots == 1 && mcs_index ==0
            continue;
        end
        %PacketLengthType: 0 for subslots, 1 for slots
        %PacketLength: min is 1, max is 16 according to Table 6.2.1-2a in part 4
        if mod(packet_length_subslots, numerology.N_SLOT_u_subslot) % need to operate in subslots
            mac_meta_tx.PacketLengthType = 0;
            mac_meta_tx.PacketLength = double(packet_length_subslots);
        else
            mac_meta_tx.PacketLengthType = 1;

            mac_meta_tx.PacketLength = double(packet_length_subslots) / numerology.N_SLOT_u_subslot;

        end

        if mac_meta_tx.PacketLength > 16
            error('Invalid packet length');
        end


        % The length of packet transmission in subslot or slots. Packet length is signalled
        % numerical value plus one subslot or slot. The length of the subslot is 5 OFDM
        % symbols as defined in ETSI TS 103 636-3 [3].

        % Beamforming

        mac_meta_tx.tm_mode_0_to_11 = 0; 	% Table PHY 7.2-1, mode determines MIMO modes closed/open loop mode, values range from 0 to 11
        mac_meta_tx.codebook_index = 0;                     % 6.3.4, any value other than 0 makes packet beamformed, throws error if out of bound (depends on tm_mode_0_to_11)



        %Modulation and coding
        mac_meta_tx.mcs_index = mcs_index;                  % Table A-1 in part 3, values range from 0 to 11
        mac_meta_tx.Z = 2048;

        % Details of how channel modeling works
        % The oversampling setup is broken HARD, setting this to anything but 2
        % breaks stuff
        mac_meta_tx.oversampling = 2;                    	% By how much do we oversample our ofdm packet compared to critical sampling (insert zeros at specturm edges before IFFT)?

        mac_meta_tx.PLCF_type = 1;                          % Type 1 is 40 bits, Type 2 is 80 bits

        mac_meta_tx.rv = 0;                                 % HARQ version, values range from 0, 1, 2 to 3 (right HARQ retransmission order is 0 2 3 1)
        mac_meta_tx.network_id = de2bi(1e6,32,'left-msb');  % 7.6.6 must be given as a 32 bit vector with network_id(1) being the MSB, network_id must be known for scrambler on PHY
        mac_meta_tx.stf_version = 1;
        mac_meta_tx.stf_cover_sec_enable = 1;

        phy_4_5 = lib_util.run_chapter_4_5(1, mac_meta_tx);
        N_TB_bits = phy_4_5.N_TB_bits;
        N_CB = phy_4_5.N_CB;

        bits = randi([0,1], N_TB_bits, 1);
        bits = lteCRCEncode(bits,'24A');
        if N_CB > 1
            c_r = lib_6_generic_procedures.Code_block_segmentation_Z_2048(bits);
        else
            c_r = {bits};
        end

        %assert(validation_tbs(val_idx) == N_TB_bits);
        cb_len = cellfun(@length, c_r)

        fprintf("packet_len: %d, mcs: %d, total_bits:%d, tbs: %d, cbs: %d, codeblocks:%d \n",packet_length_subslots*5,  mcs_index, ...
		phy_4_5.n_total_bits, N_TB_bits, round(mean(cb_len)), N_CB);
        val_idx = val_idx+ 1;

        d_turbo_tx = lteTurboEncode(c_r);
        %disp('d_turbo at TX after encoding')
        %cellfun(@length,d_turbo_tx)

        %chs.Modulation = phy_4_5.mcs.modulation;

        %f = lteRateMatchTurbo(d_turbo_tx, phy_4_5.n_total_bits, redundancy_version, chs);
        %disp('bits after rate-matching');
        %size(f)
        %cbsbuffers=d_turbo_tx;
        %d_turbo = lteRateRecoverTurbo(f, N_TB_bits, redundancy_version, chs, cbsbuffers);
        %disp('d_turbo rate recovery at RX')
        %cellfun(@length,d_turbo)
        %break
    end
    %break
end


