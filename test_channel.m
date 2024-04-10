close all;
mac_meta_tx.u = 1;                                  % mu = 1, 2, 4 or 8
mac_meta_tx.b = 1;                                  % beta = 1, 2, 4, 8, 12 or 16
mac_meta_tx.PacketLengthType = 0;                   % 0 for subslots, 1 for slots
mac_meta_tx.PacketLength = 2;                       % min is 1, max is 16 according to Table 6.2.1-2a in part 4
mac_meta_tx.tm_mode_0_to_11 = 0;                 	% Table 7.2-1, mode determines wether transmission is closed loop or not, values range from 0 to 11
mac_meta_tx.mcs_index = 0;
mac_meta_tx.Z = 6144;                               % 5.3 -> so far only Z=6144 fully supported, 2048 only at TX, RX missing (Matlab has no option for Z=2048 in LTE toolbox)
mac_meta_tx.oversampling = 2;                    	% By how much do we oversample our ofdm packet compared to critical sampling (insert zeros at specturm edges before IFFT)?
mac_meta_tx.codebook_index = 0;                     % 6.3.4, any value other than 0 makes packet beamformed, throws error if out of bound (depends on tm_mode_0_to_11)
mac_meta_tx.PLCF_type = 2;                          % Type 1 is 40 bits, Type 2 is 80 bits
mac_meta_tx.rv = 0;                                 % HARQ version, values range from 0, 1, 2 to 3 (right HARQ retransmission order is 0 2 3 1)
mac_meta_tx.network_id = de2bi(1e6,32,'left-msb');  % 7.6.6 must be given as a 32 bit vector with network_id(1) being the MSB, network_id must be known for scrambler on PHY
mac_meta_tx.stf_mode = 0;        % STF shwitch between base and update sequences 0 - legacy STF 2020, 1 - current STF 2023
mac_meta_tx.stf_cover_sec_mode = 0;  % 1 if cover sequence should be enabled


% temporary restrictions
if mac_meta_tx.Z ~= 6144
    error('Z must be 6144.');
end    
    
% create tx
verbose = 0;
tx = dect_tx(verbose, mac_meta_tx);


ch                      = lib_rf_channel.rf_channel();
ch.verbose              = 5;
ch.verbose_cp           = tx.phy_4_5.numerology.N_b_CP*tx.mac_meta.oversampling;
ch.samp_rate        	= tx.phy_4_5.numerology.B_u_b_DFT*tx.mac_meta.oversampling;
ch.amp                  = 1.0;
ch.noise                = true;
ch.snr_db             	= 5;
ch.spectrum_occupied    = tx.phy_4_5.n_spectrum_occupied/tx.mac_meta.oversampling;
ch.N_TX                	= 1;
ch.N_RX               	= 1;
ch.awgn_random_source   = 'global';
ch.awgn_randomstream 	= RandStream('mt19937ar','Seed', randi(1e9,[1 1]));
ch.sto                = 10; 
ch.cfo               	= 100; 
ch.err_phase         	= 0.0; % TODO: figure out what this means
ch.ir_interpolation = true;

%path_delays = [10e-7, 11e-7 , 13e-7, 15e-7];
%path_gains = [1.0, 0.7, 0.2, 0.05];
[path_delays, path_gains] = lib_rf_channel.rt_loader('traced_complex.mat');
ch.init_fading_channel(path_gains,path_delays);
