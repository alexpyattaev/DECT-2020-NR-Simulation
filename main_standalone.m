% Help Andrey
% Verify the block sizes!
% Add SCM support for fading channel modeling
% Model SNR sweep in fading channel and re-running AWGN

%Use real PCC maybe?

function results = main_standalone(matfile)
verbose = 0;
params = load(matfile);
params
% NEVER run without STO.
if params.sto_shift == 0
    params.sto_shift = 150
end
%params.cfo_shift = 0.0;
close all;
rng(params.seed);

num_packets = params.num_packets;
%num_packets = 1;
%params.snr_dB = 15;
params.stf_seq_version = 1;  %0=2020 or 1=2023
params.stf_cover_sec_enable = 1; % this is not optional in current spec
%warning('off');

% Current configuration is for 1.7MHz channel
% these variables need to be set before creating tx and rx
mac_meta_tx.u = 1;                                  % mu = 1, 2, 4 or 8
% Choice of u sets the subcarrier spacing
mac_meta_tx.b = 1;                                  % beta = 1, 2, 4, 8, 12 or 16
numerology = lib_4_physical_layer_principles.numerologies(mac_meta_tx.u,mac_meta_tx.b);


if  params.packet_length_subslots <=0
    error('Invalid packet length')
end
%PacketLengthType: 0 for subslots, 1 for slots
%PacketLength: min is 1, max is 16 according to Table 6.2.1-2a in part 4
if mod(params.packet_length_subslots, numerology.N_SLOT_u_subslot) % need to operate in subslots
    mac_meta_tx.PacketLengthType = 0;
    mac_meta_tx.PacketLength = double(params.packet_length_subslots); 
else
    mac_meta_tx.PacketLengthType = 1;
    
    mac_meta_tx.PacketLength = double(params.packet_length_subslots) / numerology.N_SLOT_u_subslot; 
  
end

  if mac_meta_tx.PacketLength > 16 
        error('Invalid packet length');
    end
                   
                      
% The length of packet transmission in subslot or slots. Packet length is signalled
% numerical value plus one subslot or slot. The length of the subslot is 5 OFDM
% symbols as defined in ETSI TS 103 636-3 [3].

% Beamforming
tm_mode = lib_7_Transmission_modes.transmission_modes(params.tm_mode_0_to_11)
mac_meta_tx.tm_mode_0_to_11 = params.tm_mode_0_to_11; 	% Table PHY 7.2-1, mode determines MIMO modes closed/open loop mode, values range from 0 to 11
mac_meta_tx.codebook_index = 0;                     % 6.3.4, any value other than 0 makes packet beamformed, throws error if out of bound (depends on tm_mode_0_to_11)



%Modulation and coding
mac_meta_tx.mcs_index = params.mcs_index;                  % Table A-1 in part 3, values range from 0 to 11
mac_meta_tx.Z = 6144;                               % 5.3 -> so far only Z=6144 fully supported, 2048 only at TX, RX missing (Matlab has no option for Z=2048 in LTE toolbox)
% Z is turbo coder block size. Spec allows also for 2048. 

% Details of how channel modeling works
% The oversampling setup is broken HARD, setting this to anything but 2
% breaks stuff
mac_meta_tx.oversampling = 2;                    	% By how much do we oversample our ofdm packet compared to critical sampling (insert zeros at specturm edges before IFFT)?

mac_meta_tx.PLCF_type = 2;                          % Type 1 is 40 bits, Type 2 is 80 bits
%TODO: WTF IS PLCF???

mac_meta_tx.rv = 0;                                 % HARQ version, values range from 0, 1, 2 to 3 (right HARQ retransmission order is 0 2 3 1)
mac_meta_tx.network_id = de2bi(1e6,32,'left-msb');  % 7.6.6 must be given as a 32 bit vector with network_id(1) being the MSB, network_id must be known for scrambler on PHY
mac_meta_tx.stf_version = params.stf_seq_version;        % STF switch between base and update sequences 0 - legacy STF 2020, 1 - current STF 2023
mac_meta_tx.stf_cover_sec_enable = params.stf_cover_sec_enable;  % 1 if cover sequence should be enabled


% temporary restrictions
if mac_meta_tx.Z ~= 6144
    error('Z must be 6144.');
end    
    
% create tx
%verbose = 5;
tx = dect_tx(verbose, mac_meta_tx);
%verbose = 0;
% additional rx configuration
mac_meta_rx = mac_meta_tx;
mac_meta_rx.N_RX = tm_mode.N_SS; % Match number of RX antennas to transmission mode

% synchronization based on STF
mac_meta_rx.synchronization.stf.active = params.stf_sync; % default false
if mac_meta_rx.synchronization.stf.active == true
    % STO (detection, coarse peak search, fine peak search)
    mac_meta_rx.synchronization.stf.sto_config = lib_rx.sync_STO_param(mac_meta_tx.u, mac_meta_tx.b, mac_meta_tx.oversampling);    
    mac_meta_rx.synchronization.stf.sto_config.threshold.value = 0.7;
    %mac_meta_rx.synchronization.stf.sto_config.coarse_peak.threshold  = 0.07;
    
    % CFO (fractional, integer)
    mac_meta_rx.synchronization.stf.cfo_config = lib_rx.sync_CFO_param(mac_meta_tx.u);
    mac_meta_rx.synchronization.stf.cfo_config.active_fractional = false;
    mac_meta_rx.synchronization.stf.cfo_config.active_integer = false;
end
    
% synchronization based on DRS (residual CFO)
mac_meta_rx.synchronization.drs.cfo_config.active_residual = params.drs_sync; % default false

% channel estimation
mac_meta_rx.active_ch_estim_type = params.channel_estimation;

% channel equalization
mac_meta_rx.active_equalization_detection = true; % Enables equalization for data symbols. Should be true.

% create rx
rx = dect_rx(verbose, mac_meta_rx);



% how many bits does tx need?
N_TB_bits = tx.phy_4_5.N_TB_bits;



%warning('off');                      
% Counters for stats
power_cnt = 0;
energy_tx = 0;
energy_rx = 0;
n_packets_error = 0;        
n_bits_PCC_sent = 0;
n_bits_PCC_error = 0;
n_packets_PCC_error = 0;

n_bits_PDC_sent = 0;
n_bits_PDC_error = 0;
n_packets_PDC_error = 0;
                
% extract variables for wiener weights
if strcmp(rx.mac_meta.active_ch_estim_type,'wiener') == true
    physical_resource_mapping_DRS_cell  = tx.phy_4_5.physical_resource_mapping_DRS_cell;
    N_b_DFT                             = tx.phy_4_5.numerology.N_b_DFT;
    N_PACKET_symb                       = tx.phy_4_5.N_PACKET_symb;
    N_b_CP                              = tx.phy_4_5.numerology.N_b_CP;
    samp_rate                           = tx.phy_4_5.numerology.B_u_b_DFT;
    noise_estim_wiener               	= 1/10^(params.snr_dB/10);
    f_d_hertz_wiener                    = 20;
    tau_rms_sec_wiener                  = 363e-9;

    % the wiener filter for each snr
    rx.wiener = lib_rx.channel_estimation_wiener_weights(physical_resource_mapping_DRS_cell,...
                                                            N_b_DFT,...
                                                            N_PACKET_symb,...
                                                            N_b_CP,...
                                                            samp_rate,...
                                                            noise_estim_wiener, f_d_hertz_wiener, tau_rms_sec_wiener);
end

% create channel
ch                      = lib_rf_channel.rf_channel();
ch.verbose              = verbose;
ch.verbose_cp           = tx.phy_4_5.numerology.N_b_CP*tx.mac_meta.oversampling;

ch.snr_db             	= params.snr_dB;
ch.spectrum_occupied    = tx.phy_4_5.n_spectrum_occupied/tx.mac_meta.oversampling;
ch.N_TX                	= tx.phy_4_5.tm_mode.N_TX;
ch.N_RX               	= rx.mac_meta.N_RX;
ch.awgn_random_source   = 'global';
ch.awgn_randomstream 	= RandStream('mt19937ar','Seed', randi(1e9,[1 1]));
ch.sto                  = params.sto_shift; 
ch.cfo               	= params.cfo_shift; 

ch.samp_rate        	= tx.phy_4_5.numerology.B_u_b_DFT*tx.mac_meta.oversampling;
ch.r_max_doppler     	= params.doppler_Hz;                            % 1.946 19.458

ch.ir_interpolation      = false;
ch.r_seed = double(params.seed);

%TODO: use RT channel
% USE comm.Ray for loading model rays
if params.channel_type == "awgn"
    path_delays = [100e-9];
    path_gains = [1.0];
elseif params.channel_type == "rayleigh"
    path_delays = [100e-9];
    path_gains = [1.0];
elseif params.channel_type=="deterministic"
    [path_delays, path_gains] = lib_rf_channel.rt_loader(params.channel_trace_file);
end


ch.init_fading_channel(path_gains, path_delays);

% Loop for number of packets to be simulated
for j=1:1:num_packets

    % give rx handles so it can debug
    rx.tx_handle = tx;
    rx.ch_handle = ch;
    
    % generate random PCC bits
    if tx.mac_meta.PLCF_type == 1
        PCC_user_bits = randi([0 1], 40, 1);
    elseif tx.mac_meta.PLCF_type == 2
        PCC_user_bits = randi([0 1], 80, 1);
    end

    % generate bits
    PDC_user_bits = randi([0 1], N_TB_bits, 1);       
    
    % harq abort conditions
    pcc_decoded_successfully = false;
    pdc_decoded_successfully = false;

    for z=0:1:params.max_harq_retransmissions

        % there is a specific order for the redundancy version
        if mod(z,4) == 0
            tx.mac_meta.rv = 0; % initial transmission
            rx.mac_meta.rv = 0; % initial transmission
        elseif mod(z,4) == 1
            tx.mac_meta.rv = 2;
            rx.mac_meta.rv = 2;
        elseif mod(z,4) == 2
            tx.mac_meta.rv = 3;
            rx.mac_meta.rv = 3;
        elseif mod(z,4) == 3
            tx.mac_meta.rv = 1;
            rx.mac_meta.rv = 1;
        end

        % let tx create the packet
        samples_antenna_tx =  tx.generate_packet(PCC_user_bits, PDC_user_bits);

        % pass samples through channel    
        samples_antenna_rx = ch.pass_samples(samples_antenna_tx);
        % make next channel impulse response independent from this one
        ch.reset_random_rayleigh_rician();
        
        % measure powers
        power_cnt = power_cnt + 1;
        energy_tx = energy_tx + sum(sum(abs(samples_antenna_tx).^2))/size(samples_antenna_tx,1);
        energy_rx = energy_rx + sum(sum(abs(samples_antenna_rx).^2))/size(samples_antenna_rx,1);

        % Now let rx decode the frame.
        % Rx can do so because it's mac_meta is the exact same.
        [PCC_user_bits_recovered, PDC_user_bits_recovered] = rx.demod_decode_packet(samples_antenna_rx);
        
        % measure the BER uncoded
        
        n_bits_PCC_sent = n_bits_PCC_sent + numel(tx.packet_data.pcc_enc_dbg.d);
        n_bits_PCC_error = n_bits_PCC_error + sum(abs(double(tx.packet_data.pcc_enc_dbg.d) - double(rx.packet_data.pcc_dec_dbg.d_hard)));                
        
        n_bits_PDC_sent = n_bits_PDC_sent + numel(tx.packet_data.pdc_enc_dbg.d);
        n_bits_PDC_error = n_bits_PDC_error + sum(abs(double(tx.packet_data.pdc_enc_dbg.d) - double(rx.packet_data.pdc_dec_dbg.d_hard)));
        
        % we might be done
        if numel(PCC_user_bits_recovered) ~= 0
            pcc_decoded_successfully = true;
        end                

        % we might be done
        if numel(PDC_user_bits_recovered) ~= 0
            pdc_decoded_successfully = true;
        end
        
        % we continue sending retransmissions as long as not both were decoded correctly
        if pcc_decoded_successfully == true && pdc_decoded_successfully == true
            break;
        end
    end

    % delete harq buffer
    rx.harq_buf_40 = [];
    rx.harq_buf_80 = [];
    rx.harq_buf = [];
    
    
    % check if frame was decoded correctly
    if ~pdc_decoded_successfully || ~pcc_decoded_successfully
        n_packets_error= n_packets_error+1;

        if pcc_decoded_successfully == false
            n_packets_PCC_error = n_packets_PCC_error + 1;
        end            
            
        if pdc_decoded_successfully == false
            n_packets_PDC_error = n_packets_PDC_error + 1;
        end
    
    end
end



energy_tx_global = energy_tx/power_cnt;
energy_rx_global = energy_rx/power_cnt;




% save results

results.n_bits_PCC_sent = n_bits_PCC_sent;
results.n_bits_PCC_error = n_bits_PCC_error;
results.n_packets_sent = num_packets;
results.n_packets_PCC_error = n_packets_PCC_error; 
results.n_packets_error = n_packets_error; 
results.n_bits_PDC_sent = n_bits_PDC_sent;
results.n_bits_PDC_error = n_bits_PDC_error;
results.n_packets_PDC_error = n_packets_PDC_error;
results.n_packets_error = n_packets_error;
results.energy_rx= energy_rx_global;
results.energy_tx= energy_tx_global;
save(strcat("results_", matfile), "results","params");

end
