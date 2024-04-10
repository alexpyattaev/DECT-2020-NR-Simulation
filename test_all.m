% A simple wrapper to run the simulator without cluster harness

X.num_packets=10;
X.max_harq_retransmissions=0;
X.channel_type= 'awgn';
X.channel_trace_file= 'traced_complex.mat';
X.mcs_index=3;
X.snr_dB=15;
X.sto_shift=150;
X.cfo_shift= 0;
X.doppler_Hz= 0;
X.packet_length_subslots= 2;
X.stf_sync= 1;
X.stf_seq_mode= 1;
X.stf_cover_sec_mode= 0;
X.drs_sync= 1;
X.channel_estimation= 'wiener';
X.tm_mode_0_to_11= 0;
X.seed= 1;

save('test.mat','-struct', 'X');
main_standalone('test.mat');
res = load('results_test.mat');

res.results
