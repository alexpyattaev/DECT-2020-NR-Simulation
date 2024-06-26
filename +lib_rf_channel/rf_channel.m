classdef rf_channel < handle
    
    properties
        
        verbose;
        verbose_cp;             % Cyclic prefix (of any reference length) in samples (oversampling must be included).
                                % Only used when verbose is set to a value above 0 and channel is rayleigh/rician.
                                % Very useful for debugging to see if channel response is longer than cyclic prefix (or any reference length).
                                % That the channel response is longer than the CP can be missed easily without this debugging function.
        
        % for all channel types              
        snr_db;                 % signal-to-noise-ratio in dB
        spectrum_occupied;      % How much of the total spectrum do we occupy? This information is required for adding noise.
        N_TX;                   % number of transmit antennas
        N_RX;                   % number of receive antennas

        % awgn_random_source = 'global'
        %
        %   Noise is generated by using the global rng.
        %   By controlling the global rng, we can control the regeneration of noise.
        %   By not altering the global rng, we get new noise everytime.
        %
        % awgn_random_source = 'local'
        %
        %   Noise is generated by using a local random number stream 'awgn_randomstream', which has to be initialized externally.
        %   A call of reset_random_awgn() resets the awgn_randomstream and the same values of noise are regenerated.
        %   This funcionality is convinient for switching between independent channels with only one class instance:
        %       The comm.MIMOchannel uses the global rng, which is controlled externally.
        %       By setting the global rng to a value A, we generate a channel A and jump in time.
        %       By setting the global rng to a value B, we generate a channel B and jump in time.
        %       By resetting to A again, we can go back to channel A and jump in time.
        %       The noise, however, uses this local stream and is completely independent from channel A and B.
        %
        awgn_random_source;     % can be global or local
        awgn_randomstream;      % if awgn_random_source=local, this is the random stream for the awgn function         
        
        % for deterministic channel only (prefix d_)
        sto;                  % symbol timing offset im samples, must be >= 0
        cfo;                  % carrier frequency offset in Hertz
        

        
        r_seed;                 % the seed for the matlab Channel object
        
        
        samp_rate;            % system sample rate in Samples/s
        r_max_doppler;      	% max doppler in Hertz        
             
        ir_interpolation;        % do we use interpolation of the PDP?
       
        
        % variables set internally
        r_matlab_MIMO_obj;      % reference to the matlab object -> will be set in init function
        
        r_appendix;          	% appended samples by rayleigh channel -> will be set in init function
        
        % tools for debugging
        samples_antenna_rx_no_noise;        % channel output without white noise
        samples_antenna_rx_only_noise;      % only white noise channel ouput without the user signal
    end
    
    methods (Static = true, Access = public)
       
       % initialization is done outside of this function
       function obj = rf_channel()
           
            obj.verbose                 = 0;
            obj.verbose_cp          = [];            
                      
            
            obj.snr_db                  = [];
            obj.spectrum_occupied       = [];
            obj.N_TX                    = [];
            obj.N_RX                    = []; 
            
            obj.awgn_random_source      = [];
            obj.awgn_randomstream       = [];
                        
            obj.sto                   = [];
            obj.cfo                   = [];
            
            obj.samples_antenna_rx_no_noise = [];
            obj.samples_antenna_rx_only_noise = [];
       end
    end
    
    methods (Static = false, Access = public)
        
        function init_fading_channel(obj, path_gains, path_delays)     
            % Path gains should be in linear scale, path delays in seconds
            Ts = 1/obj.samp_rate;
                   
            % sort by delays, in case they are out of order
            [path_delays, I] = sort(path_delays);
            path_gains = path_gains(I);
            if obj.verbose > 1
                figure();               
                plot(path_delays, path_gains,'ro','LineWidth',3);                
            end

            if obj.ir_interpolation == true
                
                % how often can we sample the power delay profile at our system sample rate?
                n_points = ceil(path_delays(end) / Ts);
                path_delays_res = (0:1:n_points)*Ts;

                % empty container
                path_gains_res = zeros(size(path_delays_res));

                % assign each given path delay from the profiles above to the closest sampling point
                for qq = 1:1:numel(path_delays)
                   
                    [~,idx] = min(abs(path_delays_res - path_delays(qq)));

                    path_gains_res(idx) = path_gains_res(idx) + path_gains(qq);
                end

              

                % sanity check
                if abs(sum(path_gains_res) / sum(path_gains) - 1) > 1e-6
                    error('Power sums are not equal.');
                end
                path_gains = path_gains_res;     
                path_delays = path_delays_res;
            end
            n_points = numel(path_delays);
            % plot Path gains in dB and linear
            if obj.verbose > 1
                hold on;                          
                plot(path_delays, path_gains,'bx','LineWidth',3);
                title('Path Gains after interpolation');
                xlabel('Time');
                ylabel('Path Gain');
                legend('RAW', 'Interpolation');    
                s = max(path_delays) - min(path_delays);
                xlim([min(path_delays)-s/10, max(path_delays)+s/10]);
                ylim([-0.1, max(path_gains)*1.1]);
                grid on
            end                
           
            % normalize (Matlab also normalizes)
        
                path_gains = path_gains/sqrt(sum(path_gains.^2));
         
            
            
            % sanity checks
            if numel(path_delays) ~= numel(unique(path_delays))
                error('There are equal path delays.');
            end            
            if sum(path_delays < 0) > 0
                error('There are negaive path delays.');
            end
            if sum(isnan(path_gains)) > 0
                error('There are NaNs.');
            end
            if sum(~isfinite(path_gains)) > 0
                error('There are infinity values.');
            end
            
            % show some channel properties
            if obj.verbose > 0

                mean_tau_weighted = sum(pathDelays.*avgPathGains_linear)/sum(avgPathGains_linear);
                rms_tau = sqrt(sum(((pathDelays-mean_tau_weighted).^2).*avgPathGains_linear)/sum(avgPathGains_linear));

                disp('##### Channel Properties ######');
                fprintf('Sampling Rate: %f MS/s\n', obj.samp_rate/1e6);
                fprintf('Sampling time Ts: %f ns\n', Ts/1e-9);
                fprintf('Largest delay: %f ns\n', path_delays(end)/1e-9);
                fprintf('PDP sampling points: %d\n', n_points);
                fprintf('CP length: %f ns\n', obj.verbose_cp*Ts/1e-9);
                fprintf('CP number of samples: %d\n', obj.verbose_cp);
                fprintf('CP / latest delay: %f\n', obj.verbose_cp*Ts / path_delays(end));
                fprintf('Tau mean: %f ns\n', mean(path_delays)/1e-9);
                fprintf('Tau mean weighted: %f ns\n', mean_tau_weighted/1e-9);
                fprintf('Tau rms weighted: %f ns\n', rms_tau/1e-9);
                fprintf('Coherence Bandwidth: %f MHz\n', 1/rms_tau/1e6);
                fprintf('Occupied Bandwidth: %f MHz\n', obj.samp_rate*obj.spectrum_occupied/1e6);
            end
            
            % determine how many samples the channel is adding
            obj.r_appendix = ceil(max(path_delays)/Ts);

            % create channel object
            
                
              
                
                


        obj.r_matlab_MIMO_obj = comm.MIMOChannel(   'SampleRate', obj.samp_rate, ...
                                                    'PathDelays', path_delays, ...
                                                    'AveragePathGains', path_gains, ...
                                                    'NormalizePathGains', true,...
                                                    'FadingDistribution', 'Rayleigh',... % KFactor, DirectPathDopplerShift, DirectPathInitialPhase
                                                    'MaximumDopplerShift', obj.r_max_doppler, ... %'DopplerSpectrum', 'Jakes',...
                                                    'SpatialCorrelationSpecification', 'Separate Tx Rx',... %'NumTransmitAntennas', obj.N_TX, 'NumReceiveAntennas', obj.N_RX,...
                                                    'TransmitCorrelationMatrix', eye(obj.N_TX),...
                                                    'ReceiveCorrelationMatrix', eye(obj.N_RX),... % SpatialCorrelationMatrix
                                                    'AntennaSelection', 'off',...
                                                    'NormalizeChannelOutputs', false,...
                                                    'FadingTechnique', 'Sum of sinusoids', ...
                                                    'NumSinusoids', 48,...
                                                    'InitialTimeSource', 'Input Port',... % InitialTime
                                                    'RandomStream', 'mt19937ar with seed',...
                                                    'Seed', obj.r_seed,...
                                                    'PathGainsOutputPort',0,...
                                                    'Visualization', 'off'); %AntennaPairsToDisplay, PathsForDopplerDisplay, SamplesToDisplay

                
            
            
        end
        
        function [samples_antenna_rx] = pass_samples(obj, samples_antenna_tx)
            % sanity check
            [N_TX_samples, N_TX_check] = size(samples_antenna_tx);  
            if obj.N_TX ~= N_TX_check
                error('Expect %d tx antennas, received %d.', obj.N_TX, N_TX_check);
            end                     
            if obj.sto < 0
                error('STO must be positive or zero, is %d.', obj.sto);
            end

            % Normalize the TX signal power
            for i=1:1:obj.N_TX
                tx_sig_real_power = std(samples_antenna_tx(:,i));
                if tx_sig_real_power > 1.1 || tx_sig_real_power < 0.7
                    error('Problematic TX signal power out of reasonable bounds!')
                end
                
                samples_antenna_tx(:,i) = samples_antenna_tx(:,i) /  tx_sig_real_power;                    
            end
            
            % add sto samples on all sides and thereby desynchronize signal
            sto_i =real(floor(obj.sto));           
                
            samples_antenna_rx = [zeros(sto_i, obj.N_TX); samples_antenna_tx; zeros(sto_i, obj.N_TX)];
            %Save where the actual signal is for our AGC fake
            valid_signal_slice = (sto_i+1:(N_TX_samples+sto_i));
            err_phase = (obj.sto - sto_i) * 2*pi;
                            
            % add error phase
            if err_phase ~= 0
                samples_antenna_rx = exp(1i*err_phase)*samples_antenna_rx;
            end

            % this can be very time consuming
            %channel_time_in_seconds = 0; %TODO - make this autoincrement with every packet!!!
            %samples_antenna_rx = obj.r_matlab_MIMO_obj(samples_antenna_rx, channel_time_in_seconds);
            
            % the received signal is downconverted, which creates a phase coherent CFO usually much larger than Doppler
            if obj.cfo ~= 0

                n_samples_per_antenna = numel(samples_antenna_rx(:,1));

                time_base = 0:(n_samples_per_antenna-1);
                time_base = time_base';

                for i=1:1:obj.N_RX
                    samples_antenna_rx(:,i) = samples_antenna_rx(:,i).*exp(1i*2*pi*obj.cfo*time_base);
                end
            end
                
            % Perform gain adaptation at each RX (as that is something RX would do anyway)                
            for i=1:1:obj.N_RX
                rx_sig_real_power = std(samples_antenna_rx(valid_signal_slice,i));
                if rx_sig_real_power > 2.1 || rx_sig_real_power < 0.3
                    error('Problematic RX signal power out of reasonable bounds!')
                end                
                samples_antenna_rx(:,i) = samples_antenna_rx(:,i) /  rx_sig_real_power;                    
            end
            
            % save RX signal without noise
            obj.samples_antenna_rx_no_noise = samples_antenna_rx;
            
            % Final step is to add noise to each rx antenna.
            % The SNR with value snr_db refers to the SNR of a single subcarrier, not the total signal power relative to the total noise power.
            if obj.snr_db < 50.0 % if snr is over 40dB we might just as well have perfect signal
                if strcmp(obj.awgn_random_source,'global') == true
                    for i=1:1:obj.N_RX
                        samples_antenna_rx(:,i) = awgn(samples_antenna_rx(:,i), obj.snr_db);
                        %samples_antenna_rx(:,i) = samples_antenna_rx(:,i) + randn(size(samples_antenna_rx(:,i))) / db2mag(6);
                    end
                elseif strcmp(obj.awgn_random_source,'local') == true
                    for i=1:1:obj.N_RX                        
                        samples_antenna_rx(:,i) = awgn(samples_antenna_rx(:,i), obj.snr_db, 0.0, obj.awgn_randomstream);
                    end
                end
            end
            
            % save only the noise
            obj.samples_antenna_rx_only_noise = samples_antenna_rx - obj.samples_antenna_rx_no_noise;
        end
        
        function [] = reset_random_rayleigh_rician(obj)            
            obj.r_matlab_MIMO_obj.reset();            
        end

        function [] = reset_random_awgn(obj)
          	if strcmp(obj.awgn_random_source, 'local') == true
                obj.awgn_randomstream.reset();
            end
        end
    end
end
