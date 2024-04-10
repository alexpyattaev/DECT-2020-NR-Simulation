function [samples_antenna_sto_cfo, STO_CFO_report] = sync_STF(  verbose,...
                                                                u,...
                                                                N_b_DFT,...
                                                                samples_antenna,...
                                                                STF_templates,...
                                                                n_packet_samples,...
                                                                oversampling,...
                                                                sto_config,...
                                                                cfo_config,...
                                                                stf_version,...
                                                                stf_cover_sec_enable)
    %TODO: kill legacy STF stuff for good
    %% precalculation of parameters that are known at the receiver a-priori
    %verbose=5
    STF_template_td = zeros(4, length(STF_templates.time_domain{1}));
    switch stf_version
        case 0
            N_stf = 4;
        case 1
            N_stf= 1;
        otherwise
            error('Unsupported STF version')
    end

    % Fetch the STF templates 
    for t = 1:N_stf        
    
        %if stf_cover_sec_mode == 1
        %    [STF_template_td(t,:)] = lib_6_generic_procedures.STF_signal_cover_sequence( ...
        %         STF_templates.time_domain{t}, u, oversampling);
        %end
        STF_template_td(t,:) = STF_templates.time_domain{t};
        % normalize power of stf
        STF_template_td(t,:) = STF_template_td(t,:)/rms(STF_template_td(t,:));
    end

    % decompose the shape of the data
    [n_samples_ch, N_RX] = size(samples_antenna);
    stf_len_samples = length(STF_template_td(1,:));
    stf_plot_samples = 64;
    if verbose && N_RX == 1
        figure('PaperSize',[14,8]);
        % subplot(5,1,1);
        % plot(real(samples_antenna(150:150+stf_plot_samples)))
        % hold on;
        % plot(imag(samples_antenna(150:150+stf_plot_samples)))
        % ylim([-2,2]);

        for i = 1:N_stf
            subplot(N_stf,1,i);            
            
            plot(real(STF_template_td(i,1:stf_plot_samples)))
            hold on;
            plot(imag(STF_template_td(i,1:stf_plot_samples)))
            plot(abs(STF_template_td(i,1:stf_plot_samples)),'-k', 'LineWidth',2);
            ylim([-1.5,1.7]);
            ylabel(sprintf('STF %d', i));
            xline(32, ':k');  
        end
        xlabel('Sample number');
    end

    
    if stf_version~=1
        error('The following code only works with new STF version!')
    end
    % STF pattern
    L = sto_config.n_pattern;

    % what is the length of the STF with oversampling included?
    n_samples_STF_os = sto_config.n_samples_STF_b_os;

    % sanity check
    if mod(n_samples_STF_os, L) ~= 0
        error("STF length not a multiple of L.");
    end

    % what is the length of a single pattern?
    M = n_samples_STF_os/L;

    %% detection by coarse metric threshold crossing

    % For each RX antenna, we go over all samples and calculate a metric until we cross a predefined threshold.
    % We need the sample index of this theshold crossings.
    coarse_stf_start_idx = zeros(N_RX, 1);
    if verbose > 1
        figure()
        hold on;
        title('Coarse sync peak search');
    end
    for i_rx=1:1:N_RX
        for i_stf = 1:N_stf
            CV = abs(conv(samples_antenna(:,i_rx), conj(flip(STF_template_td(i_stf,:))), 'same'));
            
            
            if verbose > 1
                plot(abs(CV));                
            end

           [maxv, maxi] = max(CV);
            if (maxv > stf_len_samples / 2) && maxi >= (stf_len_samples/2 - 1)
                coarse_stf_start_idx(i_rx) = maxi - (stf_len_samples/2 - 1);
                if verbose > 1
                    xline(maxi, 'r');  
                end
                break
            end
        end        
    end

    % 
    %     if verbose > 1
    %             figure()
    %             hold on
    %     end
    %     % Look for signal in any of the attached antennae
    %     for i=1:1:N_RX
    % 
    %         % get all samples of this particular antenna
    %         samples_antenna_single = samples_antenna(:,i);
    % 
    %         % container for detection metric
    %         metric = zeros(n_samples_ch - 2*n_samples_STF_os, 1);
    % 
    %         % the step must be small enough to find every STF at every SNR
    %         for k = 1 : sto_config.threshold.step : numel(metric)
    % 
    %             % determine the metric at this samples index
    %             metric(k) = lib_rx.sync_STO_coarse_metric(...
    %                 samples_antenna_single(k : k + n_samples_STF_os - 1), M, L, sto_config.minimum_power_threshold, u, stf_cover_sec_mode);
    % 
    %             % the threshold must be small enough to find every STF at every SNR
    %             if metric(k) >= sto_config.threshold.value
    % 
    %                 % save where we have crossed the threshold, this is an absolute index
    %                 coarse_threshold_crossing_idx(i) = k;
    %                 fprintf("Found STF around index %d",k);
    %                 % abort the search for this antenna
    %                 break;
    %             end
    %         end
    %         % if we have not made a selection based on threshold, just pick
    %         % highest value. FIXME this is not how RX works!!!
    %         if coarse_threshold_crossing_idx(i) == 0
    %             [~,coarse_threshold_crossing_idx(i)] = max(metric);
    %         end
    %         if verbose > 1
    %             plot(abs(metric))
    %             xline(coarse_threshold_crossing_idx(i), 'r');            
    %             title(append('Coarse signal detection Threshold Metric for antenna ', num2str(i)));            
    %         end
    % 
    %     end
    % 
    %     % keep the antennas at which we found a preamble
    %     coarse_threshold_crossing_idx = coarse_threshold_crossing_idx(coarse_threshold_crossing_idx > 0);
    % 
    %     % In a real receiver, we will start the rest of the synchronization algorithm once we hit a preamble at ANY antenna.
    %     % This corresponds in our case to finding the smallest index across all antennas.
    %     % If we didn't find any preamble at all, assume index 1, which will lead to a incorrect decoding if STO is sufficiently large.
    %     if isempty(coarse_threshold_crossing_idx) == true
    %         fprintf('Could not complete coarse synchronization!');
    %         coarse_threshold_crossing_idx = 1;
    %     else
    %         coarse_threshold_crossing_idx = min(coarse_threshold_crossing_idx);
    %     end
    % 
    % 
    % %% after detection, extract the range of samples what will be required for the upcoming synchronization steps
    % 
    % req_start = coarse_threshold_crossing_idx + 1;
    % search_window_len = sto_config.coarse_peak.search_length + sto_config.fine.search_area + n_samples_STF_os;
    % req_end = req_start + search_window_len;
    % 
    % % if we overreach, set index to 1 which will very likely to an erroneous packet (if the STO is suffiently large)
    % if req_end > n_samples_ch
    %     req_start = 1;
    %     req_end = req_start + search_window_len;
    % end
    % 
    % samples_antenna_required = samples_antenna(req_start:req_end, :);
    % 
    % %% coarse peak search
    % 
    % % We now know we most likely found an STF, now we need to find the peak of the coarse detection metric.
    % % We start searching for the peak right after the index threshold_cross_idx.
    % % Peaks are determined for each rx antenna individually.
    % 
    % % how strong is the detected packet?
    % coarse_peak_height = zeros(N_RX, 1);
    % 
    % % where was the packet detected? this is now an index relative to samples_antenna_required
    % coarse_peak_idx = ones(N_RX, 1);
    % 
    % for i=1:1:N_RX
    % 
    %     % get samples of this particular antenna
    %     samples_antenna_single = samples_antenna_required(:,i);
    % 
    %     % container for detection metric
    %     metric = zeros(sto_config.coarse_peak.search_length, 1);
    % 
    %     % we continue the search with a much smaller step size
    %     for m = 1 : sto_config.coarse_peak.step : sto_config.coarse_peak.search_length
    % 
    %         % We found the packet too late and we leave the rest of the matrix at zero.
    %         % This can't happen in a real receiver, but here we have a limited amount of samples.
    %         % Nevertheless, this should barely ever happen.
    %         if m + n_samples_STF_os-1 > n_samples_ch
    %             break;
    %         end
    % 
    %         % determine the metric at this samples index
    %         metric(m) = lib_rx.sync_STO_coarse_metric(...
    %             samples_antenna_single(m:m+n_samples_STF_os-1), M, L, sto_config.minimum_power_threshold, u, stf_cover_sec_mode);
    %     end
    % 
    %     % warning: if sto_config.coarse_peak.step > 1, we apply the moving average across zeros!
    % 
    %     % <modified> TODO THIS IS PROBABLY REDUNDANT
    %     % we apply some filtering to smooth the coarse detection
    %     if stf_cover_sec_mode == 0
    %         metric_mm = movmean(metric, sto_config.coarse_peak.movmean);
    %     elseif stf_cover_sec_mode == 1
    %         metric_mm = metric;
    %     end
    %     % </modified>
    % 
    %     % finally we search for the maximum
    %     [cph,cpi] = max(metric_mm);
    % 
    %     % save values, these values are relative to samples_antenna_required
    %     coarse_peak_height(i) = cph;
    %     coarse_peak_idx(i) = cpi;
    % 
    %     % debugging
    %     if  verbose > 1
    %         figure()
    %         clf()
    %         hold on
    % 
    %         subplot(2,1,1)
    %         plot(abs(metric));
    %         xline(cpi, 'r');
    %         str = append('Coarse Synchronization Metric for antenna ', num2str(i));
    %         title(str);
    %         ylim([0 1.2]);
    % 
    %         subplot(2,1,2)
    %         plot(abs(metric_mm));
    %         xline(cpi, 'r');
    %         str = append('Coarse Synchronization Metric after movmean for antenna ', num2str(i));
    %         title(str);
    %         ylim([0 1.2]);
    %         xline(cpi, 'r');
    %     end
    % end
    % 
    % 
    % 
    % 
    % %% coarse peak selection ( this is only relevant for mimo case)
    % 
    % % we have found a peak for each rx antenna, now we have to pick the optimal peak index
    % 
    % % we consider only peaks above a threshold
    % idx_peak_high = coarse_peak_height >= sto_config.coarse_peak.threshold;
    % 
    % % remove peaks that are too low
    % coarse_peak_height = coarse_peak_height(idx_peak_high);
    % coarse_peak_idx    = coarse_peak_idx(idx_peak_high);
    % 
    % % it can actually happen that we have no peak above the threshold, e.g. for a false alarm
    % if isempty(coarse_peak_height) 
    %     fprintf("Failed to lock onto coarse peak during peak selection stage")
    %     coarse_sync_idx = 1;
    % else
    %     % the sync point are the peak indices weighted by the height
    %     normalization = coarse_peak_height/sum(coarse_peak_height);
    %     coarse_sync_idx = sum( coarse_peak_idx .* normalization);
    %     coarse_sync_idx = floor(coarse_sync_idx);
    % 
    %     % make sure every the indices are within range
    %     if coarse_sync_idx <= 0
    %         coarse_sync_idx = 1;
    %     end
    %     if coarse_sync_idx >= sto_config.coarse_peak.search_length
    %         coarse_sync_idx = sto_config.coarse_peak.search_length;
    %     end
    %end
    %%%%%%%===========
    % Disabling useless code from Max
    %% fractional cfo calculation

    % create a container with the coarsly synchronized STFs only
    samples_antenna_stf_at_coarse_peak = zeros(n_samples_STF_os, N_RX);
    for i=1:1:N_RX
        %samples_antenna_stf_at_coarse_peak(:,i) = samples_antenna_required(coarse_sync_idx : coarse_sync_idx + n_samples_STF_os - 1, i);
        samples_antenna_stf_at_coarse_peak(:,i) = samples_antenna(coarse_stf_start_idx:coarse_stf_start_idx+n_samples_STF_os-1, i);
    end
    if verbose
        figure()
        plot(abs(samples_antenna_stf_at_coarse_peak(:,1)))
    end

    if cfo_config.active_fractional == true
        % use the STFs and determine a fractional CFO
        [samples_antenna_stf_at_coarse_peak,...
            fractional_cfo_report] = lib_rx.sync_CFO_fractional(...
            verbose, samples_antenna_stf_at_coarse_peak, n_samples_STF_os, u, cfo_config, stf_cover_sec_enable);
    else        
        fractional_cfo_report = 0;
    end

    %% integer cfo calculation

    if cfo_config.active_integer == true
        integer_cfo_report = lib_rx.sync_CFO_integer(   samples_antenna_stf_at_coarse_peak,...
                                                        STF_templates.freq_domain,...
                                                        N_b_DFT,...
                                                        oversampling,...
                                                        sto_config,...
                                                        cfo_config,...
                                                        u,...
                                                        stf_cover_sec_enable);
    else
        integer_cfo_report = 0;
    end

    %% perform crosscorrelation with STF templates for N_eff_TX and fine synchonization pointer determination

    % we have a coarse sync point, we have also determined the fractional + integer CFO

    % first create a time base which we will use to correct the fractional + integer CFO
    time_base = 0:1:(2*sto_config.fine.search_area + n_samples_STF_os - 1);
    time_base = time_base';
    n_search_area = length(time_base);

    % N_eff_TX = {1,2,4,8}, for each we will check one STF, search range can be limited by radio device class
    metric_maxima = zeros(4, N_RX);
    metric_maxima_index = zeros(4, N_RX);

    % every antenna has the same search range
    search_area_min = max(1, coarse_stf_start_idx - sto_config.fine.search_area);
    search_area_max = search_area_min + n_search_area-1;

    for i=1:1:N_RX
        
        % get the samples of this particular antenna around STF location
        samples_antenna_single = samples_antenna(search_area_min:search_area_max,i);
        % derotate samples in search range
        total_cfo = fractional_cfo_report + integer_cfo_report * 1/(N_b_DFT*oversampling);
        samples_antenna_single = samples_antenna_single.*exp(1i*2*pi*(-total_cfo)*time_base);
        % try every possible stf type
        for j=1:1:N_stf
 
            % perform a cross correlation between the samples and the stf template
            metric = abs(lib_rx.sync_xcorr(samples_antenna_single, STF_template_td(j, :)));

            % debugging
            if verbose > 2
                figure()
                plot(abs(metric));
                if stf_cover_sec_enable
                    s = 'with cover sequence';
                else
                    s = 'no cover sequence';
                end
                str = append('Fine sync for ant ', num2str(i), ' stf# ', num2str(j), ', ', s);
                ylabel('Metric')
                xlabel('Subsample offset')
                title(str);
            end

            % this maximum is local in the correlation window!
            [~, max_idx] = max(abs(metric));

            % save the maximum
            metric_maxima(j,i) = abs(metric(max_idx));

            % this index is converted to global
            max_idx = max_idx + search_area_min - 1;

            % save global maximum
            metric_maxima_index(j,i) = max_idx;
        end
    end
    
    %% find best fitting N_eff_TX

    % average peaks across antennas
    metric_maxima_sum = sum(metric_maxima, 2);
    
    % find the best fitting stf
    [~, j_best] = max(metric_maxima_sum);
    
    % how many effective antennas do we have?
    switch j_best
        case 1
            N_eff_TX_report = 1;
        case 2
            N_eff_TX_report = 2;
        case 3
            N_eff_TX_report = 4;
        case 4
            N_eff_TX_report = 8;
    end

    %% find best fitting fine synchronization point

    metric_maxima_fine = metric_maxima(j_best, :);

    max_idx_fine = metric_maxima_index(j_best, :);

    % weight a detection point by the height of the coarse metric of the respective antenna
    if sum(metric_maxima_fine) > 0
        normalization = metric_maxima_fine/sum(metric_maxima_fine);
    else
        normalization = 1/numel(metric_maxima_fine);
    end
    max_idx_fine = sum(max_idx_fine .* normalization);
    max_idx_fine = floor(max_idx_fine);

    % minimum is 1
    if max_idx_fine <= 0
        max_idx_fine = 1;
    end

    

    %% extract packet starting at the fine synchronization point and correct the CFO for the entire packet
    
    % packet is longer when oversampling
    n_packet_samples = n_packet_samples*oversampling;

    % we need a new time base with the length of one frame
    time_base = 0:1:(n_packet_samples - 1);
    time_base = time_base';
    
    % output container
    samples_antenna_sto_cfo = zeros(n_packet_samples, N_RX);

    % final index of packet
    max_idx_fine_plus_packet = max_idx_fine + n_packet_samples - 1;

    % make sure we don't overreach
    if max_idx_fine_plus_packet > n_samples_ch
        max_idx_fine = n_samples_ch - n_packet_samples + 1;
        max_idx_fine_plus_packet = n_samples_ch;
    end

    for i=1:1:N_RX
        % derotate samples
        total_cfo = fractional_cfo_report + integer_cfo_report * 1/(N_b_DFT*oversampling);
        samples_antenna_sto_cfo(:,i) = samples_antenna(max_idx_fine : max_idx_fine_plus_packet, i) .* exp(1i*2*pi*(-total_cfo)*time_base);
    end

    %% create output report
    STO_CFO_report.max_idx_coarse   = coarse_stf_start_idx;
    STO_CFO_report.fractional_cfo   = fractional_cfo_report;
    STO_CFO_report.integer_cfo      = integer_cfo_report;
    STO_CFO_report.N_eff_TX         = N_eff_TX_report;
    STO_CFO_report.max_idx_fine     = max_idx_fine;
    %error('boo')
end

