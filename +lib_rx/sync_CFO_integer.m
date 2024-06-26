function [cfo_report] = sync_CFO_integer(   samples_antenna_sto,...
                                            STO_templates_freq_domain,...
                                            N_b_DFT,...
                                            oversampling,...
                                            sto_config,...
                                            cfo_config, ...
                                            u, ...
                                            stf_cover_sec_mode)

    %% <modified>
    if stf_cover_sec_mode == 1
        [samples_antenna_sto] = lib_6_generic_procedures.STF_signal_cover_sequence(samples_antenna_sto,u,oversampling);
    end
    % <modified>
    
    %% transform from time domain samples into frequency domain

    % remove cp
    samples_antenna_sto_no_cp = samples_antenna_sto(sto_config.n_samples_STF_cp_only_b_os + 1 : end, :);

    % transform into frequency domain
    samples_antenna_sto_freq_os = fft(samples_antenna_sto_no_cp);
    samples_antenna_sto_freq_os = fftshift(samples_antenna_sto_freq_os, 1);

    % %% <modified>
    % figure()
    % plot(abs(samples_antenna_sto_freq_os(:,1)))
    % figure()
    % plot(abs(samples_antenna_sto_freq_os(:,2)))
    % % <modified>

    %% at this point, we don't know N_eff_TX yet, so we can use any STF template as we have to perform a simple power detection

    STF_values_freq_domain = STO_templates_freq_domain{1};

    % Matlab saves mirrored version compared to DECT-2020 NR standard
    STF_values_freq_domain = flipud(STF_values_freq_domain);

    %% we need to expand the STF template if oversampling is used

    % zeros at spectrum edge
    n_oversampling_zeros = (oversampling-1)*N_b_DFT/2;
    STF_values_freq_domain = [zeros(n_oversampling_zeros,1); STF_values_freq_domain; zeros(n_oversampling_zeros,1)];

    %% now we have to test for each antenna, which integer CFO makes most sense

    % one metric value per possible integer CFO per RX antenna
    N_RX = size(samples_antenna_sto, 2);
    metric = zeros(numel(cfo_config.integer.candidate_values), N_RX);

    % go over each rx antenna
    for i=1:1:N_RX

        % we only need the STF symbol
        STF_this_antenna = samples_antenna_sto_freq_os(:,i);

        % try for each possible frequency offset
        idx = 1;
        for cfo_candidate = cfo_config.integer.candidate_values

            % Shift spectrum:
            %
            % Within the rf channel, the cfo is applied this way:
            %
            %       samples_out.*exp(1i*2*pi*cfo*time_base);
            %
            % When cfo is negative, the spectrum is shifted towards low indices.
            % For instance, when DC lies at index 33 with cfo=0, it will be at index 32 when cfo=-1*subcarrier_spacing.
            %
            % when cfo_candidate is negative, assume a negative integer CFO, so we have to push our STF towards high indices to compensate the CFO.
            % when cfo_candidate is positive, assume a positive integer CFO, so we have to push our STF towards low indices to compensate the CFO.
            STF_this_antenna_shifted = circshift(STF_this_antenna, -cfo_candidate);

            % calculate metrix
            metric(idx, i) = sum(abs(STF_values_freq_domain .* conj(STF_this_antenna_shifted)));

            % % <modified>
            % figure()
            % plot(abs(STF_values_freq_domain).*100)
            % hold on
            % plot(abs(STF_this_antenna_shifted))
            % hold off
            % % </modified>

            idx = idx+1;
        end
    end

    % %% <modified>
    % figure()
    % plot(metric(:,1))
    % figure()
    % plot(metric(:,2))
    % % <modified>
    
    %% combine the results from all antennas and determine the most likely integer CFO

    metric_abs = abs(metric);
    metric = sum(metric_abs,2);

    % %% <modified>
    % figure()
    % plot(metric)
    % % <modified>

    % find index of best fitting candidate
    [~, CFO_report_integer_index] = max(metric);
    %CFO_report_integer_index = 5;

    % extract corresponding CFO value
    cfo_report = cfo_config.integer.candidate_values(CFO_report_integer_index);
end

