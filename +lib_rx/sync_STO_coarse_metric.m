function metric = sync_STO_coarse_metric(window, M, L, detection_E_rms_threshold, u, stf_cover_sec_mode)

    % separate into repetitive pattern
    window = reshape(window, M, L);

    % power for normalization
    E = sum(sum(window.*conj(window)));

    % we have to have a minimal power to assume there is a frame present at
    % all
    E_mean = sqrt(E/(M*L));
    if E_mean < detection_E_rms_threshold
        metric = 0;
        return;
    end

    
    % lookup cover sequence
    c_u = lib_6_generic_procedures.STF_cover_sequence(u);
     
    
    % correlation from pattern to pattern
    P = 0;
    for j=1:1:L-1
        % <modified>
        if stf_cover_sec_mode == 1
            P = P + sum((window(:,j).*c_u(j)).*conj(window(:,j+1).*c_u(j+1)));
        elseif stf_cover_sec_mode == 0
            P = P + sum(window(:,j).*conj(window(:,j+1)));
        end
        % </modified>
    end

    metric = (  (L/(L-1)) * abs(P)/E  )^2;

    % We assume the worst case:
    % We can't detect NaN and inf, instead we get maximum metric.
    if isnan(metric) || isinf(metric)
        metric = 1;
    end
end

