% Returns the cover sequence for STF given u (mu) parameter. Errors for
% invalid u.
% According to chapter 6.3.7
function [cov_seq]= STF_cover_sequence(u)
    switch  u
        case 1
            cov_seq = [1; -1; 1; 1; -1; -1; -1];
        case {2,4,8}
            cov_seq = [1; -1; 1; 1; -1; -1; -1; -1; -1];
        otherwise
            error("Invalid value for u");
    end

end