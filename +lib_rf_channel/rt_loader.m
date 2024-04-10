function [delays,powers] = rt_loader(fname)

data = load(fname);
delays = data.delays;
powers = data.powers;
% remove overall propagation time as it is not helpful for LLS
delays = delays - min(delays);
% normalize powers as well for same reason
powers = powers / max(powers);
end