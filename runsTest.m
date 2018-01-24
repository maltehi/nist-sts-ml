function [results] = runsTest(bitStream, n)

if nargin < 2
    n = size(bitStream,1);
end

% Truncate and make 0/1
bitStream = bitStream(1:n,:);
if ~isempty(find(bitStream < 0, 1))
    bitStream = (bitStream > 0);
end

% Proportion of ones (pi)
pi_r = sum(bitStream,1) / n;
estCrit = abs(pi_r - 0.5) > 2/sqrt(n);
% Zero all streams for which Pi condition is not fulfilled
bitStream(:,estCrit) = 0;

% Compute V_n
results.V = sum(xor(bitStream(1:end-1,:), bitStream(2:end,:)),1) + 1;

% Compute P-value
results.p_value = erfc(abs(results.V - 2*n * pi_r .* (1-pi_r)) ...
                       ./ (2*sqrt(2*n) * pi_r .* (1-pi_r)));

% Calculate pass ratio
results.pass_ratio = length(find(results.p_value >= 0.01))/length(results.p_value);

end

