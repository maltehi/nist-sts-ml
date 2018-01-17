% bitStream: Bit stream(s) under evaluation. Column vector(s)
% n: Length of bit stream to be evaluated (default: length of input stream)
function [ results ] = frequencyTest( bitStream, n )

if nargin < 2
    n = size(bitStream,1);
end

% Truncate and make +/-1
bitStream = bitStream(1:n,:);
if isempty(find(bitStream < 0, 1))
    bitStream = 2*bitStream-1;
end

results.sum = sum(bitStream,1);
results.s_obs = abs(results.sum)/sqrt(n);
results.f = results.s_obs/sqrt(2);
results.p_value = erfc(results.f);

end

