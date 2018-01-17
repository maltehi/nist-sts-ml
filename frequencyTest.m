% bitStream: Bit stream(s) under evaluation. Column vector(s)
% n: Length of bit stream to be evaluated (default: length of input stream)
function [ results ] = frequencyTest( bitStream, n )

if nargin < 2
    n = size(bitStream,1);
end

results.sum = sum(2*double(bitStream(1:n,:))-1,1);
results.s_obs = abs(results.sum)/sqrt(n);
results.f = results.s_obs/sqrt(2);
results.p_value = erfc(results.f);

end

