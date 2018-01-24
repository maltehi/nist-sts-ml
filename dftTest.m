function [results] = dftTest(bitStream, n)
if nargin < 2
    n = size(bitStream,1);
end

% Truncate and make +/-1
bitStream = bitStream(1:n,:);
if isempty(find(bitStream < 0, 1))
    bitStream = 2*bitStream-1;
end

% FFT
S = fft(bitStream,n,1);

% Modulus of FFT sequence
M = abs(S(1:n/2,:));

% 95% bound
results.T = sqrt(-log(0.05)*n);

% Expected count
results.N_0 = 0.95*n / 2;

% Actual count
results.N_1 = sum(M < results.T,1);

% d
results.d = (results.N_1 - results.N_0) / sqrt(n*0.95*0.05/4);

% Statistics
results.p_value = erfc(abs(results.d)/sqrt(2));
results.pass_ratio = length(find(results.p_value >= 0.01)) / length(results.p_value);

end

