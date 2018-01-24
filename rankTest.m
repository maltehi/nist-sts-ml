function [results] = rankTest(bitStream, n)
if nargin < 2
    n = size(bitStream,1);
end
nStreams = size(bitStream,2);

M = 32; Q = 32;
N = floor(n / (M*Q));

% Truncate and make 0/1 (bit representation might be irrelevant here?)
bitStream = bitStream(1:M*Q*N,:);
if ~isempty(find(bitStream < 0, 1))
    bitStream = double(bitStream > 0);
end

% Not enough elements
if ~N
    results.p_value = zeros(1,nStreams);
    results.pass_ratio = 0;
    return
end

% Probabilities
% r = 32; i = 0:(r-1);
% p_r = pow2(r*(Q+M-r)-M*Q) * prod((1-pow2(i-Q))*(1-pow2(i-M))./(1-pow2(i-r)));
p_32 = prod(1-pow2((0:31)-Q));
r = 31; i = 0:(r-1);
p_31 = pow2(r*(Q+M-r)-M*Q) .* prod((1-pow2(i-Q)).*(1-pow2(i-M)) ./ (1-pow2(i-r)));
p_30 = 1 - (p_32 + p_31);

% Rank evaluation
R = zeros(N,nStreams);
bitMat = reshape(bitStream,M,Q,N,nStreams);    
for i = 1:nStreams
    for j = 1:N
        R(j,i) = rank(bitMat(:,:,j,i));
    end
end
F_32 = sum(R == 32,1);
F_31 = sum(R == 31,1);
F_30 = N - (F_32 + F_31);

% Statistics
results.bits_truncated = n - M*Q*N;
results.chi_squared = ((F_32 - N*p_32).^2 / p_32 ...
                      + (F_31 - N*p_31).^2 / p_31 ...
                      + (F_30 - N*p_30).^2 / p_30) / N;		
results.p_value = exp(-results.chi_squared/2);
results.pass_ratio = length(find(results.p_value >= 0.01)) / length(results.p_value);

end

