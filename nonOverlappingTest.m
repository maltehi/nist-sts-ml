function [results] = nonOverlappingTest(bitStream, m, n, N, limitTemplates)
if nargin < 3
    n = size(bitStream,1);
end
if nargin < 4
    N = 8;
end
if nargin < 5
    limitTemplates = true;
end
nStreams = size(bitStream,2);

M = fix(n/N);

% Truncate and make 0/1 (bit representation might be irrelevant here?)
bitStream = reshape(bitStream(1:M*N,:),M,N,[]);
if ~isempty(find(bitStream < 0, 1))
    bitStream = double(bitStream > 0);
end

% Number of aperiodic templates for m = 1..21
if limitTemplates
    % In the C-code, the number of templates is limited to 148 by a constant.
    numOfTemplates = [0, 2, 4, 6, 12, 20, 40, 74, 148*ones(1,13)].';
else
    numOfTemplates = [0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116, ...
    			      2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152].';
end
              
% Read patterns from template file into columns of sequences
templateFile = fopen(sprintf('templates\\template%d',m),'r');
sequences = fscanf(templateFile,'%d',[m numOfTemplates(m)]);
fclose(templateFile);

% Theoretical mean and variance
mu = (M - m + 1) / pow2(m);
varWj = M * (pow2(-m) - (2*m-1)*pow2(-2*m));
              
% Compute probabilities
pi_r = zeros(6,1);
pi_r(1:5) = exp(-mu+(1:5)*log(mu)-gammaln(2:6));
% Add probability for i = 0 to first entry
pi_r(1) = pi_r(1) + exp(-mu);
% Fill up
pi_r(6) = 1 - sum(pi_r);

Wj = zeros(N,numOfTemplates(m),nStreams);

% Statistics
results.bits_truncated = n - N*M;
results.chi_squared = zeros(numOfTemplates(m),nStreams);
results.p_value = zeros(numOfTemplates(m),nStreams);

for kk = 1:nStreams
    for jj = 1:numOfTemplates(m)
        for  i = 1:N % Loop blocks
            W_obs = 0;
%             skip = 0;
            for j = 0:M-m % Loop bits of substring % !! Start from 0 here!
                match = 1;
%                 if skip
%                     skip = skip-1;
%                     continue;
%                 end
                for k = 1:m  % Loop bits of template. Does not jump blocks, but is fast and delivers equal results
                    if sequences(k,jj) ~= bitStream(j+k,i,kk)
                        match = 0;
                        break;
                    end
                end
                if match
                    W_obs = W_obs + 1;
%                     skip = m-1;
                end
            end
            Wj(i,jj,kk) = W_obs;
        end
        results.chi_squared(jj,kk) = sum((Wj(:,jj,kk) - mu).^2 / varWj, 1);
        results.p_value(jj,kk) = gammainc(results.chi_squared(jj,kk)/2, N/2, 'upper');
    end
end

results.pass_ratio = numel(find(results.p_value > 0.01)) / numel(results.p_value);

end