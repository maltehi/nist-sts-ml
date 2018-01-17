function [ results ] = cumulativeSumsTest( bitStream, n, direction )

if nargin < 3
    direction = {'forward', 'reverse'};
end
if nargin < 2
    n = size(bitStream,1);
end

% Truncate and make +/-1
bitStream = bitStream(1:n,:);
if isempty(find(bitStream < 0, 1))
    bitStream = 2*bitStream-1;
end

% Iterate directions
for i = 1:length(direction)
    S = cumsum(bitStream,1,direction{i});
    results.(direction{i}).z = max(abs(S),[],1);
    results.(direction{i}).p_value = sums(results.(direction{i}).z, n);
end

    % Calculate p_value
    function p_val = sums(z, nSeq)
        sum1 = 0;
        for k = fix((fix(-nSeq./z)+1)/4):fix((fix(nSeq./z)-1)/4)
             sum1 = sum1 + ndtr(((4*k+1)*z)/sqrt(nSeq));
             sum1 = sum1 - ndtr(((4*k-1)*z)/sqrt(nSeq));
        end
        sum2 = 0;
        for k = fix((fix(-nSeq./z)-3)/4):fix((fix(nSeq./z)-1)/4)
             sum2 = sum2 + ndtr(((4*k+3)*z)/sqrt(nSeq));
             sum2 = sum2 - ndtr(((4*k+1)*z)/sqrt(nSeq));
        end
        p_val = 1.0 - sum1 + sum2;
    end    
end

