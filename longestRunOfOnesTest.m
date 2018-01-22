function [results] = longestRunOfOnesTest(bitStream, n)

if nargin < 2
    n = size(bitStream,1);
end

if n < 128
    results.p_value = zeros(1,n);
    fprintf("\t\t\t  LONGEST RUNS OF ONES TEST\n");
    fprintf("\t\t---------------------------------------------\n");
    fprintf("\t\t   n=%d is too short\n", n);
    return
elseif n < 6272
%     V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
%     pi[0] = 0.21484375;
%     pi[1] = 0.3671875;
%     pi[2] = 0.23046875;
%     pi[3] = 0.1875;
    K = 3;
    M = 8;
    V = 1:4;
    pi_i = [0.21484375, 0.3671875, 0.23046875, 0.1875];
elseif n < 750000
%     V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
%     pi[0] = 0.1174035788;
%     pi[1] = 0.242955959;
%     pi[2] = 0.249363483;
%     pi[3] = 0.17517706;
%     pi[4] = 0.102701071;
%     pi[5] = 0.112398847;
    K = 5;
    M = 128;
    V = 4:9;
    pi_i = [0.1174035788, 0.242955959, 0.249363483, 0.17517706, 0.102701071, 0.112398847];
else 
%     V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
%     pi[0] = 0.0882;
%     pi[1] = 0.2092;
%     pi[2] = 0.2483;
%     pi[3] = 0.1933;
%     pi[4] = 0.1208;
%     pi[5] = 0.0675;
%     pi[6] = 0.0727;
    K = 6;
    M = 10000;
    V = 10:16;
    pi_i = [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727];
end

N = fix(n/M);

% Truncate and make 0/1
bitStream = bitStream(1:M*N,:);
if ~isempty(find(bitStream < 0, 1))
    bitStream = (bitStream > 0);
end

% Blocks in columns
bitStream = reshape(bitStream, M, N,[]);

% Find longest run per block
for k = 1:size(bitStream,3)    
    nu = zeros(1,length(V));
    for i = 1:N
        run = 0; v_n_obs = 0;
        for j = 1:M
            if bitStream(j, i, k) == 1
                run = run + 1;
                v_n_obs = max(v_n_obs, run);
            else
                run = 0;
            end
        end
        % Register run length in corresponding bin
        bin = discretize(v_n_obs,[-Inf, V(2:end), Inf]);
        nu(bin) = nu(bin) + 1;
    end

    % Compute chi_squared
    results.chi_squared(k) = sum((nu - N*pi_i).^2 ./ (N*pi_i));
end

% Compute p-value
results.p_value = gammainc(results.chi_squared/2, K/2, 'upper');

% Calculate pass ratio
results.pass_ratio = length(find(results.p_value > 0.01))/length(results.p_value);

end

