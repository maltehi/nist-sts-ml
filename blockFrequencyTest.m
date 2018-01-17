% bitStream: Bit stream(s) under evaluation. Column vector(s)
% M: Block length
% n: Length of bit stream to be evaluated (default: length of input stream)
function [ results ] = blockFrequencyTest( bitStream, M, n )

if nargin < 3
    n = size(bitStream,1);
end

N = floor(n/M); % Number of blocks
blocks = double(reshape(bitStream(1:M*N,:),M,N,[])); % Blocks in lines

results.sum = reshape(sum((sum(blocks,1)/M - 0.5).^2,2),1,[]);
results.chi_squared = 4 * M * results.sum;

% Definition of (cephes_)igamc: see
% https://github.com/jeremybarnes/cephes/blob/60f27df395b8322c2da22c83751a2366b82d50d1/cprob/igam.c
% and compare Matlab documentation gammainc
results.p_value = gammainc(results.chi_squared/2,N/2,'upper');

end

