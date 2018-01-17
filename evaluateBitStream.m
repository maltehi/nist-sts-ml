% Statistical testing of a pseudorandom bit streams according to
% NIST Statistical Test Suite
% https://csrc.nist.gov/Projects/Random-Bit-Generation/Documentation-and-Software
%
% Inputs:
%   bitStream:  n x m Matrix. Bit stream(s) under evaluation in columns.
%   opt:        Parameter structure
%
% Output:
%   stats:      structure containing the results returned by every test in
%               sub-structures
function [stats] = evaluateBitStream(bitStream, opt)

stats = struct;
if ~isfield(opt,'n')
    opt.n = size(bitStream,1);
end

% Frequency Test
if isfield(opt,'freq') && opt.freq.active
    stats.frequency = frequencyTest(bitStream, opt.n);
    
    % If frequency test fails, other tests are likely to fail, too
    failedIdx = find(stats.frequency.p_value < 0.01);
    if ~isempty(failedIdx)
        s = input(sprintf(['Frequency test resulted in non-randomness of ' ...
                           '%d sequence(s). Continue? (y/N) '], length(failedIdx)),'s');
        if ~strcmpi(s,'y')
            return
        end
    end
end

% Block Frequency Test
if isfield(opt,'blockFreq') && opt.blockFreq.active
    if ~isfield(opt,'M')
        error('Block frequency test not possible: M is not defined');
    end
    stats.blockFrequency = blockFrequencyTest(bitStream, opt.M, opt.n);
end

% Cumulative Sums Test
if isfield(opt,'cumSums') && opt.cumSums.active
    stats.cumulativeSums = cumulativeSumsTest(bitStream, opt.n);
end

% Runs Test
if isfield(opt,'runs') && opt.runs.active
    stats.runs = runsTest(bitStream, opt.n);
end

end

