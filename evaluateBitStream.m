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
if nargin < 2
    opt = struct;
    opt.all = true;
end
if ~isfield(opt,'all')
    opt.all = false;
end
if ~isfield(opt,'alpha')
    opt.alpha = 0.01;
end
if ~isfield(opt,'n')
    opt.n = size(bitStream,1);
end
if ~isfield(opt,'m')
    opt.m = round(log2(opt.n)/2);
end
if ~isfield(opt,'M')
    opt.M = round(sqrt(opt.n));
end  
if ~isfield(opt,'N')
    opt.N = 8;
end
nStreams = size(bitStream,2);

p_c = 1-opt.alpha;
stats.confidenceInterval = [p_c - 3*sqrt(p_c*opt.alpha/nStreams), ...
                            p_c + 3*sqrt(p_c*opt.alpha/nStreams)];

% Frequency Test
if opt.all || (isfield(opt,'freq') && opt.freq.active)
    disp('====================')
    disp('   Frequency test')
    disp('====================')
    fprintf('\n')
    stats.frequency = frequencyTest(bitStream, opt.n);
    stats.frequency.totalPass = assess(stats.frequency.pass_ratio, ...
                                       stats.confidenceInterval);
    
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
if opt.all || (isfield(opt,'blockFreq') && opt.blockFreq.active)
    fprintf('\n')
    fprintf('\n')
    disp('==========================')
    disp('   Block Frequency test')
    disp('==========================')
    fprintf('\n')
    if ~isfield(opt,'M')
        error('Block frequency test not possible: M is not defined');
    end
    stats.blockFrequency = blockFrequencyTest(bitStream, opt.M, opt.n);    
    stats.blockFrequency.totalPass = assess(stats.blockFrequency.pass_ratio, ...
                                            stats.confidenceInterval);
end

% Cumulative Sums Test
if opt.all || (isfield(opt,'cumSums') && opt.cumSums.active)
    fprintf('\n')
    fprintf('\n')
    disp('==========================')
    disp('   Cumulative Sums test')
    disp('==========================')
    fprintf('\n')
    stats.cumulativeSums = cumulativeSumsTest(bitStream, opt.n);
    if isfield(stats.cumulativeSums,'forward')
        fprintf('Forward: ')
        stats.cumulativeSums.forward.totalPass = assess(stats.cumulativeSums.forward.pass_ratio, ...
                                                        stats.confidenceInterval);
    end
    if isfield(stats.cumulativeSums,'reverse')
        fprintf('Reverse: ')
        stats.cumulativeSums.reverse.totalPass = assess(stats.cumulativeSums.reverse.pass_ratio, ...
                                                        stats.confidenceInterval);
    end
end

% Runs Test
if opt.all || (isfield(opt,'runs') && opt.runs.active)
    fprintf('\n')
    fprintf('\n')
    disp('===============')
    disp('   Runs test')
    disp('===============')
    fprintf('\n')
    stats.runs = runsTest(bitStream, opt.n);
    stats.runs.totalPass = assess(stats.runs.pass_ratio, ...
                                  stats.confidenceInterval);
end

% Longest Run of Ones Test
if opt.all || (isfield(opt,'longRuns') && opt.longRuns.active)    
    fprintf('\n')
    fprintf('\n')
    disp('==============================')
    disp('   Longest Run of Ones test')
    disp('==============================')
    fprintf('\n')
    stats.longestRunOfOnes = longestRunOfOnesTest(bitStream, opt.n);
    stats.longestRunOfOnes.totalPass = assess(stats.longestRunOfOnes.pass_ratio, ...
                                              stats.confidenceInterval);
end

% Rank Test
if opt.all || (isfield(opt,'rank') && opt.rank.active)
    fprintf('\n')
    fprintf('\n')
    disp('===============')
    disp('   Rank test')
    disp('===============')
    fprintf('\n')
    stats.rank = rankTest(bitStream, opt.n);
    stats.rank.totalPass = assess(stats.rank.pass_ratio, ...
                                  stats.confidenceInterval);                              
end

% DFT Test
if opt.all || (isfield(opt,'dft') && opt.dft.active)
    fprintf('\n')
    fprintf('\n')
    disp('=====================================')
    disp('   Discrete Fourier Transform test')
    disp('=====================================')
    fprintf('\n')
    stats.dft = dftTest(bitStream, opt.n);
    stats.dft.totalPass = assess(stats.dft.pass_ratio, ...
                                 stats.confidenceInterval);
end

% Non Overlapping Templates Test
if opt.all || (isfield(opt,'nonOverlap') && opt.nonOverlap.active)
    fprintf('\n')
    fprintf('\n')
    disp('====================================')
    disp('   Non Overlapping Templates test   ')
    disp('====================================')
    fprintf('\n')
    stats.nonOverlap = nonOverlappingTest(bitStream, opt.m, opt.n, opt.N);
    stats.nonOverlap.totalPass = assess(stats.nonOverlap.pass_ratio, ...
                                        stats.confidenceInterval);
end

% % Overlapping Templates Test
% if opt.all || (isfield(opt,'overlap') && opt.overlap.active)
%     stats.overlap = overlappingTest(bitStream, opt.n);
%     stats.overlap.totalPass = assess(stats.overlap.pass_ratio, ...
%                                      stats.confidenceInterval);
% end

    function pass = assess(pass_ratio, confInt)
        pass = pass_ratio > confInt(1) && pass_ratio < confInt(2);
        if pass
            fprintf('Passed! pass_ratio: %d\n',pass_ratio)
        else
            fprintf('Not passed! pass_ratio: %d\n',pass_ratio)
        end
    end

end

