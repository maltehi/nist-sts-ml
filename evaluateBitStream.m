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

p_c = 1-opt.alpha;
stats.metric.alpha = opt.alpha;
stats.metric.nStreams = size(bitStream,2);
stats.metric.confidenceInterval = min([p_c - 3*sqrt(p_c*opt.alpha/stats.metric.nStreams), ...
                                       p_c + 3*sqrt(p_c*opt.alpha/stats.metric.nStreams)], ...
                                      1);
                           
fprintf('alpha: %d\n', stats.metric.alpha)
fprintf('nStreams: %d\n', stats.metric.nStreams)
fprintf('confidence interval: [%d, %d]\n', stats.metric.confidenceInterval(1),...
        stats.metric.confidenceInterval(2))

% Frequency Test
if opt.all || (isfield(opt,'freq') && opt.freq.active)
    fprintf('\n')
    disp('====================')
    disp('   Frequency test')
    disp('====================')
    fprintf('\n')
    stats.frequency = frequencyTest(bitStream, opt.n);
    stats.frequency.final = assess(stats.frequency, stats.metric);
    
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
    disp('==========================')
    disp('   Block Frequency test')
    disp('==========================')
    fprintf('\n')
    if ~isfield(opt,'M')
        error('Block frequency test not possible: M is not defined');
    end
    stats.blockFrequency = blockFrequencyTest(bitStream, opt.M, opt.n);    
    stats.blockFrequency.final = assess(stats.blockFrequency, stats.metric);
end

% Cumulative Sums Test
if opt.all || (isfield(opt,'cumSums') && opt.cumSums.active)
    fprintf('\n')
    disp('==========================')
    disp('   Cumulative Sums test')
    disp('==========================')
    fprintf('\n')
    stats.cumulativeSums = cumulativeSumsTest(bitStream, opt.n);
    if isfield(stats.cumulativeSums,'forward')
        fprintf('Forward: ')
        stats.cumulativeSums.forward.final = assess(stats.cumulativeSums.forward, stats.metric);
    end
    if isfield(stats.cumulativeSums,'reverse')
        fprintf('Reverse: ')
        stats.cumulativeSums.reverse.final = assess(stats.cumulativeSums.reverse, stats.metric);
    end
end

% Runs Test
if opt.all || (isfield(opt,'runs') && opt.runs.active)
    fprintf('\n')
    disp('===============')
    disp('   Runs test')
    disp('===============')
    fprintf('\n')
    stats.runs = runsTest(bitStream, opt.n);
    stats.runs.final = assess(stats.runs, stats.metric);
end

% Longest Run of Ones Test
if opt.all || (isfield(opt,'longRuns') && opt.longRuns.active)    
    fprintf('\n')
    disp('==============================')
    disp('   Longest Run of Ones test')
    disp('==============================')
    fprintf('\n')
    stats.longestRunOfOnes = longestRunOfOnesTest(bitStream, opt.n);
    stats.longestRunOfOnes.final = assess(stats.longestRunOfOnes, stats.metric);
end

% Rank Test
if opt.all || (isfield(opt,'rank') && opt.rank.active)
    fprintf('\n')
    disp('===============')
    disp('   Rank test')
    disp('===============')
    fprintf('\n')
    stats.rank = rankTest(bitStream, opt.n);
    stats.rank.final = assess(stats.rank, stats.metric);
end

% DFT Test
if opt.all || (isfield(opt,'dft') && opt.dft.active)
    fprintf('\n')
    disp('=====================================')
    disp('   Discrete Fourier Transform test')
    disp('=====================================')
    fprintf('\n')
    stats.dft = dftTest(bitStream, opt.n);
    stats.dft.final = assess(stats.dft, stats.metric);
end

% Non Overlapping Templates Test
if opt.all || (isfield(opt,'nonOverlap') && opt.nonOverlap.active)
    fprintf('\n')
    disp('====================================')
    disp('   Non Overlapping Templates test   ')
    disp('====================================')
    fprintf('\n')
    stats.nonOverlap = nonOverlappingTest(bitStream, opt.m, opt.n, opt.N);
    stats.nonOverlap.final = assess(stats.nonOverlap, stats.metric);
end

% % Overlapping Templates Test
% if opt.all || (isfield(opt,'overlap') && opt.overlap.active)
%     stats.overlap = overlappingTest(bitStream, opt.n);
%     stats.overlap.final = assess(stats.overlap.pass_ratio, ...
%                                      stats.confidenceInterval);
% end

    function outStats = assess(testStats, metricStats)
        testStats.pass_ratio = numel(find(testStats.p_value >= metricStats.alpha)) ...
                               / numel(testStats.p_value); 
        outStats.proportion = testStats.pass_ratio >= metricStats.confidenceInterval(1) ...
                                     && testStats.pass_ratio <= metricStats.confidenceInterval(2);
        
        fprintf('Proportion of passed sequences: ')
        if outStats.proportion
            fprintf('PASSED\n');
        else
            fprintf('NOT PASSED\n');
        end 
        fprintf('pass_ratio: %d\n\n', testStats.pass_ratio)
        
%         outStats.F = discretize(testStats.p_value(:).', 0:0.1:1);
%         outStats.chi_squared = sum((outStats.F - metricStats.nStreams/10).^2 ...
%                                           / (metricStats.nStreams/10));
%         outStats.p_value_T = gammainc(outStats.chi_squared/2, 9/2, 'upper');        
%         outStats.fitting = outStats.p_value_T >= 0.0001;
%         
%         fprintf('Uniform distribution of P-values: ')
%         if outStats.fitting
%             fprintf('PASSED\n');
%         else
%             fprintf('NOT PASSED\n');
%         end 
%         fprintf('p_value_T: %d\n\n', outStats.p_value_T)
    end

end

