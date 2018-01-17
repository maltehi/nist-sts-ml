function [ stats ] = evaluateBitStream( bitStream, opt )
% Statistical testing of a pseudorandom bit stream according to
% NIST Statistical Test Suite
% https://csrc.nist.gov/Projects/Random-Bit-Generation/Documentation-and-Software

stats = struct;

if opt.frequencyTest.active
    stats.frequency = frequencyTest(bitStream, opt.n);
    
    % If frequency test fails, other tests are likely to fail, too
    failedIdx = find(stats.frequency.p_value < 0.01);
    if ~isempty(failedIdx)
        s = input(sprintf('Frequency test resulted in non-randomness of %d sequence(s). Continue? (y/N) ', length(failedIdx)),'s');
        if ~strcmpi(s,'y')
            return
        end
    end
    [p_val_min, p_idx_min] = min(stats.frequency.p_value);
    fprintf('Frequency test: min(p_value) = %d, sequence %d\n', p_val_min, p_idx_min);
end
if opt.blockFrequencyTest.active
    stats.blockFrequency = blockFrequencyTest(bitStream, opt.M, opt.n);
    [p_val_min, p_idx_min] = min(stats.blockFrequency.p_value);
    fprintf('Block frequency test: min(p_value) = %d, sequence %d\n', p_val_min, p_idx_min);
end

end

