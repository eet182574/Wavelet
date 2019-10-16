function WaveletCoefficents = ComputeWaveletCoefficents(LPFCoefficents)
    
    
    N = length(LPFCoefficents);
    WaveletCoefficents = zeros(1,N);
    for i = 1:N
        WaveletCoefficents(i) = (-1)^(i-1) * LPFCoefficents(N-i+1);
    end
    
end