%% Quantize the coefficents to b bits per coefficent.
function [ScaledCoefficents,a] = Quantize(Level_1Decomposition)

cmax = max(Level_1Decomposition(:));

min_value = min(Level_1Decomposition(:));

b = input('Enter Number of Bits: ');

if abs(min_value) > cmax
    
    a = abs(min_value) / (2^(b-1));
    
else
    
    a = abs(cmax) / (2^(b-1)-1);
    
end

ScaledCoefficents = floor(Level_1Decomposition / a);

end







