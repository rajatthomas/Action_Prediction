function [l,h] = indexfinder(array, band)

[~,I] = sort(abs(array-band(1)));
l = I(1);


[~,I] = sort(abs(array-band(2)));
h = I(1);

end