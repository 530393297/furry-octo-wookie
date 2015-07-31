function [ r ] = copper_smith( f, N, h )

k = length(f) - 1;
X = ceil(2 ^ (-1/2) * (h * k) ^ (-1 / (h * k - 1))) * ...
    N ^ ((h - 1) / (h * k - 1)) - 1;
r = zeros(1, length(f));

l = lll(generate_poly(f, N, h, k, X), 0.75);

for i = 1:length(l)
    r(i) = l(1,i) / X ^ (i - 1);
end

r = roots(fliplr(r));

end

