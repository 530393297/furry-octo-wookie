function [ poly_mat ] = generate_poly( p, N, h, k, X )

poly_mat = zeros(h * k);
n = N;
x = [1, 0];

for i = 1:h * k
    v = floor((i - 1) / k);
    u = (i - 1) - (k * v);
    q = conv(conv(poly_pow(n, h - 1 - v), poly_pow(x, u)), poly_pow(p, v));
    q = padarray(q, [0 ((h * k) - length(q))], 'pre');
    
    for j = 1:h * k
        poly_mat(i, j) = q(h * k - j + 1) * X ^ (j - 1);
    end
end

end

