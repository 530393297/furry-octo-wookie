function [ f2 ] = poly_pow(f1, x)

if x == 0
    f2 = 1;
else
    y = 1;
    while x > 1
        if ~mod(x,2)
            f1 = conv(f1, f1);
            x = x / 2;
        else
            y = conv(f1, y);
            f1 = conv(f1, f1);
            x = (x - 1) / 2;
        end
    end
    f2 = conv(y, f1);
end

end