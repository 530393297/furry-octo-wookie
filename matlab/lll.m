function [ b ] = lll( b, delta )

B = gram_schmidt(b);

for i = 2:length(b)
    
    for j = i- 1:-1:1
        c = round(dot(b(i,:), B(j,:)) / dot(B(j,:), B(j,:)));
        b(i,:) = b(i,:) - (c * b(j,:));
    end
end
for i = 1:length(b) - 1
    l = ((dot(b(i + 1, :), B(i, :)) / dot(B(i,:), B(i,:))) * B(i,:)) ...
        + B(i + 1, :);
    if delta * dot(B(i,:), B(i,:)) > dot(l, l)
        b([i i + 1], :) = b([i + 1 i], :);
        b = lll(b, delta);
        break;
    end
end

end

