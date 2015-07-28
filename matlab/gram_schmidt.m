function [ Q ] = gram_schmidt( A )

Q = A;
for i = 1:length(A)
    for j = 1:i - 1
        Q(i, :) = Q(i, :) - (dot(A(i,:), Q(j,:)) / dot(Q(j,:), Q(j,:))) ...
            * Q(j, :);
    end
end

end

