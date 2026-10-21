function distance = normalizedDistance(X,lower,upper)
%NORMALIZEDDISTANCE Pairwise Euclidean distance after bound normalization.
    lower = reshape(double(lower),1,[]);
    upper = reshape(double(upper),1,[]);
    X = double(X);
    assert(size(X,2) == numel(lower) && all(upper > lower), ...
        'EqualNiche:Bounds','Invalid decision matrix or bounds.');
    Z = (X-lower)./(upper-lower);
    distance = pdist2(Z,Z,'euclidean');
end
