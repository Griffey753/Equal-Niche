function [nicheId,centers,sizes] = temporaryNiches(X,F,rank,lower,upper,radius,centerOrder)
%TEMPORARYNICHES Rank/coordinate ordered, greedy exclusive radius partition.
    if nargin < 7
        centerOrder = 'rank_lexicographic';
    end
    n = size(X,1);
    Z = (double(X)-reshape(double(lower),1,[]))./ ...
        (reshape(double(upper),1,[])-reshape(double(lower),1,[]));
    key = [double(rank(:)),Z,double(F),(1:n)'];
    [~,order] = sortrows(key,1:size(key,2));
    distance = pdist2(Z,Z,'euclidean');
    if strcmp(centerOrder,'rank_stratified_farthest_first')
        canonical = order;
        order = zeros(n,1);
        next = 0;
        ranks = unique(rank(canonical));
        for r = reshape(ranks,1,[])
            remaining = canonical(rank(canonical) == r);
            while ~isempty(remaining)
                if next == 0
                    position = 1;
                else
                    nearest = min(distance(remaining,order(1:next)),[],2);
                    [~,position] = max(nearest);
                end
                next = next+1;
                order(next) = remaining(position);
                remaining(position) = [];
            end
        end
    else
        assert(strcmp(centerOrder,'rank_lexicographic'), ...
            'EqualNiche:CenterOrder','Unknown center order.');
    end
    nicheId = zeros(n,1);
    centers = zeros(0,1);
    sizes = zeros(0,1);
    for p = 1:n
        center = order(p);
        if nicheId(center) ~= 0
            continue
        end
        members = find(nicheId == 0 & distance(center,:)' <= radius);
        localId = numel(centers)+1;
        nicheId(members) = localId;
        centers(end+1,1) = center; %#ok<AGROW>
        sizes(end+1,1) = numel(members); %#ok<AGROW>
    end
    assert(all(nicheId > 0) && sum(sizes) == n && all(sizes > 0), ...
        'EqualNiche:NichePartition','Exclusive partition audit failed.');
end
