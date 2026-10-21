function [chosen,criticalFront,violation] = selectByFronts(X,rank,target,lower,upper,k)
%SELECTBYFRONTS Preserve complete better fronts and truncate only critical.
    n = size(X,1);
    assert(target >= 0 && target <= n,'EqualNiche:SelectionSize','Invalid target size.');
    chosen = zeros(0,1);
    criticalFront = NaN;
    for front = unique(rank(:))'
        members = find(rank == front);
        if numel(chosen)+numel(members) <= target
            chosen = [chosen;members]; %#ok<AGROW>
            if numel(chosen) == target
                break
            end
        else
            criticalFront = front;
            chosen = truncateCritical(X,chosen,members,target-numel(chosen),lower,upper,k);
            break
        end
    end
    assert(numel(chosen) == target,'EqualNiche:SelectionSize','Front fill failed.');
    picked = false(n,1); picked(chosen) = true;
    violation = any(picked) && any(~picked) && ...
        any(rank(~picked) < max(rank(picked)));
    assert(~violation,'EqualNiche:FrontPriority','Better front was displaced.');
end

function chosen = truncateCritical(X,protected,critical,slots,lower,upper,k)
    Z = (double(X)-reshape(double(lower),1,[]))./ ...
        (reshape(double(upper),1,[])-reshape(double(lower),1,[]));
    activeOrder = [protected;critical];
    n = size(X,1);
    active = false(n,1); active(activeOrder) = true;
    removable = false(n,1); removable(critical) = true;
    distance = pdist2(Z,Z,'euclidean');
    distance(1:n+1:end) = inf;
    [~,neighborOrder] = sort(distance,2,'ascend');
    effectiveK = min(max(k,1),max(numel(activeOrder)-1,1));
    nearest = zeros(n,effectiveK);
    crowding = inf(n,1);
    refresh(activeOrder);
    remaining = numel(critical);
    while remaining > slots
        candidates = activeOrder(active(activeOrder) & removable(activeOrder));
        minimum = min(crowding(candidates));
        tied = candidates(abs(crowding(candidates)-minimum) <= ...
            1e-8+1e-5*abs(minimum));
        tied = tied(:);
        [~,order] = sortrows([Z(tied,:),-tied]);
        removed = tied(order(end));
        active(removed) = false;
        remaining = remaining-1;
        affected = find(active & any(nearest == removed,2));
        refresh(affected);
    end
    chosen = activeOrder(active(activeOrder));

    function refresh(indices)
        for ii = reshape(indices,1,[])
            ordered = neighborOrder(ii,:);
            neighbors = ordered(active(ordered));
            neighbors = neighbors(1:min(effectiveK,numel(neighbors)));
            nearest(ii,:) = 0;
            nearest(ii,1:numel(neighbors)) = neighbors;
            crowding(ii) = sum(distance(ii,neighbors));
        end
    end
end
