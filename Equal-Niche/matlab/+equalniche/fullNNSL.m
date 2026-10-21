function [rank,diagnostic] = fullNNSL(X,F,lower,upper,radius,alpha)
%FULLNNSL Fixed-radius, local-dominance front sorting and alpha update.
% This is the shared sorter for the population and the exact-deduplicated
% archive. It follows src/equal_niche/full_nnsl_sorting.py.
    X = double(X); F = double(F);
    n = size(X,1);
    assert(n == size(F,1) && all(isfinite(F(:))), ...
        'EqualNiche:Objectives','Decision/objective matrices are invalid.');
    distance = equalniche.normalizedDistance(X,lower,upper);
    neighbor = distance <= radius;
    neighbor(1:n+1:end) = false;
    dominates = false(n,n);
    for j = 1:size(F,2)
        component = F(:,j) <= F(:,j)';
        if j == 1
            weak = component;
            strict = F(:,j) < F(:,j)';
        else
            weak = weak & component;
            strict = strict | (F(:,j) < F(:,j)');
        end
    end
    dominates = weak & strict & neighbor;
    rank = zeros(n,1);
    remaining = sum(dominates,1)';
    front = 0;
    while any(rank == 0)
        current = find(rank == 0 & remaining == 0);
        assert(~isempty(current),'EqualNiche:SortCycle', ...
            'Local dominance graph has no available front.');
        front = front+1;
        rank(current) = front;
        remaining = remaining-sum(dominates(current,:),1)';
    end
    first = find(rank == 1);
    sums = sum(F(first,:),2);
    assert(all(sums >= 0),'EqualNiche:Alpha','Alpha update requires nonnegative objective sums.');
    denominator = max(min(sums),eps);
    demoted = first(sums/denominator > alpha);
    rank(demoted) = 2;
    diagnostic = struct('candidateCount',n,'dominanceEdges',nnz(dominates), ...
        'rawFrontCount',front,'adjustedFrontCount',max(rank), ...
        'firstFrontBeforeAlpha',numel(first), ...
        'firstFrontAfterAlpha',sum(rank == 1), ...
        'alphaDemotedCount',numel(demoted), ...
        'meanNeighborCount',mean(sum(neighbor,2)), ...
        'isolatedRatio',mean(sum(neighbor,2) == 0));
end
