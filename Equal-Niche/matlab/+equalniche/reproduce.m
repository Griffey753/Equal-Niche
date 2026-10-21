function [children,parent1,parent2] = reproduce(X,nicheId,lower,upper,count,parent1Mode,strength)
%REPRODUCE Parent 1 sampling, global-uniform Parent 2, SBX and mutation.
% Equal-niche: P(member j in niche i)=1/(K*|N_i|).
% Individual-uniform control: P(member j)=1/n. No other operator changes.
    if nargin < 6
        parent1Mode = 'EQUAL_NICHE';
    end
    if nargin < 7
        strength = 1;
    end
    X = double(X);
    [n,d] = size(X);
    K = max(nicheId);
    sizes = accumarray(nicheId(:),1,[K,1]);
    switch char(parent1Mode)
        case 'EQUAL_NICHE'
            probability = 1./(K*sizes(nicheId(:)));
        case 'INDIVIDUAL_UNIFORM'
            probability = ones(n,1)/n;
        case 'INTERPOLATED'
            assert(isscalar(strength) && strength >= 0 && strength <= 1, ...
                'EqualNiche:Strength','Interpolation strength must be in [0,1].');
            probability = (1-strength)/n + strength./(K*sizes(nicheId(:)));
        otherwise
            error('EqualNiche:Parent1Mode','Unknown Parent 1 mode: %s',char(parent1Mode));
    end
    assert(abs(sum(probability)-1) < 1e-12, ...
        'EqualNiche:Probability','Parent 1 probabilities do not sum to one.');
    cumulative = cumsum(probability); cumulative(end) = 1;
    draws = rand(count,1);
    parent1 = 1+sum(draws >= cumulative',2);
    parent2 = randi(n,count,1);
    p1 = X(parent1,:); p2 = X(parent2,:);
    children = p1;
    crossoverProbability = 1;
    crossoverEta = 20;
    mutationProbability = 1;
    mutationEta = 20;
    for i = 1:count
        if rand <= crossoverProbability
            mu = rand(1,d);
            beta = zeros(1,d);
            left = mu <= 0.5;
            beta(left) = (2*mu(left)).^(1/(crossoverEta+1));
            beta(~left) = (2-2*mu(~left)).^(-1/(crossoverEta+1));
            signs = randi([0,1],1,d);
            beta(signs == 1) = -beta(signs == 1);
            keep = rand(1,d);
            beta(keep < 0.5) = 1;
            if rand < 0.5
                children(i,:) = 0.5*((1+beta).*p1(i,:)+(1-beta).*p2(i,:));
            else
                children(i,:) = 0.5*((1-beta).*p1(i,:)+(1+beta).*p2(i,:));
            end
        end
    end
    lower = reshape(double(lower),1,[]);
    upper = reshape(double(upper),1,[]);
    children = min(max(children,lower),upper);
    site = rand(count,d) < mutationProbability/d;
    mu = rand(count,d);
    range = upper-lower;
    for j = 1:d
        low = site(:,j) & mu(:,j) <= 0.5;
        if any(low)
            delta = (children(low,j)-lower(j))/range(j);
            term = (2*mu(low,j)+(1-2*mu(low,j)).* ...
                (1-delta).^(mutationEta+1)).^(1/(mutationEta+1))-1;
            children(low,j) = children(low,j)+range(j)*term;
        end
        high = site(:,j) & mu(:,j) > 0.5;
        if any(high)
            delta = (upper(j)-children(high,j))/range(j);
            term = 1-(2*(1-mu(high,j))+2*(mu(high,j)-0.5).* ...
                (1-delta).^(mutationEta+1)).^(1/(mutationEta+1));
            children(high,j) = children(high,j)+range(j)*term;
        end
    end
    children = min(max(children,lower),upper);
end
