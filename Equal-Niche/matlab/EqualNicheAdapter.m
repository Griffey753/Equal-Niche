classdef EqualNicheAdapter < ALGORITHM
%EQUALNICHEADAPTER MATLAB/PlatEMO v4.9 port of the formal Equal-Niche solver.
% Population and Archive use the same fullNNSL/selectByFronts functions.
    properties
        finalPopulation
        nativeFinalSet
        generationAudit
        parent1Mode = 'EQUAL_NICHE'
        neighborhoodRadius = 1/3
        interpolationStrength = 1
        centerOrder = 'rank_lexicographic'
        captureSampling = false
        samplingAudit
    end
    methods
        function main(Algorithm,Problem)
            radius = Algorithm.neighborhoodRadius;
            alpha = 5;
            crowdingNeighbors = 2;
            Population = Problem.Initialization();
            Archive = Population;
            lower = Problem.lower;
            upper = Problem.upper;
            audit = zeros(0,5);
            sampling = cell(0,1);
            while true
                try
                    Algorithm.NotTerminated(Population);
                catch err
                    if strcmp(err.identifier,'PlatEMO:Termination')
                        Algorithm.finalPopulation = Population;
                        Algorithm.nativeFinalSet = Archive;
                        Algorithm.generationAudit = audit;
                        Algorithm.samplingAudit = sampling;
                    end
                    rethrow(err);
                end
                X = Population.decs;
                F = Population.objs;
                [rank,~] = equalniche.fullNNSL(X,F,lower,upper,radius,alpha);
                [nicheId,~,sizes] = equalniche.temporaryNiches( ...
                    X,F,rank,lower,upper,radius,Algorithm.centerOrder);
                count = min(Problem.N,floor(Problem.maxFE-Problem.FE));
                [childX,parent1,~] = equalniche.reproduce( ...
                    X,nicheId,lower,upper,count,Algorithm.parent1Mode, ...
                    Algorithm.interpolationStrength);
                Offspring = Problem.Evaluation(childX);
                unionPopulation = [Population,Offspring];
                unionX = unionPopulation.decs;
                unionF = unionPopulation.objs;
                [unionRank,~] = equalniche.fullNNSL( ...
                    unionX,unionF,lower,upper,radius,alpha);
                [selected,~,~] = equalniche.selectByFronts( ...
                    unionX,unionRank,Problem.N,lower,upper,crowdingNeighbors);
                Population = unionPopulation(selected);
                rawArchive = [Archive,unionPopulation];
                [~,first] = unique(rawArchive.decs,'rows','stable');
                candidates = rawArchive(first);
                [archiveRank,~] = equalniche.fullNNSL( ...
                    candidates.decs,candidates.objs,lower,upper,radius,alpha);
                [archiveSelected,~,~] = equalniche.selectByFronts( ...
                    candidates.decs,archiveRank,min(Problem.N,numel(candidates)), ...
                    lower,upper,crowdingNeighbors);
                Archive = candidates(archiveSelected);
                parentCounts = accumarray(nicheId(parent1),1,[numel(sizes),1]);
                assert(sum(parentCounts) == count && numel(Population) == Problem.N, ...
                    'EqualNiche:GenerationAudit','Population or Parent 1 audit failed.');
                audit(end+1,:) = [Problem.FE,count,numel(sizes), ...
                    min(sizes),max(sizes)]; %#ok<AGROW>
                if Algorithm.captureSampling
                    % Observation only: no draws or other RNG-state changes.
                    switch Algorithm.parent1Mode
                        case 'EQUAL_NICHE'
                            expected = ones(numel(sizes),1)*count/numel(sizes);
                        case 'INDIVIDUAL_UNIFORM'
                            expected = count*sizes/size(X,1);
                        otherwise
                            lambda = Algorithm.interpolationStrength;
                            expected = count*((1-lambda)*sizes/size(X,1) + ...
                                lambda/numel(sizes));
                    end
                    sampling{end+1,1} = [sizes,parentCounts,expected]; %#ok<AGROW>
                end
            end
        end
    end
end
