function [finalX,finalF,archiveX,archiveF] = run_equal_niche_example()
%RUN_EQUAL_NICHE_EXAMPLE Small standalone PlatEMO MMF1 demonstration.
% Add PlatEMO v4.9 (including subdirectories) and this folder to the path.
% This example does not reproduce the paper's shared benchmark protocol.
    assert(exist('ALGORITHM','class') == 8 && exist('MMF1','class') == 8, ...
        'EqualNiche:Dependency','Add the PlatEMO v4.9 folders to the MATLAB path.');
    rng(1,'twister');
    Problem = MMF1('N',40,'maxFE',400);
    Algorithm = EqualNicheAdapter('save',0,'outputFcn',@silentOutput);
    Algorithm.Solve(Problem);
    finalX = Algorithm.finalPopulation.decs;
    finalF = Algorithm.finalPopulation.objs;
    archiveX = Algorithm.nativeFinalSet.decs;
    archiveF = Algorithm.nativeFinalSet.objs;
    assert(Problem.FE == 400 && size(finalX,1) == 40);
    fprintf('Equal-Niche example PASS: %d evaluations, %d final individuals.\n', ...
        Problem.FE,size(finalX,1));
end

function silentOutput(~,~)
end
