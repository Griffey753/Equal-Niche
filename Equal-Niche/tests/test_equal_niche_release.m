function test_equal_niche_release()
%TEST_EQUAL_NICHE_RELEASE Fixed-input tests without paper benchmark files.
    X = [0.05,0.05;0.1,0.05;0.9,0.9];
    F = [1,2;2,1;1.5,1.5];
    lower = [0,0]; upper = [1,1];
    [id,centers,sizes] = equalniche.temporaryNiches( ...
        X,F,ones(3,1),lower,upper,1/3);
    assert(isequal(id,[1;1;2]) && isequal(sizes,[2;1]));
    assert(isequal(centers,[1;3]));
    rng(75,'twister');
    expectedDraws = rand(256,1);
    expectedParent1 = 1+sum(expectedDraws >= [0.25,0.5,1],2);
    rng(75,'twister');
    [children,parent1,parent2] = equalniche.reproduce( ...
        X,id,lower,upper,256,'EQUAL_NICHE');
    assert(isequal(parent1,expectedParent1));
    assert(all(parent2 >= 1 & parent2 <= 3));
    assert(all(children(:) >= 0 & children(:) <= 1));
    rng(75,'twister');
    [childrenAgain,parentAgain] = equalniche.reproduce( ...
        X,id,lower,upper,256,'EQUAL_NICHE');
    assert(isequal(children,childrenAgain) && isequal(parent1,parentAgain));
    [rank,~] = equalniche.fullNNSL( ...
        [0,0;0.1,0.1;1,1],[1,1;2,2;0.5,0.5],lower,upper,1/3,5);
    assert(isequal(rank,[1;2;1]));
    [selected,~,violation] = equalniche.selectByFronts( ...
        X,[1;2;1],2,lower,upper,2);
    assert(isequal(selected,[1;3]) && ~violation);
    fprintf('Equal-Niche release tests PASS.\n');
end
