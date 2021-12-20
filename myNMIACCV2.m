function [res_mean,res_std,res1,res2,res3,res4]= myNMIACCV2(U,Y,numclass)

stream = RandStream.getGlobalStream;
reset(stream);
U_normalized = U ./ repmat(sqrt(sum(U.^2, 2)), 1,numclass);
maxIter = 50;
res1 = zeros(maxIter,1);
res2 = zeros(maxIter,1);
res3 = zeros(maxIter,1);
res4 = zeros(maxIter,1);
for it = 1 : maxIter
    indx = litekmeans(U_normalized,numclass, 'MaxIter',250, 'Replicates',1);
    indx = indx(:);
    [newIndx] = bestMap(Y,indx);
    res1(it) = mean(Y==newIndx);
    res2(it) = MutualInfo(Y,newIndx);
    res3(it) = purFuc(Y,newIndx);
    res4(it) = adjrandindex(Y,newIndx);
end
res_mean(1) = max(res1);
res_mean(2) = max(res2);
res_mean(3) = max(res3);
res_mean(4) = max(res4);
res_std(1) = std(res1);
res_std(2) = std(res2);
res_std(3) = std(res3);
res_std(4) = std(res4);