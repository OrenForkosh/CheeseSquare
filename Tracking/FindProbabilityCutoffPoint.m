function c = FindProbabilityCutoffPoint(input, alpha)
nsteps = 2;
c = max(input)+1;
for i=1:nsteps
    [h,x] = hist(input(input<c), 1000);
    h = h/sum(h);
    k=find(cumsum(h) > alpha, 1);
    c = x(k);
end
