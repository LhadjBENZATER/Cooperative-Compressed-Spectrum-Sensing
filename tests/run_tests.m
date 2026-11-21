function run_tests
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'examples'));
m=support_metrics([1 1 0 0],[1 0 1 0]);assert(max(abs(m-[0.5 0.5 0.5]))<1e-12);
m=support_metrics([1 1 0 0],[1 1 0 0]);assert(isequal(m,[1 0 1]));
m=support_metrics([1 1 0 0],[0 0 0 0]);assert(isequal(m,[0 0 0]));
A=eye(12);Y=zeros(12,4);Y([2 7],:)=1;V=(ones(4)-eye(4))/3;
assert(norm(Y*V-Y,'fro')<1e-12);
assert(norm(somp_baseline(A,Y,2)-Y,'fro')<1e-12);
rng(7);W=randn(10000,4);assert(max(abs(var(W*V)./var(W)-1/3))<0.03);
fprintf('All numerical checks passed.\n');
end
