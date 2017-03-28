% Finite Element Method Solver

%% parameters
% shishkin type mesh
% nPerAxis=2^8;
n={[nPerAxis/4;nPerAxis/4+1];[nPerAxis/4;nPerAxis/4];[nPerAxis/4;nPerAxis/4];[nPerAxis/4+1;nPerAxis/4]};
w={[0.90];[0.1];[0.9];[0.1]};
% uniform mesh
% nPerAxis=2^8;
% n={[nPerAxis/2];[nPerAxis/2];[nPerAxis/2];[nPerAxis/2]};
% w={[];[];[];[]};

meshType='LshapeSegUniform';
basis='Linear';


%% numerical solution    
% the following depends on dFmt, f(x) and n
% get the coefficient matrices S, C, M and vecf
%meshWidth=min(0.49,epsilon/b*2.5*log(n));
tic;
mesh0=makeMesh(meshType,n,w);
disp(['Time to makeMesh: ',num2str(toc)]);tic;
[S,~,~,~,vecf,id2fun,fun2id]=getCoeffs2D(mesh0,basis,1);
disp(['Time to getCoeffs: ',num2str(toc)]);tic;

% solve
u=S\vecf;
disp(['Time to solve linear system: ',num2str(toc)]);tic;
% save;

%% interpolant
eigenState=1;
% get node coordinates
N=mesh0.Nnodes;
xList=mesh0.nodes.x;yList=mesh0.nodes.y;

% prepare solution data
Ninner=length(u);
tmp_u=[u(:,eigenState);zeros(N-Ninner,1)];
tmp_x=[xList(id2fun(1:Ninner));xList(fun2id==0)];
tmp_y=[yList(id2fun(1:Ninner));yList(fun2id==0)];

% interpolant
numSol=scatteredInterpolant(tmp_x,tmp_y,tmp_u);

disp(['Time to interpolant: ',num2str(toc)]);
