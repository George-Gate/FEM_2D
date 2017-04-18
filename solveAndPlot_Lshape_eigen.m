% Finite Element Method Solver

%% parameters
% shishkin type mesh
% nPerAxis=2^5;
% n={[nPerAxis/4;nPerAxis/4];[nPerAxis/4;nPerAxis/4];[nPerAxis/4;nPerAxis/4];[nPerAxis/4;nPerAxis/4]};
% w={[0.9];[0.1];[0.9];[0.1]};
% uniform mesh
nPerAxis=2^6;
n={[nPerAxis/2];[nPerAxis/2];[nPerAxis/2];[nPerAxis/2]};
w={[];[];[];[]};

meshType='LshapeSegUniform';
basis='BiLinear';


%% numerical solution    
% the following depends on dFmt, f(x) and n
% get the coefficient matrices S, C, M and vecf
%meshWidth=min(0.49,epsilon/b*2.5*log(n));
mesh0=makeMesh(meshType,n,w);
[S,~,~,M,~,id2fun,fun2id]=getCoeffs2D(mesh0,basis);


% solve
opts.issym=1;
opts.p=400;
[u,D]=eigs(S,M,20,'sa',opts);
% job=batch('tic;[u,D]=eigs(S,M,20,''sa'',opts);toc;');
% wait(job);
% load(job);
% delete(job);
% 
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

% plot
figure();
nx=200;ny=200;
[x,y]=meshgrid( linspace(min(xList),max(xList),nx) , linspace(min(yList),max(yList),ny) );
% contour(x,y,imag(numSol(x,y)));hold on;
% normalize maximum to 1
valMax=max(max(numSol(x,y)));
valMin=min(min(numSol(x,y)));
if abs(valMax)<abs(valMin)
    valMax=valMin;
end
surf(x,y,(numSol(x,y)/valMax),numSol(x,y),'lineStyle','none');hold on;
colormap(gca,'jet');
% caxis([0 1]);
W=1;set(gca,'xlim',[-W,W],'ylim',[-W,W]);

% refine plot
title([num2str(eigenState),'th eigen state.  Eigen Value=',num2str(D(eigenState,eigenState))]);
xlabel('$$x$$','interpreter','latex');
ylabel('$$y$$','interpreter','latex');
% set(ax(1),'fontsize',12,'ylim',[0,1],'Ycolor','black');
% set(ax(2),'fontsize',12,'Ycolor','black');