% Finite Element Method Solver

%% parameters
b=[1;-0.5];
c=1;
epsilon=1e-1;
f=@(x,y)1;

% shishkin type mesh
nPerAxis=2^7;
n={[nPerAxis/4;nPerAxis/4];[nPerAxis/4;nPerAxis/4];[nPerAxis/4;nPerAxis/4];[nPerAxis/4;nPerAxis/4]};
w={[0.9];[0.1];[0.9];[0.1]};
% uniform mesh
% nPerAxis=2^8;
% n={[nPerAxis/2];[nPerAxis/2];[nPerAxis/2];[nPerAxis/2]};
% w={[];[];[];[]};

meshType='LshapeSegUniform';
basis='BiLinear';


%% numerical solution    
% the following depends on dFmt, f(x) and n
% get the coefficient matrices S, C, M and vecf
%meshWidth=min(0.49,epsilon/b*2.5*log(n));
tic;
mesh0=makeMesh(meshType,n,w);
disp(['Time to makeMesh: ',num2str(toc)]);tic;
[S,Cx,Cy,M,vecf,id2fun,fun2id]=getCoeffs2D(mesh0,basis,f);
disp(['Time to getCoeffs: ',num2str(toc)]);tic;

% solve
H=epsilon*S+b(1)*Cx+b(2)*Cy+c*M;
u=H\vecf;
disp(['Time to solve linear system: ',num2str(toc)]);tic;
% save;

%% interpolant
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

%% plot
plotSol=@(x,y)numSol(x,y);
figure();
nx=200;ny=200;
[x,y]=meshgrid( linspace(-1,1,nx) , linspace(-1,1,ny) );
% contour(x,y,imag(numSol(x,y)));hold on;
% normalize maximum to 1
valMax=max(max(plotSol(x,y)));
valMin=min(min(plotSol(x,y)));
if abs(valMax)<abs(valMin)
    valMax=valMin;
end
surf(x,y,(plotSol(x,y)/valMax),plotSol(x,y),'lineStyle','none');hold on;
colormap(gca,'jet');
% caxis([0 1]);
W=1;set(gca,'xlim',[-W,W],'ylim',[-W,W]);

% refine plot
title(['nPerAxis=',num2str(nPerAxis)]);
xlabel('$$x$$','interpreter','latex');
ylabel('$$y$$','interpreter','latex');
% set(ax(1),'fontsize',12,'ylim',[0,1],'Ycolor','black');
% set(ax(2),'fontsize',12,'Ycolor','black');