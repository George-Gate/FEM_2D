% Finite Element Method Solver

%% parameters
b=[1;-0.5];
c=1;
epsilon=1e-1;
f=@(x,y)x.^2+y+1;
fUb=[2,1];     % the order of f corresponding to x and y


% shishkin type mesh
nPerAxis=4;
n={[nPerAxis/4;nPerAxis/4];[nPerAxis/4;nPerAxis/4];[nPerAxis/4;nPerAxis/4];[nPerAxis/4;nPerAxis/4]};
w={[0.9];[0.1];[0.9];[0.1]};
% uniform mesh
% nPerAxis=2;
% n={nPerAxis/2;nPerAxis/2;nPerAxis/2;nPerAxis/2};
% w={[];[];[];[]};
% 
K=40;
meshType='LshapeSegUniform';
basis='Lobatto';


%% numerical solution    
% the following depends on dFmt, f(x) and n
% get the coefficient matrices S, C, M and vecf
%meshWidth=min(0.49,epsilon/b*2.5*log(n));
tic;
mesh0=makeMesh(meshType,n,w);
disp(['Time to makeMesh: ',num2str(toc)]);tic;
[S,Cx,Cy,M,vecf,No2fun,fun2No]=getCoeffs2D(mesh0,basis,f,K,'fOrderUpperBound',fUb);
disp(['Time to getCoeffs: ',num2str(toc)]);tic;

% solve
H=epsilon*S+b(1)*Cx+b(2)*Cy+c*M;
u=H\vecf;
disp(['Time to solve linear system: ',num2str(toc)]);tic;
% save;

%% interpolant
[numSol,xSample,ySample]=combineSolution(u,200,mesh0,fun2No,'Lobatto',K);
xSample=(xSample+1)/2;   % map range to [0,1]x[0,1]
ySample=(ySample+1)/2;
disp(['Time to interpolant: ',num2str(toc)]);

%% plot
figure();
axis;hold on;
for i=1:mesh0.Nsurfaces
    surf(xSample*mesh0.surfaces.hx(i)+mesh0.surfaces.x(1,i),ySample*mesh0.surfaces.hy(i)+mesh0.surfaces.y(1,i),numSol{i},'lineStyle','none');
end

% refine plot
title(['K=',num2str(K)]);
xlabel('$$x$$','interpreter','latex');
ylabel('$$y$$','interpreter','latex');
% set('fontsize',12);

