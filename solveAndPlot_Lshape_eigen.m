% Finite Element Method Solver

%% parameters
epsilon=1;
n={[2^8+1];[2^8];[2^8];[2^8]};
w={[];[];[];[]};
meshType='LshapeSegUniform';
basis='Linear';


%% numerical solution    
% the following depends on dFmt, f(x) and n
% get the coefficient matrices S, C, M and vecf
%meshWidth=min(0.49,epsilon/b*2.5*log(n));
mesh0=makeMesh(meshType,n,w);
[S,Cx,Cy,M,vecf,id2fun,fun2id]=getCoeffs2D(mesh0,basis);

% the following depends on n, epsilon, b and c
H=epsilon*S;

% solve
% u=H\vecf;
job=batch('[u,D]=eigs(H,100,''sa'')');
wait(job);
load(job);
delete(job);

save;

%% interpolant
% get node coordinates
N=mesh0.Nnodes;
xList=mesh0.nodes.x;yList=mesh0.nodes.y;

% prepare solution data
Ninner=length(u);
tmp_u=[u(:,5);zeros(N-Ninner,1)];
tmp_x=[xList(id2fun(1:Ninner));xList(fun2id==0)];
tmp_y=[yList(id2fun(1:Ninner));yList(fun2id==0)];

% interpolant
numSol=scatteredInterpolant(tmp_x,tmp_y,tmp_u);

%% plot
figure();
nx=200;ny=200;
[x,y]=meshgrid( linspace(min(xList),max(xList),nx) , linspace(min(yList),max(yList),ny) );
% contour(x,y,imag(numSol(x,y)));hold on;
mesh(x,y,(numSol(x,y)));hold on;

% refine plot
% title(['$$N=',num2str(N),'\quad \varepsilon=$$',num2str(epsilon),'$$\quad b=',num2str(b),'\quad c=',num2str(c),'\quad f(x)=x^k, k=',num2str(k),'$$  dFmt=',dFmt],'interpreter','latex');
xlabel('$$x$$','interpreter','latex');
ylabel('$$y$$','interpreter','latex');
% set(ax(1),'fontsize',12,'ylim',[0,1],'Ycolor','black');
% set(ax(2),'fontsize',12,'Ycolor','black');