% Finite Element Method Solver

%% parameters
b=[1;0];
c=100;
k=1;
f=@(x)x.^k;
epsilon=1e-5;
n={[2^7+1];[2^7]};
w={[];[]};
meshType='boxSegUniform';
basis='Linear';


%% numerical solution    
% the following depends on dFmt, f(x) and n
% get the coefficient matrices S, C, M and vecf
%meshWidth=min(0.49,epsilon/b*2.5*log(n));
mesh0=makeMesh(meshType,n,w);
[S,Cx,Cy,M,vecf,id2fun,fun2id]=getCoeffs2D(mesh0,basis,@(x,y)x.^2+1);

% the following depends on n, epsilon, b and c
H=epsilon*S+b(1)*Cx+b(2)*Cy+c*M;

% solve
u=H\vecf;
% job=batch('[u,D]=eigs(H,5,''sm'')');
% wait(job);
% load(job);
% delete(job);


%% interpolant
% get node coordinates
N=mesh0.Nnodes;
xList=zeros(N,1);yList=zeros(N,1);
for i=1:N
    xList(i)=mesh0.nodes(i).x;
    yList(i)=mesh0.nodes(i).y;
end

% prepare solution data
Ninner=length(u);
tmp_u=[u(:,1);zeros(N-Ninner,1)];
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