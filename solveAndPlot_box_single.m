% Finite Element Method Solver

%% parameters
b=[1;-0.5];
c=1;
% f=@(x,y)x.^2+y+1;
f=1;
epsilon=1e-3;

n={[2^6;2^6];[2^6;2^6]};
tau(1)=min(0.5,epsilon./abs(b(1))*2.5*log(n{1}(1)+1));
tau(2)=min(0.5,epsilon./abs(b(2))*2.5*log(n{2}(1)+1));
w={[1-tau(1)];[tau(2)]};
meshType='boxSegUniform';
basis='BiLinear';


%% numerical solution    
% the following depends on dFmt, f(x) and n
% get the coefficient matrices S, C, M and vecf
%meshWidth=min(0.49,epsilon/b*2.5*log(n));
mesh0=makeMesh(meshType,n,w);
[S,Cx,Cy,M,vecf,id2fun,fun2id]=getCoeffs2D(mesh0,basis,f);

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
xList=mesh0.nodes.x;yList=mesh0.nodes.y;

% prepare solution data
Ninner=length(u);
tmp_u=[u(:,1);zeros(N-Ninner,1)];
tmp_x=[xList(id2fun(1:Ninner));xList(fun2id==0)];
tmp_y=[yList(id2fun(1:Ninner));yList(fun2id==0)];

% interpolant
numSol=scatteredInterpolant(tmp_x,tmp_y,tmp_u);

%% plot
figure();
nxy=[200;200];
coList=cell(2,1);
for i=1:2
    coList{i}=0;
    len=length(w{i});
    for j=1:len
        tmp=linspace(coList{i}(end),coList{i}(end)+w{i}(j),nxy(i)+1)';
        coList{i}=[coList{i};tmp(2:end)];
    end
    tmp=linspace(coList{i}(end),1,nxy(i)+1)';
    coList{i}=[coList{i};tmp(2:end)];
end
[x,y]=meshgrid( coList{1} , coList{2} );
% contour(x,y,imag(numSol(x,y)));hold on;
mesh(x,y,(numSol(x,y)));hold on;

% refine plot
% title(['$$N=',num2str(N),'\quad \varepsilon=$$',num2str(epsilon),'$$\quad b=',num2str(b),'\quad c=',num2str(c),'\quad f(x)=x^k, k=',num2str(k),'$$  dFmt=',dFmt],'interpreter','latex');
xlabel('$$x$$','interpreter','latex');
ylabel('$$y$$','interpreter','latex');
% set(ax(1),'fontsize',12,'ylim',[0,1],'Ycolor','black');
% set(ax(2),'fontsize',12,'Ycolor','black');