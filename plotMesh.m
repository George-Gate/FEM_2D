figure();
displayText=0;
displayRelation=0;
%% plot node
idOnB=find(mesh0.nodes.onBoundary == 1);
idOffB=find(mesh0.nodes.onBoundary == 0);
x=mesh0.nodes.x;
y=mesh0.nodes.y;
scatter(x(idOnB ),y(idOnB ),'c^');hold on;
scatter(x(idOffB),y(idOffB),'bo');

if displayText
    % Disp the id of node
    text(x,y,num2str((1:mesh0.Nnodes)'),'VerticalAlignment','bottom','color','r');
    if displayRelation
        % Disp nodes.s
        text(x,y,num2str(mesh0.nodes.s(1,:)','\\qquad %i'),'HorizontalAlignment','left','VerticalAlignment','bottom','color','r','interpreter','latex');
        text(x,y,num2str(mesh0.nodes.s(2,:)','%i \\quad'),'HorizontalAlignment','right','VerticalAlignment','bottom','color','r','interpreter','latex');
        text(x,y,num2str(mesh0.nodes.s(3,:)','%i \\quad'),'HorizontalAlignment','right','VerticalAlignment','top','color','r','interpreter','latex');
        text(x,y,num2str(mesh0.nodes.s(4,:)','\\quad %i'),'HorizontalAlignment','left','VerticalAlignment','top','color','r','interpreter','latex');
    end
end

%% plot edge
idOnB=find(mesh0.edges.onBoundary == 1);
idOffB=find(mesh0.edges.onBoundary == 0);
clear x y;
x(1,:)=mesh0.nodes.x(mesh0.edges.n(1,:))';y(1,:)=mesh0.nodes.y(mesh0.edges.n(1,:))';
x(2,:)=mesh0.nodes.x(mesh0.edges.n(2,:))';y(2,:)=mesh0.nodes.y(mesh0.edges.n(2,:))';
line(x(:,idOnB ),y(:,idOnB ),'color','c');
line(x(:,idOffB),y(:,idOffB),'color','b');

if displayText
    if displayRelation
        % Disp the id of edge and adjacent surfaces
        text(sum(x,1)'/2,sum(y,1)'/2,num2str([mesh0.edges.s(1,:)',(1:mesh0.Nedges)',mesh0.edges.s(2,:)']),'color','b','HorizontalAlignment','center');
    else
        % Disp the id of edge
        text(sum(x,1)'/2,sum(y,1)'/2,num2str((1:mesh0.Nedges)'),'color','b','HorizontalAlignment','center');
    end
end



%% plot surface

if displayText
    if displayRelation
        for i=1:mesh0.Nsurfaces
            x=sum(mesh0.surfaces.x(:,i))/2;
            y=sum(mesh0.surfaces.y(:,i))/2;
            text(x,y,num2str([mesh0.surfaces.n(2,i),mesh0.surfaces.e(1,i),mesh0.surfaces.n(1,i);  ...
                              mesh0.surfaces.e(2,i),        i            ,mesh0.surfaces.e(4,i);       ...
                              mesh0.surfaces.n(3,i),mesh0.surfaces.e(3,i),mesh0.surfaces.n(4,i)]),'color','k','HorizontalAlignment','center');
        end
    else
        x=mesh0.surfaces.x;
        y=mesh0.surfaces.y;
        text(sum(x,1)'/2,sum(y,1)'/2,num2str((1:mesh0.Nsurfaces)'),'color','k','HorizontalAlignment','center');
    end
end

clear x y;
%%
set(gca,'xTick',[0 1],'xTickMode','manual');
set(gca,'yTick',[0 1],'yTickMode','manual');
xlabel('x');ylabel('y')