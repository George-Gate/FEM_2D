figure();
displayText=1;
displayRelation=0;
%% plot node
posInner=zeros(mesh0.Nnodes,3);
posBorder=zeros(mesh0.Nnodes,3);
tI=1;tB=1;
for i=1:mesh0.Nnodes
    if (mesh0.nodes(i).onBoundary)
        color='c^';
    else
        color='bo';
    end
    p=mesh0.nodes(i);
    x=p.x;
    y=p.y;
    scatter(x,y,color);hold on;
    if displayText
        % Disp the id of node
        text(x,y,num2str(i),'VerticalAlignment','bottom','color','r');
        if displayRelation
            % Disp nodes{i}.s
            text(x,y,['     ',num2str(p.s(1))],'HorizontalAlignment','left','VerticalAlignment','bottom','color','r');
            text(x,y,[num2str(p.s(2)),'     '],'HorizontalAlignment','right','VerticalAlignment','bottom','color','r');
            text(x,y,[num2str(p.s(3)),'     '],'HorizontalAlignment','right','VerticalAlignment','top','color','r');
            text(x,y,['     ',num2str(p.s(4))],'HorizontalAlignment','left','VerticalAlignment','top','color','r');
        end
    end
end


%% plot edge
for i=1:mesh0.Nedges
    if (mesh0.edges(i).onBoundary)
        color='c';
    else
        color='b';
    end
    n=mesh0.edges(i).n;
    x=[mesh0.nodes(n(1)).x,mesh0.nodes(n(2)).x];
    y=[mesh0.nodes(n(1)).y,mesh0.nodes(n(2)).y];
    line(x,y,'color',color);
    if displayText
        if displayRelation
            % Disp the id of edge and adjacent surfaces
            text(sum(x)/2,sum(y)/2,num2str([mesh0.edges(i).s(1),i,mesh0.edges(i).s(2)]),'color','b','HorizontalAlignment','center');
        else
            % Disp the id of edge
            text(sum(x)/2,sum(y)/2,num2str(i),'color','b','HorizontalAlignment','center');
        end
    end
end

%% plot surface
if displayText
    for i=1:mesh0.Nsurfaces
        s=mesh0.surfaces(i);
        x=sum(s.x)/2;
        y=sum(s.y)/2;
        if displayRelation
            text(x,y,num2str([s.n(2),s.e(1),s.n(1);s.e(2),i,s.e(4);s.n(3),s.e(3),s.n(4)]),'color','k');
        else
            text(x,y,num2str(i),'color','k');
        end
    end
end

%%
set(gca,'xTick',[0 1],'xTickMode','manual');
set(gca,'yTick',[0 1],'yTickMode','manual');
xlabel('x');ylabel('y')