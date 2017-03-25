figure();
displayText=0;
%% plot node
posInner=zeros(mesh.Nnodes,3);
posBorder=zeros(mesh.Nnodes,3);
tI=1;tB=1;
for i=1:mesh.Nnodes
    if (mesh.nodes(i).onBoundary)
        color='r^';
    else
        color='bo';
    end
    p=mesh.nodes(i);
    x=p.x;
    y=p.y;
    scatter(x,y,color);hold on;
    if displayText
        text(x,y,num2str(i),'VerticalAlignment','bottom');
        text(x,y,['     ',num2str(p.s(1))],'HorizontalAlignment','left','VerticalAlignment','bottom','color','r');
        text(x,y,[num2str(p.s(2)),'     '],'HorizontalAlignment','right','VerticalAlignment','bottom','color','r');
        text(x,y,[num2str(p.s(3)),'     '],'HorizontalAlignment','right','VerticalAlignment','top','color','r');
        text(x,y,['     ',num2str(p.s(4))],'HorizontalAlignment','left','VerticalAlignment','top','color','r');
    end
end


%% plot edge
for i=1:mesh.Nedges
    if (mesh.edges(i).onBoundary)
        color='r';
    else
        color='b';
    end
    n=mesh.edges(i).n;
    x=[mesh.nodes(n(1)).x,mesh.nodes(n(2)).x];
    y=[mesh.nodes(n(1)).y,mesh.nodes(n(2)).y];
    line(x,y,'color',color);
    if displayText
        text(sum(x)/2,sum(y)/2,[num2str([mesh.edges(i).s(1),i,mesh.edges(i).s(2)])],'color','b');
    end
end

%% plot surface
if displayText
    for i=1:mesh.Nsurfaces
        s=mesh.surfaces(i);
        x=sum(s.x)/2;
        y=sum(s.y)/2;
        text(x,y,num2str([s.n(2),s.e(1),s.n(1);s.e(2),i,s.e(4);s.n(3),s.e(3),s.n(4)]),'color','g');
    end
end