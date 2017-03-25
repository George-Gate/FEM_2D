function mesh = makeMesh_box( xList, yList)
%% A sub function called by function makeMesh()
% Generate a rectangle mesh for a box domain. 
% Input xList, yList: coordinate division for x, y axis
        Nx=length(xList);
        Ny=length(yList);
        
        Nnodes=Nx*Ny;
        Nedges=Nx*(Ny-1)+Ny*(Nx-1);
        Nsurfaces=(Nx-1)*(Ny-1);
        
        % struct def
        typ_nodes=struct('x',0,...
                         'y',0,...
                         's',[0;0;0;0],...
                         'onBoundary',0);
        typ_edges=struct('n',[0;0],...
                         's',[0;0],...
                         'onBoundary',0);
        typ_surfaces=struct('n',[0;0;0;0],...
                            'e',[0;0;0;0],...
                            'x',[0;0],...
                            'y',[0;0],...
                            'hx',0.0,...
                            'hy',0);
        
        mesh.Nnodes=Nx*Ny;
        mesh.Nedges=Nx*(Ny-1)+Ny*(Nx-1);
        mesh.Nsurfaces=(Nx-1)*(Ny-1);
        mesh.nodes=repmat(typ_nodes,Nnodes,1);
        mesh.edges=repmat(typ_edges,Nedges,1);
        mesh.surfaces=repmat(typ_surfaces,Nsurfaces,1);
        
        tNode=typ_nodes;
        tEdge=typ_edges;
        tSurface=typ_surfaces;
        %---------- set nodes -------------------
        top=1;
        % init all nodes
        for j=1:Ny  % row num
            for i=1:Nx  % col num
                mesh.nodes(top).x=xList(i);
                mesh.nodes(top).y=yList(j);
                mesh.nodes(top).s=[(j-1)*(Nx-1)+i;  (j-1)*(Nx-1)+i-1;  (j-2)*(Nx-1)+i-1;   (j-2)*(Nx-1)+i];
                mesh.nodes(top).onBoundary=0;
                top=top+1;
            end
        end
        if top~=mesh.Nnodes+1
            error('Nnodes not match.');
        end
        % reset boundary nodes
        for i=1:Nx
            mesh.nodes(i).onBoundary=1;
            mesh.nodes(i).s([3,4])=0;
            mesh.nodes((Ny-1)*Nx+i).onBoundary=1;
            mesh.nodes((Ny-1)*Nx+i).s([1,2])=0;
        end
        for j=1:Ny
            mesh.nodes((j-1)*Nx+1).onBoundary=1;
            mesh.nodes((j-1)*Nx+1).s([2,3])=0;
            mesh.nodes(j*Nx).onBoundary=1;
            mesh.nodes(j*Nx).s([1,4])=0;
        end
        %--------- set edges ----------------------
        top=1;
        % edges parallel to x axis
        for j=1:Ny        % row num
            for i=1:Nx-1
                tEdge.n=[(j-1)*Nx+i;        (j-1)*Nx+i+1];
                tNode=mesh.nodes(tEdge.n(1));
                tEdge.s=[tNode.s(1); tNode.s(4)];
                tEdge.onBoundary=tNode.onBoundary * mesh.nodes(tEdge.n(2)).onBoundary;
                mesh.edges(top)=tEdge;
                top=top+1;
            end
        end
        % edges parallel to y axis
        for j=1:Ny-1    % row num
            for i=1:Nx
                tEdge.n=[(j-1)*Nx+i;        j*Nx+i];
                tNode=mesh.nodes((j-1)*Nx+i);
                tEdge.s=[tNode.s(2); tNode.s(1)];
                tEdge.onBoundary=tNode.onBoundary * mesh.nodes(tEdge.n(2)).onBoundary;
                mesh.edges(top)=tEdge;
                top=top+1;
            end
        end
        if top~=mesh.Nedges+1
            error('Nedges not match.');
        end
        %---------- set surfaces -----------------------
        offset=Ny*(Nx-1);
        top=1; 
        for j=1:Ny-1     % row num
            for i=1:Nx-1
                tSurface.n=[j*Nx+i+1; j*Nx+i; (j-1)*Nx+i; (j-1)*Nx+i+1];
                tSurface.e=[top+Nx-1;
                            offset+(j-1)*Nx+i   ;
                            top   ;
                            offset+(j-1)*Nx+i+1  ];
                tSurface.x=[mesh.nodes(tSurface.n(3)).x;mesh.nodes(tSurface.n(1)).x];
                tSurface.y=[mesh.nodes(tSurface.n(3)).y;mesh.nodes(tSurface.n(1)).y];
                tSurface.hx=tSurface.x(2)-tSurface.x(1);
                tSurface.hy=tSurface.y(2)-tSurface.y(1);
                mesh.surfaces(top)=tSurface;
                top=top+1;
            end
        end
        if top~=mesh.Nsurfaces+1
            error('Nsurfaces not match.');
        end
end

