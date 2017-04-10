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
        mesh.Nnodes=Nnodes;
        mesh.Nedges=Nedges;
        mesh.Nsurfaces=Nsurfaces;
        mesh.nodes=struct('x',zeros(Nnodes,1),...
                          'y',zeros(Nnodes,1),...
                          's',zeros(4,Nnodes),...
                          'onBoundary',zeros(Nnodes,1));
        mesh.edges=struct('n',zeros(2,Nedges),...
                          's',zeros(2,Nedges),...
                          'onBoundary',zeros(Nedges,1));
        mesh.surfaces=struct('n',zeros(4,Nsurfaces),...
                             'e',zeros(4,Nsurfaces),...
                             'x',zeros(2,Nsurfaces),...
                             'y',zeros(2,Nsurfaces),...
                             'hx',zeros(Nsurfaces,1),...
                             'hy',zeros(Nsurfaces,1));
        
        %---------- set nodes -------------------
        top=1;
        % init all nodes
        for j=1:Ny  % row num
            for i=1:Nx  % col num
                mesh.nodes.x(top)=xList(i);
                mesh.nodes.y(top)=yList(j);
                mesh.nodes.s(:,top)=[(j-1)*(Nx-1)+i;  (j-1)*(Nx-1)+i-1;  (j-2)*(Nx-1)+i-1;   (j-2)*(Nx-1)+i];
                %mesh.nodes.onBoundary(top)=0;  % not necessary
                top=top+1;
            end
        end
        if top~=mesh.Nnodes+1
            error('Nnodes not match.');
        end
        % reset boundary nodes
        for i=1:Nx
            mesh.nodes.onBoundary(i)=1;
            mesh.nodes.s([3,4],i)=0;
            mesh.nodes.onBoundary((Ny-1)*Nx+i)=1;
            mesh.nodes.s([1,2],(Ny-1)*Nx+i)=0;
        end
        for j=1:Ny
            mesh.nodes.onBoundary((j-1)*Nx+1)=1;
            mesh.nodes.s([2,3],   (j-1)*Nx+1)=0;
            mesh.nodes.onBoundary(j*Nx)=1;
            mesh.nodes.s([1,4],   j*Nx)=0;
        end
        %--------- set edges ----------------------
        top=1;
        % edges that parallel to x axis
        for j=1:Ny        % row num
            for i=1:Nx-1
                n1=(j-1)*Nx+i;  n2=(j-1)*Nx+i+1;
                mesh.edges.n(:,top)=[n1;  n2];
                mesh.edges.s(:,top)=[mesh.nodes.s(1,n1); mesh.nodes.s(4,n1)];
                mesh.edges.onBoundary(top)= (mesh.edges.s(1,top)==0 || mesh.edges.s(2,top)==0);
                top=top+1;
            end
        end
        % edges that parallel to y axis
        for j=1:Ny-1    % row num
            for i=1:Nx
                n1=(j-1)*Nx+i;  n2=j*Nx+i;
                mesh.edges.n(:,top)=[n1;   n2];
                mesh.edges.s(:,top)=[mesh.nodes.s(2,n1); mesh.nodes.s(1,n1)];
                mesh.edges.onBoundary(top)= (mesh.edges.s(1,top)==0 || mesh.edges.s(2,top)==0);
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
                mesh.surfaces.n(:,top)=[j*Nx+i+1; j*Nx+i; (j-1)*Nx+i; (j-1)*Nx+i+1];
                mesh.surfaces.e(:,top)=[top+Nx-1;
                                        offset+(j-1)*Nx+i   ;
                                        top   ;
                                        offset+(j-1)*Nx+i+1  ];
                mesh.surfaces.x(:,top)=mesh.nodes.x( mesh.surfaces.n([3;1],top) );
                mesh.surfaces.y(:,top)=mesh.nodes.y( mesh.surfaces.n([3;1],top) );
                mesh.surfaces.hx(top)=mesh.surfaces.x(2,top)-mesh.surfaces.x(1,top);
                mesh.surfaces.hy(top)=mesh.surfaces.y(2,top)-mesh.surfaces.y(1,top);
                top=top+1;
            end
        end
        if top~=mesh.Nsurfaces+1
            error('Nsurfaces not match.');
        end
end

