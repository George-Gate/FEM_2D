%function [ No2fun, fun2No, Nbasis, hxList, hyList, xList, yList, Svec, Cxvec, Cyvec, Mvec ] = getCoeffs2D_Lobatto( mesh,M )
%Generate the coeff matrices for Lobatto basis
%   Called by function getCoeffs2D()
% [Def] id: the index of an object (either node, edge or surface) in the mesh structure
%       No: the numbering of a basis in the linear equation system. 
    meshEps=eps;
%---------------- fun2No and No2fun -----------------------------
    % initialize fun2No and No2fun
    Nnodal=mesh.Nnodes-sum(mesh.nodes.onBoundary);
    Nedge=(mesh.Nedges-sum(mesh.edges.onBoundary))*M;
    Nface=mesh.Nsurfaces*M*M;
    Nbasis=Nnodal+Nedge+Nface;
    No2fun.name=strings(Nbasis,1);  % possible values: nodal, edge, face
    No2fun.objid=zeros(Nbasis,1);   % to specify which node/edge/surface the basis is bound for.
    No2fun.subid=zeros(Nbasis,2);   % for nodal mode: 1, for edge mode: [1,2~M+1] or [2~M+1,1], for face mode: [2~M+1,2~M+1];
                                    % .subid(:,1): subid for x dirl .subid(:,2): subid for y dir
    fun2No.nodal=zeros(mesh.Nnodes,1);   % the No. for nodal mode bases
    fun2No.edge=zeros(M+1,mesh.Nedges);   % the No. for edge mode bases
    fun2No.face=zeros(M+1,M+1,mesh.Nsurfaces);  % the No. for face mode bases, fun2No.face(Xsubid,Ysubid,objid)
    
    % set fun2No and No2fun
    % nodal mode
    No2fun.name(1:Nnodal)='nodal';
    No2fun.objid(1:Nnodal)=find(mesh.nodes.onBoundary==0);
    No2fun.subid(1:Nnodal,1:2)=1;
    fun2No.nodal(No2fun.objid(1:Nnodal))=(1:Nnodal)';
    % edge mode
    No2fun.name(Nnodal+1:Nnodal+Nedge)='edge';
    offset=Nnodal+1;
    for i=find(mesh.edges.onBoundary==0)'
        n=mesh.edges.n(1:2,i);
        if (abs( mesh.nodes.x(n(1)) - mesh.nodes.x(n(2)) )<=meshEps)
            % edge is along y axis
            for k=2:M+1
                No2fun.objid(offset)=i;
                No2fun.subid(offset,1)=1;
                No2fun.subid(offset,2)=k;
                fun2No.edge(k,i)=offset;
                offset=offset+1;
            end
        else
            % edge is along x axis
            for k=2:M+1
                No2fun.objid(offset)=i;
                No2fun.subid(offset,1)=k;
                No2fun.subid(offset,2)=1;
                fun2No.edge(k,i)=offset;
                offset=offset+1;
            end
        end
    end
    % face mode
    No2fun.name(offset:offset+Nface-1)='face';
    xid=repmat(2:M+1,M,1);xid=xid(:);
    yid=repmat(2:M+1,1,M)';
    for i=1:mesh.Nsurfaces
        No2fun.objid(offset:offset+M*M-1)=i;
        No2fun.subid(offset:offset+M*M-1,1)=xid;  % = 2 2 ... 2 3 3 ... 3 4 4 ... 4 5...
        No2fun.subid(offset:offset+M*M-1,2)=yid;  % = 2 3 4 ... 2 3 4 ... 2 3 4 ...
        fun2No.face(2:M+1,2:M+1,i)=reshape(offset:offset+M*M-1,M,M)';
        offset=offset+M*M;
    end
    
%--------------------- Read Mesh ----------------------------------------
    % extract node from mesh
    Nobj=mesh.Nnodes+mesh.Nedges+mesh.Nsurfaces;
    adjNodeList=zeros(4,2,Nobj);    % [For node] adjNode(:,1): the id of four corner nodes, adjNode(:,2): four adjacent nodes, same numbering with edges.                   
                                    % [For edge] adjNode(:,1): the id of four corner nodes, adjNode(1:2,2): the id of edge's end points. When the edge is 
                                    %                          horizontal, the numbering is the same as we used in mesh. As for the edge is vertical, we should 
                                    %                          keep the relative position the same.
                                    % [For surface] adjNode(:,1): the id of four corner nodes
    adjEdgeList=zeros(6,2,Nobj);    % [For node/surface] adjEdge(:,1): the id of vertical(along y axis) edges, adjEdge(:,2): horizontal edges
                                    % [For edge] adjEdge(:,1): the id of two edges that patalell to the edge; adjEdge(:,1): four edges that perpendicular to it.
                                    % numbering:  1  4            1 
                                    %   __ __     __ __     __    __                2  4              1  
                                    % 1|__|__|3  |__|__|  1|__|2 |__|    __ __      __ __      __     __
                                    % 4|__|__|6  |__|__|  3|__|4 |__|  1|__|__|2   |__|__|   1|__|2  |__| 
                                    %             3   6           2                 1  3              2 
    adjFaceList=zeros(4,Nobj);      % 4~1 adjacent surface for the node/edge/surface, the same numbering as we used in mesh. 
    hxList=zeros(2,Nobj);hyList=zeros(2,Nobj);
    xList=zeros(Nobj,1);  % [For node] the (x,y) for that node; [For edge] the (x,y) for n(1) of that edge
    yList=zeros(Nobj,1);  % [For surface] the (x,y) for n(3) of that surface
    
    % for node
    for i=find(mesh.nodes.onBoundary==0)'
        Nid=i;
        face=mesh.nodes.s(:,i);
        adjFaceList(1:4,Nid)=face;
        for j=1:4
            adjNodeList(j,1,Nid)=mesh.surfaces.n(j         , face(j));
            adjNodeList(j,2,Nid)=mesh.surfaces.n(mod(j,4)+1, face(j));
        end
        % vertical edges
        adjEdgeList(1:6,1,Nid)=[mesh.surfaces.e([2,4],face(2)) ; mesh.surfaces.e(4,face(1)) ; mesh.surfaces.e([2,4],face(3)) ; mesh.surfaces.e(4,face(4))];
        % horizontal edges
        adjEdgeList(1:6,2,Nid)=[mesh.surfaces.e([1,3],face(2)) ; mesh.surfaces.e(3,face(3)) ; mesh.surfaces.e([1,3],face(1)) ; mesh.surfaces.e(3,face(4))];
        hxList(:,Nid)=[mesh.surfaces.hx(face(3));  mesh.surfaces.hx(face(1))];
        hyList(:,Nid)=[mesh.surfaces.hy(face(3));  mesh.surfaces.hy(face(1))];
        xList(Nid)=mesh.nodes.x(i);yList(Nid)=mesh.nodes.y(i);
    end
    
    % for edge 
    for i=find(mesh.edges.onBoundary==0)'
        Eid=mesh.Nnodes+i;
        face=mesh.edges.s(:,i);
        adjFaceList(1:2,Eid)=face;
        adjNodeList(1:4,1,Eid)=[mesh.surfaces.n(1:2, face(1));mesh.surfaces.n(3:4, face(2))];
        adjNodeList(1:2,2,Eid)=mesh.edges.n(:,i);
        if abs(diff(mesh.nodes.x(adjNodeList(1:2,2,Eid))))<meshEps
            % vertical edge
            % para. edges
            adjEdgeList(1:2,1,Eid)=[mesh.surfaces.e(2,face(1)) ; mesh.surfaces.e(4,face(2))];
            % perp. edges
            adjEdgeList(1:4,2,Eid)=[mesh.surfaces.e([3,1],face(1)) ; mesh.surfaces.e([3,1],face(2))];
            hxList(:,Eid)=[mesh.surfaces.hx(face(1));  mesh.surfaces.hx(face(2))];
            hyList(1,Eid)=mesh.surfaces.hy(face(1));
        else
            % horizontal edge
            % para. edges
            adjEdgeList(1:2,1,Eid)=[mesh.surfaces.e(1,face(1)) ; mesh.surfaces.e(3,face(2))];
            % perp. edges
            adjEdgeList(1:4,2,Eid)=[mesh.surfaces.e([2,4],face(1)) ; mesh.surfaces.e([2,4],face(2))];
            hxList(1,Eid)=mesh.surfaces.hx(face(1));
            hyList(:,Eid)=[mesh.surfaces.hy(face(2));  mesh.surfaces.hy(face(1))];
        end
        xList(Eid)=mesh.nodes.x(adjNodeList(1,2,Eid));yList(Eid)=mesh.nodes.y(adjNodeList(1,2,Eid));
    end
    
    % for surface
    for i=1:mesh.Nsurfaces
        Sid=mesh.Nnodes+mesh.Nedges+i;
        adjNodeList(1:4,1,Sid)=mesh.surfaces.n(:,i);
        % vertical edges
        adjEdgeList(1:2,1,Sid)=mesh.surfaces.e([2,4],i);
        % horizontal edges
        adjEdgeList(1:2,2,Sid)=mesh.surfaces.e([1,3],i);
        hxList(1,Sid)=mesh.surfaces.hx(i);
        hyList(1,Sid)=mesh.surfaces.hy(i);
        xList(Sid)=mesh.nodes.x(adjNodeList(3,1,Sid));yList(Sid)=mesh.nodes.y(adjNodeList(3,1,Sid));
    end
    
%-------------------- The following do not read mesh any more -----------------------
    % generate coeff vectors
    Svec=zeros(100*Nbasis,3);
    Cxvec=zeros(100*Nbasis,3);
    Cyvec=zeros(100*Nbasis,3);
    Mvec=zeros(100*Nbasis,3);
    topS=1;topM=1;topCx=1;topCy=1;
    for iNo=1:Nbasis
        % consider the line for ( basis(iNo) , u )
        switch No2fun.name(iNo)
            case 'nodal'
                Nid=No2fun.objid(iNo);
                No_adjNode=fun2No.nodal(adjNodeList(:,:,Nid));
                adjEdge=adjEdgeList(:,:,Nid);
                adjFace=adjFaceList(:,Nid);
                hx=hxList(1:2,Nid);
                hy=hyList(1:2,Nid);
                % ========== construct M matrix ===============    (A(x)B(y),C(x)D(y))=(A,C)(B,D)
                Nline=9+12*2+4*4;
                Mvec(topM:topM+Nline-1,1)=iNo;  
                % overlap with 9 adjacent nodal basis
                Mvec(topM:topM+8,2)=[No_adjNode(:,1);No_adjNode(:,2);iNo];
                Mvec(topM:topM+8,3)=[hx([2;1;1;2])/6;sum(hx)/3;hx(1)/6  ;sum(hx)/3;hx(2)/6  ;sum(hx)/3]...  % x part of the integral, (A,C)
                                  .*[hy([2;2;1;1])/6;hy(2)/6  ;sum(hy)/3;hy(1)/6  ;sum(hy)/3;sum(hy)/3];    % y part of the integral, (B,D)
                % overlap with 12 adjacent edge basis group
                Mvec(topM+9:topM+9+11,2)=fun2No.edge(2,[adjEdge(1:6,1);adjEdge(1:6,2)])';  % m=2
                Mvec(topM+9:topM+9+11,3)=[hx(1)/6;sum(hx)/3;hx(2)/6;hx(1)/6;sum(hx)/3;hx(2)/6;        -hx([1;1;1;2;2;2])/3    ]...                     % x part of the integral 
                                       .*[              -hy([2;2;2;1;1;1])/3                 ;hy(2)/6;sum(hy)/3;hy(1)/6;hy(2)/6;sum(hy)/3;hy(1)/6];    % y part of the integral 
                Mvec(topM+9+12:topM+9+23,2)=fun2No.edge(3,[adjEdge(1:6,1);adjEdge(1:6,2)])'; % m=3
                Mvec(topM+9+12:topM+9+23,3)=[hx(1)/6;sum(hx)/3;hx(2)/6;hx(1)/6;sum(hx)/3;hx(2)/6;    -hx([1;1;1])/15      ;    hx([2;2;2])/15       ]...  % x part of the integral 
                                          .*[     hy([2;2;2])/15      ;   -hy([1;1;1])/15       ;hy(2)/6;sum(hy)/3;hy(1)/6;hy(2)/6;sum(hy)/3;hy(1)/6];    % y part of the integral 
                % overlap with 4 adjacent face basis group
                Mvec(topM+9+24:topM+9+24+15,2)=[reshape(fun2No.face(2,2,adjFace(1:4)),4,1);reshape(fun2No.face(2,3,adjFace(1:4)),4,1);...
                                                reshape(fun2No.face(3,2,adjFace(1:4)),4,1);reshape(fun2No.face(3,3,adjFace(1:4)),4,1)];  
                Mvec(topM+9+24:topM+9+24+15,3)=[  -hx([2;1;1;2;2;1;1;2])/3                       ;  hx([2;1;1;2;2;1;1;2])/15.*[1;-1;-1;1;1;-1;-1;1] ]...  % x part of the integral 
                                             .*[  -hy([2;2;1;1])/3; hy([2;2;1;1])/15.*[1;1;-1;-1];  -hy([2;2;1;1])/3; hy([2;2;1;1])/15.*[1;1;-1;-1] ];    % y part of the integral
                topM=topM+Nline;
                % =========== construct S matrix ===============
                Nline=9+24;
                Svec(topS:topS+Nline-1,1)=iNo;
                % overlap with 9 adjacent nodal basis
                Svec(topS:topS+8,2)=[No_adjNode(:,1);No_adjNode(:,2);iNo];
                Svec(topS:topS+8,3)=[-1./hx([2;1;1;2]); sum(1./hx); -1/hx(1)  ; sum(1./hx); -1/hx(2)  ; sum(1./hx) ]...  % (A',C')(B,D)
                                  .*[  hy([2;2;1;1])/6;  hy(2)/6  ;sum(hy)/3  ;  hy(1)/6  ;sum(hy)/3  ; sum(hy)/3]...
                                   +[  hx([2;1;1;2])/6;  sum(hx)/3; hx(1)/6   ;  sum(hx)/3;  hx(2)/6  ; sum(hx)/3]...  % (A,C)(B',D')
                                  .*[-1./hy([2;2;1;1]); -1/hy(2)  ; sum(1./hy); -1/hy(1)  ; sum(1./hy); sum(1./hy) ];
                % overlap with 12 adjacent edge basis group
                Svec(topS+9:topS+9+11,2)=fun2No.edge(2,[adjEdge(1:6,1);adjEdge(1:6,2)])';  % m=2
                Svec(topS+9:topS+9+11,3)=[-1/hx(1);sum(1./hx);-1/hx(2);-1/hx(1);sum(1./hx);-1/hx(2);        -hx([1;1;1;2;2;2])/3    ]...              % (A',C')(B,D) for vertical edge
                                       .*[              -hy([2;2;2;1;1;1])/3          ;-1/hy(2);sum(1./hy);-1/hy(1);-1/hy(2);sum(1./hy);-1/hy(1)];    % and (A,C)(B',D') for horizontal
                Svec(topS+9+12:topS+9+23,2)=fun2No.edge(3,[adjEdge(1:6,1);adjEdge(1:6,2)])'; % m=3
                Svec(topS+9+12:topS+9+23,3)=[-1/hx(1);sum(1./hx);-1/hx(2);-1/hx(1);sum(1./hx);-1/hx(2);    -hx([1;1;1])/15      ;    hx([2;2;2])/15       ]...   % (A',C')(B,D) for vertical edge
                                          .*[     hy([2;2;2])/15      ;   -hy([1;1;1])/15       ;-1/hy(2);sum(1./hy);-1/hy(1);-1/hy(2);sum(1./hy);-1/hy(1)];     % and (A,C)(B',D') for horizontal
                topS=topS+Nline;
                % =========== construct Cx matrix ==============
                Nline=6+14+8;
                Cxvec(topCx:topCx+Nline-1,1)=iNo;
                % overlap with 6 adjacent nodal basis
                Cxvec(topCx:topCx+5,2)=[No_adjNode(:,1);No_adjNode([2;4],2)];
                Cxvec(topCx:topCx+5,3)=[ [1;-1;-1;1     ;-1          ;1]/2 ]...          % (A,C')
                                     .*[ hy([2;2;1;1])/6; sum(hy)/3  ;  sum(hy)/3  ];    % (B,D)
                % overlap with 12 adjacent edge basis group
                Cxvec(topCx+6:topCx+6+9,2)=fun2No.edge(2,[adjEdge([1;3;4;6],1);adjEdge(1:6,2)])';  % m=2
                Cxvec(topCx+6:topCx+6+9,3)=[   [-1;1;1;-1]/2 ;                  [1;1;1;-1;-1;-1]/3 ]...               % (A,C')
                                         .*[ -hy([2;2;1;1])/3; hy(2)/6;sum(hy)/3;hy(1)/6;hy(2)/6;sum(hy)/3;hy(1)/6 ]; % (B,D)
                Cxvec(topCx+6+10:topCx+6+13,2)=fun2No.edge(3,adjEdge([1;3;4;6],1))'; % m=3
                Cxvec(topCx+6+10:topCx+6+13,3)=[   [-1;1;1;-1]/2                  ]...  % (A,C')
                                             .*[ hy([2;2;1;1])/15.*[1;1;-1;-1]    ];    % (B,D)
                % overlap with 4 adjacent face basis group
                Cxvec(topCx+6+14:topCx+6+14+7,2)=[reshape(fun2No.face(2,2,adjFace(1:4)),4,1);reshape(fun2No.face(2,3,adjFace(1:4)),4,1)];
                Cxvec(topCx+6+14:topCx+6+14+7,3)=[   [-1;1;1;-1]/3;   [-1;1;1;-1]/3              ]...  % (A,C')
                                               .*[-hy([2;2;1;1])/3;hy([2;2;1;1])/15.*[1;1;-1;-1] ];    % (B,D)
                topCx=topCx+Nline;
                % ============= construct Cy matrix ===============
                Nline=6+14+8;
                Cyvec(topCy:topCy+Nline-1,1)=iNo;
                % overlap with 6 adjacent nodal basis
                Cxvec(topCx:topCx+5,2)=[No_adjNode(:,1);No_adjNode([1;3],2)];
                Cxvec(topCx:topCx+5,3)=[ hx([2;1;1;2])/6; sum(hx)/3  ;  sum(hx)/3 ]... % (A,C)
                                     .*[ [1;1;-1;-1     ;1           ;-1]/2  ];        % (B,D')
                % overlap with 12 adjacent edge basis group
                Cxvec(topCx+6:topCx+6+9,2)=fun2No.edge(2,[adjEdge([1;3;4;6],2);adjEdge(1:6,1)])';  % m=2
                Cxvec(topCx+6:topCx+6+9,3)=[ -hy([1;1;2;2])/3; hx(1)/6;sum(hx)/3;hx(2)/6;hx(1)/6;sum(hx)/3;hx(2)/6]...  % (A,C)
                                         .*[   [1;-1;1;-1]/2 ;                  [-1;-1;-1;1;1;1]/3  ];                  % (B,D')
                Cxvec(topCx+6+10:topCx+6+13,2)=fun2No.edge(3,adjEdge([1;3;4;6],2))'; % m=3
                Cxvec(topCx+6+10:topCx+6+13,3)=[ hx([1;1;2;2])/15.*[-1;-1;1;1]  ]...  % (A,C)
                                             .*[ [1;-1;1;-1]/2    ];                  % (B,D')
                % overlap with 4 adjacent face basis group
                Cxvec(topCx+6+14:topCx+6+14+7,2)=[reshape(fun2No.face(2,2,adjFace(1:4)),4,1);reshape(fun2No.face(2,3,adjFace(1:4)),4,1)];
                Cxvec(topCx+6+14:topCx+6+14+7,3)=[-hx([2;1;1;2])/3;hx([2;1;1;2])/15.*[1;-1;-1;1]  ]...  % (A,C)
                                               .*[   [-1;-1;1;1]/3;   [-1;-1;1;1]/3               ];    % (B,D')
                topCy=topCy+Nline;
            case 'edge'
                Eid=No2fun.objid(iNo);
                No_adjNode=[fun2No.nodal(adjNodeList(:,1,mesh.Nnodes+Eid)),[fun2No.nodal(adjNodeList(1:2,2,mesh.Nnodes+Eid));0;0]];  % ´¦ÀíÁã
                if abs(diff(mesh.nodes.x(mesh.edges.n(:,Eid))))<=meshEps
                    % for vertical edge
                    hp=hyList(1,mesh.Nnodes+Eid);
                    ho=hxList(1:2,mesh.Nnodes+Eid);
                else                                   %        hp
                    % for horizontal edge              %        __
                    hp=hxList(1,mesh.Nnodes+Eid);      % ho(1) |__|
                    ho=hyList([2,1],mesh.Nnodes+Eid);  % ho(2) |__|
                end
                
            case 'face'
                Sid=No2fun.objid(iNo);
                No_adjNode=fun2No.nodal(adjNodeList(:,1,mesh.Nnodes+mesh.Nedges+Sid));
                
            otherwise
                error('Unknow basis name!');
        end
    end
    % delete lines that have Xvec(:,2)=0, those are boundary elements.

    
%end

