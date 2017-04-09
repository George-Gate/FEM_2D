function [ No2fun, fun2No, Nbasis, hxList, hyList, xList, yList, Svec, Cxvec, Cyvec, Mvec ] = getCoeffs2D_Lobatto( mesh,K )
%Generate the coeff matrices for Lobatto basis
%   Called by function getCoeffs2D()
%        K: cutOff of the Lobatto basis
% [Def] id: the index of an object (either node, edge or surface) in the mesh structure
%       No: the numbering of a basis in the linear equation system. 
    meshEps=eps;
%---------------- fun2No and No2fun -----------------------------
    % initialize fun2No and No2fun
    Nnodal=mesh.Nnodes-sum(mesh.nodes.onBoundary);
    Nedge=(mesh.Nedges-sum(mesh.edges.onBoundary))*K;
    Nface=mesh.Nsurfaces*K*K;
    Nbasis=Nnodal+Nedge+Nface;
    No2fun.name=char(zeros(Nbasis,1));  % possible values: 'n' for nodal; 'e' for edge; 'f' for face
    No2fun.objid=zeros(Nbasis,1);   % to specify which node/edge/surface the basis is bound for.
    No2fun.subid=zeros(Nbasis,2);   % for nodal mode: 1, for edge mode: [1,2~K+1] or [2~K+1,1], for face mode: [2~K+1,2~K+1];
                                    % .subid(:,1): subid for x dirl .subid(:,2): subid for y dir
    fun2No.nodal=zeros(mesh.Nnodes,1);   % the No. for nodal mode bases
    fun2No.edge=zeros(K+1,mesh.Nedges);   % the No. for edge mode bases
    fun2No.face=zeros(K+1,K+1,mesh.Nsurfaces);  % the No. for face mode bases, fun2No.face(Xsubid,Ysubid,objid)
    
    % set fun2No and No2fun
    % nodal mode
    No2fun.name(1:Nnodal)='n';
    No2fun.objid(1:Nnodal)=find(mesh.nodes.onBoundary==0);
    No2fun.subid(1:Nnodal,1:2)=1;
    fun2No.nodal(No2fun.objid(1:Nnodal))=(1:Nnodal)';
    % edge mode
    No2fun.name(Nnodal+1:Nnodal+Nedge)='e';
    offset=Nnodal+1;
    for i=1:mesh.Nedges
        if mesh.edges.onBoundary(i)==1
            continue;
        end
        n=mesh.edges.n(1:2,i);
        if (abs( mesh.nodes.x(n(1)) - mesh.nodes.x(n(2)) )<=meshEps)
            % edge is along y axis
            for k=2:K+1
                No2fun.objid(offset)=i;
                No2fun.subid(offset,1)=1;
                No2fun.subid(offset,2)=k;
                fun2No.edge(k,i)=offset;
                offset=offset+1;
            end
        else
            % edge is along x axis
            for k=2:K+1
                No2fun.objid(offset)=i;
                No2fun.subid(offset,1)=k;
                No2fun.subid(offset,2)=1;
                fun2No.edge(k,i)=offset;
                offset=offset+1;
            end
        end
    end
    % face mode
    No2fun.name(offset:offset+Nface-1)='f';
    xid=repmat(2:K+1,K,1);xid=xid(:);
    yid=repmat(2:K+1,1,K)';
    for i=1:mesh.Nsurfaces
        No2fun.objid(offset:offset+K*K-1)=i;
        No2fun.subid(offset:offset+K*K-1,1)=xid;  % = 2 2 ... 2 3 3 ... 3 4 4 ... 4 5...
        No2fun.subid(offset:offset+K*K-1,2)=yid;  % = 2 3 4 ... 2 3 4 ... 2 3 4 ...
        fun2No.face(2:K+1,2:K+1,i)=reshape(offset:offset+K*K-1,K,K)';
        offset=offset+K*K;
    end
    
%--------------------- Read Mesh ----------------------------------------
    % extract node from mesh
    Nobj=mesh.Nnodes+mesh.Nedges+mesh.Nsurfaces;
    adjNodeList=zeros(4,2,Nobj);    % [For node] adjNode(:,1): the id of four corner nodes, adjNode(:,2): four adjacent nodes, same numbering with edges.                   
                                    % [For edge] adjNode(:,1): the id of four corner nodes, adjNode(1:2,2): the id of edge's end points. When the edge is 
                                    %                          horizontal, the numbering is the same as we used in mesh. As for the edge is vertical, we make 
                                    %                          a reflection about y=x
                                    % [For surface] adjNode(:,1): the id of four corner nodes
    adjEdgeList=zeros(6,2,Nobj);    % [For node/surface] adjEdge(:,1): the id of vertical(along y axis) edges, adjEdge(:,2): horizontal edges
                                    % [For edge] adjEdge(:,1): the id of two edges that paralell to the edge; adjEdge(:,2): four edges that perpendicular to it.
                                    % numbering:  1  4            1 
                                    %   __ __     __ __     __    __                4  2              1  
                                    % 1|__|__|3  |__|__|  1|__|2 |__|    __ __      __ __      __     __
                                    % 4|__|__|6  |__|__|  3|__|4 |__|  2|__|__|1   |__|__|   1|__|2  |__| 
                                    %             3   6           2                 3  1              2 
    adjFaceList=zeros(4,Nobj);      % 4~1 adjacent surface for the node/edge/surface, the same numbering as we used in mesh. 
    hxList=zeros(2,Nobj);hyList=zeros(2,Nobj);
    xList=zeros(Nobj,1);  % [For node] the (x,y) for that node; [For edge] the (x,y) for n(1) of that edge
    yList=zeros(Nobj,1);  % [For surface] the (x,y) for n(3) of that surface
    edgeDir=zeros(mesh.Nedges,1);  % the direction of an edge. 1: vertical; 2: horizontal; 0: boundary
    
    % for node
    for i=1:mesh.Nnodes
        if mesh.nodes.onBoundary(i)==1
            continue;
        end
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
    for i=1:mesh.Nedges
        if mesh.edges.onBoundary(i)==1
            continue;
        end
        Eid=mesh.Nnodes+i;
        face=mesh.edges.s(:,i);
        adjFaceList(1:2,Eid)=face;
        adjNodeList(1:2,2,Eid)=mesh.edges.n(:,i);
        if abs(diff(mesh.nodes.x(adjNodeList(1:2,2,Eid))))<meshEps
            % vertical edge
            edgeDir(i)=1;
            adjNodeList(1:4,1,Eid)=[mesh.surfaces.n([1;4], face(2));mesh.surfaces.n([3;2], face(1))];
            % para. edges
            adjEdgeList(1:2,1,Eid)=[mesh.surfaces.e(4,face(2)) ; mesh.surfaces.e(2,face(1))];
            % ortho. edges
            adjEdgeList(1:4,2,Eid)=[mesh.surfaces.e([3;1],face(2)) ; mesh.surfaces.e([3;1],face(1))];
            hxList(:,Eid)=[mesh.surfaces.hx(face(1));  mesh.surfaces.hx(face(2))];
            hyList(1,Eid)=mesh.surfaces.hy(face(1));
        else
            % horizontal edge
            edgeDir(i)=2;
            adjNodeList(1:4,1,Eid)=[mesh.surfaces.n([1;2], face(1));mesh.surfaces.n([3;4], face(2))];
            % para. edges
            adjEdgeList(1:2,1,Eid)=[mesh.surfaces.e(1,face(1)) ; mesh.surfaces.e(3,face(2))];
            % ortho. edges
            adjEdgeList(1:4,2,Eid)=[mesh.surfaces.e([2;4],face(1)) ; mesh.surfaces.e([2;4],face(2))];
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
    
%-------------------- The following do not read mesh any more except mesh.Nnodes and mesh.Nedges -----------------------
    % generate coeff vectors
    nVertical=sum(edgeDir(edgeDir==1));
    nHorizontal=sum(edgeDir(edgeDir==2))/2;
    Svec=zeros(2+33*Nnodal+13*Nedge+5*Nface,3);
    Cxvec=zeros(2+28*Nnodal+(14*nHorizontal+12*nVertical)*K+6*Nface,3);
    Cyvec=zeros(2+28*Nnodal+(12*nHorizontal+14*nVertical)*K+6*Nface,3);
    Mvec=zeros(2+49*Nnodal+21*Nedge+9*Nface,3);
    topS=1;topM=1;topCx=1;topCy=1;
    for iNo=1:Nbasis
        % consider the line for ( basis(iNo) , u )
        switch No2fun.name(iNo)
            case 'n'
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
            case 'e'
                Eid=No2fun.objid(iNo);
                subid=max(No2fun.subid(iNo,1:2));
                No_adjNode=[fun2No.nodal(adjNodeList(:,1,mesh.Nnodes+Eid)),[fun2No.nodal(adjNodeList(1:2,2,mesh.Nnodes+Eid));0;0]];
                adjEdge=adjEdgeList(1:4,:,mesh.Nnodes+Eid);
                adjFace=adjFaceList(1:2,mesh.Nnodes+Eid);
                if edgeDir(Eid)==1
                    % for vertical edge
                    hp=hyList(1,mesh.Nnodes+Eid);
                    ho=hxList(1:2,mesh.Nnodes+Eid);
                else                                   %        hp
                    % for horizontal edge              %        __
                    hp=hxList(1,mesh.Nnodes+Eid);      % ho(2) |__|
                    ho=hyList([1,2],mesh.Nnodes+Eid);  % ho(1) |__|
                end
                
                pos2No_node=reshape(No_adjNode([3;5;2;4;6;1]),3,2);   % use pos2No_node(pos_o+2,pos_p+1) to get nodeNo.
                pos2id_edgep=[adjEdge(2,1);Eid;adjEdge(1,1)];         % use pos2id_edgep(pos_o+2) to get edgeNo for para. edge
                pos2id_edgeo=[adjEdge(3,2),adjEdge(4,2);adjEdge(1,2),adjEdge(2,2)];  % use pos2id_edgeo(pos_o+2,pos_p+1) to get edgeNo for ortho. edge
                %------------------- enumerate bases on the orthogonal direction -------------------
                for  subid_o=1:min(3,K+1)
                    for pos_o=-1:1    % enumerate bases on the orthogonal direction
                        if subid_o>1 && pos_o==1
                            continue;
                        end
                        % calc phiphi, phidphi, dphidphi first
                        if pos_o==0 && subid_o==1
                            range=(1:2)';
                        elseif pos_o==-1
                            range=1;
                        else
                            range=2;
                        end
                        phiphi_o=phiphi(1,subid_o,pos_o,ho(range));
                        dphidphi_o=dphidphi(1,subid_o,pos_o,ho(range));
                        phidphi_o=phidphi(1,subid_o,pos_o);
                        %------------- enumerate bases on the parallel direction -------------
%                         if subid<4
%                             subid_p_List=1:subid+2;
%                         else
%                             subid_p_List=subid-2:subid+2;
%                         end
                        for subid_p=max(1,subid-2):min(subid+2,K+1)  
                            for pos_p=0:1    % enumerate bases on the parallel direction
                                if subid_p>1 && pos_p==1
                                    continue;
                                end
                                % calc jNo
                                if subid_o>1 && subid_p>1
                                    % face basis
                                    if edgeDir(Eid)==1  % if edge Eid is vertical 
                                        jNo=fun2No.face(subid_o,subid_p,adjFace(pos_o+2));
                                    else
                                        jNo=fun2No.face(subid_p,subid_o,adjFace(1-pos_o));
                                    end  
                                elseif subid_o==1 && subid_p==1
                                    % nodal basis
                                    jNo=pos2No_node(pos_o+2,pos_p+1);
                                else
                                    % edge basis
                                    if subid_o==1
                                        jNo=fun2No.face(subid_p,pos2id_edgep(pos_o+2));
                                    else
                                        jNo=fun2No.face(subid_o,pos2id_edgeo(pos_o+2,pos_p+1));
                                    end
                                end
                                phiphi_p=phiphi(subid,subid_p,pos_p,hp);
                                % set Mvec
                                Mvec(topM,1)=iNo;
                                Mvec(topM,2)=jNo;
                                Mvec(topM,3)=phiphi_o*phiphi_p;
                                if abs(Mvec(topM,3))>0; topM=topM+1; end
                                % set Svec
                                Svec(topS,1)=iNo;
                                Svec(topS,2)=jNo;
                                Svec(topS,3)=dphidphi_o*phiphi_p+phiphi_o*dphidphi(subid,subid_p,pos_p,hp);
                                if abs(Svec(topS,3))>0; topS=topS+1; end
                                % set Cxvec, Cyvec
                                Cxvec(topCx,1)=iNo;Cyvec(topCy,1)=iNo;
                                Cxvec(topCx,2)=jNo;Cyvec(topCy,2)=jNo;
                                if edgeDir(Eid)==1
                                    % for vertical edge
                                    Cxvec(topCx,3)=phidphi_o*phiphi_p;
                                    Cyvec(topCy,3)=phiphi_o*phidphi(subid,subid_p,pos_p);
                                else
                                    % for horizontal edge
                                    Cxvec(topCx,3)=phidphi(subid,subid_p,pos_p)*phiphi_o;
                                    Cyvec(topCy,3)=phiphi_p*phidphi_o;
                                end
                                if abs(Cxvec(topCx,3))>0; topCx=topCx+1; end
                                if abs(Cyvec(topCy,3))>0; topCy=topCy+1; end
                            end
                        end
                    end
                end
            case 'f'
                Sid=No2fun.objid(iNo);
                subid=No2fun.subid(iNo,:)';
                No_adjNode=fun2No.nodal(adjNodeList(:,1,mesh.Nnodes+mesh.Nedges+Sid));
                adjEdge=adjEdgeList(1:2,:,mesh.Nnodes+mesh.Nedges+Sid);
                hx=hxList(1,mesh.Nnodes+mesh.Nedges+Sid); hy=hyList(1,mesh.Nnodes+mesh.Nedges+Sid);
                
                pos2No_node=reshape(No_adjNode([3;4;2;1]),2,2);   % use pos2No_node(pos_x+1,pos_y+1) to get node No.
                %------------- enumerate bases on x dir -------------
                for subid_x=max(1,subid(1)-2):min(subid(1)+2,K+1)
                    for pos_x=0:1
                        if (subid_x>1 && pos_x==1)
                            continue;
                        end
                        % calc phiphi_x,phidphi_x, dphidphi_x first
                        phiphi_x=phiphi(subid(1),subid_x,pos_x,hx);
                        phidphi_x=phidphi(subid(1),subid_x,pos_x);
                        dphidphi_x=dphidphi(subid(1),subid_x,pos_x,hx);
                        %------------- enumerate bases on y dir -------------
                        for subid_y=max(1,subid(2)-2):min(subid(2)+2,K+1)
                            for pos_y=0:1
                                if (subid_y>1 && pos_y==1)
                                    continue;
                                end
                                % calc jNo
                                if subid_x>1 && subid_y>1
                                    % face basis
                                    jNo=fun2No.face(subid_x,subid_y,Sid);
                                elseif subid_x==1 && subid_y>1
                                    % vertical edge basis
                                    jNo=fun2No.edge(subid_y,adjEdge(pos_x+1,1));
                                elseif subid_x>1 && subid_y==1
                                    % horizontal edge basis
                                    jNo=fun2No.edge(subid_x,adjEdge(2-pos_y,2));
                                else
                                    % nodal basis
                                    jNo=pos2No_node(pos_x+1,pos_y+1);
                                end
                                phiphi_y=phiphi(subid(2),subid_y,pos_y,hy);
                                % set Mvec
                                Mvec(topM,1)=iNo;
                                Mvec(topM,2)=jNo;
                                Mvec(topM,3)=phiphi_x*phiphi_y;
                                if abs(Mvec(topM,3))>0; topM=topM+1; end
                                % set Svec
                                Svec(topS,1)=iNo;
                                Svec(topS,2)=jNo;
                                Svec(topS,3)=dphidphi_x*phiphi_y+phiphi_x*dphidphi(subid(2),subid_y,pos_y,hy);
                                if abs(Svec(topS,3))>0; topS=topS+1; end
                                % set Cxvec
                                Cxvec(topCx,1)=iNo;
                                Cxvec(topCx,2)=jNo;
                                Cxvec(topCx,3)=phidphi_x*phiphi_y;
                                if abs(Cxvec(topCx,3))>0; topCx=topCx+1; end
                                % set Cyvec
                                Cyvec(topCy,1)=iNo;
                                Cyvec(topCy,2)=jNo;
                                Cyvec(topCy,3)=phiphi_x*phidphi(subid(2),subid_y,pos_y);
                                if abs(Cyvec(topCy,3))>0; topCy=topCy+1; end
                            end
                        end
                    end
                end
                
            otherwise
                error('Unknow basis name!');
        end
    end
    % delete lines that have Xvec(:,2)=0, those are boundary elements.
    Mvec= Mvec (Mvec (1:topM-1,2)>0.5 ,1:3);
    Svec= Svec (Svec (1:topS-1,2)>0.5 ,1:3);
    Cxvec=Cxvec(Cxvec(1:topCx-1,2)>0.5 ,1:3);
    Cyvec=Cyvec(Cyvec(1:topCy-1,2)>0.5 ,1:3);
    
end


function result=phiphi(l,m,pos,h)
% calc the integral (phi_l,phi_m), pos is the position of phi_m relative to phi_l.
% pos = -1 or 0 or +1
    result=0;
    if l>1 && m>1
        if m==l
            result=2*h(1)*(1/(2*l+1)+1/(2*l-3))/(2*l-1)/(2*l-1);
        elseif m==l+2
            result=-2*h(1)/(2*l-1)/(2*l+3)/(2*l+1);
        elseif m==l-2
            result=-2*h(1)/(2*l-1)/(2*l-5)/(2*l-3);
        end
    elseif l==1 && m==1
        result=sum(h)/3/(1+abs(pos));
    elseif l==1 && m<4
        if m==2
            result=-h(1)/3;
        elseif m==3
            result=h(1)/15*(2*pos+1);
        end
    elseif l<4 && m==1
        if l==2
            result=-h(1)/3;
        elseif l==3
            result=-h(1)/15*(2*pos-1);
        end
    end
end

function result=phidphi(l,m,pos)
% calc the integral (phi_l,phi'_m), pos is the position of phi_m relative to phi_l.
% pos = -1 or 0 or +1
    result=0;
    if l>1 && m>1
        if m==l+1
            result=2/(2*l-1)/(2*l+1);
        elseif m==l-1
            result=-2/(2*l-1)/(2*l-3);
        end
    elseif l==1 && m==2
        result=-(2*pos+1)/3;
    elseif l==2 && m==1
        result=-(2*pos-1)/3;
    elseif l==1 && m==1
        result=pos/2;
    end
end

function result=dphidphi(l,m,pos,h)
% calc the integral (phi'_l,phi'_m), pos is the position of phi_m relative to phi_l.
% pos = -1 or 0 or +1
    result=0;
    if l>1 && m>1
        if m==l
            result=2/h(1)/(2*l-1);
        end
    elseif l==1 && m==1
        if pos==0
            result=sum(1./h);
        else
            result=-1/h(1);
        end
    end
end
