function mesh = makeMesh( type,n,varargin )
% Generate mesh on a given area.
%   Minimum Input: type, n
%   Supported mesh type and area:
%     boxUniform - Uniform mesh on [0,1] x [0,1], hx can be different from hy. Use a 2x1 vector as n to set hx and hy
%                  seperately.
%     boxSegUniform - Segmented uniform mesh on [0,1] x [0,1]. Something very similar to shishkin mesh. 
%                     n: 2 x 1 cell where n{1} and n{2} is a vector, giving the point number on each interval.
%                     w: 2 x 1 cell where length(w{i})=length(n{i})-1, giving the width of each interval. 
%                        Input should ensure that sum(w{i})<1, the width of last interval is set to 1-sum(w{i}).
%           [Usage] mesh = makeMesh( 'boxSegUniform', n, w );
%
%    LshapeSegUniform - Segmented uniform mesh on an L shape domain compose of [-1,0] x [-1,1] and [0,1] x [0,1]. It can
%                       be viewed as the union of two boxSegUniform mesh.
%                     n: 4 x 1 cell where n{i} is a vector, giving the point number on each interval.
%                     w: 4 x 1 cell where length(w{i})=length(n{i})-1, giving the width of each interval. 
%                        Input should ensure that sum(w{i})<1, the width of last interval is set to 1-sum(w{i}).
%                  i=1,2 describe the x axis division of [-1,0] and [0,1]. 
%                  i=3,4 describe the y axis division of [-1,0] and [0,1].
%           [Usage] mesh = makeMesh( 'LshapeSegUniform', n, w );
%
%  [About the return struct mesh]
%  mesh
%     .Nnodes, .Nesges, .Nsurfaces: the number of nodes, edges and surfaces in the mesh
%     .nodes: node information
%          .x(i)             x coordinate of node i
%          .y(i)             y coordinate of node i
%          .s(1:4,i)         the id of adjacent surfaces of node i
%          .onBoundary(i)    to indicate whether node i is on the boundary
%     .edges: edge information
%          .n(1:2,i)         the id of endpoint nodes of edge i
%          .s(1:2,i)         the id of adjacent surfaces of edge i
%          .onBoundary(i)    to indicate whether edge i is on the boundary
%     .surfaces: surface information
%          .n(1:4,i)         the id of corner nodes of surface i
%          .e(1:4,i)         the id of adjacent edge of surface i
%          .x(1:2,i)         x coordinates of surface i
%          .y(1:2,i)         y coordinates of surface i
%          .hx(i)            width of surface i along x axis
%          .hy(i)            width of surface i along y axis
use_mex=1;

switch type
    case 'boxSegUniform'
        % check input args
        w=varargin{1};
        if length(n{1})-length(w{1})~=1 || length(n{2})-length(w{2})~=1
            error('Input does not satisfy length(w{i})=length(n{i})-1.');
        end
        if sum(w{1})>=1 || sum(w{2})>=1
            error('Input does not satisfy sum(w{i}<1)');
        end
        if min(min(w{1}),min(w{2})) <= 0
            error('w shoud be larger than 0.');
        end
        if min(min(n{1}),min(n{2}))<1
            error('n should be larger than 0.');
        end
        % x, y coordinates
        coList=cell(2,1);
        for i=1:2
            coList{i}=0;
            len=length(w{i});
            for j=1:len
                tmp=linspace(coList{i}(end),coList{i}(end)+w{i}(j),n{i}(j)+1)';
                coList{i}=[coList{i};tmp(2:end)];
            end
            tmp=linspace(coList{i}(end),1,n{i}(len+1)+1)';
            coList{i}=[coList{i};tmp(2:end)];
        end
        clear tmp len;
        % make mesh
        xList=coList{1};
        yList=coList{2};
        if use_mex
            mesh=makeMesh_box_mex(xList,yList);
        else
            mesh=makeMesh_box(xList,yList);
        end
        
    case 'LshapeSegUniform'
        % Make mesh for [-1,1] x [0,1] and [-1,0] x [-1,0] by calling makeMesh_box
        % And then combine then into one mesh
        % check input args
        w=varargin{1};
        if length(n{1})-length(w{1})~=1 || length(n{2})-length(w{2})~=1 ...
                                        || length(n{3})-length(w{3})~=1 || length(n{4})-length(w{4})~=1
            error('Input does not satisfy length(w{i})=length(n{i})-1.');
        end
        if sum(w{1})>=1 || sum(w{2})>=1 || sum(w{2})>=1 || sum(w{3})>=1
            error('Input does not satisfy sum(w{i}<1)');
        end
        if min( [min(w{1}),min(w{2}),min(w{3}),min(w{4})] ) <= 0
            error('w shoud be larger than 0.');
        end
        if min([min(n{1}),min(n{2}),min(n{3}),min(n{4})])<1
            error('n should be larger than 0.');
        end
        
        % x, y coordinates
        coList=cell(4,1);
        for i=1:4
            coList{i}=0;
            len=length(w{i});
            for j=1:len
                tmp=linspace(coList{i}(end),coList{i}(end)+w{i}(j),n{i}(j)+1)';
                coList{i}=[coList{i};tmp(2:end)];
            end
            tmp=linspace(coList{i}(end),1,n{i}(len+1)+1)';
            coList{i}=[coList{i};tmp(2:end)];
        end
        clear tmp len;
        
        % make mesh
        zoom=1;
        xList1=zoom*(coList{1}-1);
        xList2=zoom*coList{2};
        yList1=zoom*(coList{3}-1);
        yList2=zoom*coList{4};
        
        if use_mex
            mesh0=makeMesh_box_mex([xList1(1:end-1);xList2],yList2);
            mesh1=makeMesh_box_mex(xList1,yList1);
        else
            mesh0=makeMesh_box([xList1(1:end-1);xList2],yList2);
            mesh1=makeMesh_box(xList1,yList1);
        end
        
        %------------ combine mesh ------------------------
        %********* delete redundant nodes and edges ********
        % delete node 1~length(xList1) in mesh0
        % delete edge 1~length(xList1)-1 in mesh0
        len=length(xList1);
        mesh0.Nnodes=mesh0.Nnodes-len;
        mesh0.Nedges=mesh0.Nedges-(len-1);
        mesh0.nodes.x=trimMatrix(mesh0.nodes.x,len,1);
        mesh0.nodes.y=trimMatrix(mesh0.nodes.y,len,1);
        mesh0.nodes.s=trimMatrix(mesh0.nodes.s,len,2);
        mesh0.nodes.onBoundary=trimMatrix(mesh0.nodes.onBoundary,len,1);
        mesh0.edges.n=trimMatrix(mesh0.edges.n,len-1,2);
        mesh0.edges.s=trimMatrix(mesh0.edges.s,len-1,2);
        mesh0.edges.onBoundary=trimMatrix(mesh0.edges.onBoundary,len-1,1);
        % shift all node No. and edge No. in mesh0
        mesh0.edges.n(mesh0.edges.n>0)=mesh0.edges.n(mesh0.edges.n>0)-len;
        mesh0.surfaces.n(mesh0.surfaces.n>0)=mesh0.surfaces.n(mesh0.surfaces.n>0)-len;
        mesh0.surfaces.e(mesh0.surfaces.e>0)=mesh0.surfaces.e(mesh0.surfaces.e>0)-(len-1);
               
        %******** shift all node, edge, surface No. in mesh1 to match with mesh0 *********
        mesh1.nodes.s(mesh1.nodes.s>0)=mesh1.nodes.s(mesh1.nodes.s>0)+mesh0.Nsurfaces;
        mesh1.edges.n(mesh1.edges.n>0)=mesh1.edges.n(mesh1.edges.n>0)+mesh0.Nnodes;
        mesh1.edges.s(mesh1.edges.s>0)=mesh1.edges.s(mesh1.edges.s>0)+mesh0.Nsurfaces;
        mesh1.surfaces.n(mesh1.surfaces.n>0)=mesh1.surfaces.n(mesh1.surfaces.n>0)+mesh0.Nnodes;
        mesh1.surfaces.e(mesh1.surfaces.e>0)=mesh1.surfaces.e(mesh1.surfaces.e>0)+mesh0.Nedges;
        
        %******** link mesh ***********
        Nnodes=mesh0.Nnodes+mesh1.Nnodes;
        Nedges=mesh0.Nedges+mesh1.Nedges;
        Nsurfaces=mesh0.Nsurfaces+mesh1.Nsurfaces;
        % fix the relative relation on the linking edge
        % mesh0
        mesh0.surfaces.n(3,1:len)=Nnodes-len+1:Nnodes;
        mesh0.surfaces.n(4,1:len-1)=Nnodes-len+2:Nnodes;
        Eid=mesh1.surfaces.e(1,end-len+2);
        mesh0.surfaces.e(3,1:len-1)=Eid:Eid+len-2;
        Eid=mesh0.surfaces.e(2,1);
        mesh0.edges.n(1,Eid:Eid+len-1)=Nnodes-len+1:Nnodes;
        mesh0.edges.n(1,1)=Nnodes;
        % mesh1 - nodes
        mesh1.nodes.s(1,end-len+1:end)=1:len;
        mesh1.nodes.s(2,end-len+2:end)=1:len-1;
        mesh1.nodes.onBoundary(end-len+2:end-1)=0;
        % mesh1 - edges
        Eid=mesh1.surfaces.e(1,end-len+2)-mesh0.Nedges;
        mesh1.edges.s(1,Eid:Eid+len-2)=1:len-1;
        mesh1.edges.onBoundary(Eid:Eid+len-2)=0;

        % link! The order of assignment matters if you want to use mex function.
        mesh.Nnodes=Nnodes; mesh.Nedges=Nedges; mesh.Nsurfaces=Nsurfaces;
        mesh.nodes.x=[mesh0.nodes.x;mesh1.nodes.x];
        mesh.nodes.y=[mesh0.nodes.y;mesh1.nodes.y];
        mesh.nodes.s=[mesh0.nodes.s,mesh1.nodes.s];
        mesh.nodes.onBoundary=[mesh0.nodes.onBoundary;mesh1.nodes.onBoundary];
        
        
        mesh.edges.n=[mesh0.edges.n,mesh1.edges.n];
        mesh.edges.s=[mesh0.edges.s,mesh1.edges.s];
        mesh.edges.onBoundary=[mesh0.edges.onBoundary;mesh1.edges.onBoundary];
        
        mesh.surfaces.n=[mesh0.surfaces.n,mesh1.surfaces.n];
        mesh.surfaces.e=[mesh0.surfaces.e,mesh1.surfaces.e];
        mesh.surfaces.x=[mesh0.surfaces.x,mesh1.surfaces.x];
        mesh.surfaces.y=[mesh0.surfaces.y,mesh1.surfaces.y];
        mesh.surfaces.hx=[mesh0.surfaces.hx;mesh1.surfaces.hx];
        mesh.surfaces.hy=[mesh0.surfaces.hy;mesh1.surfaces.hy];
        
    otherwise
        error(['Invalid mesh type: ',type]);
end

end

