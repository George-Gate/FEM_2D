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
        xList=coList{1};
        yList=coList{2};
        
        mesh=makeMesh_box_mex(xList,yList);
        
    case 'LshapeSegUniform'
        % Make mesh for [-1,1] x [0,1] and [-1,0] x [-1,0] by calling makeMesh_box
        % And then combine then into one mesh
        
    otherwise
        error(['Invalid mesh type: ',type]);
end

end

