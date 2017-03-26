function [ id2fun, fun2id, Ninner, hxList, hyList, xList, yList, Svec, Cxvec, Cyvec, Mvec ] = getCoeffs2D_Linear( mesh )
%Generate the coeff matrices for linear basis
%   Called by function getCoeffs2D()
    N=mesh.Nnodes;
    % generate id2fun and fun2id
    % id means the index of a basis function in the coeff vector, fun means the basis function name or 
    % node No. if this basis correspond to a spatial point.
    id2fun=find(mesh.nodes.onBoundary == 0);
    Ninner=length(id2fun);
    fun2id=zeros(N,1);
    fun2id(id2fun)=(1:Ninner)';
    
    % extract node info first
    adjNodeList=zeros(4,2,Ninner);  % adjNode(:,1): four corner nodes, adjNode(:,2): four adjacent nodes, same numbering with edges
    hxList=zeros(2,Ninner);hyList=zeros(2,Ninner);
    xList=mesh.nodes.x(id2fun);
    yList=mesh.nodes.y(id2fun);
    for i=1:Ninner
        Nid=id2fun(i);
        for j=1:4
            adjNodeList(j,1,i)=mesh.surfaces.n(j         , mesh.nodes.s(j,Nid));
            adjNodeList(j,2,i)=mesh.surfaces.n(mod(j,4)+1, mesh.nodes.s(j,Nid));
        end
        hxList(:,i)=[mesh.surfaces.hx(mesh.nodes.s(3,Nid));  mesh.surfaces.hx(mesh.nodes.s(1,Nid))];
        hyList(:,i)=[mesh.surfaces.hy(mesh.nodes.s(3,Nid));  mesh.surfaces.hy(mesh.nodes.s(1,Nid))];
    end
    
    
%-------------------- The following do not read mesh any more -----------------------
    % generate coeff matrices
    Svec=zeros(9*Ninner,3);
    Cxvec=zeros(6*Ninner,3);
    Cyvec=zeros(6*Ninner,3);
    Mvec=zeros(9*Ninner,3);
    for i=1:Ninner
        % set node info
        Nid=i;
        adjNode=adjNodeList(:,:,i);      % adjNode(:,1): four corner nodes, adjNode(:,2): four adjacent nodes, same numbering with edges
        % convert adjNode to id, which is the row index of a node in coeff vector.
        adjNode=fun2id(adjNode);

        hx=hxList(:,i);hy=hyList(:,i);


        % for S matrix
        top=9*(i-1);
        Svec(top+1:top+9,1)=Nid;         % row number
        % col number, central node, corner nodes and then adjacent nodes
        Svec(top+1:top+9,2)=[Nid; adjNode(:,1);  adjNode(:,2) ];
        % values
        Svec(top+1,3)=(sum(hx)/hy(1) + sum(hx)/hy(2) + sum(hy)/hx(1) + sum(hy)/hx(2))/3;
        Svec(top+2,3)=-(hx(2)/hy(2)+hy(2)/hx(2))/6;  % corner 1
        Svec(top+3,3)=-(hx(1)/hy(2)+hy(2)/hx(1))/6;
        Svec(top+4,3)=-(hx(1)/hy(1)+hy(1)/hx(1))/6;
        Svec(top+5,3)=-(hx(2)/hy(1)+hy(1)/hx(2))/6;
        Svec(top+6,3)=-(sum(hx)/hy(2))/3+(hy(2)/hx(1)+hy(2)/hx(2))/6;  % adjacent 1
        Svec(top+7,3)=-(sum(hy)/hx(1))/3+(hx(1)/hy(1)+hx(1)/hy(2))/6;
        Svec(top+8,3)=-(sum(hx)/hy(1))/3+(hy(1)/hx(1)+hy(1)/hx(2))/6;  % adjacent 3
        Svec(top+9,3)=-(sum(hy)/hx(2))/3+(hx(2)/hy(1)+hx(2)/hy(2))/6;

        % for M matrix
        top=9*(i-1);
        Mvec(top+1:top+9,1)=Nid;         % row number
        % col number, central node, corner nodes and then adjacent nodes
        Mvec(top+1:top+9,2)=[Nid; adjNode(:,1);  adjNode(:,2) ];
        % values
        Mvec(top+1,3)=sum(hx)*sum(hy)/9;
        Mvec(top+2,3)=hx(2)*hy(2)/36;  % corner 1
        Mvec(top+3,3)=hx(1)*hy(2)/36; 
        Mvec(top+4,3)=hx(1)*hy(1)/36; 
        Mvec(top+5,3)=hx(2)*hy(1)/36;
        Mvec(top+6,3)=sum(hx)*hy(2)/18;  % adjacent 1
        Mvec(top+7,3)=hx(1)*sum(hy)/18; 
        Mvec(top+8,3)=sum(hx)*hy(1)/18;  % adjacent 3
        Mvec(top+9,3)=hx(2)*sum(hy)/18; 

        % for Cx matrix
        top=6*(i-1);
        Cxvec(top+1:top+6,1)=Nid;         % row number
        % col number, corner nodes, left and right.
        Cxvec(top+1:top+6,2)=[adjNode(:,1);  adjNode([2;4],2) ];
        % values
        Cxvec(top+1,3)=+hy(2)/12;   % corner 1
        Cxvec(top+2,3)=-hy(2)/12;
        Cxvec(top+3,3)=-hy(1)/12;
        Cxvec(top+4,3)=+hy(1)/12;
        Cxvec(top+5,3)=-sum(hy)/6; % left, adjacent 2
        Cxvec(top+6,3)=+sum(hy)/6; % right, adjacent 4

        % for Cy matrix
        top=6*(i-1);
        Cyvec(top+1:top+6,1)=Nid;         % row number
        % col number, corner nodes, up and down.
        Cyvec(top+1:top+6,2)=[adjNode(:,1);  adjNode([1;3],2) ];
        % values
        Cyvec(top+1,3)=+hx(2)/12;   % corner 1
        Cyvec(top+2,3)=+hx(1)/12;
        Cyvec(top+3,3)=-hx(1)/12;
        Cyvec(top+4,3)=-hx(2)/12;
        Cyvec(top+5,3)=+sum(hx)/6; % up, adjacent 1
        Cyvec(top+6,3)=-sum(hx)/6; % down, adjacent 3
    end

    % Delete data with col=0. Those data correspond to boundary nodes, which should not be count.
    Svec=Svec(Svec(:,2)>0,:);
    Mvec=Mvec(Mvec(:,2)>0,:);
    Cxvec=Cxvec(Cxvec(:,2)>0,:);
    Cyvec=Cyvec(Cyvec(:,2)>0,:);

    
end

