function [ S,Cx,Cy,M,vecf,No2fun,fun2No ] = getCoeffs2D( mesh,basis,f,varargin )
% Generate the coefficient matrix for a given mesh and basis (2D version)
%   f: function f(x,y) or f{1}(x)*f{2}(y). Set f=0 if don't need vecf
%   mesh: a mesh structure returned by makeMesh();
%   basis: can be one of 'BiLinear' and 'Lobatto'
%   
%   No2fun, fun2No: No means the index of a basis function in the coeff vector, fun means the basis function name or 
%                   node id. if this basis correspond to a spatial point.
%
%  'fOrderUpperBound': a 2x1 or 1x2 vector, given an upper bound for the f's order of x and y. This option can acclerate the calculation of vecf and is only 
%                      avaliable whan basis='Lobatto'. The default value is [0, 0] for constant f and [inf, inf] otherwise.
%
% [Usage]
%      [ S,Cx,Cy,M,vecf,No2fun,fun2No ] = getCoeffs2D( mesh,'BiLinear',f )
%      [ S,Cx,Cy,M,vecf,No2fun,fun2No ] = getCoeffs2D( mesh,'Lobatto',f,K )   , K>=2
%
use_mex=1;

switch basis
%============================= Linear basis =============================================================
    case 'BiLinear'
        % get vectorized coeff matrix
        if use_mex
            [ No2fun, fun2No, Ninner, hxList, hyList, xList, yList, Svec, Cxvec, Cyvec, Mvec ]=getCoeffs2D_Linear_mex(mesh);
        else
            [ No2fun, fun2No, Ninner, hxList, hyList, xList, yList, Svec, Cxvec, Cyvec, Mvec ]=getCoeffs2D_Linear(mesh);
        end

        % generate sparse matrix
        S=sparse(Svec(:,1),Svec(:,2),Svec(:,3),Ninner,Ninner);
        M=sparse(Mvec(:,1),Mvec(:,2),Mvec(:,3),Ninner,Ninner);
        Cx=sparse(Cxvec(:,1),Cxvec(:,2),Cxvec(:,3),Ninner,Ninner);
        Cy=sparse(Cyvec(:,1),Cyvec(:,2),Cyvec(:,3),Ninner,Ninner);

        % check symmetry
        if (max(max(abs(S-S')))>eps) || (max(max(abs(M-M')))>eps) || (max(max(abs(Cx+Cx')))>eps) || (max(max(abs(Cy+Cy')))>eps)
            error('Symmetry test of S, Cx, Cy and M failed.');
        end

        % calc vecf=(f,phi);
        vecf=zeros(Ninner,1);
        % if f=0, do not calc vecf
        if nargin==3
            if iscell(f)
                % for the case that f(x,y) is separable
                sgns=[1;-1;1];
                parfor i=1:Ninner
                    % get node info
                    xi=xList(i);yi=yList(i);     hx=hxList(:,i);hy=hyList(:,i);

                    % integrate
                    for k=1:2      % hy(k)
                        for m=1:2   % hx(m)
                            vecf(i)=vecf(i)+hx(m)*hy(k)/16*integral( @(ksi)(1+sgns(m)*ksi).*f{1}(  hx(m)/2*(ksi-sgns(m))+xi  ) ,-1,1)...
                                                          *integral( @(eta)(1+sgns(k)*eta).*f{2}(  hy(k)/2*(eta-sgns(k))+yi  ) ,-1,1);
                        end
                    end
                end
            elseif isnumeric(f)
                % for constant f
                vecf=f*(hxList(1,:).*hyList(1,:)+hxList(1,:).*hyList(2,:)+hxList(2,:).*hyList(1,:)+hxList(2,:).*hyList(2,:))'/4;
            else
                % for the case that f(x,y) is not separable
                sgns=[1;-1;1];
                parfor i=1:Ninner
                    % get node info
                    xi=xList(i);yi=yList(i);   hx=hxList(:,i);hy=hyList(:,i);

                    % integrate
                    for k=1:2      % hy(k)
                        for m=1:2   % hx(m)
                            vecf(i)=vecf(i)+hx(m)*hy(k)/16*integral2( @(ksi,eta)(1+sgns(m)*ksi).*(1+sgns(k)*eta)...
                                                                             .*f(  hx(m)/2*(ksi-sgns(m))+xi , hy(k)/2*(eta-sgns(k))+yi ) ,-1,1,-1,1);
                        end
                    end
                end
            end
        end
%============================= Lobatto basis =============================================================
    case 'Lobatto'
        narginchk(4, 6)
        % cut off for Lobatto basis
        K=varargin{1};
        % default arguments
        if isnumeric(f)
            fUb=[0;0];
        else
            fUb=[inf;inf];
        end
        % pharse varargin
        for i=2:2:length(varargin)
            switch varargin{i}
                case 'fOrderUpperBound'
                    if length(varargin{i+1})==2
                        fUb=reshape(varargin{i+1},2,1);
                    else
                        warning('fOrderUpperBound should be a vector with len=2.');
                    end
                otherwise
                    error(['Unknow option: ',varargin{i}]);
            end
        end
        
        % get vectorized coeff matrix
        if use_mex
            [ No2fun, fun2No, Nbasis, hxList, hyList, xList, yList, edgeDir, Svec, Cxvec, Cyvec, Mvec ]=getCoeffs2D_Lobatto_mex(mesh,K);
        else
            [ No2fun, fun2No, Nbasis, hxList, hyList, xList, yList, edgeDir, Svec, Cxvec, Cyvec, Mvec ]=getCoeffs2D_Lobatto(mesh,K);
        end
        
        % generate sparse matrix
        S=sparse(Svec(:,1),Svec(:,2),Svec(:,3),Nbasis,Nbasis);
        M=sparse(Mvec(:,1),Mvec(:,2),Mvec(:,3),Nbasis,Nbasis);
        Cx=sparse(Cxvec(:,1),Cxvec(:,2),Cxvec(:,3),Nbasis,Nbasis);
        Cy=sparse(Cyvec(:,1),Cyvec(:,2),Cyvec(:,3),Nbasis,Nbasis);        
        
        % check symmetry
        if (max(max(abs(S-S')))>eps) || (max(max(abs(M-M')))>eps) || (max(max(abs(Cx+Cx')))>eps) || (max(max(abs(Cy+Cy')))>eps)
            error('Symmetry test of S, Cx, Cy and M failed.');
        end
        
        % calc vecf=(f,phi);
        vecf=zeros(Nbasis,1);
        
        
        if iscell(f)
        % for the case that f(x,y) is separable
            error('This function is not finished so far.');
            % for nodal basis
            nodalList=find(No2fun.name=='n');
            sgns=[1;-1;1];
            tmpvecf=zeros(length(nodalList),1);
            parfor i=1:length(nodalList)
                iNo=nodalList(i);
                Nid=No2fun.objid(iNo);
                % get node info
                xi=xList(Nid);yi=yList(Nid);     hx=hxList(:,Nid);hy=hyList(:,Nid);

                % integrate
                for k=1:2      % hy(k)
                    for m=1:2   % hx(m)
                        tmpvecf(i)=tmpvecf(i)+hx(m)*hy(k)/16*integral( @(ksi)(1+sgns(m)*ksi).*f{1}(  hx(m)/2*(ksi-sgns(m))+xi  ) ,-1,1)...
                                                            *integral( @(eta)(1+sgns(k)*eta).*f{2}(  hy(k)/2*(eta-sgns(k))+yi  ) ,-1,1);
                    end
                end
            end
            vecf(nodalList)=tmpvecf;
        elseif isnumeric(f)
        % for constant f
            error('This function is not finished so far.');
            if f
                % for nodal basis
                nodalList=find(No2fun.name=='n');
                Nid=No2fun.objid(nodalList);
                vecf(nodalList)=f*(hxList(1,Nid).*hyList(1,Nid)+hxList(1,Nid).*hyList(2,Nid)+hxList(2,Nid).*hyList(1,Nid)+hxList(2,Nid).*hyList(2,Nid))'/4;
            end
        else
        % for the case that f(x,y) is not separable
            integralList=repmat(struct,Nbasis,1);  % .fun; .ndim; .ub; .lb; .globalFactor; 
            basLst=zeros(Nbasis,1);    % execute vecf(basLst(1:topBL-1))=sum(integralValue(sumLst(1:topBL-1,:)),2 ) can set vecf
            sumLst=zeros(Nbasis,4);    
            topIL=1;topBL=1;
            %-------------- Generate Integral List ----------------------------
            % zero integral
            integralList(topIL).ndim=0;
            integralList(topIL).fun=0;
            integralList(topIL).globalFactor=0;
            integralList(topIL).ub=0;
            integralList(topIL).lb=0;
            topIL=topIL+1;
            
            % for nodal basis
            sgns=[1;-1;1];
            for iNo=find(No2fun.name=='n')'
                Nid=No2fun.objid(iNo);
                % get node info
                xi=xList(Nid);yi=yList(Nid);     hx=hxList(:,Nid);hy=hyList(:,Nid);
                % set integral
                basLst(topBL)=iNo;
                sumLst(topBL,1:4)=topIL:topIL+3;
                topBL=topBL+1;
                for k=1:2      % hy(k)
                    for m=1:2   % hx(m)
                        integralList(topIL).ndim=2;
                        integralList(topIL).fun=@(ksi,eta)(1+sgns(m)*ksi).*(1+sgns(k)*eta).*f(  hx(m)/2*(ksi-sgns(m))+xi , hy(k)/2*(eta-sgns(k))+yi );
                        integralList(topIL).globalFactor=hx(m)*hy(k)/16;
                        integralList(topIL).ub=[ 1; 1];
                        integralList(topIL).lb=[-1;-1];
                        topIL=topIL+1;
                    end
                end
            end
            
            % for edge basis
            sgns=[1;-1;1];
            for Eid=find(mesh.edges.onBoundary==0)'
                % get edge info
                xi=xList(Eid+mesh.Nnodes);yi=yList(Eid+mesh.Nnodes);     hx=hxList(:,Eid+mesh.Nnodes);hy=hyList(:,Eid+mesh.Nnodes);       
                if edgeDir(Eid)==1
                    % vertical edge
                    maxSubid=min(fUb(2)+2,K+1);
                else
                    maxSubid=min(fUb(1)+2,K+1);
                end
                % set integral
                for subid=2:maxSubid
                    iNo=fun2No.edge(subid,Eid);
                    basLst(topBL)=iNo;
                    sumLst(topBL,1:4)=[topIL:topIL+1,1,1];
                    topBL=topBL+1;
                    for m=1:2
                        integralList(topIL).ndim=2;
                        if edgeDir(Eid)==1  % vertical edge
                            integralList(topIL).fun=@(ksi,eta)(1+sgns(m)*ksi).*lobattoP_N(subid-1,eta).*f(  hx(m)/2*(ksi-sgns(m))+xi , hy(1)/2*(eta+1)+yi );
                            integralList(topIL).globalFactor=hx(m)*hy(1)/8;
                        elseif edgeDir(Eid)==2  % horizontal edge
                            integralList(topIL).fun=@(ksi,eta)lobattoP_N(subid-1,ksi).*(1+sgns(m)*eta).*f(  hx(1)/2*(ksi+1)+xi , hy(m)/2*(eta-sgns(m))+yi );
                            integralList(topIL).globalFactor=hx(1)*hy(m)/8;
                        end
                        integralList(topIL).ub=[ 1; 1];
                        integralList(topIL).lb=[-1;-1];
                        topIL=topIL+1;
                    end
                end
            end
            
            % for face basis
            for Sid=1:mesh.Nsurfaces
                % get edge info
                xi=xList(Sid+mesh.Nnodes+mesh.Nedges);yi=yList(Sid+mesh.Nnodes+mesh.Nedges);     
                hx=hxList(:,Sid+mesh.Nnodes+mesh.Nedges);hy=hyList(:,Sid+mesh.Nnodes+mesh.Nedges);       
                
                % set integral
                for subid_x=2:min(fUb(1)+2,K+1)
                    for subid_y=2:min(fUb(2)+2,K+1)
                        iNo=fun2No.face(subid_x,subid_y,Sid);
                        basLst(topBL)=iNo;
                        sumLst(topBL,1:4)=[topIL,1,1,1];
                        topBL=topBL+1;

                        integralList(topIL).ndim=2;
                        integralList(topIL).fun=@(ksi,eta)lobattoP_N(subid_x-1,ksi).*lobattoP_N(subid_y-1,eta).*f(  hx(1)/2*(ksi+1)+xi , hy(1)/2*(eta+1)+yi );
                        integralList(topIL).globalFactor=hx(1)*hy(1)/4;
                        integralList(topIL).ub=[ 1; 1];
                        integralList(topIL).lb=[-1;-1];
                        topIL=topIL+1;
                    end
                end
            end
            % -------------- Calc Integrals -----------------------
            integralValue=parIntegral(integralList,topIL);
            
            % -------------- Set vecf Value -----------------------
            vecf(basLst(1:topBL-1))=sum(integralValue(sumLst(1:topBL-1,:)'),1 );
            
        end
    otherwise
        error(['Unknow basis: ',basis]);
end

end

function value=parIntegral(list,topIL)
% Calc integrals in parallel
    value=zeros(topIL-1,1);
    parfor i=1:topIL-1
        ub=list(i).ub;
        lb=list(i).lb;
        if list(i).ndim==2
            value(i)=integral2(list(i).fun,lb(1),ub(1),lb(2),ub(2));
        elseif list(i).ndim==1
            value(i)=integral(list(i).fun,lb,ub);
        end
        value(i)=value(i)*list(i).globalFactor;
    end
end
