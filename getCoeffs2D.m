function [ S,Cx,Cy,M,vecf,id2fun,fun2id ] = getCoeffs2D( mesh,basis,f )
% Generate the coefficient matrix for a given mesh and basis (2D version)
%   f: function f(x,y) or f{1}(x)*f{2}(y). Set f=0 if don't need vecf
%   mesh: a mesh structure returned by makeMesh();
%   basis: can be one of 'Linear' and 'Lobatto'
%   
%   id2fun, fun2id: id means the index of a basis function in the coeff vector, fun means the basis function name or 
%                   node No. if this basis correspond to a spatial point.
%
%
switch basis
    case 'Linear'
        
        % get vectorized coeff matrix
        [ id2fun, fun2id, Ninner, hxList, hyList, xList, yList, Svec, Cxvec, Cyvec, Mvec ]=getCoeffs2D_Linear_mex(mesh);

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
                    Nid=i;  xi=xList(i);yi=yList(i);     hx=hxList(:,i);hy=hyList(:,i);

                    % integrate
                    for k=1:2      % hy(k)
                        for m=1:2   % hx(m)
                            vecf(i)=vecf(i)+hx(m)*hy(k)/16*integral( @(ksi)(1+sgns(m)*ksi).*f{1}(  hx(m)/2*(ksi-sgns(m))+xi  ) ,-1,1)...
                                                          *integral( @(eta)(1+sgns(k)*eta).*f{2}(  hy(k)/2*(eta-sgns(k))+yi  ) ,-1,1);
                        end
                    end
                end
            else
                % for the case that f(x,y) is not separable
                sgns=[1;-1;1];
                parfor i=1:Ninner
                    % get node info
                    Nid=i;    xi=xList(i);yi=yList(i);   hx=hxList(:,i);hy=hyList(:,i);

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
        
    otherwise
        error(['Unknow basis: ',basis]);
end

end

