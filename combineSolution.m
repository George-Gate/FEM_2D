function [sol,xSample,ySample] = combineSolution( u, sampleRate, mesh, fun2No, basis, K )
%Combine the solution and calc it on sampling points
%  [Usage]
%      [sol,xSample,ySample] = combineSolution( u, sampleRate, mesh, fun2No, basis, K )
%  [Input]
%        u: the solution vector
%  sampleRate: the number of sampling points per surface per axis
%    basis: which basis is used to solve the PDE
%        K: the cut off of Lobatto Basis
%  mesh, fun2No: just pass those parameters
%
%  [Output]
%    sol: An (Nsurfaces x 1) cell array. Each element is a (sampleRate x sampleRate) matrix containing 
%         the value of solution in each surface.
    switch basis
        case 'Lobatto'
            % calc base value
            [xSample,ySample]=meshgrid(linspace(-1,1,sampleRate),linspace(-1,1,sampleRate));
            xList=xSample(1,:);
            yList=ySample(:,1);
            xAxisBase=zeros(K+2,length(xList));
            xAxisBase(1,:)=(1-xList)/2;  % base @ x=-1
            xAxisBase(2,:)=(xList+1)/2;  % base @ x=1
            xAxisBase(3:K+2,:)=lobattoP_N(1:K,xList); % Lobatto base
            yAxisBase=zeros(length(yList),K+2);
            yAxisBase(:,1)=(1-yList)/2;  % base @ y=-1
            yAxisBase(:,2)=(yList+1)/2;  % base @ y=1
            yAxisBase(:,3:K+2)=lobattoP_N(1:K,yList); % Lobatto base
            baseValue=cell(K+2);
            for jx=0:K+1
                for jy=0:K+1
                    baseValue{jx+1,jy+1}=yAxisBase(:,jy+1)*xAxisBase(jx+1,:);
                end
            end
            % combine solution
            j2Nid=[3,2;4,1];    % use j2Nid(jx+1,jy+1) to get relative position of nodal basis
            
            sol=cell(mesh.Nsurfaces,1);
            for Sid=1:mesh.Nsurfaces
                summer=zeros(sampleRate);
                % read mesh
                edges=mesh.surfaces.e(:,Sid);
                nodes=mesh.surfaces.n(:,Sid);
                % sum basis
                for jx=0:K+1
                    for jy=0:K+1
                        if jx<2 && jy<2  % nodal basis
                            iNo=fun2No.nodal(nodes(j2Nid(jx+1,jy+1)));
                        elseif jx<2  % vertical edge basis
                            iNo=fun2No.edge(jy,edges(2*(jx+1)));
                        elseif jy<2  % horizontal edge basis
                            iNo=fun2No.edge(jx,edges(3-2*jy));
                        else  % face basis
                            iNo=fun2No.face(jx,jy,Sid);
                        end
                        if iNo>0 % ignore boundary basis
                            summer=summer+u(iNo)*baseValue{jx+1,jy+1};
                        end
                    end
                end
                sol{Sid}=summer;
            end
            
        otherwise
            error(['Unknow basis: ',basis]);
    end
end

