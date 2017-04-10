function [ result ] = lobatto_N( jList,xList )
%% Calculate the value of Lobatto Polynomials at given positions
%  [Usage]
%       result=lobatto_N( j, x )
%   j should be a vector and j>=1
%       if x is a col vector,
%             result=[Lo( j(1) , x ),Lo( j(2) , x ),Lo( j(3) , x ), ... ,Lo( j(end) , x )]
%       if x is a row vector,
%             result=[Lo( j(1) , x );Lo( j(2) , x );Lo( j(3) , x ); ... ;Lo( j(end) , x )]
%       if x is an array with ndim>1
%             result is a cell array with the same size of j. result{i}=Lo(j(i),x);
%             For the case that j is a scalar, result=Lo(j,x);
%   [Math]
%      lobatto_N( j, x ) = Lo_{j+1}(x) == ( L_{j+1}(x) - L_{j-1}(x) )/(2*j+1) == integral(L_j(t),-1,x)
%      where L_j(x) is Legendre Polynomials
%
    if ~isvector(jList)
        error('jList must be a vector.');
    end
    if min(jList)<1
        error('jList should larger than 0.');
    end
    maxN=max(jList)+1;
    % get Legendre Polynomials
    L=legendreP_N(0:maxN,xList);
    
    % set output
    if (isrow(xList))  % xList is row vector
        jList=reshape(jList,length(jList),1);
        result=( L(jList+2,:)-L(jList,:) )./repmat(2*jList+1,1,length(xList));
    elseif (iscolumn(xList))
        jList=reshape(jList,1,length(jList));
        result=( L(:,jList+2)-L(:,jList) )./repmat(2*jList+1,length(xList),1);
    elseif isscalar(jList)
        result=( L{jList+2}-L{jList} )/(2*jList+1);
    else % vector jList and array xList 
        result=cell(size(jList));
        for i=1:length(jList)
            result{i}=( L{jList(i)+2}-L{jList(i)} )/(2*jList(i)+1);
        end
    end
end

