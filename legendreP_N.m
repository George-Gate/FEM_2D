function [ result ] = legendreP_N( nList,xList )
%% Calculate the value of Legendre Polynomials at given position
%  [George-Gate @2017-04-10]
%  [Usage]
%       result=legendreP_N( n, x )
%   n should be a vector and n>=0
%       if x is a col vector,
%             result=[L( n(1) , x ),L( n(2) , x ),L( n(3) , x ), ... ,L( n(end) , x )]
%       if x is a row vector,
%             result=[L( n(1) , x );L( n(2) , x );L( n(3) , x ); ... ;L( n(end) , x )]
%       if x is an array with ndim>1
%             result is a cell array with the same size of n. result{i}=L(n(i),x);
%             For the case that n is a scalar, result=L(n,x);
%

    if ~isvector(nList)
        error('nList must be a vector.');
    end
    if min(nList)<0
        error('nList must be non-negative.');
    end
    maxN=max(nList);
    P=cell(maxN+1,1);
    P{1}=ones(size(xList));
    P{2}=xList;
    for n=1:maxN-1
        P{n+2}=( (2*n+1)*xList.*P{n+1}-n*P{n} )/(n+1);
    end
    % set output
    if (isrow(xList))
        result=zeros(length(nList),length(xList));
        for i=1:length(nList)
            result(i,:)=P{nList(i)+1};
        end
    elseif (iscolumn(xList))
        result=zeros(length(xList),length(nList));
        for i=1:length(nList)
            result(:,i)=P{nList(i)+1};
        end
    elseif isscalar(nList)
        result=P{nList+1};
    else
        result=cell(size(nList));
        for i=1:length(nList)
            result{i}=P{nList(i)+1};
        end
    end
end

