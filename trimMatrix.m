function [ matOut ] = trimMatrix( matIn,len,dim,dir )
% Trim a dim of 5D array.
% Examples:
%   trimMatrix( matIn,len ) = matIn(len+1:end,:,:,:,:)
%   trimMatrix( matIn,len,2 ) = matIn(:,len+1:end,:,:,:)
%   trimMatrix( matIn,len,5 ) = matIn(:,:,:,:,len+1:end)
%   trimMatrix( matIn,len,5,'first') = matIn(:,:,:,:,len+1:end)
%   trimMatrix( matIn,len,5,'last') = matIn(:,:,:,:,1:end-len)
%   trimMatrix( matIn,len,6) = matIn
%
if nargin==2
    dim=1;
    dir='first';
elseif nargin==3
    dir='first';
end

if ndims(matIn)>5
    error('Do not support array with dim higher than 5.');
end

switch dir
    case 'first'
        switch dim
            case 1
                matOut=matIn(len+1:end,:,:,:,:);
            case 2
                matOut=matIn(:,len+1:end,:,:,:);
            case 3
                matOut=matIn(:,:,len+1:end,:,:);
            case 4
                matOut=matIn(:,:,:,len+1:end,:);
            case 5
                matOut=matIn(:,:,:,:,len+1:end);
            otherwise
                matOut=matIn;
        end
    case 'last'
        switch dim
            case 1
                matOut=matIn(1:end-len,:,:,:,:);
            case 2
                matOut=matIn(:,1:end-len,:,:,:);
            case 3
                matOut=matIn(:,:,1:end-len,:,:);
            case 4
                matOut=matIn(:,:,:,1:end-len,:);
            case 5
                matOut=matIn(:,:,:,:,1:end-len);
            otherwise
                matOut=matIn;
        end
    otherwise
        error('Unknow trim direction.');
end
end

