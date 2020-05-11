function tp=KroneckerProduct(varargin)
% KroneckerProduct(X,Y,Z,...) is the Kronecker tensor product. 
% For example, for 2 by 2 matrices X,Y,Z, KroneckerProduct(X,Y,Z) returns the 
% (2*2*2) by (2*2*2) matrix
%
%  [ x11 * [Y11*Z Y12*Z          x12 * [Y11*Z Y12*Z  
%           Y21*Z Y22*Z]                Y21*Z Y22*Z] 
%                                                   
%    x21 * [Y11*Z Y12*Z          X22 * [Y11*Z Y12*Z
%           Y21*Z Y22*Z]                Y21*Z Y22*Z] ]
%
% KroneckerProduct(X,n) with n a scalar returns the n-time recursive
% KroneckerProduct of X

if nargin < 2
    error('nargin<2 in KroneckerProduct')
elseif (nargin==2) && prod(size(varargin{2})==1)
    tp=1;
    for k=1:varargin{2}
        tp=kron(tp,varargin{1});
    end
else
    tp=varargin{1};
    for k=1:nargin-1
        tp=kron(tp,varargin{k+1});
    end
end

end