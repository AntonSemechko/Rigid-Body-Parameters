function [isManifold, Enb, Eb, EA, Ne] = TriMeshIsManifold(TR)
% Check if triangular surface mesh is manifold (i.e., every non-boundary 
% edge is shared by two faces).
%
% INPUT
%   - TR  : triangular surface mesh represented as an object of 'TriRep' 
%           class, 'triangulation' class, or a cell such that TR = {Tri,V},
%           where Tri is an M-by-3 array of faces and V is an N-by-3 array
%           of vertex coordinates.
%
% OUTPUT
%   - isManifold : true if TR is manifold and false otherwise
%   - Enb        : K-by-2 non-boundary edge array
%   - Eb         : M-by-2 boundary edge array
%   - EA         : K-by-1 non-boundary edge attachment cell 
%   - Ne         : K-by-1 array containing number of elements in each EA cell
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


[Tri,V,fmt] = GetMeshData(TR);
if size(Tri,2)~=3
    error('This function is intended for TRIANGULAR surface meshes.')
end

%warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
if fmt>1
    TR = triangulation(Tri,V);
end

% Does TR have a boundary?
Enb = sort(edges(TR),2);
Eb = freeBoundary(TR);

isClosed = isempty(Eb);
if ~isClosed
    idx = ismember(Enb,sort(Eb,2),'rows');
    Enb(idx,:) = []; % non-boundary edges   
end

% Edge attachments
EA = edgeAttachments(TR,Enb);
Ne = cellfun(@length,EA);
isManifold = all(Ne==2); % all non-boundary edges connected to two faces
if isManifold && ~isClosed
    EAb = edgeAttachments(TR,Eb);
    Neb = cellfun(@length,EAb);
    isManifold = all(Neb==1); % all boundary edges connected to one face; probably redundant 
end
%warning('on','MATLAB:triangulation:PtsNotInTriWarnId')
