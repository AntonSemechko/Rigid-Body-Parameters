function RBP_cuboid_validation(abc,N,K)
% Demo illustrating that the inertia tensor of a cuboid with sides A>=B>=C
% transforms as R*T*R' when the cuboid is rotated by R and that the 
% magnitudes of the computed principal axes of inertia remain the same 
% under rotation. 
%
%   - abc   : (optional) 1-by-3 array specifying lengths of the cuboid.
%             abc=[3 2 1] is the default setting.
%   - N     : (optional) number random rotations. N=1E3 is the default
%             setting.
%   - K     : (optional) number of triangle subdivisions. K=2 is the
%             default setting. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Experimet parameters
if nargin<1 || isempty(abc)
    abc=[3 2 1];
else
    abc=sort(abs(abc(1:3)),'descend');
end

if nargin<2 || isempty(N)
    N=1E3;    
else
    N=round(abs(N));
end

if nargin<3 || isempty(K)
    K=2;
else
    K=round(abs(K));
end


% Get the mesh
fv=make_cube_mesh;                              % unit cube mesh
fv.vertices=bsxfun(@times,abc(:)',fv.vertices); % scale the sides to specified lengths

% Subdivide?
if K>0
    for k=1:K, fv=TriQuad(fv); end
end

% Visualize the mesh
figure('color','w')
patch(fv,'FaceColor',0.8*[1 1 1],'EdgeColor','k','FaceAlpha',0.8);
axis equal

xlabel('x','FontSize',15);
ylabel('y','FontSize',15);
zlabel('z','FontSize',15);
view([30 10]);

% Inertia tensor of a cubiod with sides abc aligned with the Cartesian
% basis; see the accompanying .pdf (section D)
V=prod(abc)/12;
L=V*[abc(2)^2 + abc(3)^2, abc(1)^2 + abc(3)^2, abc(1)^2 + abc(2)^2];
To=diag(L);

% Verify that the function RigidBodyParams produces correct result for the
% non-rotated primitive
RBP=RigidBodyParams(fv);
T=RBP.inertia_tensor;


fprintf(2,"\nPART 1: Verifying that the intertia tensor of the non-rotated primitive has been computed correctly\n")
fprintf('This is the expected inertia tensor of cuboid with sides (%.6f,%.6f,%.6f) <-- To: \n',abc)
disp(To)

fprintf("\nHere is what the function RigidBodyParams computed <-- T:\n")
disp(T)

Err=norm(abs(To\T - eye(3)),'fro');
fprintf('Error = |To\\T - eye(3)| = %.6E\n',Err)


% Generate random rotation vectors inside the pi-ball
r=randn(N,3);
r=bsxfun(@rdivide,r,sqrt(sum(r.^2,2))+eps); % rotation vector
t=pi*rand(N,1).^(1/3); % rotation magnitude
t=min(max(1E-5,t),pi-1E-5);

% Verify that the inertia tensor transforms as expected
w=zeros(3,3);
[Err1,Err2]=deal(zeros(N,1));
D=zeros(N,3);
Vo=fv.vertices;
L=sort(L,'descend'); 
for i=1:N

    % Map rotation vector to SO(3)
    w(1,2)=-r(i,3);
    w(1,3)= r(i,2);
    w(2,3)=-r(i,1);
    w(2,1)= r(i,3);
    w(3,1)=-r(i,2);
    w(3,2)= r(i,1);    
    Ri=expm(t(i)*w);
    
    % Apply rotation
    fv.vertices=ApplySimTform(Vo,Ri);

    % Compute rigid body parameters of the rotated cuboid 
    RBP=RigidBodyParams(fv);
    Ti=RBP.inertia_tensor;

    % Error in the computed inertia tensor
    Tref=Ri*To*Ri';    
    Err1(i)=norm(abs(Tref\Ti - eye(3)),'fro');
    
    % Error in the magnitudes of the principal axes of inertia
    Err2(i)=sum(abs(L - RBP.eigs));
    
end

fprintf(2,"\nPART 2: Verifying that the intertia tensor of object rotated by R transforms as R*T*R' and its eigenvalue remain unchanged\n")

fprintf('Maximum inertia tensor error (Frobenious norm): %.6E\n',max(Err1))
fprintf('Maximum principal axis error (L1 norm): %.6E\n\n',max(Err2))


function fv=make_cube_mesh
% Generate triangle mesh of a unit cube. 
%

% Vertices
X=[1 1; 0 1; 0 0; 1 0];
[Xd,Xu]=deal(X);
Xd(:,3)=0;
Xu(:,3)=1;
X=[Xd;Xu];

% Quad faces
F=[1 5 8 4; ...
   2 6 5 1; ...
   3 7 6 2; ...
   4 8 7 3; ...
   5 6 7 8; ...
   1 4 3 2];

% Split quads into triangles
T=[F(:,[1 2 3]) F(:,[3 4 1])];
T=reshape(T',3,[])';


fv.faces=T;
fv.vertices=X;