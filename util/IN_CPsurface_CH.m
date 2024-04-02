function IN_CPsurface_CH(IN_xtlfile,CPfile,scaler,Pmax,meshsize)
%__________________________________________________________________________________________________
%
%                  ++++++
%                +++++++++
%              +++++++++++++                DDDDDDD    FFFFFFF TTTTTTTT       CCCCCC   PPPPPPP
%             +++++++++++++++               D      D   F          T          C         P      P
%             +++++++++++++++               D       D  FFFFFF     T         C          P      P
%             +++++++++++++++               D       D  F          T    ==== C          PPPPPPP
%             +++++++++++++++               D      D   F          T          C         P
%              +++++++++++++                DDDDDDD    F          T           CCCCCC   P
%               +++++++++++
%        _____   +++++++++   _____                  The Fredrickson Group Package for
%      -       -   ++++++   -      -                  DFT-Chemical Pressure Analysis
%     -           - ++++ -           -
%     -            -    -            -                  INTERFACE NUCLEUS APPROACH
%      - ______  -  ++++  -  _____  -
%                 ++++++++                               IN CP surface function. 
%                ++++++++++                          (c) 2023 The Fredrickson Group
%               ++++++++++++            __________________________________________________________
%              +++++++++++++           |                                                          |
%             +++++++++++++++          |    Render CP interface functions for clusters of atoms   |
%             +++++++++++++++          |    in regions common to parent structures and their      |
%             +++++++++++++++          |    intergrowth children. Use these images to design      |
%              +++++++++++++           |    Modular Intermetallic Frameworks.                     |
%               +++++++++++            |__________________________________________________________|
%                +++++++++
%                  +++++                            Last modified: 5 February 2024
%
%__________________________________________________________________________________________________ 
% 
% USAGE:   IN_CPsurface_CH(IN_xtlfile,CPfile,scaler,Pmax,meshsize)
%
%   IN_xtlfile:  Name of file (minus .xtl extension) containing fractional coordinates 
%                of Interface Nucleus (IN) atoms in format  %s %f %f %f. 
%                DIFFERENT FROM figuretool xtl FORMAT - Do not include cell vectors!
%
%
%   CPfile:      Base name for CP scheme output files.
%
%
%   scaler:       1 = draws surface through atoms, 
%                >1 = draw isotropically expanded surface just outside of convex hull,
%                <1 = draw isotropically contracted surface just within convex hull.
%              
%
%   Pmax:         0 = automatic surface color scheme.
%                %f = user selected value based on largest CP interface magnitude
%
%   meshsize:    Sets fineness of surface mesh.  
%                 36 = course mesh/faster generation, 
%                 100 = fine mesh/slower generation.
% 
%
% Example: >> IN_CPsurface_CH('Mg2Zn11-CP_IN','Mg2Zn11',1.01,0.00015,100)


atomsfile=strcat(CPfile,'-geo');
cellfile=strcat(CPfile,'-cell');
IN_xtlfile=strcat(IN_xtlfile,'.xtl');


a=zeros(3,3);
[a1,a2,a3]=textread(cellfile,' %f %f %f');
a=[a1,a2,a3];
a_cell=a(1,:);
b_cell=a(2,:);
c_cell=a(3,:);

[names,xf,yf,zf] = textread(IN_xtlfile,'%s %f %f %f');

natoms = size(xf,1);
for j=1:natoms
  xcart1(j,1)=xf(j,1)*a(1,1)+yf(j,1)*a(2,1)+zf(j,1)*a(3,1);
  ycart1(j,1)=xf(j,1)*a(1,2)+yf(j,1)*a(2,2)+zf(j,1)*a(3,2);
  zcart1(j,1)=xf(j,1)*a(1,3)+yf(j,1)*a(2,3)+zf(j,1)*a(3,3);
end


xcenter = mean(xcart1);
ycenter = mean(ycart1);
zcenter = mean(zcart1);


hull=convhull(xcart1,ycart1,zcart1);
npoints=size(hull,1);
hullvertices = zeros(size(xcart1));


for j=1:npoints
     hullvertices(hull(j,1),1)=1;
     hullvertices(hull(j,2),1)=1;
     hullvertices(hull(j,3),1)=1;
end

k1=0;
for j = 1:natoms
    if(hullvertices(j,1)==1)
         k1=k1+1;
%         drawspheres(xcart1(j,1),ycart1(j,1),zcart1(j,1),40,0.45,[255 0 0]/255);
         xcart2(k1,1)=xcart1(j,1);
         ycart2(k1,1)=ycart1(j,1);
         zcart2(k1,1)=zcart1(j,1);
    end  
end

xcart3 = xcart1-xcenter;
ycart3 = ycart1-ycenter;
zcart3 = zcart1-zcenter;
%sph = alphaShape(xcart3, ycart3, zcart3,10);
%tri = alphaTriangulation(sph);

[Xj,Yj,Zj]=sphere(meshsize);
[phi,theta,sink]=cart2sph(Xj,Yj,Zj);
thetaj=pi/2*ones(size(phi))-theta;
thetaj=abs(thetaj);
sizex = size(Xj,1);
sizey = size(Xj,2);
for jx = 1:sizex
  for jy = 1:sizey
     [dx,dy,dz]=sph2cart(phi(jx,jy),theta(jx,jy),1.0);
      r_min = 100.00;
      for k=1:npoints
         x1 = xcart3(hull(k,1),1);
         y1 = ycart3(hull(k,1),1);
         z1 = zcart3(hull(k,1),1);
         x2 = xcart3(hull(k,2),1);
         y2 = ycart3(hull(k,2),1);
         z2 = zcart3(hull(k,2),1);
         x3 = xcart3(hull(k,3),1);
         y3 = ycart3(hull(k,3),1);
         z3 = zcart3(hull(k,3),1);
         v_vector = cross([x1-x2 y1-y2 z1-z2],[x3-x2 y3-y2 z3-z2]);
         v_vector = v_vector/norm(v_vector);
         const = v_vector*[x1 y1 z1]';
         r_test = const/([dx dy dz]*v_vector');
         if((r_test > 0) &&(r_test < r_min))
             r_min = r_test;
         end
      end
      r_surf(jx,jy) = r_min;
  end
end

[Xsurf,Ysurf,Zsurf] = sph2cart(phi,theta,scaler*r_surf);
CPinterface_IN(CPfile,4,Xsurf+xcenter,Ysurf+ycenter,Zsurf+zcenter,4.0,0.75,Pmax,0.5); 

function CPinterface_IN(filename,lmaxuse,surfx,surfy,surfz,dmax,zeta_factor,Pmax,deltad) 
% CPabinit(filename,templatefile,scales,res)

atomsfile=strcat(filename,'-geo');
cellfile=strcat(filename,'-cell');

[name,x,y,z]=textread(atomsfile,'%s %f %f %f');
atomnum=size(x);
atomnum=atomnum(1);

a=zeros(3,3);
[a1,a2,a3]=textread(cellfile,' %f %f %f');
a=[a1,a2,a3];
a_cell=a(1,:);
b_cell=a(2,:);
c_cell=a(3,:);
cellv=abs(cross(b_cell,c_cell)*a_cell');

basicsurf = surf(surfx,surfy,surfz,'Visible','off');
[orbname,coeff]=textread(strcat(filename,'-coeff'),'%s %f');
size_coeff=size(coeff);
size_coeff=size_coeff(1)/atomnum;
lmax=size_coeff^0.5-1;
atom_coeff=reshape(coeff,size_coeff,atomnum);
%atom_coeff(1, atomn)=scales*atom_coeff(1, atomn);
%counter=2;
%for l=1:lmax
%    for m=1:(2*l+1)
%      atom_coeff(counter, atomn)=scales*atom_coeff(counter, atomn);
%      if (l>lmaxuse)  
%         atom_coeff(counter, atomn) = 0;
%      end
%      counter=counter+1;
%    end
%end
cellrange = [-3 3 -3 3 -3 3];

sizex = size(surfx,1);
sizey = size(surfx,2);
for jx = 1:sizex
 for jy = 1:sizey
    CPsurf(jx,jy) = 0;
    sx = surfx(jx,jy);
    sy = surfy(jx,jy);
    sz = surfz(jx,jy);
    normvect(1,1) = basicsurf.VertexNormals(jx,jy,1);
    normvect(1,2) = basicsurf.VertexNormals(jx,jy,2);
    normvect(1,3) = basicsurf.VertexNormals(jx,jy,3);
    for nx = 1:size(y)
        r=[x(nx),
           y(nx),
           z(nx)];
          mu2grad_add=0;
          COHPgrad_add=0;
          for j1 = cellrange(1):cellrange(2)
          for k1 = cellrange(3):cellrange(4)
          for l1 = cellrange(5):cellrange(6)
                xnew=x(nx)+j1*a(1,1)+k1*a(2,1)+l1*a(3,1);
                ynew=y(nx)+j1*a(1,2)+k1*a(2,2)+l1*a(3,2);
                znew=z(nx)+j1*a(1,3)+k1*a(2,3)+l1*a(3,3);
                deltax = sx-xnew+normvect(1,1)*deltad;
                deltay = sy-ynew+normvect(1,2)*deltad;
                deltaz = sz-znew+normvect(1,3)*deltad;
                if(deltax*deltax+deltay*deltay+deltaz*deltaz < 1.1*dmax^2)
                  deltavect = [deltax,deltay,deltaz];
                  deltavect = deltavect/norm(deltavect);
                  if(norm(deltavect) < 0.0001) 
                      deltavect = normvect;
                  end
                  tempfound=0;
                  [phi,theta,dist]=cart2sph(deltax,deltay,deltaz);
                  if(dist < dmax) 
                      CPsurf(jx,jy) = CPsurf(jx,jy)+0.5*calcCPangle(deltax,deltay,deltaz,atom_coeff(:,nx),lmaxuse)*exp(-zeta_factor*dist)*dot(deltavect,normvect)*sign(dot(deltavect,normvect));
                  end 
                  deltax = sx-xnew-normvect(1,1)*deltad;
                  deltay = sy-ynew-normvect(1,2)*deltad;
                  deltaz = sz-znew-normvect(1,3)*deltad;
                  deltavect = [deltax,deltay,deltaz];
                  deltavect = deltavect/norm(deltavect);
                  if(norm(deltavect) < 0.0001) 
                      deltavect = normvect;
                  end
                  tempfound=0;
                  [phi,theta,dist]=cart2sph(deltax,deltay,deltaz);
                  if(dist < dmax) 
                      CPsurf(jx,jy) = CPsurf(jx,jy)+0.5*calcCPangle(deltax,deltay,deltaz,atom_coeff(:,nx),lmaxuse)*exp(-zeta_factor*dist)*dot(deltavect,normvect)*sign(dot(deltavect,normvect));
                  end
                end
          end      
          end
          end
    end
 end
end


[max(max(CPsurf)),min(min(CPsurf))]
if(Pmax==0)
   Pmax=max(abs([max(max(CPsurf)),min(min(CPsurf))])); 
end

for jx = 1:sizex
 for jy = 1:sizey
         CPnorm = CPsurf(jx,jy)/Pmax;
         if(CPnorm >= 0) 
              CPcolor(jx,jy,1:3) = CPnorm*[255 90 0]/255 + (1-CPnorm)*[1 1 1];
         end
         if(CPnorm < 0) 
              CPnorm = abs(CPnorm);
              CPcolor(jx,jy,1:3) = CPnorm*[0 90/255 1] + (1-CPnorm)*[1 1 1];
         end
         CPcolor2(jx,jy) = CPsurf(jx,jy)/Pmax;
         FAlpha(jx,jy) = (1.0*norm(CPnorm)+0.5);
         if(FAlpha(jx,jy) > 1.0)
             FAlpha(jx,jy) = 1;
         end
%         FAlpha(jx,jy) = 1;

 end
end

[max(max(max(CPcolor2))),min(min(min(CPcolor2)))]

CPcolor2;
%max(FAlpha)
%min(FAlpha)
%surf(surfx,surfy,surfz,CPcolor2,'FaceColor','interp','EdgeColor','none','AmbientStrength',1,'AlphaData',FAlpha,'FaceAlpha','interp','DiffuseStrength',0,'FaceLighting','gouraud');
%surf(surfx,surfy,surfz,'CData',CPcolor,'FaceColor','interp','EdgeColor','none','AmbientStrength',1,'AlphaData',FAlpha,'FaceAlpha','interp','AlphaDataMapping','none','DiffuseStrength',0,'FaceLighting','gouraud');
surf(surfx,surfy,surfz,CPcolor2,'FaceColor','interp','EdgeColor','none','AmbientStrength',1,'SpecularStrength',0,'FaceAlpha',.6,'DiffuseStrength',0,'FaceLighting','gouraud');
%FaceAlpha = 1 full opacity, = 0 full transparency, =0.6 default

colormap jet;
caxis([-1 1]);
%

function CP = calcCPangle(x,y,z,mu_components,lmax)
%

[phi,theta,sink]=cart2sph(x,y,z);
thetaj=pi/2*ones(size(phi))-theta;
thetaj=abs(thetaj);
rnew=0*ones(size(phi));
rnew=rnew+(1/(4*pi))^0.5*mu_components(1,1)*ones(size(phi));
counter=1;
for l=1:lmax
    m=0;
    Plm_0=2^0.5*legendre(l,cos(thetaj),'norm');
    size_Plm=size(Plm_0);
    Plm=0;
    Plm=(-1)^m*Plm_0(m+1,1);
    counter=counter+1;
    rnew=rnew+mu_components(counter,1)*Plm/((4*pi)^0.5);
    for m=1:l
         counter=counter+1;
         Plm_0=2^0.5*legendre(l,cos(thetaj),'norm');
         size_Plm=size(Plm_0);
         Plm=0;
         Plm=(-1)^m*Plm_0(m+1,1);
         rnew=rnew+mu_components(counter,1)*Plm*2^0.5.*cos(m*phi)/((4*pi)^0.5);
         counter=counter+1;
         rnew=rnew+mu_components(counter,1)*Plm*2^0.5.*sin(m*phi)/((4*pi)^0.5);
    end
end
CP = rnew;

function [Xnew, Ynew, Znew] = plot_Ylm(x,y,z,mu_components,meshsize,lmax,scaler,mode)
%
positivecolor=[0.3 0.3 0.3];
negativecolor=[0 0 0];

meshsize=round(3.8*meshsize);
[Xj,Yj,Zj]=sphere(meshsize);
[phi,theta,sink]=cart2sph(Xj,Yj,Zj);
thetaj=pi/2*ones(size(phi))-theta;
thetaj=abs(thetaj);
rnew=0*ones(size(phi));
rnew=rnew+(1/(4*pi))^0.5*mu_components(1,1)*ones(size(phi));
counter=1;
for l=1:lmax
    m=0;
    Plm_0=2^0.5*legendre(l,cos(thetaj),'norm');
    size_Plm=size(Plm_0);
    Plm=0;
    for j=1:size_Plm(1,2)
           for k=1:size_Plm(1,3)
             Plm(j,k)=(-1)^m*Plm_0(m+1,j,k);
           end
    end
    counter=counter+1;
    rnew=rnew+mu_components(counter,1)*Plm/((4*pi)^0.5);
    for m=1:l
         counter=counter+1;
         Plm_0=2^0.5*legendre(l,cos(thetaj),'norm');
         size_Plm=size(Plm_0);
         Plm=0;
         for j=1:size_Plm(1,2)
           for k=1:size_Plm(1,3)
             Plm(j,k)=(-1)^m*Plm_0(m+1,j,k);
           end
         end
         rnew=rnew+mu_components(counter,1)*Plm*2^0.5.*cos(m*phi)/((4*pi)^0.5);
         counter=counter+1;
         rnew=rnew+mu_components(counter,1)*Plm*2^0.5.*sin(m*phi)/((4*pi)^0.5);
    end
end

[Xnew,Ynew,Znew]=sph2cart(phi,theta,scaler*abs(rnew));
C=rnew;
size(C);
d=.5*ones(size(C));
alpha=ones(size(C));
alphacm=ones(size(C));
m=size(C);
for n=1:m(1)
  for o=1:m(2)

  cmH(n,o,:)=positivecolor;
  %         = RGB color code for positive CP lobes.

  cmHe(n,o,:)=[.1,.1,.1];
  cm(n,o,:)=[0.7,1.0,0.7];
  cmS(n,o,:)=[0.7,0.7,.7];
  cmF(n,o,:)=[0.8,0.7,0.7];
  if C(n,o) < 0 
           d(n,o)=-1;
           alpha(n,o)=1;
           alphacm(n,o)=1.0;
           cm(n,o,:)=[1,.8,.8];

           cmH(n,o,:)=negativecolor;
           %         = RGB color code for negative CP lobes.

           cmHe(n,o,:)=[.7,.7,.7];
           cmF(n,o,:)=[0,100/255,0];
           cmS(n,o,:)=[0,0,1];

  end

  end
end
Xnew=Xnew+ones(size(Xnew))*x;
Ynew=Ynew+ones(size(Ynew))*y;
Znew=Znew+ones(size(Znew))*z;
hold on
ec=[0,0,0];
if rnew < .4 
        ec='none';
end
gotit=0;

hold on
if(mode==0)
  g=surf(Xnew,Ynew,Znew,d,'AmbientStrength',.8,'FaceAlpha',.5,'EdgeColor','none','Linewidth',0.5,'EdgeAlpha',1,'AlphaData',alphacm,'AlphaDataMapping','none','FaceColor','interp','CData',cmH);
  colormap([.2,.20,.20;.8 .8 .8]);
  gotit=1;
end

function [Xnew,Ynew,Znew] = Ylm_surfacepoint(x,y,z,xcen,ycen,zcen, Ylm_components,lmax)
%

x = x - ones(size(x))*xcen;
y = y - ones(size(y))*ycen;
z = z - ones(size(z))*zcen;

[phi,theta,sink]=cart2sph(x,y,z);
thetaj=pi/2*ones(size(phi))-theta;
thetaj=abs(thetaj);
rnew=0*ones(size(phi));
rnew=rnew+(1/(4*pi))^0.5*Ylm_components(1,1)*ones(size(phi));
counter=1;
for l=1:lmax
    m=0;
    Plm_0=2^0.5*legendre(l,cos(thetaj),'norm');
    size_Plm=size(Plm_0);
    Plm=0;
    for j=1:size_Plm(1,2)
             Plm(j,1)=(-1)^m*Plm_0(m+1,j);
    end
    counter=counter+1;
    rnew=rnew+Ylm_components(counter,1)*Plm/((4*pi)^0.5);
    for m=1:l
         counter=counter+1;
         Plm_0=2^0.5*legendre(l,cos(thetaj),'norm');
         size_Plm=size(Plm_0);
         Plm=0;
         for j=1:size_Plm(1,2)
             Plm(j,1)=(-1)^m*Plm_0(m+1,j);
         end
         rnew=rnew+Ylm_components(counter,1)*Plm*2^0.5.*cos(m*phi)/((4*pi)^0.5);
         counter=counter+1;
         rnew=rnew+Ylm_components(counter,1)*Plm*2^0.5.*sin(m*phi)/((4*pi)^0.5);
    end
end

max(max(max(rnew)));
[Xnew,Ynew,Znew]=sph2cart(phi,theta,abs(rnew));
Xnew = Xnew + ones(size(Xnew))*xcen;
Ynew = Ynew + ones(size(Xnew))*ycen;
Znew = Znew + ones(size(Xnew))*zcen;





function drawspheres(x,y,z,meshsize,atomsize,rgbvalue)
%
J=[1,0,0,0,0,0,0,0,0];
s=1;
px=0;
py=0;
pz=0;
dx2=0;
dz2=0;
dxy=0;
dxz=0;
dyz=0;
orbsize=J*transpose(J); 
[Xj,Yj,Zj]=sphere(meshsize);
[phi,theta,sink]=cart2sph(Xj,Yj,Zj);
rs=s*ones(size(phi));
rpx=px*cos(theta).*cos(phi);
rpy=py*cos(theta).*sin(phi);
rpz=pz*sin(theta);
rdx2=1/2*dx2*cos(theta).*cos(theta).*cos(2*phi);
rdz2=1/(12)^.5*dz2*(3*sin(theta).*sin(theta)-ones(size(theta)));
rdxy=1/2*dxy*cos(theta).*cos(theta).*sin(2*phi);
rdxz=dxz*sin(theta).*cos(theta).*cos(phi);
rdyz=dyz*sin(theta).*cos(theta).*sin(phi);

rnew=atomsize*(rs+(rpx+rpy+rpz)+(rdx2+rdz2+rdxy+rdxz+rdyz));
[Xnew,Ynew,Znew]=sph2cart(phi,theta,abs(rnew));
C=rnew;
d=.5*ones(size(C));
alpha=1*ones(size(C));
alphacm=1*ones(size(C));
m=size(C);
for n=1:m(1)
  for o=1:m(2)
  cmC(n,o,:)=[0,0,0];
  cmN(n,o,:)=[0,0,1];
  cm(n,o,:)=[1,0,0];
  cmH(n,o,:)=[.9,.9,.9];
  cmS(n,o,:)=[0,0,1];
  d(n,o)=-1;
  end
end
Xnew=Xnew+ones(size(Xnew))*x;
Ynew=Ynew+ones(size(Ynew))*y;
Znew=Znew+ones(size(Znew))*z;
hold on
ec=[0,0,0];
if rnew < .4 
	ec='none';
end

g=surf(Xnew,Ynew,Znew,d,'EdgeColor','none','AmbientStrength',1,'FaceAlpha','flat','AlphaData',alpha,'AlphaDataMapping','none','FaceColor',rgbvalue,'Linewidth',.3);
%colormap([.2,.20,.20;1 1 1]);

