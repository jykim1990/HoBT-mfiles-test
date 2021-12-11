function [theta,a,L]=pcoord(cnr,cnr_con)
% theta : principal orient
% a : edge angle
% L : edge length

Ne=size(cnr_con,1);

%%% Edge angle & length %%%
a=zeros(Ne,1); % slope angle of each edge
L=zeros(Ne,1); % length of each edge
X=0;
Y=0;
for i=1:Ne % i-th edge
    x1=cnr(cnr_con(i,1),1);
    y1=cnr(cnr_con(i,1),2);        
    x2=cnr(cnr_con(i,2),1);
    y2=cnr(cnr_con(i,2),2);
    
    ai=atan2((y2-y1),(x2-x1));
    Li=sqrt((x2-x1)^2+(y2-y1)^2);
    
    a(i)=ai;
    L(i)=Li;
    
    % centroid term
    Xi=x1*Li+0.5*Li^2*cos(ai);
    Yi=y1*Li+0.5*Li^2*sin(ai);
    
    X=X+Xi;
    Y=Y+Yi;
end
cnt=[X/sum(L),Y/sum(L)]; % centroid

%%% principal axis for translation modes %%%
AA=sum(cos(a).*sin(a).*L); 
BB=sum((cos(a).^2-sin(a).^2).*L);

theta_t=0.5*atan(2*AA/BB);

%%% principal axis for rotation modes %%%
xc=cnr(cnr_con(:,1),1)-cnt(1);
yc=cnr(cnr_con(:,1),2)-cnt(2);

aa=sum(yc.^2.*L + yc.*sin(a).*L.^2 + 1/3*sin(a).^2.*L.^3);
bb=sum(xc.^2.*L + xc.*cos(a).*L.^2 + 1/3*cos(a).^2.*L.^3);
cc=sum(xc.*yc.*L + 1/2*xc.*sin(a).*L.^2 + 1/2*yc.*cos(a).*L.^2 + 1/3*sin(a).*cos(a).*L.^3);

theta_r=0.5*atan(2*cc/(bb-aa));

%%% axis determination %%%
oc=@(b) abs(sum(cos(a-b).*sin(a-b).*L)); % orthogonality check
if oc(0)<sum(L)*1e-14
    theta=0;
elseif oc(theta_r)<sum(L)*1e-14
    theta=theta_r;
else
    theta=theta_t;
end