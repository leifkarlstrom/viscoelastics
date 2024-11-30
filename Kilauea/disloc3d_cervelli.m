function [U,D,S]=disloc3d_cervelli(m,xy,nu,mu)
% function [U,D,S]=disloc3d(m,xy,nu,mu)
% m (10 x 1)
% xy (3 x nsites)
% 
%Returns the deformation at point 'x', given dislocation
%model 'm'. 'nu' specifies the Poisson's ratio and 'mu'
%specifies the shear modulus.
%
%The coordinate system is as follows: east = positive X,
%north = positive Y, and up = positive Z. Observation
%coordinates with positive Z values return a warning.
%
%The outputs are 'U', the three displacement components:
%east, north, and up (on the rows); 'D', the nine
%spatial derivatives of the displacement: Uxx, Uxy, Uxz,
%Uyx, Uyy, Uyz, Uzx, Uzy, and Uzz (on the rows); and 'S'
%the 6 independent stress components: Sxx, Sxy, Sxz, Syy,
%Syz, and Szz (on the rows). All these outputs have
%the same number of columns as 'x'.
%
%The dislocation model is specified as: length, width,
%depth, dip, strike, east, north, strike-slip, dip-slip,
%and opening. The coordinates (depth, east, and north)
%specify a point at the middle of the bottom edge of
%the fault for positive dips and the middle of the top
%edge for negative dips.
%
%Mind your units! E.g., if the lengths are given in km and the
%slips in m, the derivatives will be biased
%
% P. Cervelli ???

%Convert to Okada's parameterization

    alpha=1/(2-2*nu);
    
    dip=m(4);
    cd=cos(dip*pi/180);
    sd=sin(dip*pi/180);
    
    strike=(m(5)-90)*pi/180;
    cs=cos(strike);
    ss=sin(strike);
    
    disl1=m(8);
    disl2=m(9);
    disl3=m(10);

    depth=m(3)-0.5*m(2)*sd;
    
    al1=m(1)/2;
    al2=al1;
    aw1=m(2)/2;
    aw2=aw1;

    x=cs*(-m(6) + xy(1,:)) - ss*(-m(7) + xy(2,:));
    y=-0.5*cd*m(2) + ss*(-m(6) + xy(1,:)) + cs*(-m(7) +  xy(2,:));
    z=xy(3,:);
    
    sing=0;
    
%Check for unphysical model
 
% (m(3)-sd*m(2))<0
% m(3)
% m(1) <=0
% m(2) <=0
% m(3) < 0
     if (m(3)-sd*m(2))<0 || m(1) <=0 || m(2) <=0 || m(3) < 0
         warning(['Unphysical model.']);
         U=zeros(3,size(xy,2));
         D=zeros(9,size(xy,2));
         S=zeros(6,size(xy,2));
         return
     end
     
%Calculate medium constants and fault-dip constants

    c0=dccon0(alpha,dip);
    u=zeros(12,size(x,2));

%Real source contribution

    d=depth+z;
    p=y.*c0.cd+d.*c0.sd;
    q=y.*c0.sd-d.*c0.cd;

    I=(x+al1).*(x-al2)<=0;
    jxi(I)=1;
    jxi(~I)=0;
    
%     if (x+al1).*(x-al2) <=0
%         jxi=1;
%     else
%         jxi=0;
%     end

    I=(p+aw1).*(p-aw2) <=0;
    jet(I)=1;
    jet(~I)=0;
    
%     if (p+aw1).*(p-aw2) <=0
%         jet=1;
%     else
%         jet=0;
%     end

%     dd1=disl1;
%     dd2=disl2;
%     dd3=disl3;

    for k=1:2
        switch k
            case 1
                et=p+aw1;
            case 2
                et=p-aw2;
        end

        for  j=1:2;
            switch j
                case 1
                    xi=x+al1;
                case 2
                    xi=x-al2;
            end

            c2=dcon2(c0.cd,et,q,c0.sd,xi);

            J=(jxi==1 & q==0 & et==0) | (jet==1 & q==0 & xi==0);
            if any(J)
                sing(J)=1;
            end

            dua=ua(c0,c2,xi,et,q,disl1,disl2,disl3);

            R=sparse([1 2 3 2 3 4 5 6 5 6 7 8 9 8 9 10 11 12 11 12], ...
                     [1 2 2 3 3 4 5 5 6 6 7 8 8 9 9 10 11 11 12 12], ...
                     [-1 -c0.cd -c0.sd c0.sd -c0.cd -1 -c0.cd -c0.sd c0.sd -c0.cd -1 -c0.cd -c0.sd c0.sd -c0.cd 1 c0.cd c0.sd -c0.sd c0.cd]);

            du=R*dua;
            
            switch j+k
                case 3
                    u=u-du;
                otherwise
                    u=u+du;
            end
        end
    end

%Image source contribution

    %zz=z;
    d=depth-z;
    p=y.*c0.cd+d.*c0.sd;
    q=y.*c0.sd-d.*c0.cd;

%     I=(p+aw1).*(p-aw2)<=0
%     
%     if (p+aw1).*(p-aw2) <= 0
%         jet=1;
%     else
%         jet=0;
%     end

    for k=1:2
        switch k
            case 1
                et=p+aw1;
            case 2
                et=p-aw2;
        end

        for  j=1:2;
            switch j
                case 1
                    xi=x+al1;
                case 2
                    xi=x-al2;
            end

            c2=dcon2(c0.cd,et,q,c0.sd,xi);
            dua=ua(c0,c2,xi,et,q,disl1,disl2,disl3);
            dub=ub(c0,c2,xi,et,q,disl1,disl2,disl3);
            duc=uc(c0,c2,xi,et,q,z,disl1,disl2,disl3);
            
%             R2=-R;
%             R2(3:3:end,:)=-R2(3:3:end,:);
%             R2(10:12,1:3)=R2(1:3,1:3)/z;
%           
%             du=-R*dua+ -R*dub + z*R2*duc;
%             du(10)=du(10)+duc(1);

            
            for i=1:3:10
                du(i,:)=dua(i,:)+dub(i,:)+z.*duc(i,:);
                du(i+1,:)=(dua(i+1,:)+dub(i+1,:)+z.*duc(i+1,:))*c0.cd-(dua(i+2,:)+dub(i+2,:)+z.*duc(i+2,:))*c0.sd;
                du(i+2,:)=(dua(i+1,:)+dub(i+1,:)-z.*duc(i+1,:))*c0.sd+(dua(i+2,:)+dub(i+2,:)-z.*duc(i+2,:))*c0.cd;
            end
            du(10,:)=du(10,:)+duc(1,:);
            du(11,:)=du(11,:)+duc(2,:)*c0.cd-duc(3,:)*c0.sd;
            du(12,:)=du(12,:)-duc(2,:)*c0.sd-duc(3,:)*c0.cd;
          
            switch j+k
                case 3
                    u=u-du;
                otherwise
                    u=u+du;
            end
        end
    end

%Rotate out of fault coordinate system

      U(1,:)=cs*u(1,:)+ss*u(2,:);
      U(2,:)=-ss*u(1,:)+cs*u(2,:);
      U(3,:)=u(3,:);

      D(1,:)=cs^2*u(4,:) + cs*ss*(u(7,:) + u(5,:)) + ss^2*u(8,:);
      D(2,:)=cs^2*u(7,:) - ss^2*u(5,:) + cs*ss*(-u(4,:) + u(8,:));
      D(3,:)=cs*u(10,:) + ss*u(11,:);
      D(4,:)=-(ss*(cs*u(4,:) + ss*u(7,:))) + cs*(cs*u(5,:) + ss*u(8,:));
      D(5,:)=ss^2*u(4,:) - cs*ss*(u(7,:) + u(5,:)) + cs^2*u(8,:);
      D(6,:)=-(ss*u(10,:)) + cs*u(11,:);
      D(7,:)=cs*u(6,:) + ss*u(9,:);
      D(8,:)=-(ss*u(6,:)) + cs*u(9,:);
      D(9,:)=u(12,:);
      
%Calculate stress

      lambda=2*mu*nu/(1-2*nu);
      theta=D(1,:)+D(5,:)+D(9,:);

      S(1,:)=lambda*theta+2*mu*D(1,:);
      S(2,:)=mu*(D(2,:)+D(4,:));
      S(3,:)=mu*(D(3,:)+D(7,:));
      S(4,:)=lambda*theta+2*mu*D(5,:);
      S(5,:)=mu*(D(6,:)+D(8,:));
      S(6,:)=lambda*theta+2*mu*D(9,:);
      
%Flag singularities

    if any(sing)
        warning('At least one solution is singular.')
        U(:,sing)=0;
        D(:,sing)=0;
        S(:,sing)=0;
    end      
      
function c0=dccon0(alpha,dip)
%Calculate medium constants and fault-dip constants

EPS=1.0e-6;

c0.alp1=(1-alpha)/2;
c0.alp2= alpha/2;
c0.alp3=(1-alpha)/alpha;
c0.alp4= 1-alpha;
c0.alp5= alpha;

c0.sd=sin(dip*pi/180);
c0.cd=cos(dip*pi/180);

if abs(c0.cd)<EPS
    c0.cd=0;   % Vertical? EMB
    if c0.sd>0
        c0.sd=1;
    end
    if c0.sd<0
        c0.sd=-1;
    end
end

c0.sdsd=c0.sd.*c0.sd;
c0.cdcd=c0.cd.*c0.cd;
c0.sdcd=c0.sd.*c0.cd;
c0.s2d=2*c0.sdcd;
c0.c2d=c0.cdcd-c0.sdsd;

function c2=dcon2(cd,et,q,sd,xi)
%Calculate station geometry constants for finite source

EPS=1.0e-6;

% c2=struct('xi2', 0, 'et2', 0, 'q2', 0, 'r', 0, 'r2', 0, 'r3', 0, 'r5', 0, ...
%           'y', 0, 'd', 0, 'tt', 0, 'alx', 0, 'ale', 0, 'x11', 0, 'y11', 0, ...
%           'x32', 0, 'y32', 0, 'ey', 0, 'ez', 0, 'fy', 0, 'fz', 0, 'gy', 0, ...
%           'gz', 0, 'hy', 0, 'hz', 0);

I=abs(xi)<EPS;
xi(I)=0;

I=abs(et)<EPS;
et(I)=0;

I=abs(q)<EPS;
q(I)=0;

c2.xi2=xi.*xi;
c2.et2=et.*et;
c2.q2=q.*q;
c2.r2=c2.xi2+c2.et2+c2.q2;
c2.r=sqrt(c2.r2);

if c2.r~=0

    c2.r3=c2.r.*c2.r2;
    c2.r5=c2.r3.*c2.r2;
    c2.y =et.*cd+q.*sd;
    c2.d =et.*sd-q.*cd;

    I=q==0;
    c2.tt(I)=0;
    c2.tt(~I)=atan(xi(~I).*et(~I)./(q(~I).*c2.r(~I)));
    
%     if q==0
%         c2.tt=0;
%     else
%         c2.tt=atan(xi.*et./(q.*c2.r));
%     end

    I=~(xi<0 & q==0 & et==0);
    
    c2.alx(~I)=-log(c2.r(~I)-xi(~I));
    c2.x11(~I)=0;
    c2.x32(~I)=0;
        
    rxi=c2.r(I)+xi(I);
    c2.alx(I)=log(rxi);
    c2.x11(I)=1./(c2.r(I).*rxi);
    c2.x32(I)=(c2.r(I)+rxi).*c2.x11(I).*c2.x11(I)./c2.r(I);
    
%     if xi<0 & q==0 & et==0
%         c2.alx=-log(c2.r-xi);
%         c2.x11=0;
%         c2.x32=0;
%     else
%         rxi=c2.r+xi;
%         c2.alx=log(rxi);
%         c2.x11=1./(c2.r.*rxi);
%         c2.x32=(c2.r+rxi).*c2.x11.*c2.x11./c2.r;
%     end

    I=~(et<0 & q==0 & xi==0);
    
    c2.ale(~I)=-log(c2.r(~I)-et(~I));
    c2.y11(~I)=0;
    c2.y32(~I)=0;
    
    ret=c2.r(I)+et(I);
    c2.ale(I)=log(ret);
    c2.y11(I)=1./(c2.r(I).*ret);
    c2.y32(I)=(c2.r(I)+ret).*c2.y11(I).*c2.y11(I)./c2.r(I);
        
%     if et<0 & q==0 & xi==0
%         c2.ale=-log(c2.r-et);
%         c2.y11=0;
%         c2.y32=0;
%     else
%         ret=c2.r+et;
%         c2.ale=log(ret);
%         c2.y11=1./(c2.r.*ret);
%         c2.y32=(c2.r+ret).*c2.y11.*c2.y11./c2.r;
%     end

    c2.ey=sd./c2.r-c2.y.*q./c2.r3;
    c2.ez=cd./c2.r+c2.d.*q./c2.r3;
    c2.fy=c2.d./c2.r3+c2.xi2.*c2.y32.*sd;
    c2.fz=c2.y./c2.r3+c2.xi2.*c2.y32.*cd;
    c2.gy=2.*c2.x11.*sd-c2.y.*q.*c2.x32;
    c2.gz=2.*c2.x11.*cd+c2.d.*q.*c2.x32;
    c2.hy=c2.d.*q.*c2.x32+xi.*q.*c2.y32.*sd;
    c2.hz=c2.y.*q.*c2.x32+xi.*q.*c2.y32.*cd;

end

function du=ua(c0,c2,xi,et,q,disl1,disl2,disl3)
%Displacement and strain at depth (part-a) due to buried finite fault in a
%semiinfinite medium

xy=xi.*c2.y11;
qx=q.*c2.x11;
qy=q.*c2.y11;

du(1,:)=disl1*(c2.tt/2+c0.alp2*xi.*qy) + ...
		disl2*(c0.alp2*q./c2.r) + ...
		disl3*(-c0.alp1*c2.ale-c0.alp2*q.*qy);
du(2,:)=disl1*(c0.alp2*q./c2.r) + ...
		disl2*(c2.tt/2+c0.alp2*et.*qx) + ...
		disl3*(-c0.alp1*c2.alx-c0.alp2*q.*qx);
du(3,:)=disl1*(c0.alp1*c2.ale-c0.alp2*q.*qy) + ...
		disl2*(c0.alp1*c2.alx-c0.alp2*q.*qx) + ...
		disl3*(c2.tt/2-c0.alp2*(et.*qx+xi.*qy));
du(4,:)=disl1*(-c0.alp1*qy-c0.alp2*c2.xi2.*q.*c2.y32) + ...
		disl2*(-c0.alp2*xi.*q./c2.r3) + ...
		disl3*(-c0.alp1*xy+c0.alp2.*xi.*c2.q2.*c2.y32);
du(5,:)=disl1*(-c0.alp2*xi.*q./c2.r3) + ...
		disl2*(-qy/2-c0.alp2*et.*q./c2.r3) + ...
		disl3*(-c0.alp1./c2.r+c0.alp2*c2.q2./c2.r3);
du(6,:)=disl1*(c0.alp1*xy+c0.alp2*xi.*c2.q2.*c2.y32) + ...
		disl2*(c0.alp1./c2.r+c0.alp2*c2.q2./c2.r3) + ...
		disl3*(-c0.alp1*qy-c0.alp2*q.*c2.q2.*c2.y32);
du(7,:)=disl1*(c0.alp1*xy.*c0.sd+c0.alp2*xi.*c2.fy+c2.d/2.*c2.x11) + ...
		disl2*(c0.alp2*c2.ey) + ...
		disl3*(-c0.alp1*(c0.cd./c2.r+qy.*c0.sd)-c0.alp2*q.*c2.fy);
du(8,:)=disl1*(c0.alp2*c2.ey) + ...
		disl2*(c0.alp1*c2.d.*c2.x11+xy/2*c0.sd+c0.alp2*et.*c2.gy) + ...
		disl3*(-c0.alp1*c2.y.*c2.x11-c0.alp2*q.*c2.gy);
du(9,:)=disl1*(c0.alp1*(c0.cd./c2.r+qy.*c0.sd)-c0.alp2*q.*c2.fy) + ...
		disl2*(c0.alp1*c2.y.*c2.x11-c0.alp2*q.*c2.gy) + ...
		disl3*(c0.alp1*(c2.d.*c2.x11+xy.*c0.sd)+c0.alp2*q.*c2.hy);
du(10,:)=disl1*(c0.alp1*xy.*c0.cd+c0.alp2*xi.*c2.fz+c2.y/2.*c2.x11) + ...
		 disl2*(c0.alp2*c2.ez) + ...
		 disl3*(c0.alp1*(c0.sd./c2.r-qy.*c0.cd)-c0.alp2*q.*c2.fz);
du(11,:)=disl1*(c0.alp2*c2.ez) + ...
		 disl2*(c0.alp1*c2.y.*c2.x11+xy/2*c0.cd+c0.alp2*et.*c2.gz) + ...
		 disl3*(c0.alp1*c2.d.*c2.x11-c0.alp2*q.*c2.gz);
du(12,:)=disl1*(-c0.alp1*(c0.sd./c2.r-qy.*c0.cd)-c0.alp2*q.*c2.fz) + ...
		 disl2*(-c0.alp1*c2.d.*c2.x11-c0.alp2*q.*c2.gz) + ...
		 disl3*(c0.alp1*(c2.y.*c2.x11+xy.*c0.cd)+c0.alp2*q.*c2.hz);
     
du=du/(2*pi);

function du=ub(c0,c2,xi,et,q,disl1,disl2,disl3)
%Displacement and strain at depth (part-b) due to buried finite fault in a
%semiinfinite medium

rd=c2.r+c2.d;
d11=1./(c2.r.*rd);
aj2=xi.*c2.y./rd.*d11;
aj5=-(c2.d+c2.y.*c2.y./rd).*d11;

if c0.cd~=0 % If cosine of dip is not zero (cosd(90) == 0, vertical)
    I=xi~=0;
    ai4(~I)=0;
    x=sqrt(c2.xi2(I)+c2.q2(I));
    ai4(I)=1/c0.cdcd*(xi(I)./rd(I)*c0.sdcd+2*atan((et(I).*(x+q(I)*c0.cd)+x.*(c2.r(I)+x)*c0.sd)./(xi(I).*(c2.r(I)+x)*c0.cd)));    
   
    ai3=(c2.y*c0.cd./rd-c2.ale+c0.sd.*log(rd))./c0.cdcd;
    ak1=xi.*(d11-c2.y11.*c0.sd)/c0.cd;
    ak3=(q.*c2.y11-c2.y.*d11)/c0.cd;
    aj3=(ak1-aj2*c0.sd)/c0.cd;
    aj6=(ak3-aj5*c0.sd)/c0.cd;
else
    disp('Check verticals ... this part untested so far')
    rd2=rd.*rd;
    ai3=(et./rd+c2.y.*q./rd2-c2.ale)/2; 
    % keyboard
    % ai4=xi*c2.y./rd2/2; % Problem for vertical? EMB
    % ak1=xi*q./rd.*d11;   
    ai4=xi.*c2.y./rd2/2; % Fixed? EMB
    ak1=xi.*q./rd.*d11;     
    % ak3= sd./rd.*(c2.xi2.*d11-1);  % EMB
    ak3= c0.sd./rd.*(c2.xi2.*d11-1);
    % keyboard
    % aj3=-xi./rd2*(c2.q2.*d11-1/2); % EMB
    % aj6=-c2.y./rd2*(c2.xi2.*d11-1/2);% EMB
    aj3=-xi./rd2.*(c2.q2.*d11-1/2);% EMB
    aj6=-c2.y./rd2.*(c2.xi2.*d11-1/2);% EMB
end

xy=xi.*c2.y11;
ai1=-xi./rd*c0.cd-ai4*c0.sd;
ai2= log(rd)+ai3*c0.sd;
ak2= 1./c2.r+ak3*c0.sd;
ak4= xy*c0.cd-ak1*c0.sd;
aj1= aj5*c0.cd-aj6*c0.sd;
aj4=-xy-aj2*c0.cd+aj3*c0.sd;

qx=q.*c2.x11;
qy=q.*c2.y11;

du(1,:)=disl1*(-xi.*qy-c2.tt-c0.alp3*ai1*c0.sd) + ...
		disl2*(-q./c2.r+c0.alp3*ai3*c0.sdcd) + ...
		disl3*(q.*qy-c0.alp3.*ai3.*c0.sdsd);
du(2,:)=disl1*(-q./c2.r+c0.alp3*c2.y./rd*c0.sd) + ...
		disl2*(-et.*qx-c2.tt-c0.alp3*xi./rd*c0.sdcd) + ...
		disl3*(q.*qx+c0.alp3.*xi./rd.*c0.sdsd);
du(3,:)=disl1*(q.*qy-c0.alp3*ai2*c0.sd) + ...
		disl2*(q.*qx+c0.alp3*ai4*c0.sdcd) + ...
		disl3*(et.*qx+xi.*qy-c2.tt-c0.alp3.*ai4.*c0.sdsd);
du(4,:)=disl1*(c2.xi2.*q.*c2.y32-c0.alp3*aj1*c0.sd) + ...
		disl2*(xi.*q./c2.r3+c0.alp3*aj4*c0.sdcd) + ...
		disl3*(-xi.*c2.q2.*c2.y32-c0.alp3.*aj4.*c0.sdsd);
du(5,:)=disl1*(xi.*q./c2.r3-c0.alp3*aj2*c0.sd) + ...
		disl2*(et.*q./c2.r3+qy+c0.alp3*aj5*c0.sdcd) + ...
		disl3*(-c2.q2./c2.r3-c0.alp3.*aj5.*c0.sdsd);
du(6,:)=disl1*(-xi.*c2.q2.*c2.y32-c0.alp3*aj3*c0.sd) + ...
		disl2*(-c2.q2./c2.r3+c0.alp3*aj6*c0.sdcd) + ...
		disl3*(q.*c2.q2.*c2.y32-c0.alp3.*aj6.*c0.sdsd);
du(7,:)=disl1*(-xi.*c2.fy-c2.d.*c2.x11+c0.alp3*(xy+aj4)*c0.sd) + ...
		disl2*(-c2.ey+c0.alp3*aj1*c0.sdcd) + ...
		disl3*(q.*c2.fy-c0.alp3.*aj1.*c0.sdsd);
du(8,:)=disl1*(-c2.ey+c0.alp3*(1./c2.r+aj5)*c0.sd) + ...
		disl2*(-et.*c2.gy-xy*c0.sd+c0.alp3*aj2*c0.sdcd) + ...
		disl3*(q.*c2.gy-c0.alp3.*aj2.*c0.sdsd);
du(9,:)=disl1*(q.*c2.fy-c0.alp3*(qy-aj6)*c0.sd) + ...
		disl2*(q.*c2.gy+c0.alp3*aj3*c0.sdcd) + ...
		disl3*(-q.*c2.hy-c0.alp3.*aj3.*c0.sdsd);
du(10,:)=disl1*(-xi.*c2.fz-c2.y.*c2.x11+c0.alp3*ak1*c0.sd) + ...
		 disl2*(-c2.ez-c0.alp3*ak3*c0.sdcd) + ...
		 disl3*(q.*c2.fz+c0.alp3.*ak3.*c0.sdsd);
du(11,:)=disl1*(-c2.ez+c0.alp3*c2.y.*d11*c0.sd) + ...
		 disl2*(-et.*c2.gz-xy*c0.cd-c0.alp3*xi.*d11*c0.sdcd) + ...
		 disl3*(q.*c2.gz+c0.alp3.*xi.*d11.*c0.sdsd);
du(12,:)=disl1*(q.*c2.fz+c0.alp3*ak2*c0.sd) + ...
		 disl2*(q.*c2.gz-c0.alp3*ak4*c0.sdcd) + ...
		 disl3*(-q.*c2.hz+c0.alp3.*ak4.*c0.sdsd);

du=du/(2*pi);

function du=uc(c0,c2,xi,et,q,z,disl1,disl2,disl3)
%Displacement and strain at depth (part-c) due to buried finite fault in a
%semiinfinite medium

c=c2.d+z;
x53=(8*c2.r2+9*c2.r.*xi+3*c2.xi2).*c2.x11.*c2.x11.*c2.x11./c2.r2;
y53=(8*c2.r2+9*c2.r.*et+3*c2.et2).*c2.y11.*c2.y11.*c2.y11./c2.r2;
h=q*c0.cd-z;
z32=c0.sd./c2.r3-h.*c2.y32;
z53=3*c0.sd./c2.r5-h.*y53;
y0=c2.y11-c2.xi2.*c2.y32;
z0=z32-c2.xi2.*z53;
ppy=c0.cd./c2.r3+q.*c2.y32*c0.sd;
ppz=c0.sd./c2.r3-q.*c2.y32*c0.cd;
qq=z.*c2.y32+z32+z0;
qqy=3*c.*c2.d./c2.r5-qq*c0.sd;
qqz=3*c.*c2.y./c2.r5-qq*c0.cd+q.*c2.y32;
xy=xi.*c2.y11;
%qx=q.*c2.x11;
qy=q.*c2.y11;
qr=3*q./c2.r5;
%cqx=c.*q.*x53;
cdr=(c+c2.d)./c2.r3;
yy0=c2.y./c2.r3-y0*c0.cd;

du(1,:)=disl1*(c0.alp4*xy*c0.cd-c0.alp5*xi.*q.*z32) + ...
		disl2*(c0.alp4*c0.cd./c2.r-qy*c0.sd-c0.alp5*c.*q./c2.r3) + ...
		disl3*(-c0.alp4*(c0.sd./c2.r+qy.*c0.cd)-c0.alp5*(z.*c2.y11-c2.q2.*z32));
du(2,:)=disl1*(c0.alp4*(c0.cd./c2.r+2*qy.*c0.sd)-c0.alp5*c.*q./c2.r3) + ...
		disl2*(c0.alp4*c2.y.*c2.x11-c0.alp5*c.*et.*q.*c2.x32) + ...
		disl3*(c0.alp4*2*xy.*c0.sd+c2.d.*c2.x11-c0.alp5*c.*(c2.x11-c2.q2.*c2.x32));
du(3,:)=disl1*(c0.alp4*qy.*c0.cd-c0.alp5*(c.*et./c2.r3-z.*c2.y11+c2.xi2.*z32)) + ...
		disl2*(-c2.d.*c2.x11-xy*c0.sd-c0.alp5*c.*(c2.x11-c2.q2.*c2.x32)) + ...
		disl3*(c0.alp4*(c2.y.*c2.x11+xy*c0.cd)+c0.alp5*q.*(c.*et.*c2.x32+xi.*z32));
du(4,:)=disl1*(c0.alp4*y0*c0.cd-c0.alp5*q.*z0) + ...
		disl2*(-c0.alp4*xi./c2.r3*c0.cd+c0.alp5*c.*xi.*qr+xi.*q.*c2.y32.*c0.sd) + ...
		disl3*(c0.alp4*xi./c2.r3*c0.sd+xi.*q.*c2.y32.*c0.cd+c0.alp5*xi.*(3*c.*et./c2.r5-2*z32-z0));
du(5,:)=disl1*(-c0.alp4*xi.*(c0.cd./c2.r3+2*q.*c2.y32*c0.sd)+c0.alp5*c.*xi.*qr) + ...
		disl2*(-c0.alp4*c2.y./c2.r3+c0.alp5*c.*et.*qr) + ...
		disl3*(c0.alp4*2*y0*c0.sd-c2.d./c2.r3+c0.alp5*c./c2.r3.*(1-3.*c2.q2./c2.r2));
du(6,:)=disl1*(-c0.alp4*xi.*q.*c2.y32*c0.cd+c0.alp5*xi.*(3*c.*et./c2.r5-qq)) + ...
		disl2*(c2.d./c2.r3-y0*c0.sd+c0.alp5*c./c2.r3.*(1-3*c2.q2./c2.r2)) + ...
		disl3*(-c0.alp4*yy0-c0.alp5*(c.*et.*qr-q.*z0));
du(7,:)=disl1*(-c0.alp4*xi.*ppy*c0.cd-c0.alp5*xi.*qqy) + ...
		disl2*(-c0.alp4*et./c2.r3+y0.*c0.sdsd-c0.alp5.*(cdr*c0.sd-c.*c2.y.*qr)) + ...
		disl3*(c0.alp4*(q./c2.r3+y0.*c0.sdcd)+c0.alp5.*(z./c2.r3.*c0.cd+c.*c2.d.*qr-q.*z0.*c0.sd));
du(8,:)=disl1*(c0.alp4*2*(c2.d./c2.r3-y0*c0.sd)*c0.sd-c2.y./c2.r3*c0.cd-c0.alp5*(cdr*c0.sd-et./c2.r3-c.*c2.y.*qr)) + ...
		disl2*(c0.alp4*(c2.x11-c2.y.*c2.y.*c2.x32)-c0.alp5.*c.*((c2.d+2*q*c0.cd).*c2.x32-c2.y.*et.*q.*x53)) + ...
		disl3*(-c0.alp4*2*xi.*ppy.*c0.sd-c2.y.*c2.d.*c2.x32+c0.alp5*c.*((c2.y+2.*q.*c0.sd).*c2.x32-c2.y.*c2.q2.*x53));
du(9,:)=disl1*(-c0.alp4*q./c2.r3+yy0*c0.sd+c0.alp5*(cdr*c0.cd+c.*c2.d.*qr-(y0*c0.cd+q.*z0)*c0.sd)) + ...
		disl2*(xi.*ppy*c0.sd+c2.y.*c2.d.*c2.x32+c0.alp5*c.*((c2.y+2*q*c0.sd).*c2.x32-c2.y.*c2.q2.*x53)) + ...
		disl3*(-c0.alp4*(xi.*ppy*c0.cd-c2.x11+c2.y.*c2.y.*c2.x32)+c0.alp5*(c.*((c2.d+2*q.*c0.cd).*c2.x32-c2.y.*et.*q.*x53)+xi.*qqy));
du(10,:)=disl1*(c0.alp4*xi.*ppz*c0.cd-c0.alp5*xi.*qqz) + ...
		 disl2*(-q./c2.r3+y0*c0.sdcd-c0.alp5*(cdr.*c0.cd+c.*c2.d.*qr)) + ...
		 disl3*(-et./c2.r3+y0.*c0.cdcd-c0.alp5*(z./c2.r3.*c0.sd-c.*c2.y.*qr-y0.*c0.sdsd+q.*z0.*c0.cd));
du(11,:)=disl1*(c0.alp4*2*(c2.y./c2.r3-y0*c0.cd)*c0.sd+c2.d./c2.r3*c0.cd-c0.alp5*(cdr*c0.cd+c.*c2.d.*qr)) + ...
		 disl2*(c0.alp4.*c2.y.*c2.d.*c2.x32-c0.alp5*c.*((c2.y-2*q*c0.sd).*c2.x32+c2.d.*et.*q.*x53)) + ...
		 disl3*(c0.alp4*2*xi.*ppz.*c0.sd-c2.x11+c2.d.*c2.d.*c2.x32-c0.alp5*c.*((c2.d-2*q.*c0.cd).*c2.x32-c2.d.*c2.q2.*x53));
du(12,:)=disl1*(yy0*c0.cd-c0.alp5*(cdr*c0.sd-c.*c2.y.*qr-y0*c0.sdsd+q.*z0*c0.cd)) + ...
		 disl2*(-xi.*ppz.*c0.sd+c2.x11-c2.d.*c2.d.*c2.x32-c0.alp5.*c.*((c2.d-2*q*c0.cd).*c2.x32-c2.d.*c2.q2.*x53)) + ...
		 disl3*(c0.alp4*(xi.*ppz*c0.cd+c2.y.*c2.d.*c2.x32)+c0.alp5*(c.*((c2.y-2.*q.*c0.sd).*c2.x32+c2.d.*et.*q.*x53)+xi.*qqz));
     
du=du/(2*pi);
    
