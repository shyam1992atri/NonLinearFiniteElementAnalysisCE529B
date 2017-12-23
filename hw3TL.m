
function TL
%******************************************************
%*                                                    *
%*                     C.E. 529B                      *
%*                                                    *
%*               TOTAL LAGRANGIAN METHOD              *
%*  1-D FINITE ELEMENT PROGRAM - Matlab - HW3TotLag.m *
%*                                                    *
%*                      H.W. #3                       *
%*                                                    *
%*     GEOMETRICALLY NONLINEAR FORMULATION            *
%*     MATERIALLY NONLINEAR BUT ELASTIC               *
%*                                                    *
%*     ONE LINEAR ELEMENT IN MODEL                    *
%*     ONLY ONE LOAD STEP CONSIDERED.                 *
%*                                                    *
%*     EITHER CONCENTRATED FORCE F OR                 *
%*     PRESSURE p, BUT NOT BOTH SIMULTANEOUSLYi       *
%*                                                    *
%*     PLEASE NOTE THAT THERE IS ONLY 1 DISPLACEMENT  *
%*     DOF. DO NOT USE MATRICES LIKE [K]2x2.          * 
%*                                                    *
%*                                                    *
%*     Output File: hw3.txt                           *
%*                                                    *
%******************************************************
%
%open output unit
fid = fopen('hw3.out','w');
%

%input data
%number of iterations
niter=15 
%
%write titles
fprintf(fid,'\n');
fprintf(fid,'Total Lagrangian Method\n');
fprintf(fid,'Bathe Problem\n');
fprintf(fid,'C=E Material Model\n');
fprintf(fid,'Total Stress Method\n');
%
%initializations
%displacement
      u=0.0 ;
      eta=0.0;
      s=0.0;
% piola kichhoff 2nd stress
      pk2=0.0 ;
% element length in reference coordinate system
      h=2. ;
% bar crossectional undeformed area
      area=4. ;
      adef=area;
% modulus of elasticity
      e=0.1346E+08 ;
%
% set applied loads ;
      iload=2 ;
      if (iload==1)
        % applied concentrated force case
        % set pressure p to zero
          f=0.1000E+09 ;
          p=0.0 ;
      end
      if (iload==2)
        % applied pressure case
        % set force f to zero
          f=0.1000E+09 ;
          p=f/area ;
          %tau=p;
          %f=0.0 ;
      end
      
      if (iload==3)
        % applied pressure case
        % set force f to zero
          f=0.15E8 ;
          p=f/area ;
          
      end
%
% Cauchy stress
      tau=0.0 ;
%
%output of headings
fprintf(fid,'  H         Area           E           F           P\n') ;
fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f\n',h,area,e,f,p) ;
fprintf(fid,'  iter   u    eta   adef   pk2    tau    res    stiff\n') ;

% fprintf(fid,'  H  Area  E  F\n') ;
% fprintf(fid,'%6.2f %6.2f %6.2f %6.2f\n',h,area,e,f) ;
% fprintf(fid,'iter   u    eta   pk2    res    stiff\n') ;
%
%main iterative loop
      for iter=1:niter
        %
        %define displacement gradients and
        %Green-Saint Venant strain (gam)
        alpha=u/h;
        %
        %define elastic modulus c 
        %if material is nonlinear c would be different
        if (iload==1 || iload==2)
            c=e ; 
        end
        if (iload==3)
            c=e/(1+alpha)^4;
            
        end
        %
        %define tangent stiffness matrix (stiff) and residual vector (res).
        %use formulation on p. 26-27 of class notes. note that there is
        %only 1 displacement variable u. Thus the stiffness matrix is
        %[K]1x1, etc.
              if (iload==1)
                % applied concentrated force case
%                 stiff=(1+alpha)^2*a*c/h;
%                 res=f-a*pk2*(1/h)*(1+alpha);
%                 eta=res/stiff;
%                 u=u+eta;
%                 gam=u/h+(0.5*(u/h)^2);
%                 pk2=gam*e;
% %                 s=(c/h+(c*u/h^2))*eta;
% %                 pk2=pk2+s;

%                 stiff=(area*c/h + area*pk2/h + u*(2+u/h));
%                 res=(f-area*pk2*(1+u/h));
%                 eta=res/stiff;
%                 s=(c/h+c*u/h^2)*eta;
%                 u=u+eta;
%                 pk2=pk2+s;
                  adef=area*(1/(1+alpha));
                  stiff=(((1+alpha)^2*(area*c/h))+area*pk2/h);
                  res=(f-(area*pk2*(1+alpha)));
                  eta=res/stiff;
                  u=u+eta;
                  gam=(u/h)+0.5*(u/h)^2;
                  pk2=gam*e;
                  tau=pk2*(1+alpha)^2;
                  
                
              end
              if (iload==2)
                % applied pressure case
                  %adef=adef*(1/(1+alpha));
%                   stiff=(((1+alpha)^2*(area*c/h))+area*pk2/h+((area*tau)/(1+alpha)^2)*(1/h));
%                   res=(adef*tau-(area*pk2*(1+alpha)));
%                   stiff=(((1+alpha)^2*(area*c/h))+area*pk2/h+((area*p)/(1+alpha)^2)*(1/h));
%                   res=(adef*p-(area*pk2*(1+alpha)));
%                   stiff=(((1+alpha)^2*(area*c/h))+area*pk2/h+((area*p)/(1+alpha)^2)*(1/h));
%                   res=(area*p-(area*pk2*(1+alpha)));
% 

% this section is working
%                 stiff=(((1+alpha)^2*(area*c/h))+area*pk2/h+((area*p)/(1+alpha)^2)*(1/h));
%                 res=(adef*p-(area*pk2*(1+alpha)));
%                   eta=res/stiff;
%                   u=u+eta;
%                   gam=(u/h)+0.5*(u/h)^2;
%                   pk2=gam*e;
%                   tau=pk2*(1+alpha)^2;
%---------------------------------------------
                adef=area*(1/(1+alpha));
           stiff=(((1+alpha)^2*(area*c/h))+area*pk2/h+((area*p)/(1+alpha)^2)*(1/h));
                res=(adef*p-(area*pk2*(1+alpha)));
                  eta=res/stiff;
                  u=u+eta;
                  gam=(u/h)+0.5*(u/h)^2;
                  pk2=gam*e;
                  tau=pk2*(1+alpha)^2;




              end
              if (iload==3)
                  
                  adef=area*(1/(1+alpha));
                  stiff=(((1+alpha)^2*(area*c/h))+area*pk2/h);
                  res=(f-(area*pk2*(1+alpha)));
                  eta=res/stiff;
                  u=u+eta;
                  gam=(u/h)+0.5*(u/h)^2;
                  pk2=gam*e/(1+alpha)^4;;
                  tau=pk2*(1+alpha)^2;
              end
        %
        %solve incremental stiffness equations (eta)
        %
        %update the solution
        %
        %Define strain (gam)
        %
        %update the stress in the model (pk2)
        %
        %compute Cauchy stress (tau)
        %
        %define deformed area (adef)
        %
        %output results
        fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',iter,u,eta,adef,pk2,tau,res,stiff) ;
        %fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f \n',iter,u,eta,pk2,res,stiff) ;
        %
        
        end
%close output unit
fclose(fid);

