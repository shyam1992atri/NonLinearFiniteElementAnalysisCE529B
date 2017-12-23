%******************************************************
%*                                                    *
%*                     C.E. 529b                      *
%*                                                    *
%*                CONTINUATION ANALYSIS         I     *
%*                 RIKS-WEMPNER METHOD                *
%*                                                    *
%*                 EXAMPLE PROBLEM                    *
%*                 1-D FORMULATION                    *
%*                                                    *
%*            SOLVING NONLINEAR EQUATION              *
%*                  2*u-u**2=t                        *
%*               u =displacement                      *
%*                   t=force                          *
%*                                                    *
%*              SOLUTION HAS LIMIT POINT              *
%*                                                    *
%*            Matlab vers. riks.wemp.m                *
%*                                                    *
%*           TO BE USED AS BASIS FOR H.W. #8          *
%*            3x3 MATRIX INVERSION PROGRAM            * 
%*             (NEEDED FOR HW 8) INCLUDED             *               
%*                                                    *
%*               Output File: hw8.out                 *
%*                                                    *
%******************************************************
%
%     DEFINITIONS OF PARAMETERS
%
%     t,told  -   psuedo-time parameter
%     u,uold  -   displacement parameter
%     s       -   arc length parameter
%     dels    -   delta s
%     tdtold  -   dt/ds
%     udtold  -   du/ds
%     u,uold  -   displacement parameter
%     ak      -   jacobian, tangent stiffness matrix (fij)
%     re      -   residual (fi)
%
%******************************************************
%
%     OPEN OUTPUT UNIT
      fid=fopen('hw8for3.out','w') ;
%
%     INPUT DATA
%
%     NUMBER OF LOAD STEPS
      nstep=17  ;
%     NUMBER OF ITERATIONS
      niter=10  ;
%     STRAIN INCREMENT
      dels=.2 ;
%     ERROR TOLERANCE
      errtol=1.E-06 ;
%
%     OUTPUT OF HEADINGS
      fprintf(fid,'RIKS-WEMPNER METHOD\n') ;
      fprintf(fid,'\n') ;
      fprintf(fid,'INPUT DATA\n') ;
      fprintf(fid,'nstep      niter        dels      errtol\n') ;
      fprintf(fid,'%6.2f %6.2f %6.2f %6.2f\n',nstep,niter,dels,errtol) ;
%
      told=0.0 ;
      u1old=0.0 ;
      u2old=0.0 ;
      tdtold=1.0 ;
      u1dtold=1.0 ;
      u2dtold=1.0 ;
      S=0.0 ;
%
%     MAIN INCREMENTAL LOOP
      for istep=1:nstep
%
      S=S+dels ;
%
      t=told ;
      u1=u1old ;
      u2=u2old ;
%
%     MAIN ITERATIVE LOOP
      for iter=1:niter
          figure(1)
          plot(u1,t,'.-')
          hold on
          figure(2)
          plot(u2,t,'.-')
          hold on
%
%     FORM JACOBIAN OR COEFICIENT MATRIX (fij)
      ak(1,1)=tdtold ;
      ak(1,2)=u1dtold ;
      ak(1,3)=u2dtold ;
      ak(2,1)=-0.3 ;
      ak(2,2)=800+12*u1^2+16*u2+u2^2 ;
      ak(2,3)=8*u1*u2;
      ak(3,1)= -1.5;
      ak(3,2)=4*u1*(4+2*u2) ;
      ak(3,3)=800+24*u2+12*u2^2+4*u1^2;
%
%     INVERT COEFFICIENT MATRIX (fij)
%     3X3 MATRIX VERSION IN SUBPROGRAM
      [ajacinv,det,c]=jacinv(ak)
%
%     FORM RESIDUAL (fi)
      re(1)= dels-tdtold*(t-told)-u1dtold*(u1-u1old)-u2dtold*(u2-u2old) ;
      re(2)= (1/(2*104^1.5))*((800+4*u1^2)*u1 +(4*u1)*(4+u2)*u2)+0.3*t;
      re(3)= (1/(2*104^1.5))*((32+24*u2+4*u2^2)*u2+2*u1*(4+2*u2)*u1)+1.5*t ;
%
%     COMPUTE ITERATIVE CORRECTION.
      for i=1:3
      sol(i)=0.0 ;
      for j=1:3
      sol(i)=sol(i)+ajacinv(i,j)*re(j) ;
      end
      end
      delt=sol(1) ;
      delu1=sol(2) ;
      delu2=sol(3) ;
%
%     UPDATE SOLUTION
      t=t+delt ;
      u1=u1+delu1 ;
      u2=u2+delu2 ;
%
%
      %fprinf(fid,'ITERATIVE VALUES\n') ;
      fprintf(fid,'iter   delt   delu1 deltu2  t   u1 u2 \n') ;
      fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n',iter,delt,delu1,delu2,t,u1,u2) ;
      fprintf(fid,'re1   re2 re3  sol1  sol2 sol3\n') ;
      fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',re(1),re(2),re(3),sol(1),sol(2),sol(3)) ;
%
%     CONVERGENCE CHECK
%     error=sqrt(re(1)**2+re(2)**2) ;
%     if (error<=errtol) 
%     end
      end
%
%     UPDATE TANGENTS
      tdtold=(t-told)/dels ;
      u1dtold=(u1-u1old)/dels ;
      u2dtold=(u2-u2old)/dels ;
      told=t ;
      u1old=u1 ;
      u2old=u2;
%
%     OUTPUT FINAL RESULTS
      %fprinf(fid,'FINAL VALUES\n') ;
      fprintf(fid,'istep    t    u1  u2  \n') ;
      fprintf(fid,'%6.2f %6.2f %6.2f %6.2f\n',istep,t,u1,u2) ;
%
      end
%
%     CLOSE OUTPUT UNIT
      fclose(fid) ;
%