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
      fid=fopen('hw8.out','w') ;
%
%     INPUT DATA
%
%     NUMBER OF LOAD STEPS
      nstep=17  ;
%     NUMBER OF ITERATIONS
      niter=10  ;
%     STRAIN INCREMENT
      dels=.25 ;
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
      uold=0.0 ;
      tdtold=1.0 ;
      udtold=1.0 ;
      S=0.0 ;
%
%     MAIN INCREMENTAL LOOP
      for istep=1:nstep
%
      S=S+dels ;
%
      t=told ;
      u=uold ;
%
%     MAIN ITERATIVE LOOP
      for iter=1:niter
%
%     FORM JACOBIAN OR COEFICIENT MATRIX (fij)
      ak(1,1)=tdtold ;
      ak(1,2)=udtold ;
      ak(2,1)=1.0 ;
      ak(2,2)=-2.0+2.0*u ;
%
%     INVERT COEFFICIENT MATRIX (fij)
%     3X3 MATRIX VERSION IN SUBPROGRAM
      det=ak(1,1)*ak(2,2)-ak(1,2)*ak(2,1) ;
      akinv(1,1)=ak(2,2)/det ;
      akinv(1,2)=-ak(1,2)/det ;
      akinv(2,1)=-ak(2,1)/det ;
      akinv(2,2)=ak(1,1)/det ;
%
%     FORM RESIDUAL (fi)
      re(1)= dels-tdtold*(t-told)-udtold*(u-uold) ;
      re(2)= -t+2*u-u^2 ;
%
%     COMPUTE ITERATIVE CORRECTION.
      for i=1:2
      sol(i)=0.0 ;
      for j=1:2
      sol(i)=sol(i)+akinv(i,j)*re(j) ;
      end
      end
      delt=sol(1) ;
      delu=sol(2) ;
%
%     UPDATE SOLUTION
      t=t+delt ;
      u=u+delu ;
%
%
      %fprinf(fid,'ITERATIVE VALUES\n') ;
      fprintf(fid,'iter   delt   delu  t   u\n') ;
      fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f\n',iter,delt,delu,t,u) ;
      fprintf(fid,'re1   re2   sol1   sol2\n') ;
      fprintf(fid,'%6.2f %6.2f %6.2f %6.2f\n',re(1),re(2),sol(1),sol(2)) ;
%
%     CONVERGENCE CHECK
%     error=sqrt(re(1)**2+re(2)**2) ;
%     if (error<=errtol) 
%     end
      end
%
%     UPDATE TANGENTS
      tdtold=(t-told)/dels ;
      udtold=(u-uold)/dels ;
      told=t ;
      uold=u ;
%
%     OUTPUT FINAL RESULTS
      %fprinf(fid,'FINAL VALUES\n') ;
      fprintf(fid,'istep    t         u\n') ;
      fprintf(fid,'%6.2f %6.2f %6.2f\n',istep,t,u) ;
%
      end
%
%     CLOSE OUTPUT UNIT
      fclose(fid) ;
%