%******************************************************
%*                                                    *
%*                     C.E. 529b                      *
%*                                                    *
%*         FINITE ELEMENT ELASTO-PLASTICITY           *
%*                 1-D FORMULATION                    *
%*                                                    *
%*                  Matlab Version                    *
%*                                                    *
%*                      H.W. #7                       *
%*                                                    *
%*           GEOMETRICALLY LINEAR FORMULATION         *
%*            MATERIALLY NONLINEAR INELASTIC          *
%*                                                    *
%*             ONE LINEAR ELEMENT IN MODEL            *
%*                                                    *
%*     Output File: hw7.out                           *
%*                                                    *
%******************************************************
      global iunit
      global e et bet eppold tsgold ysc yscold akpold mstat
      global istart istep isubin jter eppnew depp cep
      global fid ihard epult epold
%
%     OPEN OUTPUT UNIT
      fid=fopen('hw7.out','w') ;
%
%     INPUT DATA
%
%     YIELD STRESS FROM UNIAXIAL TEST.
      akap=4.0*10.^8 ; 
%     ELASTIC MODULUS.
      e=2.0*10.^11 ;
%     TANGENT MODULUS.
      et=2.0*10.^10 ;
%     PLASTIC MODULUS: 1: H=E*ET/(E-ET) ; 2: H=Ho-Ho*(ep/epult)
      ihard=2 ;
%     ULTIMATE PLASTIC STRAIN
      epult=0.02 ;
%     ISOTROPIC/KINEMATIC HARDENING PARAMETER..
      bet=1.0 ;
%     NUMBER OF LOAD SUBINCREMENTS
      nsubin=10 ;
%
%     INITIALIZE LOAD HISTORY
      initep=2 ;
      if (initep==1) 
    %     CONSTANT LOAD STEPS
    %     NUMBER OF CONSTANT LOAD STEPS
          nstep=17  ;
    %     CONSTANT LOAD INCREMENT
          delf=10.0*10^5 ;
          f(1)=0.0 ;
          for istep=1:nstep
            f(istep+1)=istep*delf
          end
      end
      if (initep==2) 
    %     KINEMATIC HARDENING TEST.
          nstep=9 ;
          bet=0.0 ;
          f(1)=0.0 ;
          f(2)=16.0*10^5 ;
          f(3)=32.0*10^5 ;
          f(4)=48.0*10^5 ;
          f(5)=-5.0*10^5 ;
          f(6)=-32.0*10^5 ;
          f(7)=-40.0*10^5 ;
          f(8)=-43.0*10^5 ;
          f(9)=16.0*10^5 ;
          f(10)=37.0*10^5 ;
      end
      if (initep==3)
    %     ISOTROPIC HARDENING TEST.
          nstep=7 ;
          bet=1.0 ;
          f(1)=0.0 ;
          f(2)=16.0*10^5 ;
          f(3)=32.0*10^5 ;
          f(4)=46.0*10^5 ;
          f(5)=-20.0*10^5 ;
          f(6)=-50.0*10^5 ;
          f(7)=30.0*10^5 ;
          f(8)=56.0*10^5 ;
      end
%
%     INITIALIZATION OF SOLUTION
%
%     INTIAL DISPLACEMENT
      u=0.0 ;
%     INTIAL RESIDUAL
      res=0.0 ;
%     INITIAL TOTAL STRAIN
      eppnew=0.0 ;
      displacement=[eppnew];
      newForce=0;
%     ELEMENT LENGTH
      al=1.0 ;
%     ELEMENT CROSSSECTIONAL AREA
      a=.01 ;
%     INITIAL STRESS 
      tsgold=0.0 ;
%     INITIAL UNIAXIAL YIELD STRESS
      akpold=akap ;
%     ERROR TOLERANCE
      errtol=.001 ;
%
%     STARTING PARAMETERS
%
%     STARTING INITIALIZATION PARAMETER for material
      istart=0 ;
%
%     OUTPUT OF HEADINGS
      fprintf(fid,'1-D ELASTO-PLASTICITY (RADIAL RETURN)\n') ;
      fprintf(fid,'\n') ;
      fprintf(fid,'INPUT DATA\n') ;
      fprintf(fid,'e         et        depp      akap      bet\n') ;
      fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f\n',e,et,depp,akpold,bet) ;
%
%     MAIN INCREMENTAL LOOP
      for istep=1:nstep
    %
    %     COMPUTE AXIAL STRAIN INCREMENT SINCE LAST STATE UPDATE.
          df=f(istep+1)-f(istep) ;
    %
    %     SUBINCREMENTAL LOAD LOOP
          for isubin=1:nsubin
        %
              if (nsubin==1)
                fnew=f(istep+1);
              end
              if (nsubin>1)
                ddf=df/nsubin ;
                fnew=f(istep)+isubin*ddf ;
              end
        %
        %     INITIAL ELASTIC PLASTIC MATERIAL MODULUS
              cep=e ;
        %
        %     MAX NUMBER OF ITERATIONS
              niter=30 ;
              resabs=10.0*10^8 ;
              jter=0 ;
        %
        %     MAIN ITERATIVE LOOP
              while jter<niter && resabs>errtol
                  jter=jter+1 ;
                  res=fnew-tsgold*a ;
                  resabs=abs(res) ;
                  ak=(a*cep)/al ;
                  du=res/ak ;
                  depp=du/al ;
                  eppnew=eppnew+depp ;
                  u=u+du ;
            %
            %     VON MISES MATERIAL STATE ROUTINE.
                  VonMisesMat ;
        %
              end
        %
              fprintf(fid,'FINAL VALUES\n') ;
              fprintf(fid,'istep      isubin     jter       u          fnew\n') ;
              fprintf(fid,'eppnew     tsg        ysc        akp        mstat\n') ;
        fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f\n',istep,isubin,jter,u,fnew) ;
        fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f\n',eppnew,tsgold,ysc,akpold,mstat) ;
        %
        displacement(end+1)=u;
        newForce(end+1)=fnew;
          end
          
      end
%
%     CLOSE OUTPUT UNIT
      plot(displacement,newForce,'-*')
      
      grid on
      fclose(fid) ;
%
%
      