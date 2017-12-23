function VonMisesMat ;
%******************************************************
%*                                                    *
%*                     C.E. 529b                      *
%*                                                    *
%*               RADIAL RETURN METHOD           I     *
%*     STATE DETERMINATION FOR VON MISES MATERIAL     *
%*                 1-D FORMULATION                    *
%*                                                    *
%*            Von Mises Elasto-plasticity             *
%*                                                    *
%******************************************************
%
%     PARAMETERS TO BE SAVED
%
%     akpold  -   yield stress
%     eppold  -   Strain at last step
%     tsgold  -   stress at last step
%     yscold  -   yield surface center
%
%     DEFINITIONS OF PARAMETERS
%
%     mstat   -   material status flag (1 elastic, 2 plastic)
%     eppnew  -   total plastic strain
%     h       -   strain hardening modulus
%     bet     -   flag for isotropic/kinematic hardening
%                  (1 isotropic, 0 kinematic, other combination)
%     akap    -   starting yield stress from uniaxial test
%
%******************************************************
      global iunit
      global e et bet eppold tsgold ysc yscold akpold mstat
      global istart istep isubin jter eppnew depp cep
      global fid ihard epult epold
%
%     NUMBER OF STRESS COMPONENTS (1-D PROBLEM).
      ncomp=1 ;
%
%     THERE IS A PROBLEM INVOLVING REAL UNLOADING. TO INSURE
%     CONVERGENCE WHEN USING LARGE LOAD STEPS WITH UNLOADING IT
%     SEEMS TO BE NEC. TO USE THE ELASTIC MODULUS E DURING THE
%     FIRST UNLOADING ITER. IF THE MODULUS ET IS USED VERY LARGE
%     CORRECTIONS WILL BE MADE WHICH WILL ARTIFICIALLY THROW
%     THE SYSTEM OUT-OF-BALANCE. THUS ALWAYS CONSIDER THE SYSTEM
%     TO BE ELASTIC IF THE YIELD COND. IS EXACTLY SAT. THEN
%     THE ELEMENT IS CONSIDERED TO BE PLASTIC ONLY IF THE COMPUTED
%     STRESS POINT IS OUTSIDE THE YIELD SURFACE.
%     DEFINE ERROR TOLERANCE ALLOWABLE IN SATISFYING THE YIELD COND.
%     TO PREVENT AN ARTIFICIAL UNLOADING CONDITION STATE PRINT
%     DUE TO NUMERICAL ERROR IN SATISFYING YIELD CONDITION WHEN IN
%     FACT STRESS PT. LIES ON THE YIELD SURFACE.
%
      erry=0.001 ;
%
%     DEFINE THE MATERIAL CONSTANT
      eee=e ;
%
%     DEFINE THE PLASTIC MODULUS.
      if (ihard==1) 
    %     DEFINE LINEAR HARDENING MODULUS. 
          h=eee*et/(eee-et) ;
          
          
      end
      if (ihard==2)
          %------###############################
              hzero=(eee*et)/(eee-et);
              h=hzero*(1-et/epult);
              hzero=(eee*et)/(eee-et);
              h=hzero*(1-ep-0.5*ep^2/epult);
              Yel=akp+h %  from the table
              
          %------###############################
    %     DEFINE NONLINEAR HARDENING MODULUS USING
    %     CURVE DESCRIPTION MODEL AS USED IN ABAQUS.
    %     SEE ATTACHED.
      end
%
%     AKAP IS THE YIELD STRESS FROM A UNIAXIAL TENSION TEST.
      akp=akpold ;
%
%     ysc IS THE YIELD STRESS FROM A UNIAXIAL TENSION TEST.
      ysc=yscold ;
%
%     THE STATEMENTS EXECUTED ON THE FIRST ENTRY TO THIS ROUTINE.
      if (istart==0) 
    %
    %     SET THE MATERIAL STATUS FLAG TO INITIAL ELASTIC VALUE.
          mstat=1 ;
    %
    %     SET THE AXIAL STRAIN AND STRESS FOR THE REFERENCE STATE TO
    %     ZERO INITIALLY.
          tsgold=0.0 ;
          eppold=0.0 ;
          yscold=0.0 ;
          akpold=akp ;
    %
    %     DEFINE INITIAL YIELD SURFACE CENTER.
          ysc=0.0 ;
    %
          istart=1 ;
      end
%
%     AKAP IS THE YIELD STRESS FROM A SIMPLE TENSION TEST.
%     BKAP IS THE KAPPA TERM IN THE VON MISES YIELD CRITERIA.
      bkp=akp/(3.0^0.5) ;
%
%     YFRAD2 IS THE YIELD FUNCTION SQUARED.
      yfrad2=bkp^2 ;
%
%     COMPUTE AXIAL TRIAL STRESS INCREMENT ASSUMING
%     ELASTIC BEHAVIOR.
      dsig=eee*depp ;
%
%     COMPUTE TOTAL AXIAL STRESS ASSUMING ELASTIC MATERIAL(TRIAL STRESS).
      tsgtrl=tsgold+dsig ;
%
%     FOR THE 1-D STRESS APPLICATIONS, COMPUTE HYDROSTATIC STRESS.
      tsigm=tsgtrl/3.0 ;
%
%     COMPUTE AXIAL DEVIATORIC TRIAL STRESS.
      dvsig=tsgtrl-tsigm ;
%
%     DEFINE DIFF. BET. TRIAL DEVIATORIC AXIAL STRESS AND
%     YIELD SURFACE CENTER STRESS.
      t=dvsig-ysc ;
%
%     DEFINE YIELD FUNCTION VALUE.
%     VON MISES YIELD FUNCTION.
%     NOTE THAT THIS YIELD FUNCTION IS DEFINED IN TERMS OF DEV. STRESS.
%     EQ. 5 BAR IN LECTURE 9 NOTES
      yldf=(3.0/4.0)*(t*t)-yfrad2 ;
%
%     OUTPUT OF PARAMETERS.
      nprnt=0 ;
      if (nprnt==1)
    %
          fprintf(fid,'INTERMEDIATE CALCULATIONS\n') ;
          fprintf(fid,'istep,isubin,jter=\n',istep,isubin,jter) ;
          fprintf(fid,'eppold    tsgold    yscold    akpold\n') ;
          fprintf(fid,'eppnew    depp      dsig      tsgtrl\n') ;
          fprintf(fid,'tsigm     dvsig     t          akp\n') ;
          fprintf(fid,'yfrad2    yldf      cep\n') ;
    fprintf(fid,'%6.2f %6.2f %6.2f %6.2f\n',eppold,tsgold,yscold,akpold);
    %      akpold) ;
    fprintf(fid,'%6.2f %6.2f %6.2f %6.2f\n',eppnew,depp,dsig,tsgtrl);
    %      tsgtrl) ;
    fprintf(fid,'%6.2f %6.2f %6.2f %6.2f\n',tsigm,dvsig,t,akp);
    %      akp) ;
    fprintf(fid,'%6.2f %6.2f %6.2f %6.2f\n',yfrad2,yldf,cep);
    %      cep) ;
    %
      end
%
%     CHECK TO SEE IF TOTAL STRESS VALUE IS INSIDE OR OUTSIDE YIELD
%     SURFACE.
      mstat=1 ;
%
%     SPECIAL PROCEDURE TO CODE STATE FOR PRINTOUT AS PLASTIC
%     EVEN WHEN SMALL ARTIFICIAL NUMERICAL UNLOADING MAY HAVE
%     OCCURRED. 
      errry=-erry*yfrad2 ;
      if (yldf>errry) 
        mstat=2 ;
      end
      errry=-errry ;
%
%     IF STRESS POINT IS ONLY AN EPSILON ABOVE THE YIELD SURFACE,
%     CONTINUE TO TREAT IT AS ELASTIC.
      if (yldf<=errry) 
    %
    %     THE STRESS POINT LIES INSIDE OR ON THE YIELD SURFACE.
    %     THE ELASTIC CASE.
    %
    %     DEFINE CURRENT TOTAL STRESS AS OLD STRESS FOR NEXT STEP
    %     WITH ELASTIC STRESS INCREMENT.
          tsg=tsgtrl ;
    %
    %     DEFINE AXIAL ELASTIC-PLASTIC MATERIAL MATRIX AS ELASTIC MODULUS.
    %     THIS IS FOR USE IN FORMING INCREMENTAL STIFF. MATRIX. THIS IS
    %     [CEP].
          cep=eee ;
          ep=epold ;
    %
      end
%
      if (yldf>errry) 
    %
    %     THE TENTATIVE STRESS POINT LIES OUTSIDE THE YIELD SURFACE - THERE
    %     IS SOME PLASTIC STRAIN.
    %
    %     DEFINE DIST. FROM CENTER OF YIELD SURFACE TO TRIAL STRESS POINT.
    %     EQ. 35 BAR IN LECTURE 9 NOTES
          dtrial=3.*(yldf+yfrad2) ;
          dtrial=sqrt(dtrial) ;
    %
    %     DEFINE DISTANCE BETWEEN CENTER OF YIELD SURFACE AND YIELD SURFACE.
    %     EQ. 36 BAR IN LECTURE 9 NOTES
          dysurf=3.*yfrad2 ;
          dysurf=sqrt(dysurf) ;
    %
    %     DEFINE INCREMENTAL PLASTIC STRETCHING INTEGRAL.
    %     EQ. 36 BAR IN LECTURE 9 NOTES
          con=h+eee ;
          dep=(dtrial-dysurf)/con ;
          ep=epold+dep ;
    %
    %     UPDATE THE AXIAL STRESS VALUES AND STORE IN AXIAL
    %     STRESS STATE STORAGE AREA. RETURN TRIAL STRESS TO
    %     YIELD SURFACE.
    %     EQ. 43 BAR IN LECTURE 9 NOTES
          atsig=abs(tsgtrl) ;
          tsg=tsgtrl-eee*dep*tsgtrl/atsig ;
    %
    %     UPDATE THE YIELD SURFACE CENTER.
    %     EQ. 44 BAR IN LECTURE 9 NOTES
          sqrt23=2./3. ;
          sqrt23=sqrt(sqrt23) ;
          ysc=yscold+sqrt23*(1.0-bet)*h*dep*tsgtrl/atsig ;
    %
    %     UPDATE THE AXIAL YIELD SURFACE RADIUS-ACTUALLY THE EFFECTIVE
    %     TENSION YIELD STRESS.
          akp=akpold+bet*h*dep ;
    %
    %     ANOTHER EQUIVALENT WAY WHEN BET EQUAL 1.,
    %     TO INSURE THAT AT NEXT STEP YIELD COND. IS EXACTLY SAT.
    %     THUS AVOIDING THE REAL UNLOADING PROBLEM.-FORCING THE ALGORITM
    %     AT THE START OF THE NEXT LOAD STEP TO USE THE ELASTIC MODULUS E.
          if (bet==1.0)
            akp=abs(tsg) ;
          end
    %
    %     DEFINE AXIAL ELASTIC-PLASTIC MATERIAL MATRIX.
    %     THIS IS FOR USE IN FORMING INCREMENTAL STIFF. MATRIX.
    %     THIS IS [CEP].
          if (ihard==1)
            cep=et ;
          end
          if (ihard==2)
              %-----------#################################
              cep=eee-eee^2/(eee-h)
%               hzero=(eee*et)/(eee-et);
%               h=hzero*(1-ep-0.5*ep^2/epult);
%               cep=akp+h %  from the table
%                             
    %
              %------##########################################
    %     UPDATE NONLINEAR HARDENING MODULUS USING
    %     CURVE DESCRIPTION MODEL AS USED IN ABAQUS.
    %     SEE ATTACHED.
    %     COMPUTE CEP.
          end
    %
      end
%
%     SAVE THE OLD VALUES.
      eppold=eppnew ;
      tsgold=tsg  ;
      yscold=ysc ;
      akpold=akp ;
      epold=ep ;
%
      nprnt=0 ;
      if (nprnt==1)
%
          fprintf(fid,'VomMisesMat-FINAL VALUES\n') ;
          fprintf(fid,'eppnew     tsg        ysc        akp        mstat\n') ;
          fprintf(fid,'%6.2f %6.2f %6.2f %6.2f %6.2f\n',eppnew,tsg,ysc,akp,mstat);
    %      akp,mstat) ;
    %
      end
%
      end

