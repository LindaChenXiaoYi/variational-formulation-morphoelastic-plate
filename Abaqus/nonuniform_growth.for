C =================================================================
C GROWTH NEO-HOOKEAN UMAT
C COMPRESSIBLE / NEARLY INCOMPRESSIBLE PENALTY FORMULATION
C
C A = F * G^(-1)
C State variables stored for post-processing:
C   SDV  1-6  : original variables
C   SDV  7-10 : F11,F13,F31,F33
C   SDV 11-14 : A11,A13,A31,A33
C   SDV 15    : PR = K*(JA-1)
C   SDV 16    : BA13
C   SDV 17    : (MU/JA)*BA13
C   SDV 18    : UMAT S13
C   SDV 19    : deviatoric part of UMAT S11
C   SDV 20    : UMAT S11
C =================================================================


C -----------------------------------------------------------------
C 1. UMAT INTERFACE
C -----------------------------------------------------------------

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     3 CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,
     4 DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     5 KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME


C -----------------------------------------------------------------
C 2. ABAQUS ARRAY DEFINITIONS
C -----------------------------------------------------------------

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),
     2 DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     3 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),
     4 COORDS(3),DROT(3,3),DFGRD0(3,3),
     5 DFGRD1(3,3),JSTEP(4)


C -----------------------------------------------------------------
C 3. MATERIAL AND GROWTH VARIABLES
C -----------------------------------------------------------------

      DOUBLE PRECISION EMOD,ENU
      DOUBLE PRECISION MU,BULK,C10,D1
      DOUBLE PRECISION LAMG,KAPPA0,XLEN,PI
      DOUBLE PRECISION KAPPA,GROWTH
      DOUBLE PRECISION ALPHA,GTARGET

      DOUBLE PRECISION ZERO,ONE,TWO,THREE,FOUR,SIX
      DOUBLE PRECISION XREF,YREF,ZREF

      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0,
     1 THREE=3.D0, FOUR=4.D0, SIX=6.D0)


C -----------------------------------------------------------------
C 4. KINEMATIC AND CONSTITUTIVE VARIABLES
C -----------------------------------------------------------------

      DOUBLE PRECISION F(3,3)
      DOUBLE PRECISION G(3,3),GINV(3,3)
      DOUBLE PRECISION AEL(3,3)

      DOUBLE PRECISION BA(3,3),BABAR(3,3)
      DOUBLE PRECISION JF,JG,JA
      DOUBLE PRECISION DETBA,SCALE
      DOUBLE PRECISION TRACEBA,PR

      DOUBLE PRECISION EG,EK,EG23

      INTEGER I,J,K


C -----------------------------------------------------------------
C 5. READ MATERIAL AND GROWTH PARAMETERS
C -----------------------------------------------------------------

      EMOD   = PROPS(1)
      ENU    = PROPS(2)
      LAMG   = PROPS(3)
      KAPPA0 = PROPS(4)
      XLEN   = PROPS(5)

      PI = FOUR*DATAN(ONE)

C     MATERIAL CONSTANTS

      MU   = EMOD/(TWO*(ONE+ENU))
      BULK = EMOD/(THREE*(ONE-TWO*ENU))

      C10 = EMOD/(FOUR*(ONE+ENU))
      D1  = SIX*(ONE-TWO*ENU)/EMOD


C -----------------------------------------------------------------
C 6. REFERENCE COORDINATES AND GROWTH TENSOR
C -----------------------------------------------------------------

C     READ REFERENCE COORDINATES STORED BY SDVINI

      XREF = STATEV(1)
      YREF = STATEV(2)
      ZREF = STATEV(3)

C     NONUNIFORM GROWTH FIELD

      KAPPA = KAPPA0*DSIN(TWO*PI*XREF/XLEN)

C     TARGET GROWTH

      GTARGET = LAMG*DEXP(ZREF*KAPPA)

C     RAMP GROWTH FROM 1 TO THE TARGET VALUE OVER TOTAL TIME 10

      ALPHA = DMIN1(ONE,DMAX1(ZERO,(TIME(1)+DTIME)/10.0D0))
      GROWTH = ONE + ALPHA*(GTARGET-ONE)

C     SAVE KAPPA

      STATEV(4) = KAPPA

C     INITIALIZE G AND G^(-1)

      DO I=1,3
         DO J=1,3
            G(I,J)    = ZERO
            GINV(I,J) = ZERO
         END DO
      END DO

      G(1,1) = GROWTH
      G(2,2) = ONE
      G(3,3) = ONE

      GINV(1,1) = ONE/GROWTH
      GINV(2,2) = ONE
      GINV(3,3) = ONE

      JG = GROWTH


C -----------------------------------------------------------------
C 7. TOTAL DEFORMATION GRADIENT
C -----------------------------------------------------------------

      DO I=1,3
         DO J=1,3
            F(I,J) = DFGRD1(I,J)
         END DO
      END DO

      JF = F(1,1)*(F(2,2)*F(3,3)-F(2,3)*F(3,2))
     1   - F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))
     2   + F(1,3)*(F(2,1)*F(3,2)-F(2,2)*F(3,1))


C -----------------------------------------------------------------
C 8. ELASTIC DEFORMATION GRADIENT
C -----------------------------------------------------------------

      DO I=1,3
         DO J=1,3
            AEL(I,J) = ZERO
            DO K=1,3
               AEL(I,J) = AEL(I,J)
     1                    + F(I,K)*GINV(K,J)
            END DO
         END DO
      END DO

      JA = JF/JG


C -----------------------------------------------------------------
C 9. ELASTIC LEFT CAUCHY-GREEN TENSOR
C -----------------------------------------------------------------

      DO I=1,3
         DO J=1,3
            BA(I,J) = ZERO
            DO K=1,3
               BA(I,J) = BA(I,J)
     1                   + AEL(I,K)*AEL(J,K)
            END DO
         END DO
      END DO

      DETBA = JA*JA
      SCALE = DETBA**(-ONE/THREE)

      DO I=1,3
         DO J=1,3
            BABAR(I,J) = SCALE*BA(I,J)
         END DO
      END DO

      TRACEBA = (BABAR(1,1)+BABAR(2,2)+BABAR(3,3))/THREE


C -----------------------------------------------------------------
C 10. CAUCHY STRESS
C -----------------------------------------------------------------

      EG = TWO*C10/JA
      PR = TWO/D1*(JA-ONE)

      STRESS(1) = EG*(BABAR(1,1)-TRACEBA) + PR
      STRESS(2) = EG*(BABAR(2,2)-TRACEBA) + PR
      STRESS(3) = EG*(BABAR(3,3)-TRACEBA) + PR

      STRESS(4) = EG*BABAR(1,2)
      STRESS(5) = EG*BABAR(1,3)
      STRESS(6) = EG*BABAR(2,3)


C -----------------------------------------------------------------
C 11. MATERIAL JACOBIAN
C -----------------------------------------------------------------

      EG23 = EG*TWO/THREE
      EK   = TWO/D1*(TWO*JA-ONE)

      DO I=1,NTENS
         DO J=1,NTENS
            DDSDDE(I,J) = ZERO
         END DO
      END DO

      DDSDDE(1,1) = EG23*(BABAR(1,1)+TRACEBA) + EK
      DDSDDE(2,2) = EG23*(BABAR(2,2)+TRACEBA) + EK
      DDSDDE(3,3) = EG23*(BABAR(3,3)+TRACEBA) + EK

      DDSDDE(1,2) = -EG23*
     1 (BABAR(1,1)+BABAR(2,2)-TRACEBA) + EK

      DDSDDE(1,3) = -EG23*
     1 (BABAR(1,1)+BABAR(3,3)-TRACEBA) + EK

      DDSDDE(2,3) = -EG23*
     1 (BABAR(2,2)+BABAR(3,3)-TRACEBA) + EK

      DDSDDE(1,4) =  EG23*BABAR(1,2)/TWO
      DDSDDE(2,4) =  EG23*BABAR(1,2)/TWO
      DDSDDE(3,4) = -EG23*BABAR(1,2)

      DDSDDE(1,5) =  EG23*BABAR(1,3)/TWO
      DDSDDE(2,5) = -EG23*BABAR(1,3)
      DDSDDE(3,5) =  EG23*BABAR(1,3)/TWO

      DDSDDE(1,6) = -EG23*BABAR(2,3)
      DDSDDE(2,6) =  EG23*BABAR(2,3)/TWO
      DDSDDE(3,6) =  EG23*BABAR(2,3)/TWO

      DDSDDE(4,4) = EG*(BABAR(1,1)+BABAR(2,2))/TWO
      DDSDDE(5,5) = EG*(BABAR(1,1)+BABAR(3,3))/TWO
      DDSDDE(6,6) = EG*(BABAR(2,2)+BABAR(3,3))/TWO

      DDSDDE(4,5) = EG*BABAR(2,3)/TWO
      DDSDDE(4,6) = EG*BABAR(1,3)/TWO
      DDSDDE(5,6) = EG*BABAR(1,2)/TWO

      DO I=1,NTENS
         DO J=1,I-1
            DDSDDE(I,J) = DDSDDE(J,I)
         END DO
      END DO


C -----------------------------------------------------------------
C 12. STATE VARIABLES, POST-PROCESSING, ENERGY, AND END
C -----------------------------------------------------------------

C     ORIGINAL STATE VARIABLES

      STATEV(5) = JA
      STATEV(6) = GROWTH

C     POST-PROCESSING STATE VARIABLES
C     Requires *Depvar = 20 in the input file.

      IF (NSTATV .GE. 20) THEN
         STATEV(7)  = F(1,1)
         STATEV(8)  = F(1,3)
         STATEV(9)  = F(3,1)
         STATEV(10) = F(3,3)

         STATEV(11) = AEL(1,1)
         STATEV(12) = AEL(1,3)
         STATEV(13) = AEL(3,1)
         STATEV(14) = AEL(3,3)

         STATEV(15) = PR
         STATEV(16) = BA(1,3)

C        Theory-like shear using the same A:
C        sigma13 = (MU/JA)*(A*A^T)_13

         STATEV(17) = (MU/JA)*BA(1,3)

C        Actual UMAT shear stress

         STATEV(18) = STRESS(5)

C        Decompose sigma11 = sigma11_dev + PR

         STATEV(19) = EG*(BABAR(1,1)-TRACEBA)
         STATEV(20) = STRESS(1)
      END IF

C     ELASTIC STRAIN ENERGY

      SSE = JG*
     1 (C10*(THREE*TRACEBA-THREE)
     2 +(ONE/D1)*(JA-ONE)*(JA-ONE))

      SPD = ZERO
      SCD = ZERO

      RPL    = ZERO
      DRPLDT = ZERO

      DO I=1,NTENS
         DDSDDT(I) = ZERO
         DRPLDE(I) = ZERO
      END DO

      RETURN
      END


C =================================================================
C AUXILIARY SUBROUTINE: SDVINI
C INITIALIZE REFERENCE COORDINATES
C =================================================================

      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION STATEV(NSTATV),COORDS(NCRDS)

      INTEGER I

      DO I=1,NSTATV
         STATEV(I) = 0.D0
      END DO

      STATEV(1) = COORDS(1)
      STATEV(2) = COORDS(2)
      STATEV(3) = COORDS(3)

      RETURN
      END
