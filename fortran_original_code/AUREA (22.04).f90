!*************************************************************************
!
!  PROGRAMA: AUREA
!
!  OBJETIVO: OTIMIZACAO DE MÍNIMO - UNIVARIÁVEL
!
!*************************************************************************

!SUBROTINA AUREA PARA A SUBROTINA ORVALP
SUBROUTINE AUREA(TOL,ALI0,ALS0,XBEST,X,Y,T,P,F,dF, NC) 
    
    IMPLICIT REAL*8(A-H,O-Z)
    COMMON /DADOS/ TC(20), PC(20), W(20)
    COMMON /MODEL/ MODELV, ID
    COMMON /MRULE/ RMIX
  
    
    DIMENSION X(20), Y(20), FL(20), FV(20), CFL(20), CFV(20),  &
               RAZAO(20)
     
     ALI0 = 0D0
     ALS0 = 1D0 
     TOL = 1D-6
    
    !Inicialização -------------------------------------------------------
    ALI = ABS(ALI0)  
    ALS = ABS(ALS0)  
    DELTA = ALS - ALI
    XI = ALI + 0.382D0*DELTA
    XS = ALS - 0.382D0*DELTA
    CALL FUNCAO(XI,FI,X,Y,T,P,F,dF,nc)
    CALL FUNCAO(XS,FS,X,Y,T,P,F,dF,nc)
    !Iterações -----------------------------------------------------------
    
    DO WHILE (DELTA >= TOL)
        
        IF(FI>=FS) THEN
            ALI = XI
            XI = XS
            FI = FS
            DELTA = ALS - ALI
            XS = ALS - 0.382D0*DELTA
            CALL FUNCAO(XS,FS,X,Y,T,P,F,dF,nc)
        ELSE
            ALS = XS
            XS = XI
            FS = FI
            DELTA = ALS - ALI
            XI = ALI + 0.382D0*DELTA
            CALL FUNCAO(XI,FI,X,Y,T,P,F,dF,nc)
        END IF
        
    END DO
    
    XBEST = (XI+XS)/2D0
   
END
    
!SUBROTINA PARA A SUBROTINA AUREA ORVALP
SUBROUTINE FUNCAO(XAA,F2,X,Y,T,P_IN,F,dF,nc) !XAA é a variável do Aurea
    
    IMPLICIT REAL*8(A-H,O-Z)
    
    INTEGER :: PHA 
    
    COMMON /DADOS/ TC(20), PC(20), W(20) 
    COMMON /MODEL/ MODELV, ID
    COMMON /MRULE/ RMIX
   
    
    DIMENSION X(20), Y(20), FL(20), FV(20), &
                CFL(20), CFV(20), &
               RAZAO(20)
    
    P = P_in - (F/dF)*XAA 

    SOMA = 0D0
    PHA = 1
     
    IF(MODELV==2) THEN
       CALL PR0011( NC, PHA, P, T, X, Y, ZL, VL, FL )
       ELSE IF(MODELV==1) THEN
       CALL PHASES( NC, PHA, P, T, X, Y, ZL, VL, FL )
    END IF
     
    PHA = 2   ! INCLUIDO CALCULO DA FASE VAPOR, 
              ! no loop do 5004 o valor de CFV(I) = 0, RAZAO(I) = Inf
     
    IF(MODELV==2) THEN
        CALL PR0011( NC, PHA, P, T, X, Y, ZV, VV, FV )
    ELSE IF(MODELV==1) THEN
        CALL PHASES( NC, PHA, P, T, X, Y, ZV, VV, FV )
    END IF 
     
    DO 5004 I=1,NC
         
        IF(Y(I) /= 0D0) THEN
            CFL(I) = FL(I)/(X(I)*P)
            CFV(I) = FV(I)/(Y(I)*P)  ! INCLUIDO
            RAZAO(I) = CFL(I)/CFV(I)
            SOMA = SOMA + Y(I)/RAZAO(I)
            !CFL1(I) = FL1(I)/(X(I)*P1)
            !RAZAO1(I) = CFL1(I)/CFV1(I)
            !SOMA1(CC) = SOMA1(CC) + Y(I)/RAZAO1(I)
        END IF
            
5004 CONTINUE

    F2 = DABS(SOMA-1D0)  ! TESTAR AO QUADRADO TBM
    !F2 = (SOMA-1D0)  ! AQUI o valor da funcao deve ser o valor ORIGINAL

END
    
! SUBROTINA AUREA PARA A SUBROTINA BOLHAT
SUBROUTINE AUREABOLHAT (TOL,ALI0,ALS0,XBEST,X,Y,T,P,F,dF,nc)
 
    IMPLICIT REAL*8(A-H,O-Z)
    
    COMMON /DADOS/ TC(20), PC(20), W(20)
    COMMON /MODEL/ MODELV, ID
    COMMON /MRULE/ RMIX


    DIMENSION X(20), Y(20), FL(20), FV(20), &
            CFL(20), CFV(20),  &
               RAZAO(20)
     ALI0 = 0D0 !INTERVALO INFERIOR DO AUREA
     ALS0 = 1D0 !INTERVALO SUPERIOR DO AUREA
     TOL = 1D-5
    
    !Inicialização -------------------------------------------------------
    ALI = ABS(ALI0)
    ALS = ABS(ALS0)
    DELTA = ALS - ALI
    XI = ALI + 0.382D0*DELTA
    XS = ALS - 0.382D0*DELTA
    CALL FUNCAOBOLHAT(XI,FI,X,Y,T,P,F,dF,nc)
    CALL FUNCAOBOLHAT(XS,FS,X,Y,T,P,F,dF,nc)
    !Iterações -----------------------------------------------------------
    
    DO WHILE (DELTA >= TOL)
        
        IF(FI>=FS) THEN
            ALI = XI
            XI = XS
            FI = FS
            DELTA = ALS - ALI
            XS = ALS - 0.382D0*DELTA
            CALL FUNCAOBOLHAT(XS,FS,X,Y,T,P,F,dF,nc)
        ELSE
            ALS = XS
            XS = XI
            FS = FI
            DELTA = ALS - ALI
            XI = ALI + 0.382D0*DELTA
            CALL FUNCAOBOLHAT(XI,FI,X,Y,T,P,F,dF,nc)
        END IF
        
    END DO
    
    XBEST = (XI+XS)/2D0

END
    
!SUBROTINA PARA A SUBROTINA AUREABOLHAT 
SUBROUTINE FUNCAOBOLHAT(XAA,F2,X,Y,T_IN,P,F,dF,nc) !XAA é a variável do Aurea
    
    IMPLICIT REAL*8(A-H,O-Z)
    
    INTEGER :: PHA 
    
    COMMON /DADOS/ TC(20), PC(20), W(20) 
    COMMON /MODEL/ MODELV, ID
    COMMON /MRULE/ RMIX
    
    DIMENSION X(20), Y(20), FL(20), FV(20), &
              CFL(20), CFV(20), RAZAO(20)
    
    T = T_IN - (F/dF)*XAA
    SOMA = 0D0
    PHA = 2
    IF(MODELV==2) THEN
       CALL PR0011( NC, PHA, P, T, X, Y, ZV, VV, FV ) 
    ELSE IF(MODELV==1) THEN
       CALL PHASES( NC, PHA, P, T, X, Y, ZV, VV, FV )
    END IF
    
    DO 4004 I=1,NC
         
        IF(X(I) /= 0D0) THEN
            CFV(I) = FV(I)/(Y(I)*P)
            RAZAO(I) = CFL(I)/CFV(I)
            SOMA = SOMA + RAZAO(I)*X(I)
          !  as linhas a seguir nao sao necessarias
            !CFV1(I) = FV1(I)/(Y(I)*P1)   
            !RAZAO1(I) = CFL1(I)/CFV1(I)
            !SOMA1(CC) = SOMA1(CC) + RAZAO1(I)*X(I)
        END IF
         
4004 CONTINUE
     
    F2 = ABS(SOMA-1D0)
    
END
    
! SUBROTINA AUREA PARA A SUBROTINA ORVALT
SUBROUTINE AUREAORVALT (TOL,ALI0,ALS0,XBEST,X,Y,T,P,F,dF,nc)
 
IMPLICIT REAL*8(A-H,O-Z)
    
COMMON /DADOS/ TC(20), PC(20), W(20)
COMMON /MODEL/ MODELV, ID
COMMON /MRULE/ RMIX


DIMENSION X(20), Y(20), FL(20), FV(20), &
            CFL(20), CFV(20),  &
               RAZAO(20)
     ALI0 = 0D0 !INTERVALO INFERIOR DO AUREA
     ALS0 = 1D0 !INTERVALO SUPERIOR DO AUREA
     TOL = 1D-5
    
    !Inicialização -------------------------------------------------------
    ALI = ALI0 !TESTAR AO QUADRADO TBM
    ALS = ALS0 !TESTAR AO QUADRADO TBM
    DELTA = ALS - ALI
    XI = ALI + 0.382D0*DELTA
    XS = ALS - 0.382D0*DELTA
    CALL FUNCAOORVALT(XI,FI,X,Y,T,P,F,dF, nc)
    CALL FUNCAOORVALT(XS,FS,X,Y,T,P,F,dF, nc)
    !Iterações -----------------------------------------------------------
    
    DO WHILE (DELTA >= TOL)
       
        IF(FI>=FS) THEN
            ALI = XI
            XI = XS
            FI = FS
            DELTA = ALS - ALI
            XS = ALS - 0.382D0*DELTA
            CALL FUNCAOORVALT(XS,FS,X,Y,T,P,F,dF,nc)
        ELSE
            ALS = XS
            XS = XI
            FS = FI
            DELTA = ALS - ALI
            XI = ALI + 0.382D0*DELTA
            CALL FUNCAOORVALT(XI,FI,X,Y,T,P,F,dF,nc)
        END IF
        
    END DO
    
    XBEST = (XI+XS)/2D0

END
    
!SUBROTINA PARA A SUBRTOINA AUREA ORVALT
SUBROUTINE FUNCAOORVALT(XAA,F2,X,Y,T_IN,P,F,dF,nc) !XAA é a variável do Aurea
    
    IMPLICIT REAL*8(A-H,O-Z)
    
    INTEGER :: PHA 
    

    COMMON /DADOS/ TC(20), PC(20), W(20) 
    COMMON /MODEL/ MODELV, ID
    COMMON /MRULE/ RMIX
    
    
    DIMENSION X(20), Y(20), FL(20), FV(20), &
                CFL(20), CFV(20), &
               RAZAO(20)
    
    T = T_IN -(F/dF)*XAA
    SOMA = 0D0
    PHA = 1
    IF(MODELV==2) THEN
       CALL PR0011( NC, PHA, P, T, X, Y, ZL, VL, FL )
    ELSE IF(MODELV==1) THEN
       CALL PHASES( NC, PHA, P, T, X, Y, ZL, VL, FL )
    END IF
     
    DO 6004 I=1,NC
         
        IF(Y(I) /= 0D0) THEN
            CFL(I) = FL(I)/(X(I)*P)
            RAZAO(I) = CFL(I)/CFV(I)
            SOMA = SOMA + Y(I)/RAZAO(I)
        END IF
         
6004 CONTINUE
     
    F2 = ABS(SOMA-1D0) !TESTAR AO QUADRADO TBM 
     
END
    
!SUBROTINA PARA CÁLCULO DA REDUÇÃO NO FLASH
SUBROUTINE AUREAFLASH(TOL,ALI0,ALS0,XBEST,Z,T,Fvap_old,YXI,P,F,dF,CC,nc)
 
IMPLICIT REAL*8(A-H,O-Z)
    
COMMON /DADOS/ TC(20), PC(20), W(20)
COMMON /MODEL/ MODELV, ID
COMMON /MRULE/ RMIX


DIMENSION  YXI(20,10000), Z(20)

    ALI0 = 0D0 !INTERVALO INFERIOR DO AUREA
    ALS0 = 1D0 !INTERVALO SUPERIOR DO AUREA
    TOL = 1D-6
    
    !Inicialização -------------------------------------------------------
    ALI = ALI0    !!!!ABS(ALI0)  nao precisa de valor absoluto
    ALS = ALS0    !!!!ABS(ALS0)  nao precisa de valor absoluto
    DELTA = ALS - ALI
    XI = ALI + 0.382D0*DELTA
    XS = ALS - 0.382D0*DELTA
    CALL FUNCAOFLASH(XI,FI,Z,T,FVAP_old,YXI,P,F,dF,CC,nc)
    CALL FUNCAOFLASH(XS,FS,Z,T,FVAP_old,YXI,P,F,dF,CC,nc)    
    
    !Iterações -----------------------------------------------------------
    
    DO WHILE (DELTA >= TOL)
        
        IF(FI>=FS) THEN
            ALI = XI
            XI = XS
            FI = FS
            DELTA = ALS - ALI
            XS = ALS - 0.382D0*DELTA
            CALL FUNCAOFLASH(XS,FS,Z,T,FVAP_old,YXI,P,F,dF,CC,nc)
        ELSE
            ALS = XS
            XS = XI
            FS = FI
            DELTA = ALS - ALI
            XI = ALI + 0.382D0*DELTA
            CALL FUNCAOFLASH(XI,FI,Z,T,FVAP_old,YXI,P,F,dF,CC,nc)
        END IF
       
    END DO
    
    XBEST = (XI+XS)/2D0
    
END
    
!SUBRTOINA PARA CÁLCULO DA AUREA FLASH
SUBROUTINE FUNCAOFLASH(XAA,F2,Z,T,FVAP_old,YXI,P,F,dF,CC,nc) !XAA é a variável do Aurea
    
    IMPLICIT REAL*8(A-H,O-Z)
    
    INTEGER :: PHA 
    
    COMMON /DADOS/ TC(20), PC(20), W(20) 
    COMMON /MODEL/ MODELV, ID
    COMMON /MRULE/ RMIX
    COMMON /DIMENSIONAMENTO/ GRAND
    
    DIMENSION  YXI(20,10000), Z(20), CEQ(20)
   
    FVAP = Fvap_old-(F/dF)*XAA

    F1 = 0D0 ; 
         
        DO I=1,NC
            F1 = F1 + Z(I)*(YXI(I,CC)-1D0)/ &
                      (1D0+(FVAP/GRAND)*(YXI(I,CC)-1D0))                      
        END DO
      
    F2=dabs(F1)  ! TESTAR AO QUADRADO TBM
    !F2 = F1    ! AQUI o valor da funcao deve ser o valor ORIGINAL
    
END