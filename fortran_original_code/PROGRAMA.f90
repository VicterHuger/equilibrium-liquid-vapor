!********************************************************************************************************************************************************
!
!     PROGRAMA: CONDENSACAO
!
!     OBJETIVO: modelagem do fenômeno de condensação retrógrada de gás natural
!
!********************************************************************************************************************************************************

PROGRAM CONDENSACAO

      IMPLICIT REAL*8(A-H,O-Z)

      COMMON /DADOS/  TC(20), PC(20), W(20)    
      COMMON /NCOMP/  NC
      COMMON /MODEL/  MODELV
      COMMON /MRULE/  RMIX
      COMMON /ABM/    ALFAF, AAT1(20), AAT2(20), AAT3(20)

      DIMENSION Z(20), Y(100,20), X(100,20), T(100), P_bol(100), P_orv(100)

      CHARACTER*18 :: FILENAME1
      
!************************************************************************** INICIALIZAÇÃO       
      WRITE(*,0001)
0001  FORMAT(/,2X,'****************************************************', &
             /,22X,               'PROGRAM: CONDENSACAO'                , &
             /,2X,'****************************************************',/)
      
      MODELV = 1  !Modelo: (1) PR / (2) PR-UNIQUAC
      RMIX = 1    !Regra de Mistura: (1) HK / (2) LCVM
      ALFAF = 1   !Função Alfa: (1) PR78 / (2) AAT
!************************************************************************** DEFINIÇÃO DO SISTEMA      
      WRITE(*, "(/, 2X, 'Digite o numero do sistema:')")
      READ*, ID
!************************************************************************** LEITURA DOS DADOS DE ENTRADA      
      CALL BDADOS(ID, NC, Z)  
!************************************************************************** CÁLCULO DO PONTO DE ORVALHO      
      TI = 173.15D0 !K
      TF = 473.15d0 !k
      Npts = 100
      WRITE(*,'(2X," T(K) ",2X,"Pb(bar)",2X,"Po(bar)")')
      DO I=1,Npts
          T(I) = TI - (I - 1)*(TI-TF)/(Npts-1)
          CALL BOLHAP(NC, Z, T(I), P_bol(I), Y(I,:))
          CALL ORVALP(NC, Z, T(I), P_orv(I), X(I,:))
          WRITE(*,'(2X,F6.2,2X,F6.2,2X,F6.2)') T(I), P_bol(I), P_orv(I)
      END DO
!************************************************************************** IMPRESSÃO DE RESULTADOS     
      IF(ID<10) THEN
          WRITE(FILENAME1,'("Curva-",I1,".DAT")') ID
      ELSE 
          WRITE(FILENAME1,'("Curva-",I2,".DAT")') ID
      END IF
      OPEN(01,FILE=FILENAME1)
          WRITE(01,'(2X," T(K) ",2X,"Pb(bar)",2X, "Po(bar)")')
          DO I=1,Npts             
              WRITE(01,'(2X,F6.2,2X,F6.2,2X,F6.2)') T(I), P_bol(I), P_orv(I)
          END DO
      CLOSE(01)
      PAUSE
END PROGRAM

!****************************************************************************
SUBROUTINE WILSON (NC, Pi, X, Y, T)
    !IMPLICIT REAL*8(A-H,O-Z) 
    !Correlação de Wilson para o cálculo da pressão de orvalho inicial (Wilson, 1968)
    INTEGER :: I, NC
    REAL*8 RAZAO(20), PC, Pi, W, TC !,TCkr(20)
    REAL*8 X(20), Y(20), T !, Tkr
    COMMON /DADOS/ TC(20), PC(20), W(20)
    !
    !Tkr=T*1.8d0 !Conversão de K -> ºR
    !    
     DO I=1,NC
         !TCkr(I)=TC(I)*1.8D0 !Conversão de K-> ºR
         RAZAO(I) = (PC(I)/Pi)(DEXP(5.373D0(1D0+W(I))*(1D0-(TC(I)/T)))) ! Ki
         X(I) = Y(I)/RAZAO(I)
         !CORRELAÇÃO DE WILSON
     END DO
     X(1:NC) = X(1:NC)/SUM(X(1:NC)) 
     !WRITE(*,'(2/,2X,"Razao:", 6X, "PC", 6X, "W:", 6x, "TC:")')

END
      
!****************************************************************************************************************************************************
!
!     SUB-ROTINA: ORVALP 
!
!     OBJETIVO: Cálculo do ponto de orvalho (T e Y conhecidos)   
!
!****************************************************************************************************************************************************

SUBROUTINE ORVALP( NC, Y, T, P_orv, X )   

     IMPLICIT REAL*8(A-H,O-Z)
     
     COMMON /DADOS/  TC(20), PC(20), W(20)    
     COMMON /MODEL/  MODELV
     
     DIMENSION X(20), Y(20), FL(20), FV(20), FL1(20), &
               FV1(20), CFL(20), CFV(20), CFL1(20), CFV1(20),  &
               RAZAO(20), RAZAO1(20), SOMA(1000), SOMA1(1000)  
     
     INTEGER :: PHA

!**************************************************************************
     !ESTIMATIVAS INICIAIS
     SOMA=0D0 ; SOMA1=0D0
     P = 300D0 !bar
     CC = 1
     NCC = NC
     DO I=1,NC
         IF(Y(I)==0D0) NCC = NCC - 1
     END DO
     DO I=1,NC
         IF(Y(I)==0D0) THEN
             X(I) = 0D0
         ELSE
             X(I) = 1D0/NCC
         END IF
     END DO
!**************************************************************************
5001 CONTINUE
     IF (CC==1000) THEN
         !PRINT*, 'Pressao de orvalho nao convergiu' 
         GO TO 5006
     END IF
     
     P1 = P + 0.001*P
     T1 = T
     
     PHA = 1
     IF(MODELV==2) THEN
        CALL PR0011( NC, PHA, P, T, X, Y, ZL, VL, FL )
     ELSE IF(MODELV==1) THEN
        CALL PHASES( NC, PHA, P, T, X, Y, ZL, VL, FL )
     END IF 
     PHA = 2
     IF(MODELV==2) THEN
        CALL PR0011( NC, PHA, P, T, X, Y, ZV, VV, FV )
        CALL PR0011( NC, PHA, P1, T1, X, Y, ZV1, VV1, FV1 )
     ELSE IF(MODELV==1) THEN
        CALL PHASES( NC, PHA, P, T, X, Y, ZV, VV, FV )
        CALL PHASES( NC, PHA, P1, T1, X, Y, ZV1, VV1, FV1 )
     END IF 
     SOMA(CC) = 0D0
     DO 5002 I=1,NC
         IF(Y(I) /= 0D0) THEN
             CFL(I) = FL(I)/(X(I)*P)
             CFV(I) = FV(I)/(Y(I)*P)
             CFV1(I) = FV1(I)/(Y(I)*P1) 
             RAZAO(I) = CFL(I)/CFV(I)
             SOMA(CC) = SOMA(CC) + Y(I)/RAZAO(I)
         END IF
5002 CONTINUE
     
     XTOTAL = 0D0
     DO 5003 I=1,NC
         IF(Y(I) /= 0D0) X(I) = Y(I)*SOMA(CC)/RAZAO(I)
         XTOTAL = XTOTAL + X(I)
5003 CONTINUE
     DO 5007 I=1,NC
         X(I) = X(I)/XTOTAL !NORMALIZAÇÃO
5007 CONTINUE
     
     IF (CC==1000) THEN
         !PRINT*, 'Pressao de orvalho nao convergiu'
         GO TO 5006
     END IF
     CC = CC + 1
     SOMA(CC) = 0D0
     PHA = 1
     IF(MODELV==2) THEN
        CALL PR0011( NC, PHA, P, T, X, Y, ZL, VL, FL )
        CALL PR0011( NC, PHA, P1, T1, X, Y, ZL1, VL1, FL1 )
     ELSE IF(MODELV==1) THEN
        CALL PHASES( NC, PHA, P, T, X, Y, ZL, VL, FL )
        CALL PHASES( NC, PHA, P1, T1, X, Y, ZL1, VL1, FL1 )
     END IF
     DO 5004 I=1,NC
         IF(Y(I) /= 0D0) THEN
             CFL(I) = FL(I)/(X(I)*P)
             RAZAO(I) = CFL(I)/CFV(I)
             SOMA(CC) = SOMA(CC) + Y(I)/RAZAO(I)
             CFL1(I) = FL1(I)/(X(I)*P1)
             RAZAO1(I) = CFL1(I)/CFV1(I)
             SOMA1(CC) = SOMA1(CC) + Y(I)/RAZAO1(I)
         END IF
5004 CONTINUE
     
     TESTE1 = DABS(SOMA(CC)-SOMA(CC-1))
     IF(TESTE1<1D-5) GO TO 5005
     GO TO 5002
5005 CONTINUE
     
     F = SOMA(CC)-1D0 
     F1 = SOMA1(CC)-1D0 
     dF = (F1 - F)/(P1 - P)
     
     IF(DABS(F)<=1D-6) THEN
        GO TO 5006
     ELSE
        P = P - F/dF
        IF(P<0D0) P = 0.05D0
        GO TO 5001
     END IF   
5006 CONTINUE
     
     P_orv = P
END

!****************************************************************************************************************************************************
!
!     SUB-ROTINA: BOLHAP 
!
!     OBJETIVO: Cálculo do ponto de bolha (T e X conhecidos)   
!
!****************************************************************************************************************************************************

SUBROUTINE BOLHAP( NC, X, T, P_bol, Y )    
                 
     IMPLICIT REAL*8(A-H,O-Z)
     
     COMMON /DADOS/  TC(20), PC(20), W(20)    
     COMMON /MODEL/  MODELV
     
     DIMENSION X(20), Y(20), FL(20), FV(20), FL1(20), &
               FV1(20), CFL(20), CFV(20), CFL1(20), CFV1(20),  &
               RAZAO(20), RAZAO1(20), SOMA(1000), SOMA1(1000)  
     
     INTEGER :: PHA
     
!**************************************************************************
     !ESTIMATIVAS INICIAIS
     SOMA=0D0 ; SOMA1=0D0
     P = 300D0 !bar
     CC = 1
     NCC = NC
     DO I=1,NC
         IF(X(I)==0D0) NCC = NCC - 1
     END DO
     DO I=1,NC
         IF(X(I)==0D0) THEN
             Y(I) = 0D0
         ELSE
             Y(I) = 1D0/NCC
         END IF
     END DO
!**************************************************************************
3001 CONTINUE
     IF (CC==1000) THEN
         !PRINT*, 'Pressao de bolha nao convergiu' 
         GO TO 3006
     END IF
     
     P1 = P + 0.001*P
     T1 = T
     
     PHA = 1
     IF(MODELV==2) THEN
        CALL PR0011( NC, PHA, P, T, X, Y, ZL, VL, FL )
        CALL PR0011( NC, PHA, P1, T1, X, Y, ZL1, VL1, FL1 )
     ELSE IF(MODELV==1) THEN
        CALL PHASES( NC, PHA, P, T, X, Y, ZL, VL, FL )
        CALL PHASES( NC, PHA, P1, T1, X, Y, ZL1, VL1, FL1 )
     END IF 
     PHA = 2
     IF(MODELV==2) THEN
        CALL PR0011( NC, PHA, P, T, X, Y, ZV, VV, FV )
     ELSE IF(MODELV==1) THEN
        CALL PHASES( NC, PHA, P, T, X, Y, ZV, VV, FV )
     END IF 
     SOMA(CC) = 0D0
     DO 3002 I=1,NC
         IF(X(I) /= 0D0) THEN
             CFL(I) = FL(I)/(X(I)*P)
             CFL1(I) = FL1(I)/(X(I)*P1)
             CFV(I) = FV(I)/(Y(I)*P)
             RAZAO(I) = CFL(I)/CFV(I)
             SOMA(CC) = SOMA(CC) + RAZAO(I)*X(I)
         END IF
3002 CONTINUE
     
     YTOTAL = 0D0
     DO 3003 I=1,NC
         IF(X(I) /= 0D0) THEN
            Y(I) = RAZAO(I)*X(I)/SOMA(CC)
         END IF
         YTOTAL = YTOTAL + Y(I)
3003 CONTINUE
     DO 3007 I=1,NC
         Y(I) = Y(I)/YTOTAL !NORMALIZAÇÃO
3007 CONTINUE
     
     IF (CC==1000) THEN
         !PRINT*, 'Pressao de bolha nao convergiu' 
         GO TO 3006
     END IF
     CC = CC + 1
     SOMA(CC) = 0D0
     PHA = 2
     IF(MODELV==2) THEN
        CALL PR0011( NC, PHA, P, T, X, Y, ZV, VV, FV )
        CALL PR0011( NC, PHA, P1, T1, X, Y, ZV1, VV1, FV1 )
     ELSE IF(MODELV==1) THEN
        CALL PHASES( NC, PHA, P, T, X, Y, ZV, VV, FV )
        CALL PHASES( NC, PHA, P1, T1, X, Y, ZV1, VV1, FV1 )
     END IF
     DO 3004 I=1,NC
         IF(X(I) /= 0D0) THEN
             CFV(I) = FV(I)/(Y(I)*P)
             RAZAO(I) = CFL(I)/CFV(I)
             SOMA(CC) = SOMA(CC) + RAZAO(I)*X(I)
             CFV1(I) = FV1(I)/(Y(I)*P1)
             RAZAO1(I) = CFL1(I)/CFV1(I)
             SOMA1(CC) = SOMA1(CC) + RAZAO1(I)*X(I)
         END IF
3004 CONTINUE
     
     TESTE1 = DABS(SOMA(CC)-SOMA(CC-1))
     IF(TESTE1<1D-5) GO TO 3005 
     GO TO 3002
3005 CONTINUE
     
     F = SOMA(CC)-1D0
     F1 = SOMA1(CC)-1D0
     dF = (F1 - F)/(P1 - P)
      
     IF(DABS(F)<=1D-6) THEN
        GO TO 3006
     ELSE
        P = P - F/dF
        IF(P<0D0) P = 0.05D0
        GO TO 3001
     END IF   
3006 CONTINUE
     
     P_bol = P
END    
