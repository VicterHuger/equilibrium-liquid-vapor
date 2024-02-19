!********************************************************************************************************************************************************
!
!     PROGRAMA: CONDENSACAO
!
!     OBJETIVO: Modelagem do fenômeno de condensação retrógrada de gás natural
!
!********************************************************************************************************************************************************

      PROGRAM CONDENSACAO
      !
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER :: MODELV, RMIX, ALFAF, VTPR, SUBPDEW, PDEW, PHA !, PSEUDO
      !
      COMMON /DADOS/  TC(20), PC(20), W(20)    
      COMMON /NCOMP/  NC
      COMMON /MODEL/  MODELV,ID
      COMMON /MRULE/  RMIX
      COMMON /ABM/ ALFAF, AAT1(20), AAT2(20), AAT3(20)
      COMMON /VTPRR/ VTPR  
      COMMON /BETAA/ BETA
      COMMON /DIMENSIONAMENTO/ GRAND
      COMMON /PORV/ WM(20), d60(20), Tb(20)
      COMMON /PDEWW/ PDEW
      !
      DIMENSION Z(20), Y(20), X(20), FG(20), Pexp(25), Vexp(25), Vdrop(25), V_error(25), c(25), CMJ_L(25), CMJ_G(25)
      !
      CHARACTER*18 :: FILENAME1
      !
!************************************************************************** INICIALIZAÇÃO       
      WRITE(*,0001)
0001  FORMAT(/,2X,'****************************************************', &
             /,22X,               'PROGRAM: CONDENSACAO'                , &
             /,2X,'****************************************************',/)
     CALL CPU_TIME(START)
      !
     MODELV = 1  !Modelo: (1) PR / (2) PR-UNIQUAC
      RMIX = 1    !Regra de Mistura: (1) HK / (2) LCVM
      ALFAF = 1   !Função Alfa: (1) PR78 / (2) AAT 
      VTPR = 1    !Translação volumétrica: (1) OFF / (2) ON
      SUBPDEW = 1 !Subrotina de Pdew: (1) ORVALP / (2) PORVALIT
      PDEW= 1     !Estimativa inicial de Pdew: (1) Wilson / (2) Elsharkaw
     !PSEUDO= 1   !Carcaterização de C7+: (1) Lee-Kesler / (2) Pedersen             !
!************************************************************************** DEFINIÇÃO DO SISTEMA      
      WRITE(*, "(/, 2X, 'Digite o numero do sistema:')")
      READ*, ID
      !
!************************************************************************** LEITURA DOS DADOS DE ENTRADA      
      CALL BDADOS(ID, NC, Z, T, P_orv, Nexp, Pexp, Vexp)  
      
!************************************************************************** CÁLCULO DO PONTO DE ORVALHO      
      !
      IF (SUBPDEW==2) THEN
        CALL PorvaIt(NC, Z, T, P_calc, X) !AAT1, AAT2, AAT3,APE1, APE2, APE3, r, q, UIJ_0, UIJ_T, DIJ, MODELV, RMIX, ALFAF, VTPR, ICORR)
      ELSEIF (SUBPDEW==1) THEN 
        CALL ORVALP(NC, Z, T, P_calc, X)
      ENDIF
      !
      ERRO = DABS(P_orv - P_calc)/P_orv
      WRITE(*,'(/,"P_exp: ",F6.2)') P_orv
      WRITE(*,'("P_cal: ",F6.2)') P_calc
      WRITE(*,'("Erro: "F6.2," %",/)') 100*ERRO
!************************************************************************** LIQUID DROPOUT     
!      !         
!      Vi = 1D0
!      RG = 83.14462D0 !(cm3.bar.K-1.mol-1)
!      PHA = 2
!      !
!      IF(MODELV==2) THEN
!         CALL PR0011( NC, PHA, P_calc, T, X, Y, ZG, VG, FG ) !Z as argument?
!      ELSE IF(MODELV==1) THEN
!         CALL PHASES( NC, PHA, P_calc, T, X, Y, ZG, VG, FG ) !Z as argument?
!      END IF 
!      !
!      TotalN = P_calc*Vi/(ZG*RG*T)
!      WRITE(*,'(2/,6x,"P(bar):", 20x, "V drop(%):", 15x, "Error(%):")')
!      !
!      DO I=1,Nexp
!          IF(Pexp(I) >= P_calc) GO TO 0002
!          P = Pexp(I)
!          IF(I==1) BETA=1D-2    
!          CALL FLASH( NC, TotalN, P, T, Z, X, Y, TotalNL, TotalNG, ZL, ZG, VL, VG )
!          !
!          !Translação volumétrica
!          !CM_L=0
!          !CM_G=0
!          !DO J = 1,NC
!          !    CALL VTPELL(J, T, c)
!          !    CMJ_L(J) = X(J)*c(J)
!          !    CMJ_G(J) = Y(J)*c(j)
!          !    CM_L = CM_L + CMJ_L(j)
!          !    CM_G = CM_G + CMJ_G(J)
!          !END DO
!          !VL = VL - CM_L   !conferir se é +/- CM  -> Vcorr = v - c
!          !VG= VG - CM_G
!         !Print*, VL, vg          
!          TotalVL = TotalNL*VL
!          TotalVG = TotalNG*VG
!          TotalV = TotalVL + TotalVG
!          OUTV = TotalV - Vi
!          OUTN = P*OUTV/(ZG*RG*T)
!          TotalNG = TotalNG - OUTN
!          TotalN = TotalNL + TotalNG
!          
!          DO J=1,NC
!              Z(J) = (X(J)*TotalNL + Y(J)*TotalNG)/TotalN
!          END DO
!
!          Vdrop(I) = (TotalVL/Vi)*100D0
!          V_error(I) = 100d0*(DABS(Vdrop(I)-Vexp(I)))/Vexp(i)
!          
!          WRITE(*,'(/,6X,F6.2,23X,F8.2,18X,F8.2)') P, Vdrop(I), V_error(i)
!          !Print*, Vdrop(i), v_error(i)
!0002      CONTINUE
!      END DO
!      
!      ERRORS=0D0  !Soma dos erros no cálculo do volume de dropout
!      
!      NERROR= NEXP!Número de pontos experimentais que foram cálculado o volume dropout
!      
!      !DO I=1,NEXP
!          
!      !  IF (V_ERROR(I)==0D0) NERROR= NERROR-1d0
!      !END DO
!      
!      !Print*,NERROR
!      
!      DO I=1,NEXP
!          
!      ERRORS=ERRORS + V_error(i)
!    
!      END DO
!      
!      ERRORM=ERRORS/NERROR !Média dos erros do VDROPOUT
!      !PRINT*, ERRORS, ERRORM
!      
!      !WRITE(*,'(/, "GRAND =", F6.2)')GRAND
!      
!      WRITE(*,'(/, "ERROR MEDIO=", F6.2)') ERRORM
!      !PRINT*, ERRORM
!!************************************************************************** IMPRESSÃO DE RESULTADOS     
!      
!!      IF(ID<10) THEN
!!          WRITE(FILENAME1,'("Resultados-",I1,".DAT")') ID
!!      ELSE 
!!          WRITE(FILENAME1,'("Resultados-",I2,".DAT")') ID
!!      END IF
!!      OPEN(40,FILE=FILENAME1)
!!      WRITE(40,'(/,"Pexp= :",f6.2)') P_orv
!!      WRITE(40,'("P_cal: ",F6.2)') P_calc
!!      WRITE(40,'("Erro: "F6.2," %",/)') 100*ERRO
!!      DO I=1,Nexp
!!          IF(Pexp(I) >= P_calc) GO TO 0002
!!          P = Pexp(I)
!!          CALL FLASH( NC, TotalN, P, T, Z, X, Y, TotalNL, TotalNG, ZL, ZG, VL, VG )
!!          TotalVL = TotalNL*VL
!!          TotalVG = TotalNG*VG
!!          TotalV = TotalVL + TotalVG
!!          OUTV = TotalV - Vi
!!          OUTN = P*OUTV/(ZG*RG*T)
!!          TotalNG = TotalNG - OUTN
!!          TotalN = TotalNL + TotalNG
!!          DO J=1,NC
!!              Z(J) = (X(J)*TotalNL + Y(J)*TotalNG)/TotalN
!!          END DO
!!          Vdrop(I) = (TotalVL/Vi)*100D0
!!          WRITE(*,'("P: ",F6.2,1X,"bar",2X,"Vdrop: ",F5.2," %")') P, Vdrop(I)
!!          WRITE(40,'("P: ",F6.2,1X,"bar",2X,"Vdrop: ",F5.2," %")') P, Vdrop(I)
!!          !PAUSE
!!0002      CONTINUE
!!      END DO
!!      CLOSE(40)
!      CALL CPU_TIME(FINISH)
!      PRINT '(/,"TIME= ",F20.3" MILISECONDS.")', (FINISH-START)*1D3
      PAUSE
!      !WRITE(*,'(2/,"***********************************************")')
!      !WRITE(*,'(2/,"Voce deseja imprimir mais resultados ?")') 
!      !WRITE(*,'(2/,"Responda com 1 para resposta afirmativa ", 2/)')
!      !READ*, RESPOSTA
!      !IF (RESPOSTA == 1 ) GO TO 2222
    END PROGRAM

SUBROUTINE ORVALP( NC, Y, T, P_orv, X )   

     IMPLICIT REAL*8(A-H,O-Z)
     
     
     COMMON /DADOS/ TC(20), PC(20), W(20)    
     COMMON /MODEL/ MODELV, ID
     COMMON /MRULE/ RMIX
     COMMON /PORV/ WM(20), d60(20), Tb(20)
     COMMON /PDEWW/ PDEW
     
     DIMENSION X(20), Y(20), FL(20), FV(20), FL1(20), &
               FV1(20), CFL(20), CFV(20), CFL1(20), CFV1(20),  &
               RAZAO(20), RAZAO1(20), SOMA(10000), SOMA1(10000)  
     
     INTEGER :: PHA, PDEW, NC, ID

!**************************************************************************
     !ESTIMATIVAS INICIAIS
     SOMA=0D0 ; SOMA1=0D0
     IF (PDEW==1) P = 300D0 !bar  Estimativa inical para o Wilson
     CC = 1 !Contador igual 1
     !DO I=1,NC
     !    PRINT*, Y(I)
     !ENDDO
     !
     IF (PDEW==2) THEN
         CALL ElsharkawPd(NC, T, Pd, Y, X)
     ELSEIF (PDEW==1) THEN
         CALL WILSON(NC, P, X, Y, T)
         !PRINT*, P
     END IF 
     IF (PDEW==2) P = 0.6d0*Pd 
    !IF (PDEW==1) THEN !CORRELAÇÃO DE WILSON 
    !    DO I=1,NC
    !     RAZAO(I) = (PC(I)/P)*DEXP(5.373D0*(1D0+W(I))*(1D0-(TC(I)/T))) ! Ki
    !     X(I) = Y(I)/RAZAO(I) 
    !     CORRELAÇÃO DE WILSON
    !    END DO
    !    X(1:NC) = X(1:NC)/SUM(X(1:NC))
    !END IF

!**************************************************************************     
5001 CONTINUE
     IF (CC==10000) THEN 
       WRITE(*,*) ' ITERATIONS EXCEEDED IN ORVALP C>10000 '
       GO TO 5006
     END IF
     
     P1 = P + 0.0001*P
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
           ! Nao precisa colocar esse teste, o coeficiente de fugacidade sera calculado e sera diferente de zero
         IF(Y(I) /= 0D0) THEN
           ! Fases corrigido fornece o coeficiente de fugacidade, nao precisa dividir por ( X(I)*P ou Y(I)*P )
             CFL(I) = FL(I)/(X(I)*P)
             CFV(I) = FV(I)/(Y(I)*P)
             CFV1(I) = FV1(I)/(Y(I)*P1) 
             RAZAO(I) = CFL(I)/CFV(I)
             SOMA(CC) = SOMA(CC) + Y(I)/RAZAO(I)
         END IF
         
5002 CONTINUE
     
5020 XTOTAL = 0D0
     
     DO 5003 I=1,NC
           ! Nao precisa colocar esse teste, a razao sera calculada e sera diferente de zero
!         IF(Y(I) /= 0D0) X(I) = Y(I)*SOMA(CC)/RAZAO(I)
           ! Lembrar que : X(I) = Y(I)/RAZAO(I) e a normalizacao divide cada termo pelo somatorio
         IF(Y(I) /= 0D0) X(I) = Y(I)/RAZAO(I)/SOMA(CC)
!         XTOTAL = XTOTAL + X(I)
5003 CONTINUE
     
!     DO 5007 I=1,NC
!         X(I) = X(I)/XTOTAL !NORMALIZAÇÃO
!5007 CONTINUE
     
     IF (CC==10000) GO TO 5006
    
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
           ! Fases corrigido fornece o coeficiente de fugacidade, nao precisa dividir por ( X(I)*P ou Y(I)*P )         
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
     GO TO 5020   !!5002
5005 CONTINUE
     
     F = SOMA(CC)-1D0 
     F1 = SOMA1(CC)-1D0 
     dF = (F1 - F)/(P1 - P)
     
     IF(DABS(F)<=1D-6) THEN
        GO TO 5006
     ELSE
        !
        DO 5030 I=1,NC
        ! Nao precisa colocar esse teste, a razao sera calculada e sera diferente de zero
        ! Lembrar que : X(I) = Y(I)/RAZAO(I) e a normalizacao divide cada termo pelo somatorio
         IF(Y(I) /= 0D0) X(I) = Y(I)/RAZAO(I)/SOMA(CC)   !!! INCLUIDO
5030    CONTINUE
        !
        CALL AUREA(TOL,ALI0,ALS0,Reduc,X,Y,T,P,F,dF, NC)!MODIFICAÇÃO AUREA
        P = P - (F/dF)*REDUC
        IF(P<0D0) P = 0.05D0
        GO TO 5001
     END IF   

5006 CONTINUE
    
     P_orv = P
!PRINT*, P_ORV
    
END    
    
SUBROUTINE PorvaIt(NC, Y, T, P_orv, X)
!IMPLICIT REAL*8(A-H,O-Z) 
!PorvaIt(NC, T, P_orv, TC, PC, W, WM, d60, TB, Z, Y, X, AAT1, AAT2, AAT3, APE1, APE2, APE3, r, q, UIJ_0, UIJ_T, DIJ, MODELV, RMIX, ALFAF, VTPR, ICORR) 
INTEGER  I, ITER, J, Niter, NiterM, IterM, PHA, iflagPd, ID, KK !Ndiv
INTEGER  PDEW, NC
REAL*8 Y(20), T, P_orv, X(20)
REAL*8 Soma(1000), Soma1(1000)
REAL*8 F !, XTOTAL, Test1, F1, dF, Reduc
REAL*8 Sumxl! cm, ii, Tri, Pri, P1 
REAL*8 Tolx! TolF, trry 
REAL*8 P,Pi, dP! FF, FF1, P_est, PCm, Fr1, Fr2, derF, Pcor
REAL*8 Pd ! T1,
REAL*8 Dx ! ALFA, Alfa1, Alfa2, Alfa3, 
REAL*8 ZV, VV, ZL, VL ! ZV1, VV1, ZL1, zVL1
REAL*8 ckWls(20), ckWlsn(20)
REAL*8 FV(20), FL(20), X1(20) ! FV1(NC), FL1(NC)
REAL*8 CEQ(NC) , CFL(NC), CFV(NC), CFL1(NC), CFV1(NC), CEQ1(NC)
!
COMMON /DADOS/  TC(20), PC(20), W(20)    
COMMON /PORV/ WM(20), d60(20), Tb(20)
COMMON /MODEL/  MODELV,ID
COMMON /PDEWW/ PDEW
!
! DIMENSION :: X(20), X1(20), CEQ(20)
IF (PDEW==2) THEN
    CALL ElsharkawPd(NC, T, Pd, Y, X)
    !PRINT*, PD
ELSEIF (PDEW==1) THEN
    P = 300D0 !bar  Estimativa inical para o Wilson
    CALL WILSON(NC, P, X, Y, T)
END IF 
!
IF (PDEW==2) P = 0.6d0*Pd  ! usa 60% do valor previsto como valor inicial, desvio relativo da correlação é aprox esse valor
Niter = 0
NiterM = 300
IterM = 500
Tolx = 1D-6
dP = 10d0
Iter = 1
DO While (Niter < NiterM)
!
  Iter = 0
!
  PHA = 2
  If (MODELV .EQ. 2) Then
    !Call PR0011(NC, Pha, P, T, X, Y, zL, VL, fl)
     Call PR0011(NC, PHA, P, T, X, Y, ZV, VV, FV)
  ElseIf (MODELV .EQ. 1) Then
    Call PHASES(NC, PHA, P, T, X, Y ,ZV, VV, FV)
  End IF
!  
  iflagPd = 0
!  
  DO KK= 1, IterM   !DO While (iflagPd = 0)
!    troquei a linha anterior porque parece nao aceita dos WHILE aninhados  <====
    Iter = Iter + 1
!   
    PHA = 1
!    
    If (MODELV .EQ. 2) Then
!      Call PR0011(NC, Pha, P, T, X, Y, zL, VL, fl)
      Call PR0011(NC, PHA, P, T, X, Y, ZL, VL, FL)
    ElseIf (MODELV .EQ. 1) Then
      Call PHASES(NC, PHA, P, T, X, Y, ZL, VL, FL)
    Endif
!    
    Sumxl = 0.0d0
    DO I = 1,NC
          CEQ(I) = FL(I) / FV(I)
          X(I) = Y(I) / CEQ(I)
          Sumxl = Sumxl + X(I)
    ENDDO
!    
    DO I = 1,NC
      X(I) = X(I) / Sumxl
    ENDDO
!    
    Dx = 0.0d0
    X1=0.0D0
    CALL XYFUN(NC, Dx, X, X1)
    DO I = 1,NC
      If (Iter > 50) Then
         X1(I) = 0.5d0 * (X(I) + X1(I))
      Else
         X1(I) = X(I)
      EndIF
    ENDDO
!    
    If (Dx < Tolx) iflagPd = 1
!    
    If (Iter > IterM .OR. iflagPd == 1) GOTO 5 
    ! Nao convergiu, ou deveria avisar
!    
  ENDDO !While (iflagPd = 0)
!  
5  CONTINUE  
!  
   !Se Sumxl < 1, ENTAO, P  > Psat (Valor de P esta acima da P de saturacao) 
  If ((Sumxl - 1.0000001d0) > 0) Then ! 0.0000001
    P = P + 10.0d0
  Else
    If (dP < 0.01D0) Then
      F = (Sumxl - 1.0D0)
!      
      If (Abs(F) < TolF) Then
!        Convergencia atingida
       GoTo 10   !Exit Sub
      End If
    EndIF
!   
    P = P - 0.75D0 * dP
    dP = dP / 10D0
    P = P + dP
  EndIF
  Niter = Niter + 1
!
  ENDDO ! While (Niter < NiterM)
!
10 CONTINUE
P_orv = P
!
RETURN
END
!==================================================================================================================
!==================================================================================================================
SUBROUTINE ElsharkawPd(NC, T, P_orv, Y, X)
!IMPLICIT REAL*8(A-H,O-Z) 
!
! Predicting the dew point pressure for gas condensate reservoirs: empirical models and equations of state
! Adel m.Elsharkawy
!https://doi.org/10.1016/S0378-3812(01)00724-5
!==================================================================================================================
!==================================================================================================================
INTEGER :: I, NC
!
REAL*8 Trk, TF, TCi, T
REAL*8 TC, PC, W, Y(20), X(20)
REAL*8 WM, d60, Tb
REAL*8 fH2S, fN2, fCO2, fC1, fC2, fC3, fiC4
REAL*8 fC4, fiC5, fC5, fC6, fC7p
REAL*8 cSG, cMW, Pd, PdPsi ! A, B, C
REAL*8 A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18
REAL*8 PCm, TCm, TBm, Bet, Wmtd60, Argm, Pk
REAL*8 Tri, Sumxl, Sum3 !Pri, P1
REAL*8 ALFA(3), P_orv
REAL*8 ckWls(20)
! 
COMMON /PORV/ WM(20), d60(20), Tb(20) 
COMMON /DADOS/ TC(20), PC(20), W(20)
!
A0 = 4268.85D0
A1 = 0.094056D0
A2 = -7157.87D0
A3 = -4540.58D0
A4 = -4663.55D0
A5 = -1357.56D0
A6 = -7776.1D0
A7 = -9967.99D0
A8 = -4257.1D0
A9 = -1417.1D0
A10 = 691.5298D0
A11 = 40660.36D0
A12 = 205.26D0
A13 = -7260.32D0
A14 = -352.413D0
A15 = -114.519D0
A16 = 8.133D0
A17 = 94.916D0
A18 = 238.252D0
!
Trk = 1.8D0 * T     !Conversão Kelvin -> Rankine
TF = Trk - 459.6D0  !Conversão Rankine -> Farenheit
cSG = d60(NC) ! spgr C7+
cMW = WM(NC)  ! MW C7+
! 
DO I = 1,NC
!  PRINT*, Y(I)
  TCi = TC(I)
!Y(I) = Z(I)
!
  If (Abs(TCi - 373.53D0) < 0.01D0) Then ! H2S
     fH2S = Y(I)
  ElseIf (Abs(TCi - 126.6D0) < 0.01D0) Then ! N2
     fN2 = Y(I)
  ElseIf (Abs(TCi - 304.21D0) < 0.01D0) Then ! CO2
     fCO2 = Y(I)
  ElseIf (Abs(TCi - 190.56D0) < 0.01D0) Then ! C1
     fC1 = Y(I)
  ElseIf (Abs(TCi - 305.32D0) < 0.01D0) Then ! C2
     fC2 = Y(I)
  ElseIf (Abs(TCi - 369.83D0) < 0.01D0) Then ! C3
     fC3 = Y(I)
  ElseIf (Abs(TCi - 408.14D0) < 0.01D0) Then ! iC4
     fiC4 = Y(I)
  ElseIf (Abs(TCi - 425.12D0) < 0.01D0) Then ! nC4
     fC4 = Y(I)
  ElseIf (Abs(TCi - 460.43D0) < 0.01D0) Then ! iC5
     fiC5 = Y(I)
  ElseIf (Abs(TCi - 469.7D0) < 0.01D0) Then ! nC5
     fC5 = Y(I)
  ElseIf (Abs(TCi - 507.6D0) < 0.01D0) Then ! nC6
     fC6 = Y(I)
  ElseIf (Abs(TCi - 507.6D0) > 0.01D0) Then !C7+
     fC7p = Y(I)
     !PRINT*, Y(I)
     !PRINT*, TCi
  End If
  !
Enddo
!
PdPsi = A0 + A1 * TF + (A2 * fH2S + A3 * fC3 + A4 * fC2 + A5 * fC1 + A6 * fC2 + A7 * fC3 + A8 * fC4 + &
          A9 * fC5 + A10 * fC6 + A11 * fC7p + A12 * cMW + A13 * cSG) + &
         A14 * (fC7p * cMW) + A15 * (cMW / cSG) + A16 * (fC7p * cMW / cSG) + &
         A17 * (fC7p / (fC1 + fC2)) + A18 * (fC7p / (fC3 + fC4 + fC5 + fC6))
!
!PRINT*, fH2S, FN2,FCO2,FC1,FC2,FC3,FIC4,FC4,FIC5,FC5,FC6,FC7p
!
Pd = PdPsi / 14.5038D0 !Conversão psia -> bar
!
PCm = 0D0
TCm = 0.0D0
TBm = 0.0D0
DO I = 1,NC
   !Y(I) = Z(I)
   PCm = PCm + Y(I) * PC(I)
   TCm = TCm + Y(I) * TC(I)
   TBm = TBm + Y(I) * TB(I)
   If (Abs(TC(I) - 126.6D0) < 0.01D0) Then  ! H2S
     TBm = TBm - Y(I) * TB(I)  ! Eq(16) Aghamiri et al (2018)
   ElseIf (TC(I) > 508D0) Then  ! > nC6 -  PSEUDO
     TBm = TBm - Y(I) * TB(I)  ! Eq(16) Aghamiri et al (2018)
   Endif
   !
ENDDO
!
Bet = (TBm / T) ** (T / TCm)  ! Eq(16) Aghamiri et al (2018)
Wmtd60 = WM(NC) * d60(NC)
Argm = Wmtd60 / (Trk - 460.0D0)
!
ALFA(1) = 6124.3049D0
ALFA(2) = -2753.2538D0
ALFA(3) = 415.42049D0
Sum3 = 0.0D0
DO I = 1,3
  Sum3 = Sum3 + ALFA(I) * (Argm) ** (I)
ENDDO
!!Sum3 = Sum3 + 6124.3049 * (Argm) - 2753.2538 * (Argm) ** 2 + 415.42049 * (Argm) ** 3
Pk = (-2381.8542D0) + 46.341487D0 * Wmtd60 + Sum3 ! Eq(19) Aghamiri et al (2018) - Psia
Pk = Pk / 14.5038D0 !Conversão psia -> bar
!
Sumxl = 0.0D0
DO I = 1,NC
 Tri = T / TC(I)
 ckWls(I) = (Pk / Pd) ** (1D0 + 2D0 * Bet) * ((PC(I) / Pd) ** Bet) * Exp(5.3727D0 * Bet * (1.0D0 + W(I)) * (1.0D0 - 1.0D0 / Tri))
 X(I) = Y(I) / ckWls(I)
 Sumxl = Sumxl + X(I)
ENDDO
DO I = 1,NC
   X(I) = X(I) / Sumxl
ENDDO
!
P_orv = Pd
!Print*, P_orv
RETURN
END
!=============================================================================================================================================
!=============================================================================================================================================
SUBROUTINE XYFUN(NC, Dx, X, X1)
!IMPLICIT REAL*8(A-H,O-Z) 
!
INTEGER I, NC
Real*8 Dx, X(20), X1(20)
!
Dx = 0.0D0
Do I = 1,NC
   Dx = Dx + Abs(X(I) - X1(I))
ENDDO
!
RETURN
END   
!=============================================================================================================================================
!=============================================================================================================================================
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
         RAZAO(I) = (PC(I)/Pi)*(DEXP(5.373D0*(1D0+W(I))*(1D0-(TC(I)/T)))) ! Ki
         X(I) = Y(I)/RAZAO(I)
         !CORRELAÇÃO DE WILSON
     END DO
     X(1:NC) = X(1:NC)/SUM(X(1:NC)) 
     !WRITE(*,'(2/,2X,"Razao:", 6X, "PC", 6X, "W:", 6x, "TC:")')
     !WRITE(*, *)
END    
    
!****************************************************************************************************************************************************
!
!     SUB-ROTINA: FLASH (pág. 412 do Smith - Van Ness)
!
!     OBJETIVO: Cálculo de flash  
!
!****************************************************************************************************************************************************
    
SUBROUTINE FLASH( NC, TotalN, P, T, Z, X, Y, TotalNL, TotalNV, ZL, ZV, VL, VV )

     IMPLICIT REAL*8(A-H,O-Z)
     
     
     COMMON /DADOS/ TC(20), PC(20), W(20)  
     COMMON /TTIME/ Ti_bol, Ti_orv
     COMMON /MODEL/ MODELV,ID
     COMMON /MRULE/ RMIX
     COMMON /DIMENSIONAMENTO/ GRAND
     COMMON/BETAA/BETA
     DIMENSION :: XI(20,100000), YI(20,100000), YXI(20,100000), &
                  X(20), Y(20), Z(20), FL(20), FV(20), CFL(20), CFV(20)
     DIMENSION :: F(100000), dF(100000), Fvapor(100000), Fvap(100000)
     !DIMENSION CEQ(20)  !!!!! INCLUIDO
     
     INTEGER :: PHA
     
     REAL, PARAMETER :: R = 83.1451D0   ! (bar.cm3.mol−1.K−1)
     
!**************************************************************************      
    
     GRAND=1D0
     !VARIÁVEL PARA AUMENTAR A BASE MOLAR DA FVAP (BETA)
     CC2 = 3; CC = 1 
     Fvapor(1) = 1D0*GRAND 
     Fvap(1) = 1D0*GRAND
     Fvap(2) = BETA*GRAND 
     !BETA É O ÚLTIMO FVAP CALCULADO, FAZENDO QUE O PRÓXIMO FVAP SEJA MAIS FACILMENTE CONVERGIDO
     
     DO 7001 I=1,NC
         XI(I,CC) = X(I) ; YI(I,CC) = Y(I)
         !CORRELAÇÃO DE WILSON
         YXI(I,CC) = (PC(I)/P)*DEXP(5.373D0*(1D0+W(I))*(1D0-(TC(I)/T)))
         !CEQ(I) =  (PC(I)/P)*DEXP(5.373D0*(1D0+W(I))*(1D0-(TC(I)/T)))
7001 CONTINUE

7003 CONTINUE
     
!**************************************************************************
    CC2 = CC2 - 1 
    F(CC2)=1D0
    
    DO WHILE(DABS(F(CC2))>0.00001D0)
         CC2 = CC2 + 1
         
         IF(CC2>100000) THEN 
           WRITE(*,*) ' ITERATIONS EXCEEDED IN FLASH > 100000 '  !!!!!!!!!
           ! VERIFICAR O ALGORITMO 
           GO TO 7002
         END IF
         
         F(CC2) = 0D0 ; dF(CC2) = 0D0
         
         DO I=1,NC
             F(CC2) = F(CC2) + Z(I)*(YXI(I,CC)-1D0)/ &
                      (1D0+(Fvap(CC2-1)/GRAND)*(YXI(I,CC)-1D0))
             dF(CC2) = dF(CC2) - Z(I)*((YXI(I,CC)-1D0)**2D0)/ &
                       ((1D0+(Fvap(CC2-1)/GRAND)*(YXI(I,CC)-1D0))**2D0)   
         END DO
         
         CALL AUREAFLASH(TOL,ALI0,ALS0,REDUC,Z,T,Fvap(CC2-1),YXI,P,F(CC2),dF(CC2),CC,nc)
         Fvap(CC2) = Fvap(CC2-1) - (F(CC2)/dF(CC2))*REDUC !; PRINT*, Fvap(CC2)       
         
         !!== o proprio metodo da secao aurea pode calcular o valor de BETA [0 - 1]  => Fvap2
         !CALL AUREAFLASH(TOL,ALI0,ALS0,REDUC,Z,T,Fvap2,YXI,P,F2,dF2,CC,nc)
         !Fvap(CC2) = Fvap2 !; PRINT*, Fvap(CC2)
         
     END DO
!**************************************************************************     
     IF(CC<100000) THEN
         CC = CC + 1
     ELSE 
         GO TO 7002
     END IF
    
     Fvapor(CC) = Fvap(CC2) !; PRINT*, Fvapor(CC)
     XTOTAL1 = 0D0 ; XTOTAL2 = 0D0 ; YTOTAL1 = 0D0 ; YTOTAL2 = 0D0
     F1 = 0D0  !! INCLUIDO
     
     DO I=1,NC
        ! Aqui deve ser usado o valor de Fvap2 ou apenas Fvap(CC) sem divisao por GRAND
        !XI(I,CC) = Z(I)/(1D0+Fvapor(CC)*(YXI(I,CC-1)-1D0)) 
        XI(I,CC) = Z(I)/(1D0+(Fvapor(CC)/GRAND)*(YXI(I,CC)-1D0)) 
        XTOTAL1 = XTOTAL1 + XI(I,CC)   
        
        !! A funcao a seguir e a funcao do flash, e a que deve ser monitorada
        !F1 = F1 + Z(I)*(YXI(I,CC)-1D0)/ (1D0+Fvapor(CC)*(YXI(I,CC)-1D0))    !! INCLUIDO
     END DO  
     
     DO I=1,NC
        XI(I,CC) = XI(I,CC)/XTOTAL1 !NORMALIZAÇÃO DE X
        XTOTAL2 = XTOTAL2 + DABS(XI(I,CC)-XI(I,CC-1)) 
     END DO
     
     DO I=1,NC
        YI(I,CC) = XI(I,CC)*YXI(I,CC-1)
        YTOTAL1 = YTOTAL1 + YI(I,CC)
     END DO     
     
     DO I=1,NC   
        YI(I,CC) = YI(I,CC)/YTOTAL1 !NORMALIZAÇÃO DE Y
        YTOTAL2 = YTOTAL2 + DABS(YI(I,CC)-YI(I,CC-1))
     END DO
     
     FTOTAL = DABS(Fvapor(CC)-Fvapor(CC-1))
     
     DO I=1,NC
        X(I) = XI(I,CC) ; Y(I) = YI(I,CC) 
     END DO
     
     PHA = 1
     
     IF(MODELV==2) THEN
        CALL PR0011( NC, PHA, P, T, X, Y, ZL, VL, FL )
     ELSE IF(MODELV==1) THEN
        CALL PHASES( NC, PHA, P, T, X, Y, ZL, VL, FL )
     
     END IF 
     
     PHA = 2
     
     IF(MODELV==2) THEN
        CALL PR0011( NC, PHA, P, T, X, Y, ZV, VV, FV )
     ELSE IF(MODELV==1) THEN
        CALL PHASES( NC, PHA, P, T, X, Y, ZV, VV, FV )
     END IF
     
     DO I=1,NC
         ! Fases corrigido fornece o coeficiente de fugacidade, nao precisa dividir por ( X(I)*P ou Y(I)*P )
        CFL(I) = FL(I)/(X(I)*P) !COEF. DE FUGACIDADE
        CFV(I) = FV(I)/(Y(I)*P) !COEF. DE FUGACIDADE
        YXI(I,CC) = CFL(I)/CFV(I)
        !  As 3 linhas anteriores podem ser substituidas pela linha a seguir
        !YXI(I,CC) = FL(I)/FV(I)
     END DO
!**************************************************************************     
       ! Seria melhor usar o teste com F1 anterior, equacao do Flash
     IF(XTOTAL2 < 1D-5 .AND. YTOTAL2 < 1D-5 .AND. FTOTAL < 1D-5) THEN
        GO TO 7002
     ELSE
          ! Os novos compiladores nao aceitam esse tipo de retorno
        !GO TO 7001 
        GO TO 7003  ! INCLUIDO
     END IF

7002 CONTINUE       
!**************************************************************************
     TotalNV = TotalN*Fvapor(CC)/GRAND
     !TotalNV = TotalN*Fvapor(CC)
     TotalNL = TotalN - TotalNV 
     BETA= FVAPOR(CC)/GRAND
     !BETA = FVAPOR(CC)
     
END

!****************************************************************************************************************************************************
!
!     SUB-ROTINA: BOLHAT 
!
!     OBJETIVO: Cálculo do ponto de bolha (P e X conhecidos)   
!
!****************************************************************************************************************************************************
!
!SUBROUTINE BOLHAT( NC, X, P, T_bol, Y, Tinicial )   
!
!     IMPLICIT REAL*8(A-H,O-Z)
!     
!     COMMON /DADOS/ TC(20), PC(20), W(20)    
!     COMMON /MODEL/ MODELV,ID
!     COMMON /MRULE/ RMIX
!     
!     DIMENSION X(20), Y(20), FL(20), FV(20), FL1(20), &
!               FV1(20), CFL(20), CFV(20), CFL1(20), CFV1(20),  &
!               RAZAO(20), RAZAO1(20), SOMA(10000), SOMA1(10000)  
!     
!     INTEGER :: PHA
!     
!!**************************************************************************
!     !ESTIMATIVAS INICIAIS
!     SOMA=0D0 ; SOMA1=0D0
!     T = Tinicial
!     CC = 1
!     DO I=1,NC
!         RAZAO(I) = (PC(I)/P)*DEXP(5.373D0*(1D0+W(I))*(1D0-(TC(I)/T))) !Ki
!         Y(I) = X(I)*RAZAO(I) !CORRELAÇÃO DE WILSON
!     END DO
!
!!**************************************************************************     
!4001 CONTINUE
!     IF (CC==10000) GO TO 4006
!    
!     T1 = T + 0.0000001*T
!     P1 = P
!     PHA = 1
!     
!     IF(MODELV==2) THEN
!        CALL PR0011( NC, PHA, P, T, X, Y, ZL, VL, FL )
!        CALL PR0011( NC, PHA, P1, T1, X, Y, ZL1, VL1, FL1 )
!     ELSE IF(MODELV==1) THEN
!        CALL PHASES( NC, PHA, P, T, X, Y, ZL, VL, FL )
!        CALL PHASES( NC, PHA, P1, T1, X, Y, ZL1, VL1, FL1 )
!     END IF 
!     
!     PHA = 2
!     
!     IF(MODELV==2) THEN
!        CALL PR0011( NC, PHA, P, T, X, Y, ZV, VV, FV )
!     ELSE IF(MODELV==1) THEN
!        CALL PHASES( NC, PHA, P, T, X, Y, ZV, VV, FV )
!     END IF 
!     
!     SOMA(CC) = 0D0
!     
!     DO 4002 I=1,NC
!         
!         IF(X(I) /= 0D0) THEN
!             CFL(I) = FL(I)/(X(I)*P) 
!             CFL1(I) = FL1(I)/(X(I)*P1)
!             CFV(I) = FV(I)/(Y(I)*P) 
!             RAZAO(I) = CFL(I)/CFV(I)
!             SOMA(CC) = SOMA(CC) + RAZAO(I)*X(I)
!         END IF
!         
!4002 CONTINUE
!     
!     YTOTAL = 0D0
!     
!     DO 4003 I=1,NC
!         
!         IF(X(I) /= 0D0) THEN
!            Y(I) = RAZAO(I)*X(I)/SOMA(CC)
!         END IF
!         
!         YTOTAL = YTOTAL + Y(I)
!4003 CONTINUE
!     
!     DO 4007 I=1,NC
!         Y(I) = Y(I)/YTOTAL !NORMALIZAÇÃO
!4007 CONTINUE
!     
!     IF (CC==10000) GO TO 4006
!     
!     CC = CC + 1
!     SOMA(CC) = 0D0
!     PHA = 2
!     
!     IF(MODELV==2) THEN
!        CALL PR0011( NC, PHA, P, T, X, Y, ZV, VV, FV )
!        CALL PR0011( NC, PHA, P1, T1, X, Y, ZV1, VV1, FV1 )
!     ELSE IF(MODELV==1) THEN
!        CALL PHASES( NC, PHA, P, T, X, Y, ZV, VV, FV )
!        CALL PHASES( NC, PHA, P1, T1, X, Y, ZV1, VV1, FV1 )
!     END IF
!     
!     DO 4004 I=1,NC
!         
!         IF(X(I) /= 0D0) THEN
!             CFV(I) = FV(I)/(Y(I)*P)
!             RAZAO(I) = CFL(I)/CFV(I)
!             SOMA(CC) = SOMA(CC) + RAZAO(I)*X(I)
!             CFV1(I) = FV1(I)/(Y(I)*P1)
!             RAZAO1(I) = CFL1(I)/CFV1(I)
!             SOMA1(CC) = SOMA1(CC) + RAZAO1(I)*X(I)
!         END IF
!         
!4004 CONTINUE
!     
!     TESTE1 = DABS(SOMA(CC)-SOMA(CC-1))
!     IF(TESTE1<1D-5) GO TO 4005 
!     GO TO 4002
!4005 CONTINUE
!     
!     F = SOMA(CC)-1D0
!     F1 = SOMA1(CC)-1D0
!     dF = (F1 - F)/(T1 - T)
!      
!     IF(DABS(F)<=1D-6 .AND. DABS(T-Told)<=1D-3) THEN
!        GO TO 4006
!     ELSE
!        
!        Told = T 
!        CALL AUREABOLHAT(TOL,ALI0,ALS0,REDUC,X,Y,T,P,F,dF,nc)
!        T = T - F/dF/REDUC
!        
!        IF(T<1D2 .OR. T>3D2) T = 200D0 - 10D0*RAND(0)
!        GO TO 4001
!     END IF 
!     
!4006 CONTINUE
!     
!     T_bol = T
!END    

!****************************************************************************************************************************************************
!
!     SUB-ROTINA: ORVALT 
!
!     OBJETIVO: Cálculo do ponto de orvalho (P e Y conhecidos)   
!
!****************************************************************************************************************************************************
!
!SUBROUTINE ORVALT( NC, Y, P, T_orv, X, Tinicial )   
!
!     IMPLICIT REAL*8(A-H,O-Z)
!     
!     COMMON /DADOS/ TC(20), PC(20), W(20)    
!     COMMON /MODEL/ MODELV,ID
!     COMMON /MRULE/ RMIX
!     
!     DIMENSION X(20), Y(20), FL(20), FV(20), FL1(20), &
!               FV1(20), CFL(20), CFV(20), CFL1(20), CFV1(20),  &
!               RAZAO(20), RAZAO1(20), SOMA(10000), SOMA1(10000)  
!     
!     INTEGER :: PHA
!     
!!**************************************************************************
!     !ESTIMATIVAS INICIAIS
!     SOMA=0D0 ; SOMA1=0D0
!     T = Tinicial
!     CC = 1
!     DO I=1,NC
!         RAZAO(I) = (PC(I)/P)*DEXP(5.373D0*(1D0+W(I))*(1D0-(TC(I)/T)))
!         X(I) = Y(I)/RAZAO(I) !CORRELAÇÃO DE WILSON
!     END DO
!X(1:NC) = X(1:NC)/SUM(X(1:NC))
!     
!!**************************************************************************
!    
!6001 CONTINUE
!     IF (CC==10000) GO TO 6006
!    
!     T1 = T + 0.001*T
!     P1 = P
!     PHA = 1
!     
!     IF(MODELV==2) THEN
!        CALL PR0011( NC, PHA, P, T, X, Y, ZL, VL, FL )
!     ELSE IF(MODELV==1) THEN
!        CALL PHASES( NC, PHA, P, T, X, Y, ZL, VL, FL )
!     END IF 
!     
!     PHA = 2
!     
!     IF(MODELV==2) THEN
!        CALL PR0011( NC, PHA, P, T, X, Y, ZV, VV, FV )
!        CALL PR0011( NC, PHA, P1, T1, X, Y, ZV1, VV1, FV1 )
!     ELSE IF(MODELV==1) THEN
!        CALL PHASES( NC, PHA, P, T, X, Y, ZV, VV, FV )
!        CALL PHASES( NC, PHA, P1, T1, X, Y, ZV1, VV1, FV1 )
!     END IF 
!     
!     SOMA(CC) = 0D0
!     
!     DO 6002 I=1,NC
!         
!         IF(Y(I) /= 0D0) THEN
!             CFL(I) = FL(I)/(X(I)*P)
!             CFV(I) = FV(I)/(Y(I)*P)
!             CFV1(I) = FV1(I)/(Y(I)*P1) 
!             RAZAO(I) = CFL(I)/CFV(I)
!             SOMA(CC) = SOMA(CC) + Y(I)/RAZAO(I)
!         END IF
!         
!6002 CONTINUE
!     
!     XTOTAL = 0D0
!     
!     DO 6003 I=1,NC
!         IF(Y(I) /= 0D0) X(I) = Y(I)*SOMA(CC)/RAZAO(I)
!         XTOTAL = XTOTAL + X(I)
!
!6003 CONTINUE
!     
!     DO 6007 I=1,NC
!         X(I) = X(I)/XTOTAL !NORMALIZAÇÃO
!         
!6007 CONTINUE
!     
!     IF (CC==10000) GO TO 6006
!
!     CC = CC + 1
!     SOMA(CC) = 0D0
!     PHA = 1
!     
!     IF(MODELV==2) THEN
!        CALL PR0011( NC, PHA, P, T, X, Y, ZL, VL, FL )
!        CALL PR0011( NC, PHA, P1, T1, X, Y, ZL1, VL1, FL1 )
!     ELSE IF(MODELV==1) THEN
!        CALL PHASES( NC, PHA, P, T, X, Y, ZL, VL, FL )
!        CALL PHASES( NC, PHA, P1, T1, X, Y, ZL1, VL1, FL1 )
!     END IF
!     
!     DO 6004 I=1,NC
!         
!         IF(Y(I) /= 0D0) THEN
!             CFL(I) = FL(I)/(X(I)*P)
!             RAZAO(I) = CFL(I)/CFV(I)
!             SOMA(CC) = SOMA(CC) + Y(I)/RAZAO(I)
!             CFL1(I) = FL1(I)/(X(I)*P1)
!             RAZAO1(I) = CFL1(I)/CFV1(I)
!             SOMA1(CC) = SOMA1(CC) + Y(I)/RAZAO1(I)
!         END IF
!         
!6004 CONTINUE
!     
!     TESTE1 = DABS(SOMA(CC)-SOMA(CC-1))
!     IF(TESTE1<1D-5) GO TO 6005
!     GO TO 6002
!     
!6005 CONTINUE
!     
!     F = SOMA(CC)-1D0 
!     F1 = SOMA1(CC)-1D0 
!     dF = (F1 - F)/(T1 - T)
!     
!     IF(DABS(F)<=1D-6 .AND. DABS(T-Told)<=1D-3) THEN
!        GO TO 6006
!     ELSE
!        Told = T
!        CALL AUREAORVALT(TOL,ALI0,ALS0,REDUC,X,Y,T,P,F,dF,nc)
!        T = T - (F/dF)*REDUC 
!        
!        IF(T<3D2 .OR. T>5D2) T = 350D0 - 10D0*RAND(0)
!        GO TO 6001
!        
!     END IF   
!6006 CONTINUE
!     
!     T_orv = T
!    END  

