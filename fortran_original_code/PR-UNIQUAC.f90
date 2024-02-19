!****************************************************************************
!
!  SUBROUTINE: PR0011
!
!  PURPOSE:  Cálculo de coeficiente de fugacidade pelo modelo PR-UNIQUAC-HK                   
!
!****************************************************************************

SUBROUTINE PR0011( NM, PHA, P, T, XI, YI, Z_p, V_p, Fug_p )   

    IMPLICIT REAL*8(A-H,O-Z) 
    
    COMMON /MRULE/ RMIX
    
    DIMENSION :: ZI(20),XI(20),YI(20),ZRI(20),GAMA(20),Fug_p(20), &
                 AI(20),BI(20),CI(20)
    
    INTEGER :: PHA
    
    R = 83.14462D0 !(cm3.bar.K-1.mol-1)
    
    ! VARIÁVEIS PRINCIPAIS --------------------------------------------------
    !NC    ------ Número de componentes (espécies) no sistema
    !NM    ------ Número de espécies moleculares no sistema
    !GAMA(I) ---- Coeficiente de atividade do componente I
    !T     ------ Temperatura (K)
    !P     ------ Pressão (bar)
    !ZI(I) ------ Fração molar aparente do componente I na fase p
    !ZRI(I) ----- Fração molar real do componente I na fase p
    !Fug_p(I) --- Fugacidade do componente I na fase p
    
    ! VÁRIÁVEIS AUXILIARES --------------------------------------------------
    !PHA
    
    ! INICIALIZAÇÃO ---------------------------------------------------------
    GAMA = 1D0  
    NC = NM
    IF(PHA==1) ZI = XI 
    IF(PHA==2) ZI = YI    
    ZRI = ZI

    ! CÁLCULO DOS COEFICIENTES DE ATIVIDADE --------------------------------- 
    IF(RMIX==1) THEN
        DO I=1,NM
            IF(DABS(ZI(I)-1D0)<1D-10 .OR. PHA==2) GO TO 001
        END DO
        CALL UNIQUAC(NC,T,ZRI,GAMA)  
        DO I=1,NC
            IF(ZRI(I)==0D0) GAMA(I) = 1D0
        END DO
001     CONTINUE
    END IF

    !DO I=1,NC
    !    WRITE(*,"(I2,2F16.8)") I, ZRI(I), GAMA(I)
    !END DO
    !PAUSE
    
    ! CÁLCULO DOS PARÂMETROS ENERGÉTICOS E CO-VOLUMES DAS ESPÉCIES PURAS ----
    DO 002 I=1,NM
        CALL ABPRAL( I, T, a, b )   
        AI(I) = a
        BI(I) = b 
002 CONTINUE

    ! CÁLCULO DA FUGACIDADE DAS ESPÉCIES MOLECULARES -------------------------
    IF(RMIX==1) CALL HKOKAL(NC, NM, T, P, ZRI, AI, BI, GAMA, V_p, Fug_p, PHA)
    IF(RMIX==2) THEN
        CALL FUGCHK(PHA, NM, T, P, ZRI, BI, AI, V_p, Fug_p)   
        DO 003 I=1,NC
            Fug_p(I) = Fug_p(I)*(ZI(I)*P) 
003     CONTINUE
    END IF
    
    ! CÁLCULO DA TRANSLAÇÃO VOLUMÉTRICA -------------------------------------
    CS = 0D0
    DO 004 I=1,NM
         ! Passagem de parametro CI(I), nem todo compilador entende bem
         !!!!CALL VTPELL( I, T, CI(I) ) 
         CALL VTPELL( I, T, CCI ) 
         CI(I) = CCI  !!!!
         CS = CS + ZI(I)*CI(I)
004 CONTINUE
    V_p = V_p - CS
    Z_p = V_p/(R*T/P)
    
END

!****************************************************************************
!
!  PROGRAM: UNIQUAC
!
!  PURPOSE:  Cálculo de G de excesso e coeficiente de atividade.
!
!****************************************************************************

SUBROUTINE UNIQUAC(NC,T,X,GAMA)
    
    IMPLICIT REAL*8(A-H,O-Z)
    
    COMMON /UNIQ/ Uij_0(20,20), Uij_T(20,20), r(20), q(20)
    
    DIMENSION GAMA(20), GAMAC(20), GAMAR(20), X(20), FS(20), FV(20), &
              Uij(20,20), TAU(20,20), ALNGC(20), ALNGR(20), S(20),   &
              SOMA_q(20), SOMA_r(20), SOMA_R2(20)
    
    ! VARIÁVEIS PRINCIPAIS --------------------------------------------------
    !NC    ------ Número de componentes (espécies) no sistema
    !GE    ------ Energia livre de Gibbs molar
    !GEC   ------ Termo combinatorial
    !GER   ------ Termo residual
    !GAMA(I) ---- Coeficiente de atividade da espécie I
    !GAMAC(I) --- Termo combinatorial 
    !GAMAR(I) --- Termo residual      
    !RG    ------ Constante universal dos gases
    !T     ------ Temperatura
    !X(I)  ------ Fração molar da espécie I
    !q(I)  ------ Parâmetro de superfície da espécie I
    !r(I)  ------ Parâmetro de volume da espécie I
    !Z     ------ Número de coordenação
    !FS(I) ------ Fração de superfície da espécie I
    !FV(I) ------ Fração de volume da espécie I
    !Uij(I,J) --- Parâmetro de interação binária entre espécies I e J
    
    ! VARIÁVEIS AUXILIARES --------------------------------------------------
    !ALNGC(I) 
    !ALNGR(I) 
    !S(I)
    !TAU(I,J)
    !SOMA_q(I)
    !SOMA_r(I)
    !SOMA_C1
    !SOMA_C2
    !SOMA_R1
    !SOMA_R2(I)

    ! PARÂMETROS DE INTERAÇÃO BINÁRIA ---------------------------------------
    DO I=1,NC
        DO J=1,NC
            Uij(I,J) = Uij_0(I,J) + Uij_T(I,J)*(T-298.15D0)
        END DO
    END DO

    ! CÁLCULOS --------------------------------------------------------------
    Z = 10D0
    RG = 8.314462D0 !(J/mol-K)
    
    ! CÁLCULO DA FRAÇÃO DE SUPERFÍCIE
    SOMA_q = 0D0
    DO 101 I=1,NC
        DO 100 J=1,NC
            SOMA_q(I) = SOMA_q(I) + X(J)*q(J)
100     CONTINUE
    FS(I) = X(I)*q(I)/SOMA_q(I)
101 CONTINUE
    
    ! CÁLCULO DA FRAÇÃO DE VOLUME
    SOMA_r = 0D0
    DO 103 I=1,NC
        DO 102 J=1,NC
            SOMA_r(I) = SOMA_r(I) + X(J)*r(J)
102     CONTINUE
    FV(I) = X(I)*r(I)/SOMA_r(I)
103 CONTINUE   
    
    ! CÁLCULO DO TERMO COMBINATORIAL
    SOMA_C1 = 0D0 ; SOMA_C2 = 0D0 ; GAMAC = 1D0
    DO 104 I=1,NC
        IF(X(I)<1D-15) GO TO 104
        SOMA_C1 = SOMA_C1 + X(I)*DLOG(FV(I)/X(I))
        SOMA_C2 = SOMA_C2 + X(I)*q(I)*DLOG(FV(I)/FS(I))
        ALNGC(I) = DLOG(FV(I)/X(I)) + 1D0 - FV(I)/X(I)   &
                   - (Z/2D0)*q(I)*(DLOG(FV(I)/FS(I)) + 1D0 - FV(I)/FS(I))
        GAMAC(I) = DEXP(ALNGC(I))
104 CONTINUE  
    GEC = RG*T*(SOMA_C1 - (Z/2D0)*SOMA_C2)
    
    ! CÁLCULO DE TAU(I,J)
    TAU = 0D0
    DO 106 I=1,NC
        DO 105 J=1,NC            
            TAU(I,J) = DEXP(-Uij(I,J)/T)
105     CONTINUE
106 CONTINUE    
    
    ! CÁLCULO DE S(I)
    S = 0D0
    DO 108 I=1,NC
        DO 107 J=1,NC            
            S(I) = S(I) + FS(J)*TAU(J,I)
107     CONTINUE
108 CONTINUE
    
    ! CÁLCULO DO TERMO RESIDUAL
    SOMA_R1 = 0D0 ; SOMA_R2 = 0D0
    DO 110 I=1,NC
        SOMA_R1 = SOMA_R1 + X(I)*q(I)*DLOG(S(I))
        DO 109 J=1,NC
            SOMA_R2(I) = SOMA_R2(I) + TAU(I,J)*FS(J)/S(J)
109     CONTINUE
        ALNGR(I) = q(I)*(1D0 - DLOG(S(I)) - SOMA_R2(I))
        GAMAR(I) = DEXP(ALNGR(I))
110 CONTINUE
    GER = - RG*T*SOMA_R1
    
    ! CÁLCULO DE G DE EXECESSO E COEFICIENTES DE ATIVIDADE
    GE = GEC + GER
    DO 111 I=1,NC
        GAMA(I) = GAMAC(I)*GAMAR(I)  ! CONVENÇÃO SIMÉTRICA
111 CONTINUE

END SUBROUTINE
    
!****************************************************************************
!
!  SUBROUTINE: HKOKAL
!
!  PURPOSE:  Cálculo de ELV pela equação PR via metodologia de 
!            Heidemann-Kokal (1990) (PR + UNIQUAC + Debye-Hückel)
!
!****************************************************************************

SUBROUTINE HKOKAL(NC, NM, T, P, X, A, B, GAMA, V_p, Fug, PHA) 
    
    IMPLICIT REAL*8(A-H,O-Z)
    
    DIMENSION GAMA(20), X(20), A(20), B(20), ALFA(20), QSI(20), FI(20), &
              D1(20), D2(20), D4(20), D5(20), G(20), H(20), ALNFI(20), Fug(20)
    
    DIMENSION :: G1(20), G2(20), G3(20), G4(20), G5(20), AIJ(20,20), AS(20), &
                 H10(20), B1(20), B2(20), B3(20), B4(20), ZZ(3)
    
    COMMON /RGAS  / RG
    COMMON /CC    / C1,C2
    
    INTEGER :: PHA, CCC, NR
    
    RG = 83.14462D0 !(cm3.bar.K-1.mol-1)
    
    ! VARIÁVEIS PRINCIPAIS --------------------------------------------------
    !NC    ------ Número de componentes (espécies) no sistema
    !NM    ------ Número de espécies moleculares no sistema
    !GAMA(I) ---- Coeficiente de atividade do componente I
    !RG    ------ Constante universal do gases
    !T     ------ Temperatura (K)
    !P     ------ Pressão (bar)
    !X(I)  ------ Fração molar do componente I
    !A(I)  ------ Parâmetro de energia do componente I puro
    !AM    ------ Parâmetro de energia da mistura
    !B(I)  ------ Co-volume do componente I puro
    !BM    ------ Co-volume da mistura
    !ALFA(I) ---- a/bRT para o componente I puro
    !ALFAM ------ a/bRT para a mistura
    !QSI(I) ----- b/v para o componente I puro
    !QSIM  ------ b/v para a mistura
    !FI(I) ------ Coeficiente de fugacidade do componente I na mistura
    !V_p   ------ Volume molar da fase
    !DGE   ------ Variação da energia livre de Gibbs na mistura
    !Fug(I) ----- Fugacidade do componente I na mistura
    
    ! VARIÁVEIS AUXILIARES --------------------------------------------------
    !ALFAMIN
    !C1
    !C2
    !D1(I)
    !D2(I)
    !BETA
    !DELTA
    !F1
    !F2
    !G(I)
    !SOMA1
    !H(I)
    !GM
    !HM
    !FMIN
    !QSIMIN
    !DELTAF
    !D3
    !D4(I)
    !D5(I)
    !ALNFI(I)
    !SOMA0
    !GX, HX e BX
    !CCC
    !PRINT*, 'PHA:', PHA
    ! CÁLCULO DE DGE
    DGE = 0D0
    DO I=1,NC
        DGE = DGE + RG*T*X(I)*DLOG(GAMA(I)) 
    END DO
    
    ! NORMALIZAÇÃO DAS FRAÇÕES MOLARES --------------------------------------
    SOMA0 = 0D0
    DO I=1,NM
        SOMA0 = SOMA0 + X(I)
    END DO
    DO I=1,NM
        X(I) = X(I)/SOMA0  ! Frações molares aparentes
    END DO
    
    ! MÉTODO DE HEIDEMANN-KOKAL ---------------------------------------------
    
    ! CÁLCULO DE QSI(I)
    C1 = 1D0 + DSQRT(2D0) ; C2 = 1D0 - DSQRT(2D0)
    ALFAMIN = 1.15D0*((2D0 + C1 + C2) + 2D0*DSQRT((C1+1D0)*(C2+1D0)))
    BETA = -0.030082157186D0 ; DELTA = 0.003041741701D0
    DO 001 I=1,NM
        ALFA(I) = A(I)/(B(I)*RG*T)
        D1(I) = (C1 + C2) - ALFA(I)
        D2(I) = C1*C2 + ALFA(I)
        IF(ALFA(I) >= ALFAMIN) THEN
            QSI(I) = (- D1(I) + DSQRT(D1(I)**2D0-4D0*D2(I)))/(2D0*D2(I))
        ELSE
            QSI(I) = 1D0 + BETA*ALFA(I)**2D0 + DELTA*ALFA(I)**3D0
        END IF
001 CONTINUE
    
    ! CÁLCULO DE BM
    BM = 0D0
    DO 002 I=1,NM
        BM = BM + X(I)*B(I)
002 CONTINUE
    
    ! CÁLCULO DE F
    SOMA1 = 0D0 !; DGE = 0D0
    DO 003 I=1,NM
        G(I) = (1D0/(C2-C1))*DLOG((1D0+C2*QSI(I))/(1D0+C1*QSI(I)))/QSI(I)
        SOMA1 = SOMA1 + X(I)*(ALFA(I)*G(I)*QSI(I)+DLOG((1D0-QSI(I))/QSI(I)) &
                +DLOG(B(I)/BM))
        !DGE = DGE + RG*T*X(I)*DLOG(GAMA(I)) 
003 CONTINUE
    F1 = - DGE/(RG*T) + SOMA1
    
    ! CÁLCULO DE QSIM
    HM = 0D0
    DO 004 I=1,NM
        H(I) = 1D0/(1D0-QSI(I))-ALFA(I)*QSI(I)/ &
               ((1D0+C1*QSI(I))*(1D0+C2*QSI(I)))
        HM = HM + X(I)*H(I)
004 CONTINUE 
    IF(DABS(HM)<1D-5) THEN
        FMIN = 2.4744D0
        QSIMIN = DSQRT(2D0)-1D0
    ELSE
        FMIN = 2.47D0-1.25D0*HM
        QSIMIN = 0.414D0+0.015D0*HM
    END IF   
    IF(F1<=FMIN) THEN 
        !PRINT*, 'ERRO 1 (HK)' 
        !PRINT*, 'F1 = ', F1
        !PRINT*, 'Fmin = ', FMIN
        !PAUSE
        GO TO 006
    ELSE
        QSIM = 0.9999D0 ; DELTAF = 1D0 ; CCC = 1
        DO WHILE(DABS(DELTAF)>1D-10 .AND. CCC<=100)
            GM = (1D0/(C2-C1))*DLOG((1D0+C2*QSIM)/(1D0+C1*QSIM))/QSIM
            F2 = (1D0+C1*QSIM)*(1D0+C2*QSIM)*GM*(1D0/(1D0-QSIM)-HM)+ &
                    DLOG((1D0-QSIM)/QSIM)
            DELTAF = F2 - F1
            D3 = -HM*(1D0+(C1*C2*QSIM**2D0-1D0)*GM)+GM*              &
                (-1D0+2D0*QSIM+(C1+C2+C1*C2)*QSIM**2D0)/((1D0-QSIM)**2)
            QSIM = QSIM - DELTAF/(D3/QSIM)
            CCC = CCC + 1
        END DO
        IF(QSIM < QSIMIN .OR. QSIM > 1D0) THEN 
            !PRINT*, 'ERRO 2 (HK)' 
            !PRINT*, 'QSI = ', QSIM
            !PRINT*, 'QSImin = ', QSIMIN
            !PAUSE
        END IF
    END IF
        
    ! CÁLCULO DE ALFAM
    ALFAM = (1D0/(1D0-QSIM)-HM)*(1D0+C1*QSIM)*(1D0+C2*QSIM)/QSIM
    AM = ALFAM*BM*RG*T
    
    ! CÁLCULO DO COEFICIENTE DE FUGACIDADE    
    r1 = - C1
    r2 = - C2 
    !PRINT*, r1, r2, RG
    !PRINT*, T, P
    !PRINT*, BM, AM
    A12 = - BM*(r1*r2*(BM**2) + r1*r2*BM*RG*T/P + AM/P)      ;!PRINT*, 'A10:', A12
    A11 = (BM**2)*(r1*r2+r1+r2) + RG*T*BM*(r1+r2)/P + AM/P   ;!PRINT*, 'A11:', A11
    A10 = - (BM*(r1+r2+1D0) + RG*T/P)                        ;!PRINT*, 'A12:', A10
    CALL CUBICA(A10, A11, A12, ZZ, NR)
    IF(NR<=0D0) THEN
        IF(PHA==2) V_p = ZZ(1)
        IF(PHA==1) V_p = ZZ(3)
    ELSE
        V_p = ZZ(1)
    END IF
    Z_p = V_p/(RG*T/P)
    !PRINT*, Z_p;pause
    DO 005 I=1,NM 
        !Cálculo baseado em parâmetros de componente puro e BM -----
        G1(I) = -DLOG(GAMA(I))  
        G2(I) = QSI(I)*ALFA(I)*G(I)
        G3(I) = DLOG((1D0-QSI(I))/QSI(I))
        G4(I) = DLOG(B(I)/BM)
        G5(I) = -(B(I)-BM)/BM
        D4(I) = G1(I) + G2(I) + G3(I) + G4(I) + G5(I)
        !Cálculo baseado no estado de referência P=0 ---------------
        H1 = (C1 + C2 + C1*C2)*QSIM*QSIM + 2D0*QSIM - 1D0 
        H2 = (1D0-QSIM)**2D0                         
        H3 = (1D0+C2*QSIM)/(1D0+C1*QSIM)                   
        H4 = (1D0+C2*QSIM)*(1D0+C1*QSIM)                   
        H5 = DLOG (H3)/(C2-C1)/QSIM                  
        H6 = H5*H1/H2 - ((C1*C2*QSIM*QSIM-1D0)*H5+1D0)*HM    
        H7 = HM/H6                                         
        H8 = 1D0/(1D0-QSIM)                             
        H9 = DLOG((1D0-QSIM)/QSIM)  
        H10(I) = ((D4(I)-H9)*(1D0+H7)+H5*H4*(H(I)-H8)*H7)/(QSIM*H5)
        !Cálculo baseado nas condições de mistura ------------------
        H11 = BM / V_p
        H12 = (1D0+C2*H11)/(1D0+C1*H11)                    
        H13 = DLOG (H12)/(C2-C1)  
        D5(I) = H13*H10(I) 
        B1(I) = DLOG(V_p/(V_p-BM)) 
        B2(I) = B(I)/(V_p-BM)
        B3(I) = -ALFAM*B(I)*V_p/((V_p+C1*BM)*(V_p+C2*BM))
        B4(I) = -DLOG(Z_p)
        ALNFI(I) = B1(I) + B2(I) + B3(I) + B4(I) - D5(I) 
        FI(I) = DEXP(ALNFI(I))
        Fug(I) = X(I)*P*FI(I)
005 CONTINUE
    
006 CONTINUE
    IF(F1<=FMIN) THEN
        AS=0D0 ; AM=0D0
        DO 008 I=1,NM
             DO 007 J=1,NM
                AIJ(I,J) = DSQRT(A(I)*A(J)) 
                AM = AM + X(I)*X(J)*AIJ(I,J)
                AS(I) = AS(I) + X(J)*AIJ(I,J)
    007      CONTINUE
    008 CONTINUE
        r1 = - C1
        r2 = - C2 
        A12 = - BM*(r1*r2*(BM**2) + r1*r2*BM*RG*T/P + AM/P)
        A11 = (BM**2)*(r1*r2+r1+r2) + RG*T*BM*(r1+r2)/P + AM/P
        A10 = - (BM*(r1+r2+1D0) + RG*T/P)  
        CALL CUBICA(A10, A11, A12, ZZ, NR)
        IF(NR<=0D0) THEN
            IF(PHA==2) V_p = ZZ(1)
            IF(PHA==1) V_p = ZZ(3)
        ELSE
            V_p = ZZ(1)
        END IF
        Z_p = V_p/(RG*T/P)
        a2 = P*AM/((RG*T)**2D0) 
        y2 = P*BM/(RG*T)   
        DO I=1,NM 
            Fug(I) = X(I)*P*DEXP((B(I)/BM)*(Z_p-1D0)-DLOG(Z_p-y2)-   &
                     (a2/y2)*(1D0/(C1-C2))*(2D0*AS(I)/AM - B(I)/BM)* & 
                     DLOG((Z_p+C1*y2)/(Z_p+C2*y2))) 
            FI(I) = Fug(I)/X(I)/P
            ALNFI(I) = DLOG(FI(I))
        END DO
    END IF
    
END SUBROUTINE
    
!****************************************************************************
!
!     SUB-ROTINA: ABPRAL
!
!     OBJETIVO: Cálculo dos parâmetros atrativo (a) e repulsivo (b).
!               PENG E ROBINSON 1978 / ALMEIDA ET AL. 1991
!
!**************************************************************************** 

SUBROUTINE ABPRAL(I, T, a, b)

    IMPLICIT REAL*8(A-H,O-Z)
     
    COMMON /DADOS/ TC(20), PC(20), W(20)     
    COMMON /ABM/   ALFAF, AAT1(20), AAT2(20), AAT3(20)
    
    REAL, PARAMETER :: R = 83.1451D0   ! (bar.cm3.mol−1.K−1)
    
    ohm_a = 0.457235528921D0
    ohm_b = 0.0777960739039D0
    
    IF(ALFAF==1) THEN
        IF (W(I) <= 0.491) THEN
            SO1 = 0.37464D0 + 1.54226D0*W(I) - 0.26992D0*W(I)**2D0  
        ELSE
            SO1 = 0.379642D0 + 1.487503D0*W(I) - 0.164423D0*W(I)**2D0 + 0.016666*W(I)**3D0 
        END IF
        alfa = (1D0+SO1*(1D0-DSQRT(T/TC(I))))**2D0
    ELSE    
        alfa = DEXP(AAT1(I)*(1D0-T/TC(I))*(DABS(1D0-T/TC(I)))**(AAT2(I)-1)+AAT3(I)*(TC(I)/T-1D0))
        IF(T==TC(I)) alfa = 1D0
    END IF
    
    a = ohm_a*((R**2D0)*(TC(I)**2D0)/PC(I))*alfa
    b = ohm_b*R*TC(I)/PC(I)

END  
    
!****************************************************************************
!
!     SUB-ROTINA: VTPELL
!     L.A. Pellegrini, S. Moiolia, S. Gambaa, P. Ceragioli (2012)
!     OBJETIVO: Volume Translation
!
!**************************************************************************** 

SUBROUTINE VTPELL( I, T, c )

    IMPLICIT REAL*8(A-H,O-Z)
     
    COMMON /DADOS/ TC(20), PC(20), W(20)    
    COMMON /VTR/   APE1(20), APE2(20), APE3(20)
    COMMON /VTPRR/ VTPR
    !!!!DIMENSION C(20)
    REAL, PARAMETER :: R = 83.1451D0   ! (bar.cm3.mol−1.K−1)
   
    TR = T/TC(I)
    IF (VTPR==1) THEN
        !!!!c(i)=0d0
        c=0d0
    ELSE IF (VTPR==2) THEN
     !!!!c(i) = APE1(I)*TR**2D0+APE2(I)*TR+APE3(I)
     c = APE1(I)*TR**2D0+APE2(I)*TR+APE3(I)
    END IF
    !c = APE1(I)*TR**2D0+APE2(I)*TR+APE3(I)
END  
    
!****************************************************************************
!
!     SUB-ROTINA: CUBICA
!
!     OBJETIVO: Solução da equação cúbica.
!
!**************************************************************************** 
      
SUBROUTINE CUBICA(A, B, c, X, nr)
      !
      implicit none
       real*8 ex, tres, g, r, d, q, s, z1, z2, u, z3, theta, vacos, te 
       real*8 pi, X(3), A, B, c
       integer nr
       !
       pi = 3.1415926535897932384626433832795d0
       te = 0.d0
       ex = 1.0 / 3.0
       tres = dsqrt(3.d0)
       X(1) = 0.0d0
       X(2) = 0.0d0
       X(3) = 0.0d0
       q = (3.0 * B - (A ** 2)) / 9.0
       r = (9.0 * A * B - 27.0 * c - 2.0 * A ** 3) / 54.0
       d = q ** 3 + r ** 2
       If (dAbs(d) .le. 0.00000001)  d = 0.0
      !
      If (d .eq. 0.0) Then
       nr = 0
       If (r .ge. 0) Then 
         s = dAbs(r) ** ex 
       Else 
         s = -dAbs(r) ** ex
       Endif
        !   s = Sign((Abs(r) ** ex), r)
       z1 = 2.0 * s - A / 3.0
       If (z1 .le. 0.00000001)  z1 = 0.0
       z2 = -s - A / 3.0
       If (z2 .le. 0.00000001)  z2 = 0.0
       X(1) = dmax1(z1, z2)
       X(3) = dmin1(z1, z2)
       X(2) = X(3)
      ElseIf (d .gt. 0.0) Then
       nr = 1
       u = r + dsqrt(d)
       If (u .ge. 0) Then 
         s = dAbs(u) ** ex 
       Else 
         s = -dAbs(u) ** ex
       Endif
       If (dAbs(s) .le. 0.00000001)  s = 0.0
       u = r - Sqrt(d)
       If (u .ge. 0) Then 
         te = dAbs(u) ** ex 
       Else 
         te = -dAbs(u) ** ex
       Endif
       If (dAbs(te) .le. 0.00000001)  te = 0.0
       X(1) = s + te - A / 3.0
       X(2) = 0.0d0
       X(3) = 0.0d0
       return
      ElseIf (d .lt. 0.0) Then
       nr = -1
       u = 2.0 * dSqrt(-q)
       vacos = r / dsqrt(-q ** 3)
       If (vacos .gt. 1.0)  vacos = 1.0
       theta = dAcos(vacos) / 3.0
       If (dAbs(theta) .le. 0.00000001)  theta = 0.0
       z1 = u * dCos(theta) - A / 3.0
       If (z1 .le. 0.00000001)  z1 = 0.0
       z2 = u * dCos(theta + 2.0 * Pi / 3.0) - A / 3.0
       If (z2 .le. 0.00000001)  z2 = 0.0
       z3 = u * dCos(theta + 4.0 * Pi / 3.0) - A / 3.0
       If (z3 .le. 0.00000001)  z3 = 0.0
       X(1) = dmax1(z1, dmax1(z2, z3))
       X(3) = dmin1(z1, dmin1(z2, z3))
       If (X(1) .eq. z1)  X(2) = dmax1(z2, z3)
       If (X(1) .eq. z2)  X(2) = dmax1(z1, z3)
       If (X(1) .eq. z3)  X(2) = dmax1(z1, z2)
       return
      End If
     RETURN
END  

!**************************************************************************
!
!     SUB-ROTINA: PHASES
!
!     OBJETIVO: Cálculo da fugacidade das fases via PR EoS    
!
!**************************************************************************

SUBROUTINE PHASES( NC, PHA, P, T, XI, YI, Z_p, V_p, fug_p)   

     IMPLICIT REAL*8(A-H,O-Z) 
    
     COMMON /BIP/ DIJ(20,20)
     
     DIMENSION :: AIJ(20,20), ZI(20), XI(20), YI(20), &
                  AI(20), BI(20), CI(20), AS(20), fug_p(20), ZZ(3)
     
     INTEGER :: PHA, NR
     
     REAL, PARAMETER :: R = 83.1451D0   ! (bar.cm3.mol−1.K−1)
!**************************************************************************
     DO 2001 I=1,NC 
         CALL ABPRAL( I, T, a, b )   
         AI(I) = a
         BI(I) = b
2001 CONTINUE
!**************************************************************************
     IF(PHA==1) ZI = XI 
     IF(PHA==2) ZI = YI        
     AIJ = 0D0 ; AM = 0D0 ; BM = 0D0
     DO 2003 I=1,NC
         DO 2002 J=1,NC
             AIJ(I,J) = (1D0-DIJ(I,J))*DSQRT(AI(I)*AI(J))
             AM = AM + ZI(I)*ZI(J)*AIJ(I,J)
2002     CONTINUE
         BM = BM + ZI(I)*BI(I)
2003 CONTINUE
!**************************************************************************    
     a2 = P*AM/((R*T)**2D0) 
     b2 = P*BM/(R*T)        
     sig = 1D0 + DSQRT(2D0)
     eps = 1D0 - DSQRT(2D0)
     r1 = - sig
     r2 = - eps 
     A12 = - BM*(r1*r2*(BM**2) + r1*r2*BM*R*T/P + AM/P)
     A11 = (BM**2)*(r1*r2+r1+r2) + R*T*BM*(r1+r2)/P + AM/P
     A10 = - (BM*(r1+r2+1D0) + R*T/P)  
     CALL CUBICA(A10, A11, A12, ZZ, NR)
!**************************************************************************  
     IF(NR<=0D0) THEN
         IF(PHA==2) V_p = ZZ(1)
         IF(PHA==1) V_p = ZZ(3)
     ELSE
         V_p = ZZ(1)
     END IF  
     Z_p = V_p/(R*T/P)
     AS = 0D0 
     DO 2005 I=1,NC
         DO 2004 J=1,NC
            AS(I) = AS(I) + ZI(J)*AIJ(I,J)
2004     CONTINUE
          !   não precisa obter a fugacidade, apenas o coeficiente de fugacidade de componente com composição igual a zero pode ser calculado.
          !    O efeito de Z(I) = 0 será incluído no cálculo de equilíbrio de fases ( X_i*Phil_i = Y_i*Phiv_i)
          !    pode eliminar a parcela ==  ZI(I)*P* == da linha a seguir
         fug_p(I) = ZI(I)*P*DEXP((BI(I)/BM)*(Z_p-1D0)-DLOG(Z_p-b2)-    &
                    (a2/b2)*(1D0/(sig-eps))*(2D0*AS(I)/AM - BI(I)/BM)* & 
                    DLOG((Z_p+sig*b2)/(Z_p+eps*b2)))  
2005 CONTINUE
!************************************************************************** 
     CS = 0D0
     DO 2006 I=1,NC
          ! Passagem de parametro CI(I), nem todo compilador entende bem
          !!!!CALL VTPELL( I, T, CI(I) ) 
          CALL VTPELL( I, T, CCI )
          CI(I) = CCI !!!!
          CS = CS + ZI(I)*CI(I)
2006 CONTINUE
     V_p = V_p - CS
     !!!!Z_p = V_p/(R*T/P)

END
    
!C********************************************************************
!C
!C                   PROGRAMA FUGCHK
!C
!C   CALCULO DO COEFICIENTE DE FUGACIDADE
!C
!C    VARIAVEIS:
!C       IFASE    = ESCOLHA DO TIPO DE FASE (1=VAP, 2=LIQ)
!C       NC       = NUM. DE ESPEC.
!C       T        = TEMP.
!C       P        = PRES.
!C       COMP(NC) = COMPOSICAO
!C       FUG(NC)  = COEF. DE FUGACIDADE
!C
!C********************************************************************
!C
      SUBROUTINE FUGCHK(IFASE, NC, T, P, COMP, B, A, V, FUG)
!
      IMPLICIT REAL*8(A-H,O-Z)
!C
      DIMENSION COMP(NC) ,FUG(NC), B(NC), A(NC), ALFF(NC), GAM(NC)
!C
      COMMON /RGAS  / RGAS
      COMMON /ITDM2 / ITDM2
!C
      RGAS = 8.314D01
!C
      ITDM1  = 0
      M01    =    1 + ITDM1
      M02    =  M01 + NC
      M03    =  M02 + NC
      M04    =  M03 + NC
      M05    =  M04 + NC
      M06    =  M05 + NC
      M07    =  M06 + NC
      M08    =  M07 + NC
      M09    =  M08 + NC
      M10    =  M09 + NC
      M11    =  M10 + NC
      M12    =  M11 + NC
      M13    =  M12 + NC
      M14    =  M13 + NC
      M15    =  M14 + NC
      M16    =  M15 + NC
      ITDM2  =  M16
!C
!C IFASE = 2   VAP
!C
      IF (IFASE.EQ.2) ITYP=-1
      IF (IFASE.EQ.1) ITYP=1
!C
      CALL MIX(NC, P, T, COMP, B, A, ALFF, GAM, BM, AM)
      CALL FATORZ(ITYP, P, T, BM, AM, AA, BB, Z)
      CALL FUGCOF(NC, T, P, B, BM, ALFF, GAM, Z, FUG)
!C
      V = Z*RGAS*T/P
!C
      RETURN
      END
!C
!C--------------------------------------------------------------------
!C
!C     SUBROTINA: MIX
!C
!C     FINALIDADE: Calcular os coeficientes da equação cúbica em Z
!C                 e partes que vão interessar ao cálculo do
!C                 coeficiente de fugacidade.
!C
!C     REFERÊNCIAS: Boukouvalas et al., FPE, 92, (1994), 75-106
!C
!C     AUTORA: G.M.N. Costa
!C
!C     MODIFICAÇÕES: S.A.B. Vieira de Melo (MAI.97)
!C                   F.L.P.PESSOA (NOV.97)
!C
!C----------------------------------------------------------------------
!C
      SUBROUTINE MIX(NC, P, TT, ZZ, B, A, ALFF, GAM, BM, AM)
!C
      IMPLICIT REAL*8(A-H,O-Z)
!C
      COMMON /RGAS  / RGAS
      COMMON /CC    / CLAMB,AV,AMM,ALF
!C
      DIMENSION B(NC),A(NC),ALFF(NC),ZZ(NC),GAM(NC)
!C
!C     Cálculo de B da mistura
!C
      BM = 0.D0
      DO 6000 I=1,NC
         BM = BM + ZZ(I)*B(I)
6000  CONTINUE
!C
      ALF1=0.D0
      BBLOG=0.D0
      DO 100 I=1,NC
!C        CAlculo do somatOrio de X(I) e ALFA(I)
         ALFF(I)=A(I)/B(I)/RGAS/TT
         ALF1=ALF1+ZZ(I)*ALFF(I)
!C        CALculo do somatOrio de X(I)LNB(I)
         BBLOG=BBLOG+ZZ(I)*DLOG(BM/B(I))
  100 CONTINUE
!C
!C     Valores dos parAmetros LAMBDA,AV,AM da EOS e da regra de mistura
!C
      AV = -0.623
      AMM = -0.52
      CLAMB = 0.36
!C
!C     Cálculo do valor de ALFA dado pela equaCAO (5) do artigo do
!C     BOUKOUVALAS ET AL.(1994).
!C
      ALF2=(1.D0-CLAMB)/AMM
      ALF3=(CLAMB/AV)+ALF2
      ALF4=BBLOG*ALF2
!C
!C     Chamada do UNIQUAC
!C
      CALL UNIQUAC(NC, TT, ZZ, GAM)
      GERT=0.D0
      DO 135 I=1,NC 
        GERT=GERT+ZZ(I)*DLOG(GAM(I)) 
135   CONTINUE
!C
!C     GE=GE/RT e GAM é o LN GAM(I)
!C
      ALF=ALF3*GERT+ALF1
      ALF=ALF+ALF4
!C
!C     Cálculo de A da mistura
!C
      AM=ALF*BM*RGAS*TT
!C
      RETURN
      END
!C
!C*******************************************************************
!C
!C                   PROGRAMA FATORZ
!C
!C   CALCULO DO FATOR DE COMP.
!C
!C    VARIAVEIS:
!C       ITYP  = TIPO DE FASE (-1=VAP  1=LIQ)
!C       P     = PRESSAO      .
!C       T     = TEMPERATURA
!C       BM    = PAR. B DE MIST.
!C       AM    = PARAM. A DE MIST.
!C       Z     = FATOR DE COMP.
!C
!C*******************************************************************
!C
      SUBROUTINE FATORZ(ITYP, P, T, BM, AM, AA, BB, Z)
!C
      IMPLICIT REAL*8(A-H,O-Z)
!C
      COMMON /RGAS  / RGAS
!C
      DIMENSION :: ZZ(3)
!C
      INTEGER :: PAR
!C
      r1 = - (1D0 + DSQRT(2D0)) 
      r2 = - (1D0 - DSQRT(2D0))
!C
      A12 = - BM*(r1*r2*(BM**2) + r1*r2*BM*RGAS*T/P + AM/P)
      A11 = (BM**2)*(r1*r2+r1+r2) + RGAS*T*BM*(r1+r2)/P + AM/P
      A10 = - (BM*(r1+r2+1D0) + RGAS*T/P)  
      CALL CUBICA(A10, A11, A12, ZZ, PAR)
      IF(PAR<=0D0) THEN
          IF(ITYP==-1) V = ZZ(1)
          IF(ITYP==1) V = ZZ(3)
      ELSE
          V = ZZ(1)
      END IF
!C
      Z = V/(RGAS*T/P)
!C
      RETURN
      END
!C
!C *******************************************************************
!C *                                                                 *
!C *   TITULO:                                                       *
!C *     SUBROTINA FUGCOF                                            *
!C *                                                                 *
!C *******************************************************************
!C
      SUBROUTINE FUGCOF(NC,TT,P,B,BM,ALFF,GAM,Z,FUG)
!C
      IMPLICIT REAL*8(A-H,O-Z)
!C
      COMMON /RGAS  / RGAS
      COMMON /CC    / CLAMB,AV,AMM,ALF
!C
      DIMENSION B(NC),ALFF(NC),GAM(NC),FUG(NC)
!C
!C    BOUKOUVALAS ET AL. (1994).
!C
      D1=1.D0-DSQRT(2.D0)
      D2=1.D0+DSQRT(2.D0)
      PRT=P/RGAS/TT
      BMIX=BM*PRT
      T1=Z+BMIX*D2
      T2=Z+BMIX*D1
      TA=-DLOG(T1/T2)
      D3=2.D0*DSQRT(2.D0)
      TTA=TA/D3
!C
      TTB=-DLOG(Z-BMIX)
!C
      TTC=(Z-1.D0)/BM
!C
      ALF2=(1.D0-CLAMB)/AMM
      ALF3=(CLAMB/AV)+ALF2
      DO 200 I=1,NC
         ALF5=B(I)/BM
         ALF6=-DLOG(ALF5)+ALF5-1.D0
         ALF7=ALF2*ALF6+ALFF(I)
         ADD=ALF3*DLOG(GAM(I))+ALF7
         FUG(I)=TTC*B(I)+TTB+ADD*TTA
  200 CONTINUE
      DO 10 I=1,NC
        FUG(I)=DEXP(FUG(I))
  10  CONTINUE
      
!C
      RETURN
      END