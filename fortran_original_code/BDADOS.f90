!****************************************************************************
!
!  PROGRAMA: BDADOS
!
!  OBJETIVO:  Leitura dos dados do arquivo DADOS.DAT
!
!****************************************************************************

SUBROUTINE BDADOS(ID, NC, Z, T, P_orv, Nexp, Pexp, Vexp)

    IMPLICIT REAL*8(A-H,O-Z)
    
    COMMON /DADOS/ TC(20), PC(20), W(20)
    COMMON /ABM/   ALFAF, AAT1(20), AAT2(20), AAT3(20)
    COMMON /VTR/   APE1(20), APE2(20), APE3(20)
    COMMON /UNIQ/  UIJ_0(20,20), UIJ_T(20,20), R(20), Q(20)
    COMMON /BIP/   DIJ(20,20)
    
    DIMENSION Z(20), Pexp(25), Vexp(25), PEXPLIT(25), VEXPLIT(25), PLIT(25), VLIT(25)
    
    OPEN(01,file='DADOS.DAT')
    
!***********************************************************************
!                               SISTEMAS
!***********************************************************************
    
    IF(ID == 1) THEN
        NC = 12
        READ(01,101) Z(1)
        DO I=2,NC
            READ(01,102) Z(I)
        END DO
        READ(01,102) CSG, CMW
        READ(01,103) T
        READ(01,104) P_orv
    ELSE IF(ID == 2) THEN
        NC = 12
        READ(01,201) Z(1)
        DO I=2,NC
            READ(01,202) Z(I)      
        END DO
        READ(01,202) CSG, CMW
        READ(01,203) T
        READ(01,204) P_orv
    ELSE IF(ID == 3) THEN
        NC = 11 !NC=12
        READ(01,305) Z(1) !READ(01,301) Z(1)
        DO I=2,NC
            READ(01,302) Z(I)      
        END DO
        READ(01,302) CSG, CMW
        READ(01,303) T
        READ(01,304) P_orv
    ELSE IF(ID == 4) THEN
        NC = 11 !NC=12
        READ(01,405) Z(1) !READ(01,401) Z(1)
        DO I=2,NC
            READ(01,402) Z(I)      
        END DO
        READ(01,402) CSG, CMW
        READ(01,403) T
        READ(01,404) P_orv
    ELSE IF(ID == 5) THEN
        NC=11 !NC = 12
        READ(01,505) Z(1) !READ(01,501) Z(1)
        DO I=2,NC
            READ(01,502) Z(I)      
        END DO
        READ(01,502) CSG, CMW
        READ(01,503) T
        READ(01,504) P_orv
    ELSE IF(ID == 6) THEN
        NC = 12
        READ(01,601) Z(1)
        DO I=2,NC
            READ(01,602) Z(I)      
        END DO
        READ(01,602) CSG, CMW
        READ(01,603) T
        READ(01,604) P_orv
    ELSE IF(ID == 7) THEN
        NC = 11    !NC=12
        READ (01,705) Z(1) !READ(01,701) Z(1)
        DO I=2,NC
            READ(01,702) Z(I)      
        END DO
        READ(01,702) CSG, CMW
        READ(01,703) T
        READ(01,704) P_orv
    ELSE IF(ID == 8) THEN
        NC = 11 !NC=12
        READ(01,805) Z(1) !READ(01,801) Z(1)
        DO I=2,NC
            READ(01,802) Z(I)      
        END DO
        READ(01,802) CSG, CMW
        READ(01,803) T
        READ(01,804) P_orv
    ELSE IF(ID == 9) THEN
        NC = 11 !NC=12
        READ (01,905) Z(1) !READ(01,901) Z(1)
        DO I=2,NC
            READ(01,902) Z(I)      
        END DO
        READ(01,902) CSG, CMW
        READ(01,903) T
        READ(01,904) P_orv
    ELSE IF(ID == 10) THEN
        NC = 11 !NC=12
        READ(01,1005) Z(1) !READ(01,1005) Z(1)
        DO I=2,NC
            READ(01,1002) Z(I)      
        END DO
        READ(01,1002) CSG, CMW
        READ(01,1003) T
        READ(01,1004) P_orv
    ELSE IF(ID == 11) THEN
        NC = 11
        READ(01,1105) Z(1)
        DO I=2,NC
            READ(01,1102) Z(I)      
        END DO
        READ(01,1102) CSG, CMW
        READ(01,1103) T
        READ(01,1104) P_orv
    ELSE IF(ID == 12) THEN
        NC = 12
        READ(01,1201) Z(1)
        DO I=2,NC
            READ(01,1202) Z(I)      
        END DO
        READ(01,1202) CSG, CMW
        READ(01,1203) T
        READ(01,1204) P_orv
    ELSE IF(ID == 13) THEN
        NC = 12
        READ(01,1301) Z(1)
        DO I=2,NC
            READ(01,1302) Z(I)      
        END DO
        READ(01,1302) CSG, CMW
        READ(01,1303) T
        READ(01,1304) P_orv
    ELSE IF(ID == 14) THEN
        NC = 12
        READ(01,1401) Z(1)
        DO I=2,NC
            READ(01,1402) Z(I)      
        END DO
        READ(01,1402) CSG, CMW
        READ(01,1403) T
        READ(01,1404) P_orv
    ELSE IF(ID == 15) THEN
        NC = 12
        READ(01,1501) Z(1)
        DO I=2,NC
            READ(01,1502) Z(I)      
        END DO
        READ(01,1502) CSG, CMW
        READ(01,1503) T
        READ(01,1504) P_orv
    ELSE IF(ID == 16) THEN
        NC = 12
        READ(01,1601) Z(1)
        DO I=2,NC
            READ(01,1602) Z(I)      
        END DO
        READ(01,1602) CSG, CMW
        READ(01,1603) T
        READ(01,1604) P_orv
    ELSE IF(ID == 17) THEN
        NC = 12
        READ(01,1701) Z(1)
        DO I=2,NC
            READ(01,1702) Z(I)      
        END DO
        READ(01,1702) CSG, CMW
        READ(01,1703) T
        READ(01,1704) P_orv
    END IF

!***********************************************************************
!                             PROPRIEDADES
!***********************************************************************
    IF (ID==3 .OR. ID==4 .OR. ID==5 .OR. ID==7 .OR. ID==8 .OR. ID==9 .OR. ID==10 .OR. ID==11) THEN
        READ(01,2003) TC(1), PC(1), W(1), AAT1(1), AAT2(1), AAT3(1), &
                  APE1(1), APE2(1), APE3(1), R(1), Q(1)
        
        !READ(01,2001) TC(1), PC(1), W(1), AAT1(1), AAT2(1), AAT3(1), &
                  !APE1(1), APE2(1), APE3(1), R(1), Q(1)
        DO I=2,NC
        READ(01,2002) TC(I), PC(I), W(I), AAT1(I), AAT2(I), AAT3(I), &
                      APE1(I), APE2(I), APE3(I), R(I), Q(I) 
        END DO
    ELSE  
        READ(01,2001) TC(1), PC(1), W(1), AAT1(1), AAT2(1), AAT3(1), &
                  APE1(1), APE2(1), APE3(1), R(1), Q(1)
        DO I=2,NC
        READ(01,2002) TC(I), PC(I), W(I), AAT1(I), AAT2(I), AAT3(I), &
                      APE1(I), APE2(I), APE3(I), R(I), Q(I)
        END DO
    END IF
    
    !Propriedades dos pseudo-componentes -------------------------------
    !RIAZI-DAUBERT CORRELATION - Est� na referencia Propriedades &
    !para Pseudos-componentes
    !Ck=Tb/1.8 - onde tb � a temperatura de ebuli��o e ck parametro da equa��o
    CK = 250D0
    ERRO = 1D0
    AAA=42.965D0*(CSG**4.98308D0)    
    BBB=7.78712D0*CSG 
    CCCA=0.0002097D0+0.00208476D0*CSG
    DO WHILE(ERRO > 1D-10)
        F = 42.965D0*(CK**1.26007D0)*(CSG**4.98308D0)* &
            DEXP(0.0002097D0*CK-7.78712D0*CSG+0.00208476D0*CK*CSG) - CMW
        !CK1 = CK + 0.0001D0*CK
        !F1 = 42.965D0*(CK1**1.26007D0)*(CSG**4.98308D0)* &
            !DEXP(0.0002097D0*CK1-7.78712D0*CSG+0.00208476D0*CK1*CSG) - CMW
        F1=AAA*1.26007D0*(CK**0.26007D0)*DEXP((CK*CCCA)-BBB)+AAA*(CK**1.26007D0)*(CCCA)*DEXP((CK*CCCA)-BBB)
        !CK = CK - (0.0001D0*CK)*F/(F1-F)
        CK=CK-F/F1
        !ERRO = DABS((CK-CKN)/CK)
       ERRO = DABS((CK-CKN)/CK)
        !CKN = CK
       CKN=CK
    END DO

    TC(NC) = 17.14194D0*(CK**0.81067D0)*(CSG**0.53691D0)* &
             DEXP(-0.0009314D0*CK-0.544442D0*CSG+0.00064791D0*CK*CSG) 
    
    PC(NC) = 46.35124404D5*(CK**(-0.4844D0))*(CSG**4.0846D0)* &
             DEXP(-0.008505D0*CK-4.8014D0*CSG+0.005749D0*CK*CSG)
    
    Tbr = (1.8D0*CK)/TC(NC) !Temperatura de ebuli��o reduzida
    CK2 = ((1.8D0*CK)**(1/3))/CSG !Watson Factor
    !Kesler-Lee Correlation - C�lculo do fator ac�ntrico (W)
    IF(Tbr > 0.8) THEN
        W(NC) = - 7.904D0 + 0.1352D0*CK2 - 0.007465D0*(CK2**2) + &
                8.359D0*Tbr + (1.408D0 - 0.01063D0*CK2)/Tbr
    ELSE
        W(NC) = (- DLOG(PC(NC)/14.696D0) - 5.92714D0 + 6.09648D0/Tbr +  &
                 1.28862D0*DLOG(Tbr) - 0.169347D0*(Tbr**6))/(15.2518 - &
                 15.6875D0/Tbr - 13.4721D0*DLOG(Tbr) +                  &
                 0.43577D0*(Tbr**6))  ! O coeficiente do termo elevado a sexta era -0.1693478D0.
                                      ! Na refer�ncia das propriedades para pseudos-componentes � - 0.169347D0
    END IF
    RK=5D0/9D0!0.6 bizarramente, da um erro bem minusculo !Convers�o �R -> K COM O VALOR TOTAL DA OPERA��O
    TC(NC) = TC(NC)*RK !0.555556D0 (5D0/9D0) !Convers�o �R -> K
    PC(NC) = PC(NC)*0.0689476D0 !Convers�o psia -> bar
    !-------------------------------------------------------------------
    IF (ID==3 .OR. ID==4 .OR. ID==5 .OR. ID==7 .OR. ID==8 .OR. ID==9 .OR. ID==10 .OR. ID==11) THEN
        READ(01,"(4/)")
        DO I=1,NC
            READ(01,"(11X,11F8.2)") (UIJ_0(I,J),J=1,NC)
        END DO
        READ(01,"(2/)")
        DO I=1,NC
            READ(01,"(11X,11F8.2)") (UIJ_T(I,J),J=1,NC)
        END DO
        READ(01,"(4/)")
        DO I=1,NC
            READ(01,"(11X,11F8.5)") (DIJ(I,J),J=1,NC)
        END DO
     ELSE   
        READ(01,"(3/)")
        DO I=1,NC
            READ(01,"(3X,20F8.2)") (UIJ_0(I,J),J=1,NC)
        END DO
        READ(01,"(/)")
        DO I=1,NC
            READ(01,"(3X,20F8.2)") (UIJ_T(I,J),J=1,NC)
        END DO
        READ(01,"(3/)")
        DO I=1,NC
            READ(01,"(3X,20F8.5)") (DIJ(I,J),J=1,NC)
        END DO
    END IF

    CLOSE(01)

!***********************************************************************
!                          DADOS EXPERIMENTAIS
!***********************************************************************
    
    OPEN(02,file='CURVAS.DAT')
    IF(ID == 1) THEN
        Nexp = 19
        READ(02,'(3/,F6.2,2X,F5.2)') Pexp(1), Vexp(1)
    ELSE IF(ID == 2) THEN
        Nexp = 22
        READ(02,'(24/,F6.2,2X,F5.2)') Pexp(1), Vexp(1)
    ELSE IF(ID == 3) THEN
        Nexp = 22
        READ(02,'(48/,F6.2,2X,F5.2)') Pexp(1), Vexp(1)
    ELSE IF(ID == 4) THEN
        Nexp = 17
        READ(02,'(72/,F6.2,2X,F5.2)') Pexp(1), Vexp(1)
    ELSE IF(ID == 5) THEN
        Nexp = 5
        READ(02,'(91/,F6.2,2X,F5.2)') Pexp(1), Vexp(1)
    ELSE IF(ID == 6) THEN
        Nexp = 24
        READ(02,'(98/,F6.2,2X,F5.2)') Pexp(1), Vexp(1)
    ELSE IF(ID == 7) THEN
        Nexp = 8
        READ(02,'(124/,F6.2,2X,F5.2)') Pexp(1), Vexp(1)
    ELSE IF(ID == 8) THEN
        Nexp = 9
        READ(02,'(134/,F6.2,2X,F5.2)') Pexp(1), Vexp(1)
    ELSE IF(ID == 9) THEN
        Nexp = 15
        READ(02,'(145/,F6.2,2X,F5.2)') Pexp(1), Vexp(1)
    ELSE IF(ID == 10) THEN
        Nexp = 25
        READ(02,'(162/,F6.2,2X,F5.2)') Pexp(1), Vexp(1)
    END IF
    DO I=2,Nexp
        READ(02,'(F6.2,2X,F5.2)') Pexp(I), Vexp(I)
    END DO
        
    CLOSE(02)
    
    
    
!***********************************************************************
!                   DADOS EXPERIMENTAIS DA LITERATURA
!***********************************************************************
    
    OPEN(03,file='LITERATURAEXP.DAT')
    IF(ID == 1) THEN
        Nlit = 14
        READ(03,'(4/,4x,F8.4,5X,F7.4)') Pexplit(1), Vexplit(1)
    ELSE IF(ID == 2) THEN
        Nlit = 13
        READ(03,'(21/,4x,F8.4,5X,F7.4)') Pexplit(1), Vexplit(1)
    ELSE IF(ID == 3) THEN
        Nlit = 11
        READ(03,'(37/,4x,F8.4,5X,F7.4)') Pexplit(1), Vexplit(1)
    ELSE IF(ID == 4) THEN
        Nlit = 11
        READ(03,'(51/,4x,F8.4,5X,F7.4)') Pexplit(1), Vexplit(1)
    ELSE IF(ID == 5) THEN
        Nlit = 5
        READ(03,'(65/,4x,F8.4,5X,F7.4)') Pexplit(1), Vexplit(1)
    ELSE IF(ID == 6) THEN
        Nlit = 17
        READ(03,'(73/,4x,F8.4,5X,F7.4)') Pexplit(1), Vexplit(1)
    ELSE IF(ID == 7) THEN
        Nlit = 6
        READ(03,'(93/,4x,F8.4,5X,F7.4)') Pexplit(1), Vexplit(1)
    ELSE IF(ID == 8) THEN
        Nlit = 6
        READ(03,'(102/,4x,F8.4,5X,F7.4)') Pexplit(1), Vexplit(1)
    ELSE IF(ID == 10) THEN
        Nlit = 19
        READ(03,'(111/,4x,F8.4,5X,F7.4)') Pexplit(1), Vexplit(1)
    END IF
    DO I=2,Nlit
        READ(03,'(4x,F8.4,5X,F7.4)') Pexplit(I), Vexplit(I)
    END DO
        
    CLOSE(03)
    
!***********************************************************************
!                    DADOS CALCULADOS NA LITERATURA
!***********************************************************************
    
    OPEN(04,file='LITERATURACALC.DAT')
    IF(ID == 1) THEN
        Nlit = 14
        READ(04,'(4/,4x,F8.4,5X,F7.4)') Plit(1), Vlit(1)
    ELSE IF(ID == 2) THEN
        Nlit = 13
        READ(04,'(21/,4x,F8.4,5X,F7.4)') Plit(1), Vlit(1)
    ELSE IF(ID == 3) THEN
        Nlit = 11
        READ(04,'(37/,4x,F8.4,5X,F7.4)') Plit(1), Vlit(1)
    ELSE IF(ID == 4) THEN
        Nlit = 11
        READ(04,'(51/,4x,F8.4,5X,F7.4)') Plit(1), Vlit(1)
    ELSE IF(ID == 5) THEN
        Nlit = 5
        READ(04,'(65/,4x,F8.4,5X,F7.4)') Plit(1), Vlit(1)
    ELSE IF(ID == 6) THEN
        Nlit = 17
        READ(04,'(73/,4x,F8.4,5X,F7.4)') Plit(1), Vlit(1)
    ELSE IF(ID == 7) THEN
        Nlit = 6
        READ(04,'(93/,4x,F8.4,5X,F7.4)') Plit(1), Vlit(1)
    ELSE IF(ID == 8) THEN
        Nlit = 6
        READ(04,'(102/,4x,F8.4,5X,F7.4)') Plit(1), Vlit(1)
    ELSE IF(ID == 10) THEN
        Nlit = 19
        READ(04,'(111/,4x,F8.4,5X,F7.4)') Plit(1), Vlit(1)
    END IF
    DO I=2,Nlit
        READ(04,'(4x,F8.4,5X,F7.4)') Plit(I), Vlit(I)
    END DO
        
    CLOSE(04)

!***********************************************************************
!                       ERROS DA LITERATURA
!***********************************************************************  
 ERRORLIT=0
    DO I=1,NLIT
    ERRORLIT=ERRORLIT+100D0*(DABS((VEXPLIT(I)-VLIT(I))/VEXPLIT(I)))
    !print*, errorlit
    END DO
    ERRORMEDIOLIT=ERRORLIT/NLIT
   ! WRITE(*,'(/, "ERROR LITERATURA=", F6.2)') ERRORMEDIOLIT
!***********************************************************************
!                               FORMATOS
!***********************************************************************
    
101  FORMAT(3/, 5X, f6.4)
102  FORMAT(    5X, f6.4)       
103  FORMAT(4/, 4X, f7.3)
104  FORMAT(    4X, f7.3)
201  FORMAT(3/,13X, f6.4)
202  FORMAT(   13X, f6.4)       
203  FORMAT(4/,12X, f7.3)
204  FORMAT(   12X, f7.3)
301  FORMAT(3/,21X, f6.4)
302  FORMAT(   21X, f6.4)       
303  FORMAT(4/,20X, f7.3)
304  FORMAT(   20X, f7.3)
305  FORMAT(4/,21X, F6.4) ! Novo formato para n�o pegar os Z(H2S)=0
401  FORMAT(3/,29X, f6.4)
402  FORMAT(   29X, f6.4)       
403  FORMAT(4/,28X, f7.3)
404  FORMAT(   28X, f7.3)
405  FORMAT(4/,29X, f6.4) ! Novo formato para n�o pegar Z(H2S) =0
501  FORMAT(3/,37X, f6.4)
502  FORMAT(   37X, f6.4)       
503  FORMAT(4/,36X, f7.3)
504  FORMAT(   36X, f7.3)
505  FORMAT(4/,37X, f6.4) ! Novo formato para n�o pegar Z(H2S) =0
601  FORMAT(3/,45X, f6.4)
602  FORMAT(   45X, f6.4)       
603  FORMAT(4/,44X, f7.3)
604  FORMAT(   44X, f7.3)
701  FORMAT(3/,53X, f6.4)
702  FORMAT(   53X, f6.4)       
703  FORMAT(4/,52X, f7.3)
704  FORMAT(   52X, f7.3)
705  FORMAT(4/,53X, f6.4) ! Novo formato para n�o pegar Z(H2S) =0
801  FORMAT(3/,61X, f6.4)
802  FORMAT(   61X, f6.4)       
803  FORMAT(4/,60X, f7.3)
804  FORMAT(   60X, f7.3)
805  FORMAT(4/,61X, f6.4) ! Novo formato para n�o pegar Z(H2S) =0
901  FORMAT(3/,69X, f6.4)
902  FORMAT(   69X, f6.4)       
903  FORMAT(4/,68X, f7.3)
904  FORMAT(   68X, f7.3)
905  FORMAT(4/,69X, f6.4) ! Novo formato para n�o pegar Z(H2S) =0
1001 FORMAT(3/,77X, f6.4)
1002 FORMAT(   77X, f6.4)       
1003 FORMAT(4/,76X, f7.3)
1004 FORMAT(   76X, f7.3)
1005 FORMAT(4/,77X, f6.4) ! Novo formato para n�o pegar Z(H2S) =0
1101 FORMAT(3/,85X, f6.4)
1102 FORMAT(   85X, f6.4)       
1103 FORMAT(4/,84X, f7.3)
1104 FORMAT(   84X, f7.3)
1105 FORMAT(4/,85X, f6.4) ! Novo formato para n�o pegar Z(H2S) =0
1201 FORMAT(3/,93X, f6.4)
1202 FORMAT(   93X, f6.4)       
1203 FORMAT(4/,92X, f7.3)
1204 FORMAT(   92X, f7.3)
1301 FORMAT(3/,101X, f6.4)
1302 FORMAT(   101X, f6.4)       
1303 FORMAT(4/,100X, f7.3)
1304 FORMAT(   100X, f7.3)
1401 FORMAT(3/,109X, f6.4)
1402 FORMAT(   109X, f6.4)       
1403 FORMAT(4/,108X, f7.3)
1404 FORMAT(   108X, f7.3)
1501 FORMAT(3/,117X, f6.4)
1502 FORMAT(   117X, f6.4)       
1503 FORMAT(4/,116X, f7.3)
1504 FORMAT(   116X, f7.3)
1601 FORMAT(3/,125X, f6.4)
1602 FORMAT(   125X, f6.4)       
1603 FORMAT(4/,124X, f7.3)
1604 FORMAT(   124X, f7.3)
1701 FORMAT(3/,133X, f6.4)
1702 FORMAT(   133X, f6.4)       
1703 FORMAT(4/,132X, f7.3)
1704 FORMAT(   132X, f7.3)
2001 FORMAT(4/,5X,f6.2,2X,f6.3,2X,f7.5,2X,f7.5,2X,f7.5,2X,f7.5,2X,f9.5,&
               1X,f10.5,2X,f9.5,2X,f7.5,2X,f7.5)   
2002 FORMAT(   5X,f6.2,2X,f6.3,2X,f7.5,2X,f7.5,2X,f7.5,2X,f7.5,2X,f9.5,&
               1X,f10.5,2X,f9.5,2X,f7.5,2X,f7.5) 
2003 FORMAT(5/,5X,f6.2,2X,f6.3,2X,f7.5,2X,f7.5,2X,f7.5,2X,f7.5,2X,f9.5,&
               1X,f10.5,2X,f9.5,2X,f7.5,2X,f7.5) ! PROPRIEDADES PARA NC =11 e primeiro termo e zero               
                
END