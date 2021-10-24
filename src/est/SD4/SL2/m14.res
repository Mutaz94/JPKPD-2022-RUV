Sun Oct 24 03:11:26 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/SD4/SL2/dat14.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2

$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr
$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m14.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1687.24319511529        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0249E+02 -7.3769E+01 -6.2681E+01  4.2073E+00  9.2612E+01  4.2150E+01  2.9730E+00  1.3772E+01  3.7292E+01  1.0357E+00
            -2.0984E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1699.68910982627        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  1.0158E+00  1.1454E+00  1.1430E+00  9.9068E-01  1.0536E+00  1.0236E+00  9.6173E-01  9.1964E-01  8.1682E-01  9.6181E-01
             1.1199E+00
 PARAMETER:  1.1568E-01  2.3579E-01  2.3366E-01  9.0632E-02  1.5218E-01  1.2334E-01  6.0979E-02  1.6223E-02 -1.0234E-01  6.1060E-02
             2.1326E-01
 GRADIENT:   5.3413E+00  1.0063E+01 -5.8942E+00  2.2455E+01  1.7482E+01  5.0439E+00 -2.8024E+00  3.4371E+00 -5.3438E+00 -1.2705E+01
             1.5588E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1700.86778776085        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  1.0147E+00  1.1212E+00  1.3452E+00  1.0047E+00  1.1478E+00  9.7039E-01  9.3357E-01  6.1415E-01  8.6946E-01  1.2177E+00
             1.0945E+00
 PARAMETER:  1.1459E-01  2.1436E-01  3.9654E-01  1.0464E-01  2.3789E-01  6.9943E-02  3.1261E-02 -3.8751E-01 -3.9882E-02  2.9698E-01
             1.9034E-01
 GRADIENT:   6.6682E+00  1.0436E+00 -1.7688E+00  5.7040E+00  1.8410E+01 -1.6153E+01  2.7300E+00  4.5231E-01  4.0789E-01  2.0658E+00
             8.2818E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1702.57496109491        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  1.0136E+00  1.1308E+00  1.0434E+00  9.8431E-01  1.0112E+00  1.0116E+00  9.5278E-01  3.5040E-01  8.6445E-01  1.0693E+00
             1.0615E+00
 PARAMETER:  1.1354E-01  2.2292E-01  1.4247E-01  8.4184E-02  1.1115E-01  1.1153E-01  5.1632E-02 -9.4868E-01 -4.5660E-02  1.6696E-01
             1.5964E-01
 GRADIENT:   1.4574E+00 -4.1457E-01 -5.3049E-01  2.5866E-01 -1.8724E+00  2.4361E-01 -1.4209E-01  4.7672E-01  6.7147E-01  8.9514E-01
             2.1951E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1702.75559845502        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  1.0140E+00  1.2493E+00  9.6418E-01  9.0676E-01  1.0352E+00  1.0103E+00  8.9056E-01  1.2756E-01  9.1281E-01  1.0683E+00
             1.0635E+00
 PARAMETER:  1.1387E-01  3.2255E-01  6.3524E-02  2.1198E-03  1.3458E-01  1.1027E-01 -1.5909E-02 -1.9592E+00  8.7742E-03  1.6606E-01
             1.6155E-01
 GRADIENT:   5.0247E-01 -1.4821E+00 -1.2006E+00  1.5344E-01  1.9696E+00 -5.7662E-01  4.2923E-02  5.5570E-02  3.3900E-01  7.6093E-02
             2.9109E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1702.75965801992        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      876
 NPARAMETR:  1.0140E+00  1.2878E+00  9.4581E-01  8.8250E-01  1.0441E+00  1.0117E+00  8.6972E-01  8.6474E-02  9.3030E-01  1.0706E+00
             1.0637E+00
 PARAMETER:  1.1393E-01  3.5291E-01  4.4285E-02 -2.4998E-02  1.4311E-01  1.1160E-01 -3.9585E-02 -2.3479E+00  2.7753E-02  1.6825E-01
             1.6176E-01
 GRADIENT:   2.7760E-01  3.7466E-01 -6.4630E-02  4.1291E-01 -8.0194E-02 -1.2444E-01 -4.7568E-02  2.4411E-02  1.7788E-02  2.2078E-02
             2.5716E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1702.77055534550        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1034
 NPARAMETR:  1.0145E+00  1.2886E+00  9.4520E-01  8.8168E-01  1.0444E+00  1.0121E+00  8.6961E-01  2.8971E-02  9.3085E-01  1.0708E+00
             1.0637E+00
 PARAMETER:  1.1435E-01  3.5356E-01  4.3637E-02 -2.5921E-02  1.4342E-01  1.1203E-01 -3.9710E-02 -3.4415E+00  2.8338E-02  1.6845E-01
             1.6174E-01
 GRADIENT:   1.2033E+00  5.9257E-02  6.4866E-02 -3.5430E-02 -3.5239E-02  4.5827E-02 -9.4262E-04  2.7266E-03  1.8359E-02 -3.9999E-02
            -2.3459E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1702.77074965947        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1209
 NPARAMETR:  1.0140E+00  1.2884E+00  9.4506E-01  8.8181E-01  1.0443E+00  1.0120E+00  8.6979E-01  2.0361E-02  9.3061E-01  1.0710E+00
             1.0637E+00
 PARAMETER:  1.1391E-01  3.5343E-01  4.3492E-02 -2.5775E-02  1.4333E-01  1.1191E-01 -3.9507E-02 -3.7941E+00  2.8085E-02  1.6859E-01
             1.6174E-01
 GRADIENT:   2.2617E-01  8.6384E-02  7.7876E-03  4.1048E-02  1.3571E-02 -2.2397E-03 -1.2185E-03  1.3595E-03  7.2061E-03 -1.5803E-03
            -2.0699E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1702.77220728667        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1391
 NPARAMETR:  1.0135E+00  1.2822E+00  9.4498E-01  8.8558E-01  1.0404E+00  1.0132E+00  8.7440E-01  1.0000E-02  9.2634E-01  1.0683E+00
             1.0634E+00
 PARAMETER:  1.1341E-01  3.4861E-01  4.3406E-02 -2.1511E-02  1.3964E-01  1.1311E-01 -3.4221E-02 -1.2035E+01  2.3491E-02  1.6602E-01
             1.6147E-01
 GRADIENT:  -8.4292E-01  2.9912E-01  2.8363E-01 -2.8005E-01 -8.2267E-01  4.7541E-01  2.8806E-02  0.0000E+00 -5.0781E-02  1.8676E-02
            -4.3894E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1702.77385895895        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1569
 NPARAMETR:  1.0152E+00  1.2830E+00  9.4525E-01  8.8480E-01  1.0417E+00  1.0124E+00  8.7319E-01  1.0000E-02  9.2752E-01  1.0691E+00
             1.0635E+00
 PARAMETER:  1.1509E-01  3.4923E-01  4.3693E-02 -2.2390E-02  1.4086E-01  1.1228E-01 -3.5603E-02 -1.0193E+01  2.4756E-02  1.6680E-01
             1.6157E-01
 GRADIENT:   2.8257E+00 -5.2698E-01 -2.2147E-02 -5.7259E-01 -1.7563E-02  1.4872E-01  2.0590E-03  0.0000E+00  1.1919E-02 -8.6781E-04
            -1.4795E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1702.77385895895        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1591
 NPARAMETR:  1.0152E+00  1.2830E+00  9.4525E-01  8.8480E-01  1.0417E+00  1.0124E+00  8.7319E-01  1.0000E-02  9.2752E-01  1.0691E+00
             1.0635E+00
 PARAMETER:  1.1509E-01  3.4923E-01  4.3693E-02 -2.2390E-02  1.4086E-01  1.1228E-01 -3.5603E-02 -1.0193E+01  2.4756E-02  1.6680E-01
             1.6157E-01
 GRADIENT:   2.8257E+00 -5.2698E-01 -2.2147E-02 -5.7259E-01 -1.7563E-02  1.4872E-01  2.0590E-03  0.0000E+00  1.1919E-02 -8.6781E-04
            -1.4795E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1591
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -7.5015E-04 -1.4773E-02 -2.7368E-04  4.5007E-03 -2.5635E-02
 SE:             2.9805E-02  2.0343E-02  1.1786E-04  2.2890E-02  2.3999E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7992E-01  4.6774E-01  2.0223E-02  8.4412E-01  2.8544E-01

 ETASHRINKSD(%)  1.4998E-01  3.1848E+01  9.9605E+01  2.3316E+01  1.9600E+01
 ETASHRINKVR(%)  2.9973E-01  5.3552E+01  9.9998E+01  4.1195E+01  3.5358E+01
 EBVSHRINKSD(%)  4.6217E-01  3.1360E+01  9.9625E+01  2.4280E+01  1.7278E+01
 EBVSHRINKVR(%)  9.2220E-01  5.2885E+01  9.9999E+01  4.2665E+01  3.1571E+01
 RELATIVEINF(%)  9.8561E+01  1.0546E+00  1.1613E-04  1.3828E+00  8.6851E+00
 EPSSHRINKSD(%)  4.2153E+01
 EPSSHRINKVR(%)  6.6537E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1702.7738589589537     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -967.62303239521555     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1702.774       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.28E+00  9.45E-01  8.85E-01  1.04E+00  1.01E+00  8.73E-01  1.00E-02  9.28E-01  1.07E+00  1.06E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       45.029
Stop Time:
Sun Oct 24 03:11:36 CDT 2021
