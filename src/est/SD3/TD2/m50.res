Sun Oct 24 00:45:34 CDT 2021
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
$DATA ../../../../data/SD3/TD2/dat50.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2140.05448925936        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4117E+02  1.1794E+01 -2.0314E+01  2.5164E+01 -2.9885E+01  4.3088E+01  1.5601E+01  2.1900E+01  7.4847E+00  2.9191E+01
            -1.4260E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2148.35452529361        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.9342E-01  1.0239E+00  1.2250E+00  1.0319E+00  1.1113E+00  9.9235E-01  9.1885E-01  8.6707E-01  9.9653E-01  8.6463E-01
             1.0152E+00
 PARAMETER:  9.3395E-02  1.2358E-01  3.0295E-01  1.3143E-01  2.0555E-01  9.2320E-02  1.5373E-02 -4.2630E-02  9.6527E-02 -4.5457E-02
             1.1512E-01
 GRADIENT:   3.8887E+00  2.5279E+01  1.7725E+01  8.5234E-01 -1.4529E+01 -9.3109E-01  7.6342E+00  5.9331E+00 -6.5229E+00 -1.6257E+01
            -9.6548E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2153.03200776078        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      347
 NPARAMETR:  9.9391E-01  9.3184E-01  1.2229E+00  1.0836E+00  1.0911E+00  1.0028E+00  6.7607E-01  5.6238E-01  9.7137E-01  9.8955E-01
             9.9826E-01
 PARAMETER:  9.3896E-02  2.9409E-02  3.0123E-01  1.8024E-01  1.8717E-01  1.0283E-01 -2.9146E-01 -4.7558E-01  7.0953E-02  8.9497E-02
             9.8257E-02
 GRADIENT:   7.2757E+00  1.9461E+01  6.1757E+00  4.8682E+00 -8.8859E+00  3.8591E+00 -1.0825E+00  1.6142E+00 -1.9210E+01 -2.0335E+00
            -1.2171E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2154.82567590992        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.8861E-01  7.7473E-01  1.2618E+00  1.1804E+00  1.0482E+00  9.9323E-01  6.3233E-01  4.0938E-01  9.8250E-01  1.0138E+00
             1.0153E+00
 PARAMETER:  8.8548E-02 -1.5524E-01  3.3251E-01  2.6589E-01  1.4712E-01  9.3210E-02 -3.5835E-01 -7.9312E-01  8.2340E-02  1.1367E-01
             1.1521E-01
 GRADIENT:  -8.1284E-01  5.9378E+00  3.0330E+00  6.7610E+00 -6.6640E+00  9.5186E-01  8.4783E-01  6.5662E-01  3.7227E+00  1.7823E+00
             2.1517E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2155.27256978209        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  9.8687E-01  5.9917E-01  1.2943E+00  1.2864E+00  1.0019E+00  9.8828E-01  5.2098E-01  2.2231E-01  9.0725E-01  1.0178E+00
             1.0166E+00
 PARAMETER:  8.6784E-02 -4.1221E-01  3.5795E-01  3.5188E-01  1.0189E-01  8.8215E-02 -5.5205E-01 -1.4037E+00  2.6656E-03  1.1766E-01
             1.1651E-01
 GRADIENT:   3.5796E-01  1.3127E+00  5.5165E-01  2.6205E+00 -1.3974E+00  1.1632E-02  8.8841E-02  1.0527E-01 -2.0889E-01  1.8659E-01
            -9.5098E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2155.32253423868        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  9.8640E-01  5.9427E-01  1.2937E+00  1.2879E+00  1.0019E+00  9.8815E-01  4.7475E-01  8.2901E-02  9.0968E-01  1.0265E+00
             1.0180E+00
 PARAMETER:  8.6302E-02 -4.2042E-01  3.5752E-01  3.5299E-01  1.0192E-01  8.8083E-02 -6.4497E-01 -2.3901E+00  5.3415E-03  1.2620E-01
             1.1788E-01
 GRADIENT:  -5.2736E-01  5.9882E-02 -6.9497E-01  4.6646E-01  1.0897E-01  2.6660E-02  3.3186E-02  1.7981E-02  4.9152E-01  4.9603E-01
             2.1269E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2155.33228878187        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1059             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8735E-01  5.9245E-01  1.2959E+00  1.2874E+00  1.0022E+00  9.8833E-01  4.6665E-01  1.7389E-02  9.0854E-01  1.0259E+00
             1.0182E+00
 PARAMETER:  8.7271E-02 -4.2349E-01  3.5924E-01  3.5259E-01  1.0223E-01  8.8264E-02 -6.6218E-01 -3.9519E+00  4.0869E-03  1.2559E-01
             1.1801E-01
 GRADIENT:   4.0953E+02  5.7811E+01  5.7514E+00  4.0625E+02  6.1852E+00  4.0037E+01  1.7583E+00  2.5485E-03  8.7091E+00  8.4563E-01
             1.1568E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2155.33506135465        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1245             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8703E-01  5.9515E-01  1.2945E+00  1.2856E+00  1.0031E+00  9.8830E-01  4.7035E-01  1.0000E-02  9.0936E-01  1.0253E+00
             1.0181E+00
 PARAMETER:  8.6947E-02 -4.1894E-01  3.5815E-01  3.5121E-01  1.0314E-01  8.8235E-02 -6.5427E-01 -5.9042E+00  4.9849E-03  1.2501E-01
             1.1797E-01
 GRADIENT:   4.0875E+02  5.7486E+01  5.0729E+00  4.0339E+02  7.2449E+00  4.0059E+01  1.7797E+00  0.0000E+00  8.6434E+00  7.2528E-01
             1.1560E+00

0ITERATION NO.:   38    OBJECTIVE VALUE:  -2155.33565362871        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1339
 NPARAMETR:  9.8690E-01  5.9577E-01  1.2949E+00  1.2861E+00  1.0029E+00  9.8823E-01  4.6729E-01  1.0000E-02  9.0926E-01  1.0256E+00
             1.0181E+00
 PARAMETER:  8.6814E-02 -4.1790E-01  3.5840E-01  3.5159E-01  1.0293E-01  8.8160E-02 -6.6082E-01 -5.9042E+00  4.8797E-03  1.2524E-01
             1.1798E-01
 GRADIENT:   5.8988E-01 -2.2927E-01 -9.3350E-02 -1.2883E+00  8.1952E-02  5.6112E-02  5.9321E-03  0.0000E+00 -6.8916E-02 -8.9735E-02
            -6.6183E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1339
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.9736E-04 -1.0979E-02 -2.5229E-04 -1.7751E-03 -2.3405E-02
 SE:             2.9873E-02  5.4110E-03  1.7392E-04  2.9117E-02  2.4996E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7337E-01  4.2451E-02  1.4689E-01  9.5139E-01  3.4910E-01

 ETASHRINKSD(%)  1.0000E-10  8.1873E+01  9.9417E+01  2.4558E+00  1.6261E+01
 ETASHRINKVR(%)  1.0000E-10  9.6714E+01  9.9997E+01  4.8513E+00  2.9877E+01
 EBVSHRINKSD(%)  3.4796E-01  8.2905E+01  9.9419E+01  2.5388E+00  1.3717E+01
 EBVSHRINKVR(%)  6.9471E-01  9.7077E+01  9.9997E+01  5.0131E+00  2.5552E+01
 RELATIVEINF(%)  9.8459E+01  1.9253E-01  7.6614E-04  8.1402E+00  1.0860E+01
 EPSSHRINKSD(%)  3.1441E+01
 EPSSHRINKVR(%)  5.2996E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2155.3356536287120     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1236.3971204240393     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.35
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2155.336       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  5.96E-01  1.29E+00  1.29E+00  1.00E+00  9.88E-01  4.67E-01  1.00E-02  9.09E-01  1.03E+00  1.02E+00
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       40.906
Stop Time:
Sun Oct 24 00:45:43 CDT 2021
