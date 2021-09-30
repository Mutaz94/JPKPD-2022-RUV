Wed Sep 29 12:02:29 CDT 2021
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
$DATA ../../../../data/spa/A1/dat29.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1136.81121757569        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6157E+02  9.8692E+00 -1.2278E+01  3.2469E+01  1.1104E+02  2.2813E+01 -3.2466E+01  4.3868E+00 -3.9936E+01 -5.0961E+01
            -8.8715E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1387.37271214455        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.9577E-01  9.6984E-01  1.1604E+00  1.0739E+00  1.0195E+00  1.0728E+00  1.1525E+00  8.4961E-01  1.3183E+00  1.0299E+00
             1.6289E+00
 PARAMETER:  9.5761E-02  6.9374E-02  2.4877E-01  1.7130E-01  1.1936E-01  1.7026E-01  2.4190E-01 -6.2972E-02  3.7631E-01  1.2942E-01
             5.8793E-01
 GRADIENT:   1.7860E+02  1.4754E+01 -9.5111E+00  6.6734E+01  4.2243E+01  4.0749E+01  6.8087E+00  4.4616E+00  4.4918E+01 -2.3360E+01
            -2.0381E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1408.83041125420        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.6955E-01  6.1025E-01  1.2603E+00  1.2938E+00  8.9303E-01  1.1020E+00  8.7432E-01  3.7598E-02  9.8975E-01  1.1024E+00
             1.7661E+00
 PARAMETER:  6.9080E-02 -3.9389E-01  3.3136E-01  3.5760E-01 -1.3140E-02  1.9716E-01 -3.4314E-02 -3.1808E+00  8.9695E-02  1.9746E-01
             6.6875E-01
 GRADIENT:   1.0577E+02  3.9833E+01  1.1546E+01  1.4013E+02 -8.7735E+00  5.2985E+01 -4.8044E+00  1.1498E-02 -1.8627E+01 -5.5991E+00
            -1.3566E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1420.33975169799        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.4840E-01  7.6031E-01  1.1100E+00  1.1377E+00  8.8262E-01  9.9472E-01  1.0757E+00  2.9596E-02  1.0722E+00  1.0107E+00
             2.1912E+00
 PARAMETER:  4.7025E-02 -1.7403E-01  2.0436E-01  2.2900E-01 -2.4863E-02  9.4708E-02  1.7294E-01 -3.4201E+00  1.6968E-01  1.1060E-01
             8.8446E-01
 GRADIENT:  -5.8352E+00 -9.9794E+00 -2.2241E+00 -1.9128E+01  9.9006E-01 -2.8117E+00 -8.6430E-01  9.5207E-03  1.0110E+00  3.3901E+00
             7.9033E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1423.80924856145        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      363
 NPARAMETR:  9.9376E-01  8.1760E-01  1.2153E+00  1.1441E+00  9.5427E-01  1.0121E+00  1.0812E+00  4.3503E-02  1.0903E+00  1.0489E+00
             2.2130E+00
 PARAMETER:  9.3737E-02 -1.0139E-01  2.9501E-01  2.3462E-01  5.3189E-02  1.1202E-01  1.7812E-01 -3.0349E+00  1.8647E-01  1.4770E-01
             8.9433E-01
 GRADIENT:   3.9104E+00 -5.5665E+00 -3.5999E+00 -2.3473E+01  4.1722E+00 -2.7688E+00 -4.1210E-01  1.4237E-02  6.6799E-01  1.2383E+00
             1.4666E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1427.09068890703        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      540
 NPARAMETR:  9.7593E-01  4.3667E-01  1.1701E+00  1.4048E+00  8.0447E-01  1.0223E+00  1.6597E+00  3.4797E-02  9.3626E-01  9.9242E-01
             2.1873E+00
 PARAMETER:  7.5635E-02 -7.2857E-01  2.5706E-01  4.3986E-01 -1.1758E-01  1.2210E-01  6.0662E-01 -3.2582E+00  3.4134E-02  9.2392E-02
             8.8265E-01
 GRADIENT:  -2.7286E+01  8.7945E+00  1.6269E+00  1.8463E+01 -6.1895E+00  6.3913E-01 -1.1438E-01  1.2478E-02 -8.6553E-01  7.8441E-01
             2.3389E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1429.92739420319        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      716
 NPARAMETR:  9.8529E-01  1.3220E-01  1.0808E+00  1.5765E+00  6.9047E-01  1.0166E+00  3.4182E+00  1.8664E-02  8.6924E-01  9.6678E-01
             2.1782E+00
 PARAMETER:  8.5183E-02 -1.9234E+00  1.7771E-01  5.5521E-01 -2.7038E-01  1.1650E-01  1.3291E+00 -3.8812E+00 -4.0132E-02  6.6216E-02
             8.7850E-01
 GRADIENT:   1.9807E+00  3.4983E+00  3.7108E+00  2.8520E+01 -1.0458E+01 -2.4081E-01  6.9950E-01  4.3552E-03  7.1306E-01  1.9960E+00
             6.9108E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1431.13334556237        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      891
 NPARAMETR:  9.8223E-01  2.9134E-02  1.0053E+00  1.6033E+00  6.4176E-01  1.0162E+00  5.8823E+00  1.0000E-02  8.4756E-01  9.2652E-01
             2.1587E+00
 PARAMETER:  8.2069E-02 -3.4358E+00  1.0531E-01  5.7207E-01 -3.4354E-01  1.1608E-01  1.8719E+00 -5.4279E+00 -6.5392E-02  2.3676E-02
             8.6951E-01
 GRADIENT:   1.5042E+00  2.3751E-01 -2.2478E+00 -2.1931E+00  2.4968E+00  2.2824E-03 -1.3577E-01  0.0000E+00 -1.5041E+00 -1.3232E+00
             4.0463E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1431.31876230363        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  9.8161E-01  1.0000E-02  1.0594E+00  1.6192E+00  6.6392E-01  1.0150E+00  8.8034E+00  1.0000E-02  8.4688E-01  9.5576E-01
             2.1387E+00
 PARAMETER:  8.1436E-02 -4.5295E+00  1.5774E-01  5.8192E-01 -3.0960E-01  1.1485E-01  2.2751E+00 -6.2847E+00 -6.6196E-02  5.4749E-02
             8.6019E-01
 GRADIENT:   1.9517E+00  4.2070E-02 -3.9226E+00 -7.0512E+00  7.5608E+00 -3.9093E-01 -3.5033E-02  0.0000E+00  5.1859E-02 -8.3968E-01
            -3.2260E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1431.34858808521        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1266             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8086E-01  1.0000E-02  1.0616E+00  1.6216E+00  6.6159E-01  1.0160E+00  1.1259E+01  1.0000E-02  8.4655E-01  9.5968E-01
             2.1392E+00
 PARAMETER:  8.0673E-02 -4.5583E+00  1.5974E-01  5.8342E-01 -3.1311E-01  1.1589E-01  2.5212E+00 -6.2847E+00 -6.6585E-02  5.8843E-02
             8.6042E-01
 GRADIENT:   8.8161E+01  0.0000E+00  1.2262E+00  2.3230E+02  9.0020E+00  1.0093E+01  1.5635E-01  0.0000E+00  3.8765E+00  4.9911E-01
             6.9442E+00

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1431.34895110835        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1337
 NPARAMETR:  9.8084E-01  1.0000E-02  1.0616E+00  1.6216E+00  6.6159E-01  1.0160E+00  1.1769E+01  1.0000E-02  8.4655E-01  9.5968E-01
             2.1392E+00
 PARAMETER:  8.0658E-02 -4.5583E+00  1.5974E-01  5.8339E-01 -3.1311E-01  1.1588E-01  2.5654E+00 -6.2847E+00 -6.6589E-02  5.8841E-02
             8.6042E-01
 GRADIENT:   1.1392E-01  0.0000E+00  4.7006E-02 -2.7912E+00  7.1222E-01 -6.1824E-03  1.7547E-02  0.0000E+00  1.2789E-01  6.3535E-02
             6.8405E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1337
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2117E-04 -1.1739E-04  1.2658E-05 -1.1445E-02 -2.4657E-02
 SE:             2.9340E-02  1.7895E-03  1.3080E-04  2.7816E-02  2.1821E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9670E-01  9.4769E-01  9.2291E-01  6.8075E-01  2.5848E-01

 ETASHRINKSD(%)  1.7065E+00  9.4005E+01  9.9562E+01  6.8128E+00  2.6897E+01
 ETASHRINKVR(%)  3.3839E+00  9.9641E+01  9.9998E+01  1.3162E+01  4.6560E+01
 EBVSHRINKSD(%)  1.6972E+00  9.4398E+01  9.9490E+01  6.5092E+00  2.6075E+01
 EBVSHRINKVR(%)  3.3655E+00  9.9686E+01  9.9997E+01  1.2595E+01  4.5350E+01
 RELATIVEINF(%)  9.0427E+01  9.7863E-03  1.7150E-04  3.8182E+00  2.7652E+00
 EPSSHRINKSD(%)  3.4869E+01
 EPSSHRINKVR(%)  5.7580E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1431.3489511083451     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -696.19812454460691     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1431.349       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.00E-02  1.06E+00  1.62E+00  6.62E-01  1.02E+00  1.18E+01  1.00E-02  8.47E-01  9.60E-01  2.14E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.08E+03
 
 TH 2
+        0.00E+00  2.81E+03
 
 TH 3
+        1.14E+00  0.00E+00  2.60E+02
 
 TH 4
+       -2.74E+01  0.00E+00 -3.68E+01  5.16E+02
 
 TH 5
+        1.24E+01  0.00E+00 -5.87E+02 -8.43E+01  1.48E+03
 
 TH 6
+        1.25E+00  0.00E+00  2.79E+00 -7.32E+00 -3.10E+00  1.77E+02
 
 TH 7
+        3.77E-03  0.00E+00 -1.66E-02 -5.62E-02  4.98E-02 -1.64E-03  1.53E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.99E+00  0.00E+00  1.06E+01 -7.92E+00 -4.32E+00 -4.28E-01 -1.39E-02  0.00E+00  2.15E+02
 
 TH10
+       -4.02E+00  0.00E+00  3.32E+00 -2.24E+00 -6.39E+01  7.59E-01 -3.20E-02  0.00E+00  9.26E-01  7.25E+01
 
 TH11
+       -1.21E+01  0.00E+00 -7.81E+00 -8.91E+00 -4.13E-01  3.60E+00 -4.59E-03  0.00E+00  9.18E+00  1.83E+01  5.99E+01
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.02
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       22.408
Stop Time:
Wed Sep 29 12:02:53 CDT 2021
