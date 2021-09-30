Wed Sep 29 19:30:17 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat92.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1651.51232106701        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4736E+02 -4.4820E+01 -6.3782E+01  2.8235E+01  8.5901E+01  5.4745E+01 -1.5390E+01  9.8513E+00 -9.8836E+00  1.4119E+01
             2.3033E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1657.55618198027        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.8488E-01  1.1182E+00  1.2676E+00  9.4839E-01  1.0810E+00  8.8863E-01  1.1579E+00  9.2272E-01  1.1731E+00  8.3730E-01
             9.9759E-01
 PARAMETER:  8.4767E-02  2.1170E-01  3.3709E-01  4.7015E-02  1.7789E-01 -1.8076E-02  2.4657E-01  1.9568E-02  2.5965E-01 -7.7570E-02
             9.7591E-02
 GRADIENT:   9.6940E+00 -1.4847E+01  2.0204E+01 -2.0724E+01  1.6625E+01 -3.0772E+01  2.3705E+00 -9.6939E+00  1.1043E+01 -2.2033E+01
            -1.1354E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1658.89617155201        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.8833E-01  1.3945E+00  1.1171E+00  8.0878E-01  1.1217E+00  9.1927E-01  1.0679E+00  9.4744E-01  1.3467E+00  8.7128E-01
             1.0298E+00
 PARAMETER:  8.8259E-02  4.3255E-01  2.1072E-01 -1.1223E-01  2.1487E-01  1.5827E-02  1.6570E-01  4.6004E-02  3.9767E-01 -3.7789E-02
             1.2939E-01
 GRADIENT:   1.3766E+01  3.3005E+01  2.3239E+01  1.1835E+01 -1.5725E+01 -1.6524E+01  1.0494E+01 -5.8268E+00  1.2011E+01 -1.4562E+01
             3.3297E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1662.05779289094        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.8320E-01  1.3306E+00  1.1422E+00  8.3621E-01  1.1190E+00  9.5617E-01  9.7917E-01  1.1831E+00  1.2200E+00  9.4798E-01
             1.0078E+00
 PARAMETER:  8.3053E-02  3.8563E-01  2.3297E-01 -7.8879E-02  2.1241E-01  5.5179E-02  7.8953E-02  2.6818E-01  2.9882E-01  4.6574E-02
             1.0773E-01
 GRADIENT:   6.2935E-01  6.7538E+00  3.9632E+00  3.1571E+00 -4.5270E+00 -3.3103E-01 -1.3291E+00 -3.2847E-01 -2.0648E+00 -1.7345E+00
            -6.8122E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1662.46788342658        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.8359E-01  1.4688E+00  8.1058E-01  7.3859E-01  1.0645E+00  9.5569E-01  9.8683E-01  6.8016E-01  1.2496E+00  8.8643E-01
             1.0040E+00
 PARAMETER:  8.3450E-02  4.8446E-01 -1.1001E-01 -2.0302E-01  1.6246E-01  5.4683E-02  8.6738E-02 -2.8543E-01  3.2281E-01 -2.0558E-02
             1.0400E-01
 GRADIENT:  -1.9895E+00  1.2421E+01  2.8868E+00  5.5111E+00 -7.8826E+00 -1.1899E+00 -5.2881E-01  2.4372E-01 -3.7167E+00 -1.7968E+00
            -1.6015E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1662.63055978075        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.8375E-01  1.6545E+00  6.4990E-01  6.2275E-01  1.0907E+00  9.5433E-01  8.9613E-01  4.2582E-01  1.4122E+00  8.8960E-01
             1.0044E+00
 PARAMETER:  8.3618E-02  6.0352E-01 -3.3094E-01 -3.7361E-01  1.8686E-01  5.3249E-02 -9.6682E-03 -7.5374E-01  4.4517E-01 -1.6988E-02
             1.0434E-01
 GRADIENT:  -3.4724E+00  2.2753E+01 -7.4484E-01  1.4370E+01 -8.1767E+00 -2.0874E+00 -3.4283E+00  3.5930E-01 -2.6968E+00 -1.2997E+00
            -1.5624E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1663.31094523785        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1062             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8573E-01  1.8286E+00  5.7051E-01  4.8068E-01  1.1708E+00  9.5947E-01  8.4106E-01  1.4436E-01  1.6955E+00  9.4263E-01
             1.0086E+00
 PARAMETER:  8.5627E-02  7.0355E-01 -4.6122E-01 -6.3255E-01  2.5770E-01  5.8622E-02 -7.3092E-02 -1.8354E+00  6.2797E-01  4.0923E-02
             1.0853E-01
 GRADIENT:   4.1380E+02  7.5428E+02  6.0992E+00  9.5730E+01  3.4672E+00  4.1483E+01  8.3652E+00  4.9713E-02  2.5266E+01  8.6333E-01
             8.5306E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1663.41930282807        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  9.8560E-01  1.8339E+00  5.6613E-01  4.8439E-01  1.1773E+00  9.5942E-01  8.4092E-01  4.8576E-02  1.7087E+00  9.4428E-01
             1.0083E+00
 PARAMETER:  8.5492E-02  7.0643E-01 -4.6893E-01 -6.2487E-01  2.6324E-01  5.8572E-02 -7.3256E-02 -2.9246E+00  6.3573E-01  4.2667E-02
             1.0831E-01
 GRADIENT:   1.2404E+00 -8.1926E+00  2.0780E-01 -1.2047E-01 -1.7102E+00  9.9105E-02 -5.0255E-02  4.5902E-03 -2.1218E-01 -1.0868E-01
            -4.1389E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1663.42367489745        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1427
 NPARAMETR:  9.8575E-01  1.8328E+00  5.6647E-01  4.8353E-01  1.1790E+00  9.5947E-01  8.4091E-01  1.0000E-02  1.7130E+00  9.4574E-01
             1.0084E+00
 PARAMETER:  8.5647E-02  7.0587E-01 -4.6834E-01 -6.2664E-01  2.6463E-01  5.8621E-02 -7.3273E-02 -5.0923E+00  6.3827E-01  4.4209E-02
             1.0834E-01
 GRADIENT:   1.6427E+00 -1.1341E+01 -1.6704E-01 -6.0199E-01 -2.7348E-01  1.2239E-01  3.8330E-02  0.0000E+00  5.2899E-02 -1.0720E-02
             4.7410E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1663.42367489745        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1449
 NPARAMETR:  9.8575E-01  1.8328E+00  5.6647E-01  4.8353E-01  1.1790E+00  9.5947E-01  8.4091E-01  1.0000E-02  1.7130E+00  9.4574E-01
             1.0084E+00
 PARAMETER:  8.5647E-02  7.0587E-01 -4.6834E-01 -6.2664E-01  2.6463E-01  5.8621E-02 -7.3273E-02 -5.0923E+00  6.3827E-01  4.4209E-02
             1.0834E-01
 GRADIENT:   1.6427E+00 -1.1341E+01 -1.6704E-01 -6.0199E-01 -2.7348E-01  1.2239E-01  3.8330E-02  0.0000E+00  5.2899E-02 -1.0720E-02
             4.7410E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1449
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8746E-05 -3.1586E-02 -2.4050E-04  3.0706E-02 -3.9687E-02
 SE:             2.9848E-02  2.4107E-02  9.0585E-05  2.2384E-02  2.2068E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9843E-01  1.9012E-01  7.9324E-03  1.7014E-01  7.2112E-02

 ETASHRINKSD(%)  4.2291E-03  1.9237E+01  9.9697E+01  2.5010E+01  2.6071E+01
 ETASHRINKVR(%)  8.4581E-03  3.4774E+01  9.9999E+01  4.3765E+01  4.5344E+01
 EBVSHRINKSD(%)  4.6924E-01  1.8020E+01  9.9742E+01  2.7621E+01  2.4611E+01
 EBVSHRINKVR(%)  9.3628E-01  3.2793E+01  9.9999E+01  4.7613E+01  4.3165E+01
 RELATIVEINF(%)  9.9013E+01  6.5961E+00  1.1374E-04  5.1107E+00  1.5590E+01
 EPSSHRINKSD(%)  4.4228E+01
 EPSSHRINKVR(%)  6.8895E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1663.4236748974549     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -928.27284833371675     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.15
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1663.424       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.83E+00  5.66E-01  4.84E-01  1.18E+00  9.59E-01  8.41E-01  1.00E-02  1.71E+00  9.46E-01  1.01E+00
 


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
+        1.23E+03
 
 TH 2
+       -5.21E+00  3.32E+02
 
 TH 3
+        7.71E+00  1.16E+02  2.69E+02
 
 TH 4
+       -1.40E+01  2.72E+02 -2.38E+02  8.86E+02
 
 TH 5
+       -5.21E+00 -1.50E+02 -2.31E+02  2.45E+02  5.09E+02
 
 TH 6
+       -2.22E-01 -6.41E-01  1.89E+00 -3.34E+00 -6.92E-01  2.13E+02
 
 TH 7
+        7.27E-01  7.22E+00 -7.54E+00 -1.42E+01 -1.58E+01 -5.69E-01  1.33E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.57E-01 -1.51E+01 -2.80E+01  5.61E+01  1.33E+00 -4.50E-01  1.58E+01  0.00E+00  2.62E+01
 
 TH10
+       -3.09E-03 -1.38E+01 -2.58E+01 -3.89E+00 -6.27E+01  6.75E-01  9.87E+00  0.00E+00  4.97E+00  8.80E+01
 
 TH11
+       -7.65E+00 -1.60E+01 -1.85E+01  1.18E+00 -6.56E+00  2.63E+00  8.01E+00  0.00E+00  3.87E+00  1.74E+01  2.12E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       25.514
Stop Time:
Wed Sep 29 19:30:44 CDT 2021
