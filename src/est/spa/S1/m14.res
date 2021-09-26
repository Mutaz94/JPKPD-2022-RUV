Sat Sep 25 09:43:22 CDT 2021
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
$DATA ../../../../data/spa/S1/dat14.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
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
 NO. OF SIG. FIGURES REQUIRED:            3
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1701.55521116683        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5439E+01 -1.0261E+02 -6.1168E+01 -5.9583E+01  7.5348E+01 -2.2993E+00  1.3158E+00  1.3975E+01  3.1415E+01  7.6979E+00
             1.4379E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1712.34523557372        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9585E-01  1.1184E+00  1.1371E+00  1.0072E+00  1.0492E+00  1.0012E+00  9.6199E-01  8.8105E-01  8.1193E-01  9.7271E-01
             1.0482E+00
 PARAMETER:  9.5841E-02  2.1189E-01  2.2845E-01  1.0714E-01  1.4800E-01  1.0118E-01  6.1245E-02 -2.6637E-02 -1.0834E-01  7.2336E-02
             1.4705E-01
 GRADIENT:   1.0952E+00  2.1571E+01 -5.7514E+00  3.3803E+01  1.7950E+01 -2.2441E+00 -3.8943E+00  3.1080E+00 -8.4288E+00 -7.1694E+00
             1.0220E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1713.96509301746        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0072E+00  9.8576E-01  1.0941E+00  1.0872E+00  9.8502E-01  1.0150E+00  9.7308E-01  3.7937E-01  8.7165E-01  1.0077E+00
             1.0354E+00
 PARAMETER:  1.0716E-01  8.5663E-02  1.8996E-01  1.8361E-01  8.4910E-02  1.1485E-01  7.2709E-02 -8.6924E-01 -3.7366E-02  1.0765E-01
             1.3476E-01
 GRADIENT:   3.1264E+01  5.9929E+00 -8.3015E+00  4.3701E+01  1.9030E+01  4.3887E+00 -3.0315E+00  3.2816E-01  8.1879E+00 -5.8854E+00
             9.5410E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1714.64765935547        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.9709E-01  1.0402E+00  9.7699E-01  1.0315E+00  9.5106E-01  1.0085E+00  1.0501E+00  3.5132E-01  8.3166E-01  9.7426E-01
             1.0114E+00
 PARAMETER:  9.7084E-02  1.3945E-01  7.6717E-02  1.3101E-01  4.9818E-02  1.0848E-01  1.4884E-01 -9.4606E-01 -8.4332E-02  7.3924E-02
             1.1137E-01
 GRADIENT:   6.4722E+00 -1.5141E+00 -6.6375E+00  9.3296E+00  8.4517E+00  1.1345E+00  7.4704E-01  7.7692E-01  1.9987E+00  6.7066E-01
             2.3242E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1714.67600943303        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.9407E-01  1.0621E+00  9.5753E-01  1.0138E+00  9.4799E-01  1.0061E+00  1.0291E+00  2.5996E-01  8.3484E-01  9.7415E-01
             1.0072E+00
 PARAMETER:  9.4053E-02  1.6021E-01  5.6603E-02  1.1371E-01  4.6588E-02  1.0608E-01  1.2866E-01 -1.2472E+00 -8.0518E-02  7.3812E-02
             1.0720E-01
 GRADIENT:  -2.3736E-01 -1.0946E+00 -1.3207E+00 -2.2093E-02  5.5633E-01 -1.5708E-01 -4.4178E-02  3.7129E-01  4.2736E-01  6.1173E-01
             2.8191E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1714.92153342645        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.9617E-01  1.0535E+00  9.5407E-01  1.0213E+00  9.4380E-01  1.0083E+00  1.0432E+00  1.2661E-02  8.2951E-01  9.6606E-01
             1.0097E+00
 PARAMETER:  9.6159E-02  1.5207E-01  5.2977E-02  1.2111E-01  4.2157E-02  1.0827E-01  1.4233E-01 -4.2692E+00 -8.6919E-02  6.5475E-02
             1.0968E-01
 GRADIENT:   4.2358E+00  1.6959E+00 -1.0199E-01  5.4582E+00  2.9297E+00  9.3889E-01  1.8490E-01  6.9213E-04  1.9524E-01 -1.7090E+00
             6.6926E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1715.35529488500        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      555
 NPARAMETR:  1.0155E+00  1.1656E+00  9.1389E-01  9.5568E-01  9.7458E-01  1.0200E+00  9.6133E-01  1.3370E-02  8.7434E-01  9.8028E-01
             1.0094E+00
 PARAMETER:  1.1541E-01  2.5321E-01  9.9572E-03  5.4673E-02  7.4248E-02  1.1984E-01  6.0568E-02 -4.2148E+00 -3.4286E-02  8.0078E-02
             1.0940E-01
 GRADIENT:   2.5401E+00  2.2390E+00 -8.0418E-01  4.6094E+00  1.7091E+00  6.8710E-01 -8.1678E-02  7.6496E-04 -2.9067E-01 -5.2404E-01
             9.2727E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1715.37673144906        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  1.0150E+00  1.2324E+00  8.8938E-01  9.1141E-01  9.9440E-01  1.0191E+00  9.1634E-01  2.8061E-02  9.0926E-01  9.8862E-01
             1.0097E+00
 PARAMETER:  1.1491E-01  3.0899E-01 -1.7226E-02  7.2331E-03  9.4385E-02  1.1892E-01  1.2631E-02 -3.4734E+00  4.8809E-03  8.8552E-02
             1.0970E-01
 GRADIENT:   9.3170E-01  8.0720E-01  4.5789E-02  9.3748E-01  7.3039E-03  2.0445E-01  5.5291E-02  3.4268E-03  7.5861E-02 -4.6113E-02
             5.2749E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1715.37889880332        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      905
 NPARAMETR:  1.0146E+00  1.2349E+00  8.8680E-01  9.0915E-01  9.9449E-01  1.0186E+00  9.1465E-01  1.0000E-02  9.1011E-01  9.8815E-01
             1.0095E+00
 PARAMETER:  1.1449E-01  3.1100E-01 -2.0130E-02  4.7497E-03  9.4471E-02  1.1842E-01  1.0782E-02 -4.5172E+00  5.8122E-03  8.8077E-02
             1.0950E-01
 GRADIENT:   3.7342E-04 -1.0645E-03 -1.4516E-03  7.0493E-04  3.1176E-03  7.4455E-04 -2.4293E-05  0.0000E+00  1.2554E-04 -7.5092E-04
             6.8764E-04

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1715.37889880332        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      927
 NPARAMETR:  1.0146E+00  1.2349E+00  8.8680E-01  9.0915E-01  9.9449E-01  1.0186E+00  9.1465E-01  1.0000E-02  9.1011E-01  9.8815E-01
             1.0095E+00
 PARAMETER:  1.1449E-01  3.1100E-01 -2.0130E-02  4.7497E-03  9.4471E-02  1.1842E-01  1.0782E-02 -4.5172E+00  5.8122E-03  8.8077E-02
             1.0950E-01
 GRADIENT:   3.7342E-04 -1.0645E-03 -1.4516E-03  7.0493E-04  3.1176E-03  7.4455E-04 -2.4293E-05  0.0000E+00  1.2554E-04 -7.5092E-04
             6.8764E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      927
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.5297E-04 -1.2925E-02 -3.1833E-04  4.0959E-03 -2.2956E-02
 SE:             2.9829E-02  2.0784E-02  1.3575E-04  2.3224E-02  2.3866E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9323E-01  5.3403E-01  1.9024E-02  8.6001E-01  3.3611E-01

 ETASHRINKSD(%)  7.0501E-02  3.0371E+01  9.9545E+01  2.2196E+01  2.0045E+01
 ETASHRINKVR(%)  1.4095E-01  5.1518E+01  9.9998E+01  3.9466E+01  3.6072E+01
 EBVSHRINKSD(%)  4.1506E-01  3.0062E+01  9.9572E+01  2.3013E+01  1.8114E+01
 EBVSHRINKVR(%)  8.2840E-01  5.1087E+01  9.9998E+01  4.0730E+01  3.2947E+01
 RELATIVEINF(%)  9.8749E+01  1.4084E+00  1.6060E-04  1.8690E+00  7.7925E+00
 EPSSHRINKSD(%)  4.2831E+01
 EPSSHRINKVR(%)  6.7317E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1715.3788988033243     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -980.22807223958614     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1715.379       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.23E+00  8.87E-01  9.09E-01  9.94E-01  1.02E+00  9.15E-01  1.00E-02  9.10E-01  9.88E-01  1.01E+00
 


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
+        1.03E+03
 
 TH 2
+       -7.29E+00  4.44E+02
 
 TH 3
+        1.11E+01  1.32E+02  2.87E+02
 
 TH 4
+       -1.02E+01  4.56E+02 -1.81E+02  9.66E+02
 
 TH 5
+       -3.31E+00 -2.49E+02 -3.87E+02  2.19E+02  7.75E+02
 
 TH 6
+        7.52E-01 -1.32E+00  4.26E+00 -2.37E+00 -4.59E+00  1.90E+02
 
 TH 7
+        4.05E+00  2.04E+01  2.95E+00 -1.21E+01 -1.22E+01 -2.13E+00  6.69E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.69E-01 -2.09E+01 -2.20E+01  3.27E+01  8.94E+00 -7.94E-01  3.12E+01  0.00E+00  9.38E+01
 
 TH10
+       -4.76E-01 -7.96E+00 -3.38E+01 -1.47E+01 -5.71E+01  2.63E+00  1.56E+01  0.00E+00  1.13E+01  8.95E+01
 
 TH11
+       -6.36E+00 -2.08E+01 -3.28E+01  9.48E-01  5.09E+00  1.54E+00  3.24E+00  0.00E+00  1.16E+01  2.38E+01  2.19E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.320
Stop Time:
Sat Sep 25 09:43:39 CDT 2021
