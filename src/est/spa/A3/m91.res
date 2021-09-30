Wed Sep 29 13:54:07 CDT 2021
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
$DATA ../../../../data/spa/A3/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   95.4258968522314        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.5455E+02  2.9180E+01  3.9856E+01 -1.4020E+01  1.7856E+02 -8.2850E+00 -4.7715E+01 -1.0468E+01 -9.8853E+01 -1.0530E+02
            -3.0972E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1215.65212376286        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.5704E-01  9.8807E-01  1.0304E+00  1.0813E+00  9.8656E-01  9.7664E-01  9.4444E-01  9.4101E-01  1.0618E+00  7.5920E-01
             3.3414E+00
 PARAMETER:  5.6087E-02  8.8001E-02  1.2994E-01  1.7813E-01  8.6472E-02  7.6363E-02  4.2832E-02  3.9200E-02  1.5993E-01 -1.7549E-01
             1.3064E+00
 GRADIENT:   4.4591E+01 -1.9500E+01 -1.8382E+01 -8.5642E+00  4.4407E+01 -2.9167E+00  4.6935E+00  4.5616E+00 -2.7308E+00  1.2096E+01
            -8.2932E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1231.71844519905        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.5272E-01  7.7096E-01  5.6203E-01  1.1914E+00  5.8632E-01  1.0509E+00  3.9441E-01  4.5591E-01  1.1664E+00  3.8636E-01
             3.4147E+00
 PARAMETER:  5.1563E-02 -1.6012E-01 -4.7621E-01  2.7509E-01 -4.3389E-01  1.4968E-01 -8.3038E-01 -6.8547E-01  2.5393E-01 -8.5099E-01
             1.3281E+00
 GRADIENT:   1.7796E+01  2.3432E+01 -2.3779E+00  4.5455E+01  6.7356E+00  1.8638E+01 -1.1743E+00  3.1466E+00  5.2159E-01  6.9877E+00
            -4.0405E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1237.80652024844        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.4177E-01  5.9083E-01  4.6414E-01  1.2733E+00  4.7496E-01  9.7188E-01  3.8889E-01  2.1864E-01  1.1559E+00  2.2721E-01
             3.6532E+00
 PARAMETER:  4.0010E-02 -4.2623E-01 -6.6757E-01  3.4159E-01 -6.4452E-01  7.1481E-02 -8.4445E-01 -1.4203E+00  2.4487E-01 -1.3819E+00
             1.3956E+00
 GRADIENT:  -2.1471E+01  1.1367E+01 -2.0477E+01  6.2726E+01  3.5304E+01 -9.5714E+00 -3.1810E-01  1.1502E+00  1.1599E+01  2.8033E+00
             1.9820E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1240.81637123826        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      367
 NPARAMETR:  9.5785E-01  5.0436E-01  4.5903E-01  1.2853E+00  4.3422E-01  9.8678E-01  6.5582E-01  2.0114E-01  1.0842E+00  1.2770E-01
             3.7588E+00
 PARAMETER:  5.6936E-02 -5.8447E-01 -6.7864E-01  3.5101E-01 -7.3420E-01  8.6688E-02 -3.2187E-01 -1.5038E+00  1.8084E-01 -1.9580E+00
             1.4241E+00
 GRADIENT:  -1.4962E+01  1.8535E+01  2.2378E+01  1.4315E+01 -3.9571E+01 -4.3308E+00 -5.1499E-01  8.4400E-01  7.3598E+00  6.6020E-01
             3.0313E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1244.38674390752        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      542
 NPARAMETR:  9.5552E-01  3.4376E-01  4.3644E-01  1.3368E+00  3.9008E-01  9.9721E-01  1.1633E+00  3.4820E-02  9.9743E-01  4.3284E-02
             3.5733E+00
 PARAMETER:  5.4503E-02 -9.6781E-01 -7.2910E-01  3.9027E-01 -8.4140E-01  9.7205E-02  2.5130E-01 -3.2576E+00  9.7428E-02 -3.0400E+00
             1.3735E+00
 GRADIENT:  -6.8496E+00  6.0792E+00  1.2635E+01  1.7274E+00 -1.8870E+01 -1.3942E+00 -3.2851E-01  1.9110E-02 -8.6271E-02  3.6290E-02
            -9.2096E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1245.41640546794        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      717
 NPARAMETR:  9.5229E-01  1.6849E-01  4.2632E-01  1.4085E+00  3.5942E-01  9.9763E-01  2.0410E+00  1.0000E-02  9.4845E-01  1.0000E-02
             3.5994E+00
 PARAMETER:  5.1115E-02 -1.6809E+00 -7.5256E-01  4.4254E-01 -9.2326E-01  9.7626E-02  8.1343E-01 -8.4890E+00  4.7074E-02 -5.2146E+00
             1.3808E+00
 GRADIENT:  -1.5695E+00  8.9394E-01  7.3370E+00  4.9371E+00 -1.2010E+01  1.0420E-02 -4.3609E-01  0.0000E+00 -1.2199E+00  0.0000E+00
            -2.9485E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1248.11507913732        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      896
 NPARAMETR:  9.4758E-01  9.9719E-02  4.3360E-01  1.4549E+00  3.5486E-01  9.9180E-01  3.9994E+00  1.0000E-02  9.2599E-01  1.0000E-02
             3.6104E+00
 PARAMETER:  4.6151E-02 -2.2054E+00 -7.3564E-01  4.7491E-01 -9.3602E-01  9.1768E-02  1.4862E+00 -1.2951E+01  2.3111E-02 -6.9643E+00
             1.3838E+00
 GRADIENT:  -5.7830E+00  6.2322E+00  1.3416E+01  2.1507E+01 -2.4936E+01 -2.6025E+00  5.0340E+00  0.0000E+00  1.3351E+00  0.0000E+00
             3.3317E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1249.25371427580        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  9.4452E-01  5.6900E-02  4.0610E-01  1.4449E+00  3.3854E-01  9.9344E-01  5.3096E+00  1.0000E-02  9.1936E-01  1.0000E-02
             3.5634E+00
 PARAMETER:  4.2925E-02 -2.7665E+00 -8.0114E-01  4.6806E-01 -9.8310E-01  9.3420E-02  1.7695E+00 -1.7654E+01  1.5925E-02 -9.1619E+00
             1.3707E+00
 GRADIENT:  -4.9907E+00  1.6948E+00 -6.2583E+00  1.4199E+01  3.8519E+00 -1.9139E+00  1.2528E+00  0.0000E+00 -2.0762E+00  0.0000E+00
            -3.3047E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1249.52354821664        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1250
 NPARAMETR:  9.4568E-01  4.4040E-02  3.8635E-01  1.4214E+00  3.2508E-01  9.9881E-01  5.9293E+00  1.0000E-02  9.3829E-01  1.0000E-02
             3.5488E+00
 PARAMETER:  4.4150E-02 -3.0227E+00 -8.5100E-01  4.5168E-01 -1.0237E+00  9.8809E-02  1.8799E+00 -1.9713E+01  3.6303E-02 -1.0260E+01
             1.3666E+00
 GRADIENT:   2.0308E+00  8.4115E+00 -6.2991E+00 -1.4764E+01  3.9100E+00 -2.7689E+00  1.5653E+01  0.0000E+00 -5.2191E-01  0.0000E+00
            -1.0580E+01

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1249.52354821664        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1272
 NPARAMETR:  9.4568E-01  4.4040E-02  3.8635E-01  1.4214E+00  3.2508E-01  9.9881E-01  5.9293E+00  1.0000E-02  9.3829E-01  1.0000E-02
             3.5488E+00
 PARAMETER:  4.4150E-02 -3.0227E+00 -8.5100E-01  4.5168E-01 -1.0237E+00  9.8809E-02  1.8799E+00 -1.9713E+01  3.6303E-02 -1.0260E+01
             1.3666E+00
 GRADIENT:   2.0308E+00  8.4115E+00 -6.2991E+00 -1.4764E+01  3.9100E+00 -2.7689E+00  1.5653E+01  0.0000E+00 -5.2191E-01  0.0000E+00
            -1.0580E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1272
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.9432E-04  6.7793E-03  1.4812E-04 -1.7165E-02  3.6604E-05
 SE:             2.8404E-02  7.0778E-03  2.3454E-04  2.5426E-02  3.4218E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8892E-01  3.3815E-01  5.2768E-01  4.9961E-01  9.1481E-01

 ETASHRINKSD(%)  4.8442E+00  7.6289E+01  9.9214E+01  1.4820E+01  9.8854E+01
 ETASHRINKVR(%)  9.4537E+00  9.4378E+01  9.9994E+01  2.7443E+01  9.9987E+01
 EBVSHRINKSD(%)  4.1275E+00  8.3938E+01  9.9156E+01  1.3826E+01  9.8841E+01
 EBVSHRINKVR(%)  8.0846E+00  9.7420E+01  9.9993E+01  2.5741E+01  9.9987E+01
 RELATIVEINF(%)  8.9184E+01  1.1170E+00  2.5427E-04  2.5074E+01  4.6180E-04
 EPSSHRINKSD(%)  2.2073E+01
 EPSSHRINKVR(%)  3.9275E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1249.5235482166372     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -514.37272165289903     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1249.524       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.46E-01  4.40E-02  3.86E-01  1.42E+00  3.25E-01  9.99E-01  5.93E+00  1.00E-02  9.38E-01  1.00E-02  3.55E+00
 


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
+        2.14E+05
 
 TH 2
+       -3.65E+02  4.18E+05
 
 TH 3
+       -2.39E+01  1.35E+05  2.19E+04
 
 TH 4
+       -2.18E+01 -6.52E+04 -4.64E+02  4.93E+03
 
 TH 5
+        1.92E+02 -1.34E+05 -7.10E+03 -4.30E+01  3.03E+04
 
 TH 6
+       -2.01E+05  3.00E+02 -7.05E+01 -5.89E+01  1.12E+02  1.91E+05
 
 TH 7
+       -2.91E+00  3.67E+03  7.72E+01  2.03E+01 -1.57E+03  3.59E+00  5.87E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.15E+05  6.50E+02  5.13E+01 -1.86E+01  8.78E+01 -1.49E+02  7.86E+00  0.00E+00  2.17E+05
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.61E+01 -8.28E+03 -7.20E+01 -2.31E+01  9.66E+01 -1.61E+00 -9.76E+01  0.00E+00  6.69E+00  0.00E+00  1.11E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.737
Stop Time:
Wed Sep 29 13:54:32 CDT 2021
