Sat Sep 25 13:48:08 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat79.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1719.70819664378        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.3404E+00 -6.9845E+01 -4.1065E+01 -4.4645E+01  2.9515E+01  2.5241E+01 -3.8593E+00  1.1749E+01  1.7236E+01  2.0150E+01
             5.0678E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1726.52234643636        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0328E+00  1.1476E+00  1.1503E+00  9.6782E-01  1.0886E+00  9.1044E-01  1.0753E+00  9.1214E-01  8.6573E-01  8.5135E-01
             1.0354E+00
 PARAMETER:  1.3228E-01  2.3770E-01  2.4006E-01  6.7292E-02  1.8491E-01  6.1719E-03  1.7256E-01  8.0338E-03 -4.4183E-02 -6.0928E-02
             1.3483E-01
 GRADIENT:   7.7432E+01  2.7895E+01  1.6346E+01  1.8275E+01  9.8028E+00 -9.4733E+00 -4.6232E+00 -4.6372E+00 -6.4103E+00 -1.3485E+01
             1.6714E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1727.43766649876        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0272E+00  1.0416E+00  1.0261E+00  1.0230E+00  9.9169E-01  9.2050E-01  1.2278E+00  6.8909E-01  8.2985E-01  7.5344E-01
             1.0166E+00
 PARAMETER:  1.2685E-01  1.4077E-01  1.2577E-01  1.2270E-01  9.1659E-02  1.7160E-02  3.0523E-01 -2.7239E-01 -8.6509E-02 -1.8311E-01
             1.1646E-01
 GRADIENT:   6.2149E+01  1.7847E+01  1.6387E+01  2.0198E+01  1.3424E+01 -4.7992E+00  5.6004E-01 -4.1271E+00 -2.2967E+00 -1.8415E+01
            -4.6231E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1729.51821310746        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0090E+00  1.0449E+00  8.2323E-01  9.9779E-01  9.0416E-01  9.3626E-01  1.2238E+00  4.0146E-01  8.4657E-01  7.8451E-01
             9.9839E-01
 PARAMETER:  1.0897E-01  1.4392E-01 -9.4520E-02  9.7787E-02 -7.4372E-04  3.4138E-02  3.0199E-01 -8.1265E-01 -6.6559E-02 -1.4269E-01
             9.8385E-02
 GRADIENT:   6.6201E+00  2.9315E-02 -7.3754E+00  1.1454E+01  1.3144E+01  8.2238E-01  2.5198E+00  8.0229E-01  2.9970E+00 -6.6700E-01
            -8.8761E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1729.53792171642        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0071E+00  1.0496E+00  7.9629E-01  9.8952E-01  8.8877E-01  9.3607E-01  1.2092E+00  3.4144E-01  8.3945E-01  7.7843E-01
             9.9983E-01
 PARAMETER:  1.0705E-01  1.4837E-01 -1.2780E-01  8.9462E-02 -1.7916E-02  3.3940E-02  2.8997E-01 -9.7459E-01 -7.5013E-02 -1.5047E-01
             9.9828E-02
 GRADIENT:   4.4080E-02 -1.1274E+00 -4.0332E+00  3.7120E+00  4.6886E+00  5.0812E-01  1.0269E+00  6.8634E-01  1.1667E+00  6.9010E-01
             1.0885E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1729.56451211333        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0062E+00  1.0593E+00  7.6869E-01  9.7903E-01  8.7595E-01  9.3557E-01  1.1922E+00  2.4442E-01  8.3731E-01  7.7112E-01
             1.0006E+00
 PARAMETER:  1.0614E-01  1.5764E-01 -1.6306E-01  7.8812E-02 -3.2448E-02  3.3401E-02  2.7579E-01 -1.3089E+00 -7.7561E-02 -1.5991E-01
             1.0059E-01
 GRADIENT:  -3.5272E+00 -1.3688E+00 -8.1749E-01 -2.0422E+00 -1.8010E+00  9.6757E-02 -2.2364E-01  3.8540E-01 -2.9628E-01  1.2744E+00
             6.5005E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1729.81527409498        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.0077E+00  1.0677E+00  7.6102E-01  9.7464E-01  8.7710E-01  9.3556E-01  1.1908E+00  3.3356E-02  8.4581E-01  7.6802E-01
             9.9888E-01
 PARAMETER:  1.0763E-01  1.6552E-01 -1.7309E-01  7.4308E-02 -3.1136E-02  3.3390E-02  2.7459E-01 -3.3005E+00 -6.7456E-02 -1.6394E-01
             9.8882E-02
 GRADIENT:   8.6925E-01  4.8906E-01  1.4545E-01  1.2340E+00  1.2539E+00  1.6619E-01  3.3355E-01  6.0757E-03  1.4001E-01 -5.1958E-01
            -3.0141E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1730.33523940956        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  1.0284E+00  1.1793E+00  7.4003E-01  9.1489E-01  9.1620E-01  9.4298E-01  1.1033E+00  1.0000E-02  8.9717E-01  7.7849E-01
             1.0014E+00
 PARAMETER:  1.2805E-01  2.6492E-01 -2.0106E-01  1.1044E-02  1.2483E-02  4.1292E-02  1.9827E-01 -5.3429E+00 -8.5145E-03 -1.5040E-01
             1.0140E-01
 GRADIENT:   4.5090E+00  3.6984E+00  6.2393E-01  4.4043E+00 -2.7693E-02  6.4110E-02  2.4011E-01  0.0000E+00 -3.7686E-02 -1.0963E+00
            -2.3724E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1730.36195160334        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  1.0268E+00  1.2326E+00  7.2537E-01  8.7924E-01  9.3633E-01  9.4315E-01  1.0576E+00  1.0000E-02  9.2574E-01  7.9172E-01
             1.0014E+00
 PARAMETER:  1.2643E-01  3.0909E-01 -2.2107E-01 -2.8703E-02  3.4210E-02  4.1474E-02  1.5599E-01 -5.3097E+00  2.2839E-02 -1.3355E-01
             1.0145E-01
 GRADIENT:   1.1577E-01  1.3796E-01  7.9943E-02  5.6366E-02 -1.0894E-01  3.6886E-03 -6.0146E-03  0.0000E+00 -1.4838E-02 -2.3805E-02
            -8.9008E-03

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1730.36196690706        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      890
 NPARAMETR:  1.0267E+00  1.2338E+00  7.2491E-01  8.7837E-01  9.3680E-01  9.4315E-01  1.0566E+00  1.0000E-02  9.2650E-01  7.9197E-01
             1.0014E+00
 PARAMETER:  1.2639E-01  3.1012E-01 -2.2170E-01 -2.9684E-02  3.4720E-02  4.1475E-02  1.5506E-01 -5.3148E+00  2.3655E-02 -1.3323E-01
             1.0145E-01
 GRADIENT:   1.1888E-03  1.2326E-05  1.5219E-03 -2.1897E-03 -2.3755E-03  1.3631E-04  6.7184E-05  0.0000E+00  8.1066E-05 -1.6479E-04
             2.1465E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      890
 NO. OF SIG. DIGITS IN FINAL EST.:  4.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.4028E-05 -7.3300E-03 -3.8595E-04  2.9881E-03 -1.8503E-02
 SE:             2.9805E-02  2.3308E-02  1.6128E-04  2.3157E-02  2.1883E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9775E-01  7.5316E-01  1.6706E-02  8.9733E-01  3.9781E-01

 ETASHRINKSD(%)  1.4995E-01  2.1914E+01  9.9460E+01  2.2420E+01  2.6688E+01
 ETASHRINKVR(%)  2.9967E-01  3.9026E+01  9.9997E+01  3.9813E+01  4.6253E+01
 EBVSHRINKSD(%)  5.0119E-01  2.1466E+01  9.9502E+01  2.3062E+01  2.6026E+01
 EBVSHRINKVR(%)  9.9987E-01  3.8323E+01  9.9998E+01  4.0806E+01  4.5278E+01
 RELATIVEINF(%)  9.8745E+01  3.7124E+00  2.5548E-04  3.7037E+00  5.5966E+00
 EPSSHRINKSD(%)  4.3151E+01
 EPSSHRINKVR(%)  6.7682E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1730.3619669070615     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -995.21114034332334     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1730.362       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.23E+00  7.25E-01  8.78E-01  9.37E-01  9.43E-01  1.06E+00  1.00E-02  9.26E-01  7.92E-01  1.00E+00
 


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
+        1.17E+03
 
 TH 2
+       -6.86E+00  4.17E+02
 
 TH 3
+        1.68E+01  1.98E+02  5.80E+02
 
 TH 4
+       -1.42E+01  3.39E+02 -3.58E+02  1.00E+03
 
 TH 5
+       -5.30E+00 -3.20E+02 -6.76E+02  4.15E+02  1.10E+03
 
 TH 6
+       -1.93E-01 -1.57E+00  5.47E+00 -5.25E+00 -2.57E+00  2.18E+02
 
 TH 7
+        1.70E+00  2.89E+01 -1.59E+01 -1.59E+01  1.34E+00 -1.16E+00  7.32E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.79E+00 -2.20E+01 -3.54E+01  4.51E+01  4.91E-01 -2.55E+00  1.88E+01  0.00E+00  9.28E+01
 
 TH10
+       -1.96E+00 -1.24E+01 -5.49E+01 -2.00E+01 -6.73E+01 -6.49E-01  1.49E+01  0.00E+00  1.37E+01  1.01E+02
 
 TH11
+       -1.02E+01 -1.63E+01 -3.10E+01  3.29E+00  2.36E+00  3.28E+00  6.95E+00  0.00E+00  1.21E+01  2.72E+01  2.12E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.119
Stop Time:
Sat Sep 25 13:48:25 CDT 2021
