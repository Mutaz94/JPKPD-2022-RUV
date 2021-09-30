Thu Sep 30 03:00:12 CDT 2021
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
$DATA ../../../../data/spa1/D/dat37.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17057.9093917011        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1508E+02  3.0439E+02  9.4144E+01  6.3606E+01  1.0043E+02 -1.3415E+03 -6.1986E+02 -6.0108E+02 -1.2876E+03 -3.5295E+02
            -3.3776E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -691.705647258370        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4524E+00  1.0560E+00  9.9494E-01  1.6439E+00  1.1962E+00  2.4194E+00  1.2799E+00  9.6348E-01  1.6068E+00  1.0673E+00
             1.3953E+01
 PARAMETER:  4.7325E-01  1.5448E-01  9.4929E-02  5.9707E-01  2.7918E-01  9.8350E-01  3.4677E-01  6.2796E-02  5.7422E-01  1.6516E-01
             2.7357E+00
 GRADIENT:   2.4446E+01  1.1095E+01 -8.0441E+00  1.3078E+00 -3.8897E+00  8.8824E+01 -8.0169E-01  5.2636E+00 -8.7929E+00  3.5864E+00
             2.4248E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -727.088010916345        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.4090E+00  7.6872E-01  1.7384E+00  2.0903E+00  5.4448E+00  2.0740E+00  5.0331E+00  3.3638E-01  2.1431E+00  7.4729E+00
             1.2253E+01
 PARAMETER:  4.4289E-01 -1.6303E-01  6.5296E-01  8.3732E-01  1.7947E+00  8.2949E-01  1.7160E+00 -9.8951E-01  8.6227E-01  2.1113E+00
             2.6058E+00
 GRADIENT:   3.9602E+01  1.7389E+01  7.1231E+00  5.0545E+01 -2.5952E+00  1.7192E+01  1.4777E+01 -1.8059E-01  3.5281E+01  2.5578E+00
             1.9729E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -773.411526596378        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.9293E-01  4.3059E-01  6.8809E-01  1.4229E+00  6.1835E+00  1.7847E+00  2.4201E+00  1.8210E-02  1.1032E+00  6.0220E+00
             9.8311E+00
 PARAMETER:  9.2901E-02 -7.4259E-01 -2.7383E-01  4.5268E-01  1.9219E+00  6.7925E-01  9.8381E-01 -3.9058E+00  1.9822E-01  1.8954E+00
             2.3856E+00
 GRADIENT:  -6.7082E+01  2.7493E+01  2.0381E+01 -3.3560E+01 -1.6343E+01  1.7986E+01  3.7867E+00 -2.4209E-03 -3.2628E+01  6.6094E+00
             5.7815E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -823.838770493661        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  8.8145E-01  1.5063E-01  2.3459E-01  1.1658E+00  2.2047E+01  1.6132E+00  4.2772E-01  1.0000E-02  1.3042E+00  5.2042E+00
             7.7269E+00
 PARAMETER: -2.6182E-02 -1.7930E+00 -1.3499E+00  2.5343E-01  3.1932E+00  5.7824E-01 -7.4930E-01 -7.9023E+00  3.6561E-01  1.7495E+00
             2.1447E+00
 GRADIENT:  -3.8846E+00  1.8206E+01  1.1808E+01  6.5372E+01 -7.8098E-01 -2.5907E+01  5.1724E-01  0.0000E+00  8.4721E+00  1.8563E-01
            -1.5123E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -847.189665103224        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      397
 NPARAMETR:  7.3677E-01  4.3927E-02  6.5989E-02  6.0109E-01  1.4685E+02  1.4322E+00  1.0000E-02  1.0000E-02  1.0151E+00  5.4420E+00
             7.8600E+00
 PARAMETER: -2.0547E-01 -3.0252E+00 -2.6183E+00 -4.0901E-01  5.0894E+00  4.5922E-01 -5.1296E+00 -1.5425E+01  1.1503E-01  1.7941E+00
             2.1618E+00
 GRADIENT:   1.5414E+02  2.7366E+00 -9.8853E+01  9.3157E+01  5.3257E-02 -5.6021E+01  0.0000E+00  0.0000E+00 -3.3515E+00  2.2589E-03
            -1.7349E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -875.158764374484        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      572
 NPARAMETR:  5.4207E-01  2.1769E-02  4.1474E-02  4.1348E-01  6.5863E+02  1.4819E+00  1.0000E-02  1.0000E-02  9.0616E-01  6.3571E+00
             8.8470E+00
 PARAMETER: -5.1236E-01 -3.7273E+00 -3.0827E+00 -7.8315E-01  6.5902E+00  4.9334E-01 -8.1854E+00 -2.0021E+01  1.4637E-03  1.9496E+00
             2.2801E+00
 GRADIENT:  -1.0280E+01  1.0035E+00  1.1388E+00  5.3201E+00  3.7739E-04 -1.9194E+00  0.0000E+00  0.0000E+00 -8.2297E+00  1.2613E-06
             2.2428E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -875.900305122307        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      758
 NPARAMETR:  5.1905E-01  1.4992E-02  3.5933E-02  3.7109E-01  8.1619E+02  1.4750E+00  1.0000E-02  1.0000E-02  9.5137E-01  7.1318E+00
             8.8152E+00
 PARAMETER: -5.5576E-01 -4.1003E+00 -3.2261E+00 -8.9131E-01  6.8047E+00  4.8862E-01 -8.9098E+00 -2.1066E+01  5.0153E-02  2.0646E+00
             2.2765E+00
 GRADIENT:   1.3160E+00  1.6538E-02  3.3540E+00 -5.7741E+00  1.0676E-03 -1.9850E-01  0.0000E+00  0.0000E+00  2.0740E+00 -6.5768E-07
             1.1954E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -875.936301543121        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  5.1687E-01  1.3209E-02  3.6008E-02  3.7313E-01  2.8245E+01  1.4762E+00  1.0000E-02  1.0000E-02  9.4344E-01  1.0913E+01
             8.7954E+00
 PARAMETER: -5.5996E-01 -4.2269E+00 -3.2240E+00 -8.8584E-01  3.4409E+00  4.8944E-01 -8.9098E+00 -2.1066E+01  4.1774E-02  2.4899E+00
             2.2742E+00
 GRADIENT:  -3.7660E+00  8.5166E-03  8.7553E-01  1.5047E-01  2.7060E-02 -3.0559E-01  0.0000E+00  0.0000E+00  2.7392E-01  1.9635E-03
            -5.7275E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -875.994093267184        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1114
 NPARAMETR:  5.1926E-01  1.2490E-02  3.5898E-02  3.7275E-01  1.1847E+01  1.4776E+00  1.0000E-02  1.0000E-02  9.4072E-01  5.4264E-02
             8.8019E+00
 PARAMETER: -5.5535E-01 -4.2829E+00 -3.2271E+00 -8.8684E-01  2.5721E+00  4.9039E-01 -8.9098E+00 -2.1066E+01  3.8895E-02 -2.8139E+00
             2.2750E+00
 GRADIENT:  -5.5231E-01  1.7996E-02 -5.8244E-01  3.0254E-01  8.4213E-03 -1.1974E-01  0.0000E+00  0.0000E+00 -7.8374E-02  2.6256E-07
            -6.0222E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -875.996285992305        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1179
 NPARAMETR:  5.2100E-01  1.1870E-02  3.5377E-02  3.7248E-01  1.1843E+01  1.4789E+00  1.0000E-02  1.0000E-02  9.4041E-01  5.8120E-02
             8.7833E+00
 PARAMETER: -5.5467E-01 -4.3139E+00 -3.2276E+00 -8.8764E-01  2.5463E+00  4.9085E-01 -8.9098E+00 -2.1066E+01  3.8820E-02 -2.7542E+00
             2.2760E+00
 GRADIENT:  -3.7324E-01  1.4735E-03  2.1142E+00 -1.8171E-02 -3.5467E-03 -2.1068E-02  0.0000E+00  0.0000E+00  7.6471E-03 -4.7073E-07
             4.6095E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1179
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3243E-03  6.0603E-07  6.1172E-05 -1.9709E-02 -1.3300E-06
 SE:             2.8965E-02  1.2220E-06  2.2371E-04  2.4909E-02  6.7667E-06
 N:                     100         100         100         100         100

 P VAL.:         9.3604E-01  6.1995E-01  7.8451E-01  4.2882E-01  8.4418E-01

 ETASHRINKSD(%)  2.9644E+00  9.9996E+01  9.9251E+01  1.6550E+01  9.9977E+01
 ETASHRINKVR(%)  5.8410E+00  1.0000E+02  9.9994E+01  3.0361E+01  1.0000E+02
 EBVSHRINKSD(%)  2.8555E+00  9.9995E+01  9.9306E+01  1.6698E+01  9.9978E+01
 EBVSHRINKVR(%)  5.6294E+00  1.0000E+02  9.9995E+01  3.0608E+01  1.0000E+02
 RELATIVEINF(%)  4.3530E+00  2.1735E-08  4.7592E-05  6.4659E-01  4.3394E-07
 EPSSHRINKSD(%)  1.1920E+01
 EPSSHRINKVR(%)  2.2420E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -875.99628599230527     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       42.942247212367420     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -875.996       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.20E-01  1.21E-02  3.59E-02  3.72E-01  1.15E+01  1.48E+00  1.00E-02  1.00E-02  9.41E-01  5.76E-02  8.81E+00
 


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
+        1.76E+03
 
 TH 2
+       -5.30E+02  1.77E+03
 
 TH 3
+       -8.53E+03  3.20E+02  3.95E+05
 
 TH 4
+       -8.45E+01  2.35E+02 -4.51E+04  5.97E+03
 
 TH 5
+        3.63E-01 -9.68E-01 -4.51E+00  4.00E-01  3.60E-03
 
 TH 6
+       -2.49E-01  2.24E+01 -1.01E+02 -2.09E+01  8.66E-03  7.72E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.07E-01 -1.86E+01  7.70E+02 -1.03E+02 -2.13E-02 -2.29E+00  0.00E+00  0.00E+00  1.14E+02
 
 TH10
+       -3.81E-03 -2.94E-02 -3.91E-03  3.42E-03 -1.11E-05  1.39E-03  0.00E+00  0.00E+00  1.14E-02  7.16E-05
 
 TH11
+       -2.22E+01  5.99E+00  2.36E+02 -1.76E+01 -5.90E-03  1.00E+00  0.00E+00  0.00E+00  3.07E+00 -5.79E-05  6.43E+00
 
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
 #CPUT: Total CPU Time in Seconds,       30.230
Stop Time:
Thu Sep 30 03:00:44 CDT 2021
