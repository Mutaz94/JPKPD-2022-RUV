Sat Sep 18 13:36:05 CDT 2021
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
$DATA ../../../../data/spa/S2/dat70.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m70.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1687.84923959472        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.9018E+01 -5.8936E+00  8.6216E+00 -2.8119E+01 -3.7526E+01  1.5982E+00 -3.4102E+00  7.4468E+00  1.6357E+00  1.9441E+01
            -2.1245E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1691.93662684833        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0266E+00  1.0178E+00  1.0277E+00  1.0038E+00  1.0563E+00  9.9013E-01  1.0552E+00  9.3176E-01  9.7448E-01  8.6629E-01
             1.0831E+00
 PARAMETER:  1.2624E-01  1.1763E-01  1.2734E-01  1.0379E-01  1.5476E-01  9.0082E-02  1.5372E-01  2.9316E-02  7.4154E-02 -4.3537E-02
             1.7980E-01
 GRADIENT:   3.0715E+01 -8.0352E+00 -1.4230E+00 -1.9494E+00  2.6006E+01 -5.6159E-01 -4.9373E+00  1.0983E+00  5.8801E-01 -3.2477E+00
             5.8128E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1693.88160298489        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0200E+00  9.1067E-01  6.6108E-01  1.0544E+00  7.7099E-01  1.0071E+00  1.4188E+00  5.8422E-01  8.5799E-01  5.1783E-01
             1.0499E+00
 PARAMETER:  1.1985E-01  6.4211E-03 -3.1388E-01  1.5296E-01 -1.6007E-01  1.0706E-01  4.4982E-01 -4.3748E-01 -5.3164E-02 -5.5811E-01
             1.4868E-01
 GRADIENT:   9.6667E+00  1.8581E+01 -1.7621E+01  4.7650E+01  2.0877E+01  4.4386E+00  1.2103E+01  3.6510E+00  3.1780E+00 -2.0739E+00
            -2.9895E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1695.25225251580        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0169E+00  9.2763E-01  6.3731E-01  1.0207E+00  7.6855E-01  9.9907E-01  1.3049E+00  3.7273E-01  8.6075E-01  5.8798E-01
             1.0528E+00
 PARAMETER:  1.1673E-01  2.4879E-02 -3.5050E-01  1.2053E-01 -1.6325E-01  9.9066E-02  3.6612E-01 -8.8691E-01 -4.9952E-02 -4.3106E-01
             1.5144E-01
 GRADIENT:   2.8630E-01  1.5331E+00 -4.6249E+00  6.2418E+00  3.3617E+00  7.3903E-01  2.0471E+00  1.2190E+00  6.5432E-01  1.0530E+00
            -1.4712E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1695.31603108720        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0170E+00  9.3415E-01  6.0862E-01  1.0086E+00  7.5274E-01  9.9696E-01  1.2807E+00  2.4575E-01  8.6083E-01  5.7845E-01
             1.0538E+00
 PARAMETER:  1.1683E-01  3.1880E-02 -3.9656E-01  1.0858E-01 -1.8403E-01  9.6955E-02  3.4738E-01 -1.3034E+00 -4.9864E-02 -4.4741E-01
             1.5239E-01
 GRADIENT:  -4.3139E-01 -9.9208E-01 -5.0148E-01 -1.8070E+00 -1.5218E+00 -2.2590E-01 -2.6826E-01  4.5566E-01  1.9342E-01  9.9272E-01
             2.5447E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1695.37418817518        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0173E+00  9.4341E-01  5.8755E-01  1.0002E+00  7.4345E-01  9.9772E-01  1.2770E+00  9.9454E-02  8.6184E-01  5.6239E-01
             1.0546E+00
 PARAMETER:  1.1720E-01  4.1745E-02 -4.3179E-01  1.0016E-01 -1.9645E-01  9.7713E-02  3.4450E-01 -2.2081E+00 -4.8683E-02 -4.7556E-01
             1.5313E-01
 GRADIENT:  -8.5244E-03 -3.4788E-01 -5.3312E-01 -2.9077E-01  4.0561E-01 -2.9073E-03  1.1099E-01  6.2365E-02  5.8626E-02  1.4646E-01
             8.7882E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1696.05931543503        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  1.0389E+00  9.2765E-01  6.3901E-01  1.0221E+00  7.7172E-01  1.0059E+00  1.3020E+00  6.3446E-02  8.6342E-01  6.2343E-01
             1.0556E+00
 PARAMETER:  1.3814E-01  2.4905E-02 -3.4783E-01  1.2191E-01 -1.5914E-01  1.0587E-01  3.6391E-01 -2.6576E+00 -4.6857E-02 -3.7251E-01
             1.5414E-01
 GRADIENT:   1.5795E+00  9.4620E-02  1.7736E+00 -2.1419E+00 -2.0028E+00  1.9344E-02 -6.5737E-01  1.6396E-03  8.5336E-02 -4.0051E-01
            -1.4393E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1696.07091818736        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  1.0380E+00  9.0487E-01  6.4413E-01  1.0365E+00  7.6563E-01  1.0058E+00  1.3332E+00  5.5217E-02  8.5339E-01  6.2859E-01
             1.0555E+00
 PARAMETER:  1.3732E-01  3.1525E-05 -3.3985E-01  1.3587E-01 -1.6706E-01  1.0577E-01  3.8757E-01 -2.7965E+00 -5.8544E-02 -3.6427E-01
             1.5400E-01
 GRADIENT:   3.5362E-02 -1.3637E-02 -2.0905E-02 -5.3393E-03  4.0126E-02  1.0911E-02  2.3816E-02  7.4106E-04  7.3751E-03 -6.4521E-03
             5.6870E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1696.07121436515        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  1.0380E+00  9.0665E-01  6.4272E-01  1.0353E+00  7.6552E-01  1.0058E+00  1.3311E+00  1.9144E-02  8.5391E-01  6.2822E-01
             1.0555E+00
 PARAMETER:  1.3734E-01  1.9970E-03 -3.4205E-01  1.3469E-01 -1.6720E-01  1.0578E-01  3.8602E-01 -3.8558E+00 -5.7933E-02 -3.6486E-01
             1.5401E-01
 GRADIENT:   2.6479E-02 -1.9578E-02 -5.3099E-02  1.1421E-02  7.9300E-02  8.1006E-03  2.1133E-02  7.7800E-05  4.4470E-03 -2.8175E-03
             3.0269E-03

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1696.07124176330        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  1.0380E+00  9.0666E-01  6.4278E-01  1.0353E+00  7.6554E-01  1.0058E+00  1.3308E+00  1.0000E-02  8.5391E-01  6.2849E-01
             1.0555E+00
 PARAMETER:  1.3731E-01  2.0139E-03 -3.4195E-01  1.3470E-01 -1.6717E-01  1.0575E-01  3.8579E-01 -4.6778E+00 -5.7934E-02 -3.6444E-01
             1.5400E-01
 GRADIENT:  -1.8070E-02  8.2628E-05  1.2169E-02 -2.6643E-03 -1.8559E-02 -5.6185E-03 -1.2827E-02  0.0000E+00 -4.0922E-03  2.1863E-03
            -4.8979E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      979
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1581E-04  6.2617E-03 -5.1201E-04 -7.3031E-03 -4.8237E-03
 SE:             2.9839E-02  2.3103E-02  2.2140E-04  2.5020E-02  1.9900E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9690E-01  7.8637E-01  2.0745E-02  7.7037E-01  8.0847E-01

 ETASHRINKSD(%)  3.6742E-02  2.2603E+01  9.9258E+01  1.6181E+01  3.3332E+01
 ETASHRINKVR(%)  7.3471E-02  4.0096E+01  9.9994E+01  2.9744E+01  5.5554E+01
 EBVSHRINKSD(%)  4.6646E-01  2.1853E+01  9.9316E+01  1.6329E+01  3.3438E+01
 EBVSHRINKVR(%)  9.3075E-01  3.8931E+01  9.9995E+01  2.9992E+01  5.5695E+01
 RELATIVEINF(%)  9.8856E+01  5.6476E+00  4.2930E-04  8.2797E+00  2.8460E+00
 EPSSHRINKSD(%)  4.3007E+01
 EPSSHRINKVR(%)  6.7518E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1696.0712417633017     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -960.92041519956354     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.09
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1696.071       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  9.07E-01  6.43E-01  1.04E+00  7.66E-01  1.01E+00  1.33E+00  1.00E-02  8.54E-01  6.28E-01  1.06E+00
 


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
+        1.01E+03
 
 TH 2
+       -7.05E+00  4.81E+02
 
 TH 3
+        1.97E+01  3.92E+02  1.35E+03
 
 TH 4
+       -1.35E+01  2.77E+02 -5.42E+02  1.01E+03
 
 TH 5
+       -6.20E+00 -6.00E+02 -1.53E+03  5.71E+02  2.11E+03
 
 TH 6
+       -3.00E-01 -1.12E+00  4.13E+00 -3.67E+00 -3.96E-01  1.94E+02
 
 TH 7
+        9.16E-01  3.47E+01 -2.40E+01 -1.22E+01  4.78E+00  7.68E-03  4.51E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.91E+00 -2.38E+01 -4.49E+01  3.26E+01 -6.16E+00 -6.87E-01  1.40E+01  0.00E+00  1.40E+02
 
 TH10
+       -1.16E+00 -9.66E+00 -1.03E+02 -3.86E+01 -4.10E+01  2.80E-02  1.73E+01  0.00E+00  1.85E+01  1.11E+02
 
 TH11
+       -6.89E+00 -1.28E+01 -4.79E+01 -6.88E+00  7.79E+00  2.41E+00  6.14E+00  0.00E+00  1.52E+01  3.09E+01  1.92E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.800
Stop Time:
Sat Sep 18 13:36:23 CDT 2021
