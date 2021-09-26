Sat Sep 25 10:40:21 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat69.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1664.34358586572        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.0447E+01  2.8036E+00 -8.2065E+01  1.2252E+02  1.2400E+02  3.4049E+00 -1.1900E+01  9.4948E+00 -1.4825E+00 -4.6268E+00
            -3.9427E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1677.92881020609        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0470E+00  1.0982E+00  1.1972E+00  8.8639E-01  1.0394E+00  9.8825E-01  1.1252E+00  9.1053E-01  9.8469E-01  9.3201E-01
             1.1063E+00
 PARAMETER:  1.4591E-01  1.9364E-01  2.7995E-01 -2.0604E-02  1.3860E-01  8.8181E-02  2.1792E-01  6.2727E-03  8.4568E-02  2.9593E-02
             2.0099E-01
 GRADIENT:   1.0390E+02  4.5120E+00  1.2677E+00  4.2473E+00  2.9645E+01 -9.5769E-01  2.2250E+00 -1.2733E+00 -1.4468E+00 -1.3724E+01
             3.0679E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1679.24622811363        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0366E+00  9.9202E-01  9.9963E-01  9.2894E-01  9.1729E-01  9.9607E-01  1.2774E+00  5.2220E-01  9.5795E-01  8.7325E-01
             1.0943E+00
 PARAMETER:  1.3593E-01  9.1993E-02  9.9630E-02  2.6294E-02  1.3665E-02  9.6059E-02  3.4484E-01 -5.4970E-01  5.7040E-02 -3.5535E-02
             1.9009E-01
 GRADIENT:   7.7123E+01 -3.7155E+00  1.1673E+01 -1.7946E+01 -7.3322E+00  2.6673E+00  7.0443E+00 -5.2773E-02  4.7640E+00 -5.7772E+00
             1.7131E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1679.77203172555        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0151E+00  1.0113E+00  8.9076E-01  9.1486E-01  8.8401E-01  9.8982E-01  1.2249E+00  3.8466E-01  9.5462E-01  8.7162E-01
             1.0837E+00
 PARAMETER:  1.1500E-01  1.1123E-01 -1.5675E-02  1.1011E-02 -2.3285E-02  8.9765E-02  3.0287E-01 -8.5540E-01  5.3557E-02 -3.7405E-02
             1.8037E-01
 GRADIENT:   2.0158E+01 -4.4362E+00  1.3973E+00 -9.0990E+00 -4.6504E+00 -1.0274E-02  2.3986E+00  8.0610E-01  3.4797E+00  4.4886E-01
             8.5798E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1680.37529580422        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  1.0312E+00  8.6379E-01  8.7196E-01  1.0129E+00  8.1808E-01  1.0022E+00  1.4062E+00  2.1511E-01  8.6778E-01  8.3567E-01
             1.0828E+00
 PARAMETER:  1.3074E-01 -4.6430E-02 -3.7009E-02  1.1279E-01 -1.0079E-01  1.0218E-01  4.4086E-01 -1.4366E+00 -4.1813E-02 -7.9519E-02
             1.7953E-01
 GRADIENT:   6.7251E+00  8.5679E-01 -2.7899E+00  3.7684E+00  1.7329E+00  8.0477E-01 -1.3542E-01  2.9007E-01  4.4296E-03  7.7129E-01
            -2.6733E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1680.47184507805        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  1.0274E+00  7.9403E-01  8.6979E-01  1.0501E+00  7.9140E-01  1.0000E+00  1.5091E+00  1.0935E-01  8.4291E-01  8.2074E-01
             1.0836E+00
 PARAMETER:  1.2707E-01 -1.3063E-01 -3.9500E-02  1.4888E-01 -1.3395E-01  1.0001E-01  5.1151E-01 -2.1132E+00 -7.0894E-02 -9.7555E-02
             1.8028E-01
 GRADIENT:  -5.1785E-01 -7.0784E-01 -1.7763E-01 -6.6268E-01 -1.5344E-01  2.0037E-02 -7.9213E-03  5.9718E-02  7.1960E-02  2.2497E-01
             2.5675E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1680.51191337665        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      760
 NPARAMETR:  1.0280E+00  8.4729E-01  8.5925E-01  1.0191E+00  8.0623E-01  1.0001E+00  1.4332E+00  1.9021E-02  8.6191E-01  8.2522E-01
             1.0838E+00
 PARAMETER:  1.2766E-01 -6.5718E-02 -5.1697E-02  1.1890E-01 -1.1538E-01  1.0015E-01  4.5994E-01 -3.8622E+00 -4.8602E-02 -9.2111E-02
             1.8052E-01
 GRADIENT:  -6.6117E-02 -7.9453E-04 -2.9714E-01  1.0095E-01  1.8526E-01 -4.7508E-03  3.9605E-02  1.8499E-03  2.7023E-02  8.6182E-02
             5.1273E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1680.51277433151        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      936
 NPARAMETR:  1.0280E+00  8.4325E-01  8.6085E-01  1.0215E+00  8.0541E-01  1.0001E+00  1.4384E+00  1.0000E-02  8.6044E-01  8.2509E-01
             1.0838E+00
 PARAMETER:  1.2765E-01 -7.0488E-02 -4.9837E-02  1.2128E-01 -1.1640E-01  1.0014E-01  4.6353E-01 -4.5347E+00 -5.0311E-02 -9.2261E-02
             1.8049E-01
 GRADIENT:   2.5299E-04 -5.3975E-04  1.0625E-03 -2.4862E-03 -1.0983E-03 -1.6321E-04 -1.0620E-05  0.0000E+00  2.1293E-04 -3.4174E-04
            -3.3385E-04

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1680.51277433151        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      958
 NPARAMETR:  1.0280E+00  8.4325E-01  8.6085E-01  1.0215E+00  8.0541E-01  1.0001E+00  1.4384E+00  1.0000E-02  8.6044E-01  8.2509E-01
             1.0838E+00
 PARAMETER:  1.2765E-01 -7.0488E-02 -4.9837E-02  1.2128E-01 -1.1640E-01  1.0014E-01  4.6353E-01 -4.5347E+00 -5.0311E-02 -9.2261E-02
             1.8049E-01
 GRADIENT:   2.5299E-04 -5.3975E-04  1.0625E-03 -2.4862E-03 -1.0983E-03 -1.6321E-04 -1.0620E-05  0.0000E+00  2.1293E-04 -3.4174E-04
            -3.3385E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      958
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.9487E-04  6.4127E-03 -3.5402E-04 -1.0249E-02 -9.3080E-03
 SE:             2.9809E-02  2.1602E-02  1.6654E-04  2.3419E-02  2.3061E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9478E-01  7.6658E-01  3.3525E-02  6.6165E-01  6.8648E-01

 ETASHRINKSD(%)  1.3569E-01  2.7629E+01  9.9442E+01  2.1543E+01  2.2744E+01
 ETASHRINKVR(%)  2.7119E-01  4.7625E+01  9.9997E+01  3.8445E+01  4.0315E+01
 EBVSHRINKSD(%)  5.0698E-01  2.7958E+01  9.9470E+01  2.1279E+01  2.1454E+01
 EBVSHRINKVR(%)  1.0114E+00  4.8100E+01  9.9997E+01  3.8029E+01  3.8306E+01
 RELATIVEINF(%)  9.8527E+01  3.8693E+00  3.2711E-04  5.1626E+00  6.0766E+00
 EPSSHRINKSD(%)  4.2777E+01
 EPSSHRINKVR(%)  6.7256E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1680.5127743315063     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -945.36194776776813     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1680.513       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  8.43E-01  8.61E-01  1.02E+00  8.05E-01  1.00E+00  1.44E+00  1.00E-02  8.60E-01  8.25E-01  1.08E+00
 


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
+        1.04E+03
 
 TH 2
+       -8.43E+00  4.10E+02
 
 TH 3
+        1.65E+01  1.57E+02  4.79E+02
 
 TH 4
+       -1.33E+01  3.65E+02 -2.47E+02  9.08E+02
 
 TH 5
+       -8.33E+00 -3.67E+02 -7.43E+02  3.53E+02  1.55E+03
 
 TH 6
+        5.34E-02 -8.83E-01  3.13E+00 -3.49E+00 -2.73E+00  1.96E+02
 
 TH 7
+        1.17E+00  3.33E+01 -2.20E+00 -1.35E+01  4.03E-01 -3.82E-02  3.40E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.38E+00 -2.24E+01 -2.05E+01  2.27E+01  6.75E+00 -8.02E-01  1.13E+01  0.00E+00  1.28E+02
 
 TH10
+       -5.60E-01 -8.15E+00 -5.25E+01 -2.99E+01 -5.17E+01  1.66E+00  1.27E+01  0.00E+00  7.02E+00  1.15E+02
 
 TH11
+       -7.26E+00 -1.10E+01 -2.77E+01 -4.98E+00  1.05E+01  1.97E+00  2.79E+00  0.00E+00  9.11E+00  2.39E+01  1.91E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.627
Stop Time:
Sat Sep 25 10:40:39 CDT 2021
