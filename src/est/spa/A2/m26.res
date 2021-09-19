Sat Sep 18 09:45:10 CDT 2021
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
$DATA ../../../../data/spa/A2/dat26.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1131.61155916887        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.8703E+01 -1.8039E+01  5.5800E+01 -7.2765E+01  1.9097E+01  2.7333E+01 -2.4827E+01 -2.8985E+01 -1.7592E+01 -5.4590E+01
            -9.2176E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1426.04900766258        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0071E+00  1.0477E+00  9.6577E-01  1.0536E+00  9.9346E-01  8.6364E-01  1.0338E+00  1.0819E+00  9.4012E-01  9.8200E-01
             2.3927E+00
 PARAMETER:  1.0708E-01  1.4659E-01  6.5167E-02  1.5223E-01  9.3436E-02 -4.6597E-02  1.3324E-01  1.7869E-01  3.8253E-02  8.1838E-02
             9.7244E-01
 GRADIENT:   1.1262E+01  1.1770E+01 -3.8748E-01  1.8405E+01 -5.1631E+00 -1.7285E+01 -1.1114E+00 -4.1229E+00  1.3647E+00  6.3476E+00
             1.4716E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1427.64013569385        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0063E+00  1.0452E+00  1.1067E+00  1.0539E+00  1.0555E+00  8.9247E-01  1.0455E+00  1.4358E+00  9.2388E-01  8.2492E-01
             2.4377E+00
 PARAMETER:  1.0626E-01  1.4420E-01  2.0142E-01  1.5247E-01  1.5406E-01 -1.3757E-02  1.4450E-01  4.6175E-01  2.0829E-02 -9.2475E-02
             9.9107E-01
 GRADIENT:   1.1743E+01  5.5539E+00 -1.0218E+00  9.5432E+00 -1.4106E+00 -4.1092E+00 -1.2438E-01 -5.6625E-01  1.9327E+00 -7.2133E-01
             1.2628E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1430.32941201860        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.9318E-01  8.2472E-01  2.2157E+00  1.2126E+00  1.3737E+00  9.0134E-01  1.1674E+00  2.3943E+00  8.5683E-01  1.3093E+00
             2.2926E+00
 PARAMETER:  9.3161E-02 -9.2715E-02  8.9558E-01  2.9275E-01  4.1748E-01 -3.8751E-03  2.5476E-01  9.7309E-01 -5.4513E-02  3.6947E-01
             9.2967E-01
 GRADIENT:  -7.3951E+00 -3.8565E+00 -7.3330E+00 -2.7769E+00  1.2257E+01  3.0120E-01  6.0830E-02  2.0096E+00  5.4961E-01 -4.4409E-01
            -6.8238E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1431.51828319870        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.9629E-01  6.8902E-01  2.8048E+00  1.3115E+00  1.3712E+00  8.9730E-01  1.2715E+00  2.7655E+00  8.0243E-01  1.2727E+00
             2.3359E+00
 PARAMETER:  9.6284E-02 -2.7249E-01  1.1313E+00  3.7116E-01  4.1567E-01 -8.3696E-03  3.4021E-01  1.1172E+00 -1.2011E-01  3.4116E-01
             9.4838E-01
 GRADIENT:   4.6854E-01  1.9146E+00 -8.9738E-01  5.0245E+00  4.3588E-01  1.8087E-01 -2.1094E-01  6.8459E-02 -4.7535E-01  2.7477E-01
             6.1368E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1431.90958866213        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      416
 NPARAMETR:  9.9781E-01  6.0803E-01  3.5118E+00  1.3764E+00  1.4170E+00  8.9562E-01  1.3726E+00  3.2470E+00  7.7315E-01  1.2891E+00
             2.3407E+00
 PARAMETER:  9.7810E-02 -3.9753E-01  1.3561E+00  4.1947E-01  4.4856E-01 -1.0244E-02  4.1669E-01  1.2777E+00 -1.5729E-01  3.5397E-01
             9.5046E-01
 GRADIENT:  -3.6159E+00  3.5016E+00 -8.4888E-02  5.7479E+00 -4.4548E-01 -2.4997E-01 -3.2388E-01 -8.6472E-02 -1.0689E+00  4.8506E-02
            -5.1347E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1432.16023007532        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      591
 NPARAMETR:  9.9813E-01  4.7619E-01  4.5396E+00  1.4699E+00  1.4451E+00  8.9220E-01  1.6132E+00  3.6805E+00  7.4240E-01  1.2950E+00
             2.3489E+00
 PARAMETER:  9.8131E-02 -6.4195E-01  1.6128E+00  4.8523E-01  4.6819E-01 -1.4064E-02  5.7822E-01  1.4030E+00 -1.9787E-01  3.5850E-01
             9.5397E-01
 GRADIENT:  -2.5177E+00  3.2576E+00  2.3719E+00  5.5295E+00 -1.0705E+00 -2.5126E-01  5.0011E-01 -2.1134E+00 -2.0247E-01  3.9686E-01
             1.4876E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1433.01369863623        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      768
 NPARAMETR:  9.9618E-01  2.2071E-01  4.2373E+00  1.6336E+00  1.3866E+00  8.8690E-01  1.9481E+00  3.4968E+00  7.1566E-01  1.2828E+00
             2.3300E+00
 PARAMETER:  9.6173E-02 -1.4109E+00  1.5439E+00  5.9077E-01  4.2682E-01 -2.0028E-02  7.6687E-01  1.3518E+00 -2.3454E-01  3.4908E-01
             9.4587E-01
 GRADIENT:   2.2259E-01  1.1453E+00  6.4328E-01  8.4875E+00  3.5701E-01  4.6110E-02 -1.1884E-01 -1.1221E+00 -1.2320E+00  2.2449E-01
             3.7137E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1433.02387741732        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  9.9575E-01  2.1874E-01  4.2118E+00  1.6346E+00  1.3834E+00  8.8650E-01  2.0134E+00  3.4861E+00  7.1905E-01  1.2819E+00
             2.3273E+00
 PARAMETER:  9.5744E-02 -1.4199E+00  1.5379E+00  5.9140E-01  4.2451E-01 -2.0470E-02  7.9980E-01  1.3488E+00 -2.2983E-01  3.4833E-01
             9.4470E-01
 GRADIENT:   6.6747E+00  1.6140E+00  9.0266E-01  2.5271E+01  2.9597E-01  3.7072E-01  9.0817E-02  3.4455E-01  5.9220E-01  1.7558E-01
             1.0337E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1433.02591281088        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1115            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9597E-01  2.1843E-01  4.2100E+00  1.6346E+00  1.3834E+00  8.8683E-01  2.0018E+00  3.4860E+00  7.1862E-01  1.2818E+00
             2.3265E+00
 PARAMETER:  9.5957E-02 -1.4213E+00  1.5375E+00  5.9141E-01  4.2455E-01 -2.0099E-02  7.9403E-01  1.3488E+00 -2.3042E-01  3.4829E-01
             9.4435E-01
 GRADIENT:   7.3043E+00  1.5763E+00  8.7080E-01  2.4941E+01  3.8891E-01  5.0444E-01  6.2696E-02  3.7251E-01  3.6630E-01  1.4148E-01
             7.8039E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1433.02591281088        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1179
 NPARAMETR:  9.9596E-01  2.1845E-01  4.2105E+00  1.6347E+00  1.3834E+00  8.8683E-01  2.0018E+00  3.4856E+00  7.1861E-01  1.2814E+00
             2.3264E+00
 PARAMETER:  9.5957E-02 -1.4213E+00  1.5375E+00  5.9141E-01  4.2455E-01 -2.0099E-02  7.9403E-01  1.3488E+00 -2.3042E-01  3.4829E-01
             9.4435E-01
 GRADIENT:   1.2348E-02 -8.0914E+03 -7.4782E+03 -1.9446E+04  7.3719E-02  1.0629E-02 -3.1106E-04  8.5081E+03  1.0158E-02  8.2684E-02
             2.5209E-02
 NUMSIGDIG:         4.3         3.3         3.3         3.3         3.4         3.4         3.9         3.3         3.7         2.2
                    4.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1179
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.2359E-04 -7.4163E-03 -2.5067E-02 -1.5458E-02 -5.0766E-02
 SE:             2.9153E-02  6.5892E-03  1.3973E-02  2.6174E-02  1.7184E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9388E-01  2.6037E-01  7.2810E-02  5.5478E-01  3.1347E-03

 ETASHRINKSD(%)  2.3330E+00  7.7925E+01  5.3189E+01  1.2315E+01  4.2431E+01
 ETASHRINKVR(%)  4.6115E+00  9.5127E+01  7.8088E+01  2.3114E+01  6.6858E+01
 EBVSHRINKSD(%)  2.3696E+00  7.7963E+01  6.0897E+01  1.2690E+01  4.3847E+01
 EBVSHRINKVR(%)  4.6830E+00  9.5144E+01  8.4709E+01  2.3770E+01  6.8469E+01
 RELATIVEINF(%)  9.3036E+01  2.9623E-01  7.2144E+00  5.0719E+00  1.2049E+01
 EPSSHRINKSD(%)  3.4096E+01
 EPSSHRINKVR(%)  5.6567E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1433.0259128108812     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -697.87508624714303     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.94
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1433.026       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  2.18E-01  4.21E+00  1.63E+00  1.38E+00  8.87E-01  2.00E+00  3.49E+00  7.19E-01  1.28E+00  2.33E+00
 


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
+        1.36E+03
 
 TH 2
+        9.30E+07  2.98E+07
 
 TH 3
+       -7.25E+01 -6.97E+01  6.86E+04
 
 TH 4
+       -5.34E+02 -5.83E+01 -2.96E+02  3.08E+06
 
 TH 5
+       -2.58E+00 -1.58E+07 -6.02E+02 -4.06E+03  1.23E+02
 
 TH 6
+       -1.02E+01 -2.95E+03 -1.40E+02 -9.52E+02 -3.02E+00  2.16E+02
 
 TH 7
+       -1.65E-01 -5.83E+06 -1.51E+00 -1.12E+01  1.61E-01 -1.95E-01  3.47E-01
 
 TH 8
+        9.87E+01  9.61E+01  5.81E+01 -6.31E+05  8.18E+02  1.94E+02  2.12E+00  1.30E+05
 
 TH 9
+        6.46E+00  5.59E+07 -1.84E+01 -1.30E+02 -6.39E-01  9.17E+00  4.76E+00  2.67E+01  2.26E+02
 
 TH10
+        1.49E+00  2.08E+07  1.07E+03  7.14E+03 -2.07E+01  3.45E-01  1.10E-02 -1.47E+03 -3.81E+00  1.44E+07
 
 TH11
+       -1.59E+01  4.21E+06 -4.19E+02 -2.83E+03 -5.16E+00  4.43E+00  4.89E-01  5.77E+02  1.89E+01  1.35E+01  4.98E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.868
Stop Time:
Sat Sep 18 09:45:33 CDT 2021
