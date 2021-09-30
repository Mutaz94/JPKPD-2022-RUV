Wed Sep 29 18:48:14 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat15.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1645.45921197527        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2921E+02 -6.4879E+01 -1.2649E+01 -5.2745E+01  8.0113E+01  5.8460E+01 -1.2161E+01 -1.1012E+01 -2.1686E+01 -1.1397E+01
            -4.2345E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1659.55670332200        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0480E+00  1.0506E+00  9.6476E-01  1.0596E+00  9.5412E-01  9.2515E-01  1.0506E+00  1.0696E+00  1.0700E+00  1.0021E+00
             1.1001E+00
 PARAMETER:  1.4690E-01  1.4934E-01  6.4123E-02  1.5790E-01  5.3034E-02  2.2202E-02  1.4940E-01  1.6730E-01  1.6768E-01  1.0208E-01
             1.9540E-01
 GRADIENT:   4.3932E+00 -8.0148E+00 -9.7210E+00  9.5892E+00  1.7536E+01 -3.5793E+00 -2.6454E+00 -4.2245E+00  1.4182E+00  2.1350E+00
             1.5860E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1660.70076108138        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0505E+00  1.2477E+00  9.5694E-01  9.4841E-01  1.0110E+00  9.3860E-01  9.6959E-01  1.4437E+00  1.1465E+00  9.3793E-01
             1.0937E+00
 PARAMETER:  1.4925E-01  3.2133E-01  5.5982E-02  4.7035E-02  1.1091E-01  3.6634E-02  6.9118E-02  4.6722E-01  2.3670E-01  3.5919E-02
             1.8956E-01
 GRADIENT:   8.8159E+00  1.2997E+01 -1.9122E+00  1.8731E+01 -1.0758E+00  1.8737E+00 -5.6670E-01  9.2914E-01 -1.3745E+00 -3.3422E+00
            -1.7091E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1661.22259617895        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.0463E+00  1.3432E+00  1.0914E+00  8.7933E-01  1.1142E+00  9.3381E-01  8.3422E-01  1.7334E+00  1.2545E+00  1.0672E+00
             1.0973E+00
 PARAMETER:  1.4524E-01  3.9508E-01  1.8749E-01 -2.8593E-02  2.0818E-01  3.1513E-02 -8.1263E-02  6.5011E-01  3.2675E-01  1.6506E-01
             1.9282E-01
 GRADIENT:   2.1178E-01  4.4891E+00  1.9249E+00  3.5120E+00 -3.6225E+00  3.1506E-01 -8.5196E-01 -3.7431E-01 -8.2564E-01  1.8156E-01
            -4.9126E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1661.30344909913        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0487E+00  1.5288E+00  9.6326E-01  7.5779E-01  1.1528E+00  9.3594E-01  8.2229E-01  1.8284E+00  1.3684E+00  1.0795E+00
             1.0973E+00
 PARAMETER:  1.4755E-01  5.2448E-01  6.2568E-02 -1.7735E-01  2.4219E-01  3.3799E-02 -9.5665E-02  7.0345E-01  4.1366E-01  1.7646E-01
             1.9287E-01
 GRADIENT:   4.3367E+00  8.8364E+00  3.6007E+00  3.1335E+00 -5.9389E+00  9.0775E-01  8.9628E-01 -8.8596E-02 -9.8381E-01  8.5685E-02
            -1.0186E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1661.35383250602        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  1.0489E+00  1.7167E+00  7.8309E-01  6.3258E-01  1.1845E+00  9.3556E-01  7.9484E-01  1.8513E+00  1.5300E+00  1.0873E+00
             1.1004E+00
 PARAMETER:  1.4772E-01  6.4039E-01 -1.4451E-01 -3.5795E-01  2.6932E-01  3.3392E-02 -1.2961E-01  7.1587E-01  5.2527E-01  1.8373E-01
             1.9568E-01
 GRADIENT:   3.1524E+00  1.1054E+01  3.0659E+00  4.0318E+00 -4.5940E+00  4.3371E-01  1.3013E+00 -1.5462E-01 -2.0727E-01 -3.0865E-01
            -4.4510E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1661.42620414522        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1080
 NPARAMETR:  1.0489E+00  1.7913E+00  6.7053E-01  5.6983E-01  1.1975E+00  9.3528E-01  7.7692E-01  1.8213E+00  1.6177E+00  1.0888E+00
             1.1018E+00
 PARAMETER:  1.4776E-01  6.8293E-01 -2.9969E-01 -4.6243E-01  2.8022E-01  3.3085E-02 -1.5242E-01  6.9954E-01  5.8103E-01  1.8511E-01
             1.9693E-01
 GRADIENT:   2.7239E+00 -1.2482E+01 -1.2753E+00 -8.4460E-01  5.7298E+00  1.7011E-01  2.5557E-01  7.8357E-01  9.9447E-02 -1.5638E-01
             3.1973E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1661.45856389305        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1263
 NPARAMETR:  1.0491E+00  1.7935E+00  6.6386E-01  5.7138E-01  1.1885E+00  9.3505E-01  7.7782E-01  1.7739E+00  1.6207E+00  1.0885E+00
             1.1011E+00
 PARAMETER:  1.4789E-01  6.8416E-01 -3.0968E-01 -4.5970E-01  2.7269E-01  3.2844E-02 -1.5126E-01  6.7319E-01  5.8287E-01  1.8477E-01
             1.9633E-01
 GRADIENT:   2.9276E+00 -4.7269E+00  1.5413E-01  6.9778E-01  5.8341E-02  4.2232E-02  1.3265E-01  3.3229E-01  6.8898E-01  7.0061E-01
             1.4471E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1661.46587166778        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1443
 NPARAMETR:  1.0491E+00  1.7931E+00  6.5991E-01  5.7082E-01  1.1839E+00  9.3512E-01  7.7903E-01  1.7377E+00  1.6159E+00  1.0827E+00
             1.1010E+00
 PARAMETER:  1.4793E-01  6.8396E-01 -3.1565E-01 -4.6069E-01  2.6883E-01  3.2917E-02 -1.4971E-01  6.5258E-01  5.7987E-01  1.7945E-01
             1.9620E-01
 GRADIENT:   2.9940E+00 -3.3295E+00  1.2866E+00 -5.4075E-01 -1.8140E+00  6.4229E-02 -1.3236E-01 -1.4822E-01  1.6199E-01  3.4917E-01
            -3.9166E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1661.47493836708        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1625             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0491E+00  1.7914E+00  6.4991E-01  5.7027E-01  1.1834E+00  9.3517E-01  7.8113E-01  1.7336E+00  1.6142E+00  1.0791E+00
             1.1010E+00
 PARAMETER:  1.4797E-01  6.8300E-01 -3.3092E-01 -4.6165E-01  2.6838E-01  3.2977E-02 -1.4702E-01  6.5021E-01  5.7886E-01  1.7609E-01
             1.9626E-01
 GRADIENT:   5.3569E+02  5.8424E+02  1.4636E+00  7.8118E+01  1.5737E+01  2.7769E+01  7.6179E+00  1.3558E+00  2.1111E+01  1.3788E+00
             1.7474E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1661.47867738891        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1795
 NPARAMETR:  1.0491E+00  1.7915E+00  6.4725E-01  5.7033E-01  1.1809E+00  9.3519E-01  7.8198E-01  1.7169E+00  1.6123E+00  1.0766E+00
             1.1008E+00
 PARAMETER:  1.4742E-01  6.8386E-01 -3.3171E-01 -4.6119E-01  2.6683E-01  3.2868E-02 -1.4704E-01  6.4181E-01  5.7769E-01  1.7445E-01
             1.9608E-01
 GRADIENT:  -1.4316E+00  1.9136E+00  2.3337E-01  1.8852E-01  7.0348E-01 -5.7218E-02 -8.4388E-02  4.1021E-02  1.0967E-02  1.0758E-01
            -2.4071E-03
 NUMSIGDIG:         2.4         2.9         1.9         3.1         2.6         2.8         2.1         2.6         4.0         2.4
                    4.5

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1795
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3075E-03 -3.9355E-02 -3.2475E-02  3.0547E-02 -5.0770E-02
 SE:             2.9820E-02  2.1972E-02  1.2394E-02  2.2900E-02  2.0894E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6503E-01  7.3274E-02  8.7881E-03  1.8221E-01  1.5104E-02

 ETASHRINKSD(%)  9.7576E-02  2.6390E+01  5.8478E+01  2.3284E+01  3.0003E+01
 ETASHRINKVR(%)  1.9506E-01  4.5815E+01  8.2759E+01  4.1146E+01  5.1004E+01
 EBVSHRINKSD(%)  6.1091E-01  2.5809E+01  6.0907E+01  2.4692E+01  2.7203E+01
 EBVSHRINKVR(%)  1.2181E+00  4.4957E+01  8.4717E+01  4.3287E+01  4.7005E+01
 RELATIVEINF(%)  9.8718E+01  3.1563E+00  1.4805E+00  3.4206E+00  1.6423E+01
 EPSSHRINKSD(%)  4.5463E+01
 EPSSHRINKVR(%)  7.0257E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1661.4786773889130     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -926.32785082517478     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1661.479       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.79E+00  6.49E-01  5.71E-01  1.18E+00  9.35E-01  7.81E-01  1.72E+00  1.61E+00  1.08E+00  1.10E+00
 


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
+        1.14E+03
 
 TH 2
+       -7.63E+00  3.36E+02
 
 TH 3
+        5.65E+00  8.80E+01  2.23E+06
 
 TH 4
+       -1.49E+01  3.09E+02  1.83E+06  7.45E+02
 
 TH 5
+       -2.12E+06 -1.26E+02 -1.58E+02  1.25E+02  3.79E+02
 
 TH 6
+        1.65E+00 -1.39E+00  1.85E+00 -3.99E+00 -9.53E-01  2.22E+02
 
 TH 7
+        6.66E-01 -6.65E+00  2.04E+00 -1.10E+01 -2.86E+06 -6.81E-01  1.10E+02
 
 TH 8
+        4.23E-01 -9.50E+00 -4.31E+05 -3.53E+05  1.09E-02  4.10E-02  5.59E+00  4.69E+00
 
 TH 9
+        1.50E+00 -1.42E+01  5.16E+05  4.23E+05 -4.67E+00 -2.54E-01  1.65E+01 -9.99E+04  1.20E+05
 
 TH10
+        3.36E-01 -6.82E+00 -1.23E+01 -2.53E+00  1.75E+06  3.23E-01  1.07E+01  3.84E+00  1.33E+00  6.03E+01
 
 TH11
+       -8.81E+00 -1.29E+01 -9.80E+00 -8.43E-01 -1.52E+06  3.44E+00  8.15E+00  1.85E+00  2.83E+00  1.18E+01  1.68E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.240
Stop Time:
Wed Sep 29 18:48:48 CDT 2021
