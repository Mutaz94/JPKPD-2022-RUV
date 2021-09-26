Sat Sep 25 02:51:48 CDT 2021
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
$DATA ../../../../data/int/SL3/dat93.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      980
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      880
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
 RAW OUTPUT FILE (FILE): m93.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -438.962001087038        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.8883E+01 -9.8732E+01  3.0485E+02  9.5372E+01  1.3571E+02  1.0251E+01 -7.6768E+01 -2.2693E+02 -3.2170E+01 -1.4295E+01
            -6.3354E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2712.19914243009        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0421E+00  1.2424E+00  7.8747E-01  9.0539E-01  1.0028E+00  9.8206E-01  8.7622E-01  9.7164E-01  8.6259E-01  1.0028E+00
             2.7425E+00
 PARAMETER:  1.4125E-01  3.1703E-01 -1.3893E-01  6.1335E-04  1.0277E-01  8.1900E-02 -3.2132E-02  7.1232E-02 -4.7815E-02  1.0282E-01
             1.1089E+00
 GRADIENT:   9.7090E+01 -4.2729E-02  7.6981E+00  1.4817E+01 -9.7477E+00  3.7424E+00  1.0622E+01  6.7963E+00  5.2379E+00 -4.4458E+00
             4.3519E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2716.35631844626        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0304E+00  1.3598E+00  7.2753E-01  8.5079E-01  1.0365E+00  9.6805E-01  7.4001E-01  7.1132E-01  9.6223E-01  1.1169E+00
             2.7387E+00
 PARAMETER:  1.2997E-01  4.0736E-01 -2.1810E-01 -6.1588E-02  1.3585E-01  6.7533E-02 -2.0108E-01 -2.4063E-01  6.1503E-02  2.1056E-01
             1.1075E+00
 GRADIENT:   7.2572E+01  4.2633E+01  7.9308E+00  3.6364E+01 -3.6027E+01 -2.7617E-01  4.5289E+00  1.8564E+00  1.2722E+01  1.1412E+00
             4.2602E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2722.41223105789        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0035E+00  1.5615E+00  6.9481E-01  7.1873E-01  1.2086E+00  9.7367E-01  6.4685E-01  6.0491E-01  9.9592E-01  1.2319E+00
             2.6700E+00
 PARAMETER:  1.0350E-01  5.4564E-01 -2.6412E-01 -2.3027E-01  2.8946E-01  7.3321E-02 -3.3563E-01 -4.0267E-01  9.5911E-02  3.0859E-01
             1.0821E+00
 GRADIENT:   1.2634E+01  3.4352E+01  5.8548E+00  2.8841E+01 -3.5719E+00  3.2648E+00 -9.6069E-01 -5.7684E-01 -1.9608E-01  1.1284E+00
            -8.3192E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2726.30390269327        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.9910E-01  1.8941E+00  5.3619E-01  5.0109E-01  1.4368E+00  9.7115E-01  5.7106E-01  8.7061E-01  1.2621E+00  1.3820E+00
             2.6600E+00
 PARAMETER:  9.9103E-02  7.3877E-01 -5.2326E-01 -5.9097E-01  4.6241E-01  7.0726E-02 -4.6026E-01 -3.8560E-02  3.3278E-01  4.2352E-01
             1.0783E+00
 GRADIENT:   2.6237E+00  4.5491E+01  1.4330E+00  1.6695E+01  3.9917E-01  2.1822E+00 -3.5714E+00 -6.2331E-01 -1.5395E+00  3.7868E+00
            -1.6318E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2727.04682339767        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      422
 NPARAMETR:  1.0003E+00  1.9769E+00  4.9963E-01  4.4636E-01  1.5051E+00  9.6731E-01  5.8056E-01  1.1856E+00  1.3719E+00  1.3888E+00
             2.6567E+00
 PARAMETER:  1.0029E-01  7.8155E-01 -5.9389E-01 -7.0664E-01  5.0888E-01  6.6761E-02 -4.4377E-01  2.7028E-01  4.1622E-01  4.2844E-01
             1.0771E+00
 GRADIENT:  -6.1075E-01  2.0562E+01 -3.4476E-01  1.2838E+01  1.8281E+00  3.1998E-02  4.5502E-01  4.6060E-01  1.0075E+00 -4.6172E-01
            -1.7531E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2728.02352107151        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      597
 NPARAMETR:  1.0003E+00  2.1978E+00  3.3533E-01  2.9260E-01  1.6506E+00  9.6479E-01  5.5923E-01  9.8049E-01  1.8004E+00  1.4683E+00
             2.6507E+00
 PARAMETER:  1.0026E-01  8.8743E-01 -9.9264E-01 -1.1290E+00  6.0117E-01  6.4157E-02 -4.8120E-01  8.0296E-02  6.8800E-01  4.8412E-01
             1.0748E+00
 GRADIENT:  -1.3117E-01  1.4519E+01  7.8515E-01  1.5405E+00 -3.2278E+00 -9.1929E-01  2.1979E-01  2.7501E-01 -3.0661E-01 -1.0015E-01
            -2.1462E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2728.15992000348        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      772
 NPARAMETR:  1.0002E+00  2.2512E+00  2.7491E-01  2.5200E-01  1.6915E+00  9.6629E-01  5.5523E-01  4.4689E-01  1.9875E+00  1.4891E+00
             2.6503E+00
 PARAMETER:  1.0016E-01  9.1146E-01 -1.1913E+00 -1.2783E+00  6.2562E-01  6.5712E-02 -4.8837E-01 -7.0544E-01  7.8687E-01  4.9819E-01
             1.0747E+00
 GRADIENT:  -1.2180E-01 -2.8263E-01 -6.7653E-01  5.0119E-02 -3.2074E-01 -3.6647E-01  1.4076E-01  6.6199E-02 -4.0439E-01 -1.1770E-01
            -4.4382E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2728.18704433955        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      947
 NPARAMETR:  1.0002E+00  2.2558E+00  2.7279E-01  2.4835E-01  1.6957E+00  9.6728E-01  5.5399E-01  1.1538E-01  2.0282E+00  1.4891E+00
             2.6507E+00
 PARAMETER:  1.0019E-01  9.1350E-01 -1.1991E+00 -1.2929E+00  6.2812E-01  6.6729E-02 -4.9061E-01 -2.0596E+00  8.0713E-01  4.9819E-01
             1.0748E+00
 GRADIENT:  -4.3807E-02 -7.3214E-01 -4.4279E-02 -1.6173E-01 -1.0315E-01 -7.5608E-03 -2.6487E-02  2.6135E-03  4.4285E-04  9.6271E-03
             1.3907E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2728.18848617992        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1122
 NPARAMETR:  1.0002E+00  2.2539E+00  2.7414E-01  2.4987E-01  1.6943E+00  9.6733E-01  5.5429E-01  2.2202E-02  2.0215E+00  1.4882E+00
             2.6506E+00
 PARAMETER:  1.0021E-01  9.1268E-01 -1.1941E+00 -1.2868E+00  6.2730E-01  6.6780E-02 -4.9007E-01 -3.7076E+00  8.0382E-01  4.9757E-01
             1.0748E+00
 GRADIENT:  -6.2851E-03  3.2251E-02 -1.3217E-02  2.6793E-02  2.9770E-02  4.6673E-03  1.9751E-03  8.1712E-05  5.9923E-04 -5.6573E-04
             1.0628E-02

0ITERATION NO.:   49    OBJECTIVE VALUE:  -2728.18851848532        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1249
 NPARAMETR:  1.0002E+00  2.2542E+00  2.7403E-01  2.4966E-01  1.6945E+00  9.6731E-01  5.5426E-01  1.0000E-02  2.0226E+00  1.4883E+00
             2.6506E+00
 PARAMETER:  1.0021E-01  9.1279E-01 -1.1945E+00 -1.2877E+00  6.2740E-01  6.6762E-02 -4.9013E-01 -4.6580E+00  8.0439E-01  4.9762E-01
             1.0748E+00
 GRADIENT:   4.1777E-03 -1.2722E-02  2.9830E-03 -7.3358E-03 -7.1445E-03 -1.2196E-03 -5.4742E-04  0.0000E+00  5.0844E-04 -2.2647E-04
            -3.7242E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1249
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7723E-03 -2.7644E-02 -2.3648E-05  2.9246E-02 -2.0606E-02
 SE:             2.9411E-02  2.3678E-02  3.8981E-05  1.8343E-02  2.6073E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5195E-01  2.4302E-01  5.4407E-01  1.1085E-01  4.2933E-01

 ETASHRINKSD(%)  1.4689E+00  2.0674E+01  9.9869E+01  3.8549E+01  1.2652E+01
 ETASHRINKVR(%)  2.9161E+00  3.7075E+01  1.0000E+02  6.2237E+01  2.3703E+01
 EBVSHRINKSD(%)  1.6640E+00  1.9239E+01  9.9871E+01  4.3268E+01  1.1117E+01
 EBVSHRINKVR(%)  3.3003E+00  3.4776E+01  1.0000E+02  6.7815E+01  2.0998E+01
 RELATIVEINF(%)  9.6617E+01  9.0503E+00  6.5010E-05  4.2083E+00  3.6469E+01
 EPSSHRINKSD(%)  1.6260E+01
 EPSSHRINKVR(%)  2.9876E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          880
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1617.3318184402240     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2728.1885184853236     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1110.8567000450996     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2728.189       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.25E+00  2.74E-01  2.50E-01  1.69E+00  9.67E-01  5.54E-01  1.00E-02  2.02E+00  1.49E+00  2.65E+00
 


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
+        1.15E+03
 
 TH 2
+       -1.84E+01  4.81E+02
 
 TH 3
+        2.12E-01  9.29E+01  3.54E+02
 
 TH 4
+       -3.07E+01  5.23E+02 -3.08E+02  1.54E+03
 
 TH 5
+       -2.23E+00 -5.64E+01 -4.44E+01  1.17E+02  1.38E+02
 
 TH 6
+        5.28E+00 -5.07E+00  2.01E-01 -1.20E+01 -9.84E-01  1.96E+02
 
 TH 7
+        2.57E+00 -2.03E+01 -1.25E+01  1.26E+01 -6.68E+00  6.79E-02  2.82E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.12E-02 -6.10E+00 -1.45E+01  5.64E+01  1.67E-01 -4.70E-01  1.30E+01  0.00E+00  1.06E+01
 
 TH10
+       -7.18E-02 -4.52E+00  1.80E+01  8.33E+00 -8.51E+00  1.05E-01  5.23E+00  0.00E+00  2.51E+00  5.38E+01
 
 TH11
+       -1.51E+01 -2.05E+01 -6.53E+00 -1.58E+01  1.26E+00  2.73E+00  1.36E+01  0.00E+00  2.17E+00  4.92E+00  1.63E+02
 
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
 #CPUT: Total CPU Time in Seconds,       39.921
Stop Time:
Sat Sep 25 02:52:29 CDT 2021
