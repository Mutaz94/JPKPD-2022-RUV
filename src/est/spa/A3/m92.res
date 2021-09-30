Wed Sep 29 13:54:32 CDT 2021
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
$DATA ../../../../data/spa/A3/dat92.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   153.555992874165        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4965E+02  1.4146E+01  7.5941E+01 -7.8498E+01  1.5637E+02  4.6874E+01 -3.9642E+01 -4.4624E+01 -1.1968E+02 -1.1286E+02
            -3.2365E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -983.710126576662        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0019E+00  9.2430E-01  8.5970E-01  1.0894E+00  8.1376E-01  9.3125E-01  9.8466E-01  1.0088E+00  1.3181E+00  1.1348E+00
             1.8257E+00
 PARAMETER:  1.0190E-01  2.1283E-02 -5.1175E-02  1.8567E-01 -1.0609E-01  2.8772E-02  8.4542E-02  1.0874E-01  3.7620E-01  2.2646E-01
             7.0194E-01
 GRADIENT:   1.9178E+02  1.9851E+01  3.1439E+01  5.4881E-01 -2.6686E+00  2.5617E+00  2.5384E+00 -2.6186E+00  7.1815E+00 -6.0799E+00
            -8.1647E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1029.34893328705        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      200
 NPARAMETR:  1.0071E+00  6.8244E-01  4.1460E-01  1.2292E+00  5.0796E-01  9.3157E-01  2.0971E-01  3.4833E-01  1.4531E+00  1.0939E+00
             1.9170E+00
 PARAMETER:  1.0705E-01 -2.8208E-01 -7.8045E-01  3.0640E-01 -5.7735E-01  2.9113E-02 -1.4620E+00 -9.5461E-01  4.7373E-01  1.8972E-01
             7.5074E-01
 GRADIENT:   7.8595E+01 -3.0169E+01 -7.2694E+01  5.3653E+01  1.3953E+02 -9.0315E+00  6.7980E-02  9.5211E-01  4.0188E+01  1.2236E+01
            -6.8343E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1219.21333520546        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      382
 NPARAMETR:  9.7729E-01  6.3658E-01  4.7279E-01  1.1402E+00  5.0507E-01  8.4752E-01  9.1047E-02  6.4563E-01  1.0605E+00  7.9266E-01
             3.4347E+00
 PARAMETER:  7.7033E-02 -3.5164E-01 -6.4911E-01  2.3116E-01 -5.8305E-01 -6.5442E-02 -2.2964E+00 -3.3753E-01  1.5871E-01 -1.3236E-01
             1.3339E+00
 GRADIENT:  -3.3464E+01 -4.9464E+01 -3.0928E+01 -1.1156E+02  6.2791E+01 -1.9358E+01 -4.2730E-05  6.3378E+00 -1.0816E+01  1.4758E+01
            -5.0616E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1236.32116344506        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      557
 NPARAMETR:  1.0055E+00  4.0860E-01  4.4009E-01  1.3400E+00  4.0343E-01  9.1121E-01  8.8570E-02  3.3218E-01  9.6592E-01  4.8138E-01
             3.8066E+00
 PARAMETER:  1.0550E-01 -7.9502E-01 -7.2078E-01  3.9268E-01 -8.0776E-01  7.0213E-03 -2.3240E+00 -1.0021E+00  6.5323E-02 -6.3110E-01
             1.4367E+00
 GRADIENT:   1.1864E+01  1.3413E+01  2.1033E+01  3.4160E-01 -3.5177E+01  4.4957E+00 -1.8858E-02  1.1239E+00 -9.8698E-01  7.1617E-02
            -2.0199E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1243.36342127117        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      733
 NPARAMETR:  9.8212E-01  1.3790E-01  8.6879E-01  1.5967E+00  6.3386E-01  8.6129E-01  5.4339E-01  2.1489E-02  8.2104E-01  7.4438E-02
             3.9604E+00
 PARAMETER:  8.1958E-02 -1.8812E+00 -4.0657E-02  5.6795E-01 -3.5593E-01 -4.9326E-02 -5.0992E-01 -3.7402E+00 -9.7181E-02 -2.4978E+00
             1.4764E+00
 GRADIENT:   4.5790E+00 -8.8609E-01 -6.6315E+00 -9.0211E+00  9.6210E+00  2.0603E-01  1.7727E-02  5.3231E-03  2.2053E+00  5.2632E-02
            -2.9072E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1243.49821372568        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  9.8099E-01  1.2480E-01  9.3388E-01  1.6165E+00  6.5794E-01  8.5963E-01  1.4602E-01  1.0000E-02  7.9812E-01  2.0334E-02
             3.9843E+00
 PARAMETER:  8.0806E-02 -1.9810E+00  3.1591E-02  5.8023E-01 -3.1864E-01 -5.1251E-02 -1.8240E+00 -5.0655E+00 -1.2549E-01 -3.7955E+00
             1.4824E+00
 GRADIENT:   2.3338E+01  4.3969E-01  9.7094E-01  4.3768E+01  1.5418E+00  2.0715E+00  2.0682E-03  0.0000E+00  7.8907E-01  4.8500E-03
             1.1748E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1243.49896744913        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      963
 NPARAMETR:  9.8071E-01  1.2514E-01  9.3360E-01  1.6156E+00  6.5777E-01  8.5943E-01  6.9287E-02  1.0000E-02  7.9813E-01  1.1495E-02
             3.9829E+00
 PARAMETER:  8.0526E-02 -1.9784E+00  3.1288E-02  5.7970E-01 -3.1890E-01 -5.1491E-02 -2.5695E+00 -5.7449E+00 -1.2549E-01 -4.3658E+00
             1.4820E+00
 GRADIENT:   2.2779E+01  4.2235E-01  1.0678E+00  4.3110E+01  1.4696E+00  2.0164E+00  5.7398E-04  0.0000E+00  7.2927E-01  1.5948E-03
             1.1466E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1243.50084379410        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:     1084
 NPARAMETR:  9.8028E-01  1.2581E-01  9.3257E-01  1.6163E+00  6.5807E-01  8.5874E-01  7.6470E-02  1.0000E-02  7.9786E-01  1.2326E-02
             3.9821E+00
 PARAMETER:  8.0087E-02 -1.9730E+00  3.0193E-02  5.8011E-01 -3.1845E-01 -5.2286E-02 -2.4709E+00 -5.5887E+00 -1.2583E-01 -4.2961E+00
             1.4818E+00
 GRADIENT:  -1.0452E+00 -6.7627E-02  1.8572E-01 -2.2568E+00  4.8513E-03 -6.8698E-02  2.6859E-04  0.0000E+00 -1.1308E-01  1.4216E-03
            -5.9764E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1243.50311268553        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1267
 NPARAMETR:  9.8094E-01  1.3012E-01  9.2914E-01  1.6147E+00  6.5753E-01  8.5917E-01  3.2483E-02  1.0000E-02  8.0007E-01  1.0000E-02
             3.9849E+00
 PARAMETER:  8.0757E-02 -1.9393E+00  2.6505E-02  5.7918E-01 -3.1927E-01 -5.1788E-02 -3.3270E+00 -5.6162E+00 -1.2306E-01 -4.5471E+00
             1.4825E+00
 GRADIENT:   4.4483E-02 -1.2313E-02 -4.2211E-02 -7.5578E-01  4.4142E-02  6.3773E-03  5.2785E-05  0.0000E+00  1.2871E-02  3.7242E-05
            -8.5834E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1243.50316332928        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:     1418
 NPARAMETR:  9.8098E-01  1.3047E-01  9.2922E-01  1.6147E+00  6.5752E-01  8.5918E-01  1.0000E-02  1.0000E-02  8.0003E-01  1.0000E-02
             3.9851E+00
 PARAMETER:  8.0797E-02 -1.9366E+00  2.6590E-02  5.7914E-01 -3.1928E-01 -5.1780E-02 -1.3822E+01 -5.6162E+00 -1.2310E-01 -4.5874E+00
             1.4826E+00
 GRADIENT:   7.6276E-02  8.8566E-03  7.6605E-02 -5.1866E-01 -1.6618E-01  4.0558E-03  0.0000E+00  0.0000E+00 -1.8984E-02  0.0000E+00
            -8.1813E-02

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1243.50316443504        NO. OF FUNC. EVALS.:  68
 CUMULATIVE NO. OF FUNC. EVALS.:     1486
 NPARAMETR:  9.8103E-01  1.3021E-01  9.2827E-01  1.6132E+00  6.5836E-01  8.5924E-01  1.0000E-02  1.0000E-02  8.0084E-01  1.0000E-02
             3.9852E+00
 PARAMETER:  8.0792E-02 -1.9365E+00  2.6570E-02  5.7916E-01 -3.1925E-01 -5.1784E-02 -1.5482E+01 -5.6162E+00 -1.2309E-01 -4.5935E+00
             1.4826E+00
 GRADIENT:  -1.5855E-02  2.0161E-03  5.0165E-02  2.4725E-01 -1.3219E-01 -2.4539E-03  0.0000E+00  0.0000E+00 -1.8045E-02  0.0000E+00
            -2.7760E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1486
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.0976E-04 -3.9828E-05  1.0568E-04 -1.5721E-02 -2.0101E-05
 SE:             2.7964E-02  1.8501E-05  1.2236E-04  2.4138E-02  1.8792E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8546E-01  3.1341E-02  3.8778E-01  5.1485E-01  9.1482E-01

 ETASHRINKSD(%)  6.3184E+00  9.9938E+01  9.9590E+01  1.9133E+01  9.9370E+01
 ETASHRINKVR(%)  1.2238E+01  1.0000E+02  9.9998E+01  3.4605E+01  9.9996E+01
 EBVSHRINKSD(%)  6.0116E+00  9.9940E+01  9.9528E+01  1.8686E+01  9.9335E+01
 EBVSHRINKVR(%)  1.1662E+01  1.0000E+02  9.9998E+01  3.3881E+01  9.9996E+01
 RELATIVEINF(%)  7.6405E+01  1.0089E-06  8.5458E-05  3.5139E+00  9.4036E-05
 EPSSHRINKSD(%)  1.9333E+01
 EPSSHRINKVR(%)  3.4928E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1243.5031644350429     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -508.35233787130471     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.01
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1243.503       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.30E-01  9.29E-01  1.61E+00  6.58E-01  8.59E-01  1.00E-02  1.00E-02  8.00E-01  1.00E-02  3.99E+00
 


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
+        1.38E+03
 
 TH 2
+       -1.07E+02  2.51E+02
 
 TH 3
+       -1.50E+00  1.14E+02  2.58E+02
 
 TH 4
+       -1.08E+02  2.81E+02  1.18E+01  4.33E+02
 
 TH 5
+        5.13E+01 -3.45E+02 -5.07E+02 -1.71E+02  1.09E+03
 
 TH 6
+       -5.58E+00 -1.61E+01  1.13E+01 -2.58E+01 -4.97E+00  2.07E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.36E+01 -5.56E+01  6.03E+00 -1.96E+01  2.61E+01  4.98E+00  0.00E+00  0.00E+00  1.25E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  2.23E+00
 
 TH11
+       -2.25E+01 -1.08E+01 -2.63E+00 -8.99E+00  8.48E+00  5.08E+00  0.00E+00  0.00E+00  1.29E+01  0.00E+00  2.80E+01
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       25.338
Stop Time:
Wed Sep 29 13:54:59 CDT 2021
