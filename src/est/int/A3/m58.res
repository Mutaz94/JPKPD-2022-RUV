Sat Sep 18 01:48:23 CDT 2021
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
$DATA ../../../../data/int/A3/dat58.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
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

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m58.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -158.966636637462        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.2053E+02  1.5131E+02  1.8204E+02 -4.1167E+01  1.2470E+02  2.8471E+01 -1.3914E+02 -1.4991E+02 -3.5976E+01 -6.7504E+01
            -6.8349E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2627.36447274955        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.3719E-01  9.9137E-01  9.5207E-01  1.0407E+00  9.7414E-01  8.1129E-01  1.0381E+00  9.7095E-01  9.9946E-01  9.1778E-01
             2.8970E+00
 PARAMETER:  3.5128E-02  9.1335E-02  5.0879E-02  1.3988E-01  7.3800E-02 -1.0913E-01  1.3738E-01  7.0519E-02  9.9462E-02  1.4201E-02
             1.1637E+00
 GRADIENT:  -2.0319E+00 -3.8801E+00 -2.7825E+00 -7.0255E-01  9.9394E+00 -2.1412E+01 -3.2334E+00  8.6225E+00  2.1128E+00  8.4066E+00
             1.0599E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2633.95435948998        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.3612E-01  7.2410E-01  6.6319E-01  1.2114E+00  6.8034E-01  8.4322E-01  1.2660E+00  5.8332E-01  9.4641E-01  6.2876E-01
             2.8571E+00
 PARAMETER:  3.3993E-02 -2.2282E-01 -3.1070E-01  2.9174E-01 -2.8517E-01 -7.0529E-02  3.3584E-01 -4.3902E-01  4.4916E-02 -3.6400E-01
             1.1498E+00
 GRADIENT:  -4.8487E+00  2.0441E+01 -2.4742E+01  8.9211E+01  2.9599E+01 -9.1253E+00 -5.0223E+00  5.8192E+00 -8.4380E+00  2.0428E+00
             8.9452E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2643.68318965815        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.3149E-01  4.7206E-01  4.2030E-01  1.2517E+00  4.3093E-01  8.8781E-01  1.5452E+00  1.4988E-01  1.0169E+00  4.4966E-01
             2.6420E+00
 PARAMETER:  2.9025E-02 -6.5065E-01 -7.6678E-01  3.2451E-01 -7.4181E-01 -1.9000E-02  5.3513E-01 -1.7980E+00  1.1679E-01 -6.9925E-01
             1.0715E+00
 GRADIENT:  -3.9587E+00  2.9733E+01  2.0439E+01  6.5110E+01 -1.6564E+01  6.4176E+00  3.0377E+00 -4.1524E-01 -1.8066E+00 -6.2329E+00
            -5.5574E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2645.80232985798        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.3046E-01  3.1793E-01  2.6646E-01  1.1897E+00  2.8924E-01  9.0526E-01  1.4902E+00  4.3378E-02  1.1167E+00  4.2184E-01
             2.5066E+00
 PARAMETER:  2.7924E-02 -1.0459E+00 -1.2225E+00  2.7372E-01 -1.1405E+00  4.6205E-04  4.9890E-01 -3.0378E+00  2.1035E-01 -7.6312E-01
             1.0189E+00
 GRADIENT:   2.5573E-01  2.6275E+01  4.7581E+01  4.1594E+01 -4.2366E+01  1.0669E+01 -9.9391E+00 -1.9452E-01 -1.0284E+01 -3.2626E+01
            -1.3192E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2645.97872972593        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.3082E-01  2.6872E-01  2.1944E-01  1.1433E+00  2.4587E-01  9.0888E-01  1.4376E+00  2.0462E-02  1.1617E+00  4.4689E-01
             2.4566E+00
 PARAMETER:  2.8311E-02 -1.2141E+00 -1.4167E+00  2.3394E-01 -1.3029E+00  4.4569E-03  4.6296E-01 -3.7892E+00  2.4986E-01 -7.0545E-01
             9.9878E-01
 GRADIENT:   2.4526E+00  1.8085E+01  7.5137E+01  4.1301E+01 -8.8160E+01  1.1541E+01 -1.4663E+01 -4.6284E-02 -1.9819E+01 -4.6761E+01
            -1.5305E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2646.18981293188        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  9.3180E-01  2.3890E-01  1.9082E-01  1.0964E+00  2.2023E-01  9.0568E-01  1.4027E+00  1.1499E-02  1.1896E+00  4.7255E-01
             2.4434E+00
 PARAMETER:  2.9368E-02 -1.3317E+00 -1.5565E+00  1.9202E-01 -1.4131E+00  9.3134E-04  4.3842E-01 -4.3655E+00  2.7359E-01 -6.4961E-01
             9.9339E-01
 GRADIENT:   4.6685E+00  6.4322E+00  8.9100E+01  3.3356E+01 -1.1571E+02  1.0401E+01 -1.7477E+01 -1.2630E-02 -2.7005E+01 -5.3396E+01
            -1.4904E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2646.66586843545        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      510
 NPARAMETR:  9.3327E-01  2.0856E-01  1.6037E-01  1.0265E+00  1.9387E-01  8.9619E-01  1.3794E+00  1.0000E-02  1.2343E+00  5.1374E-01
             2.4540E+00
 PARAMETER:  3.0944E-02 -1.4676E+00 -1.7302E+00  1.2612E-01 -1.5406E+00 -9.5997E-03  4.2164E-01 -5.0371E+00  3.1053E-01 -5.6605E-01
             9.9772E-01
 GRADIENT:   6.8837E+00 -9.1265E+00  9.2550E+01  1.8152E+01 -1.3162E+02  7.3319E+00 -1.8705E+01  0.0000E+00 -3.2435E+01 -5.4064E+01
            -1.2182E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2647.48145118997        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      580
 NPARAMETR:  9.3472E-01  1.8288E-01  1.3181E-01  9.3870E-01  1.7067E-01  8.8123E-01  1.3950E+00  1.0000E-02  1.3173E+00  5.7481E-01
             2.5014E+00
 PARAMETER:  3.2493E-02 -1.5989E+00 -1.9264E+00  3.6746E-02 -1.6680E+00 -2.6440E-02  4.3292E-01 -5.6460E+00  3.7556E-01 -4.5371E-01
             1.0169E+00
 GRADIENT:   7.0629E+00 -1.9655E+01  7.1765E+01 -2.2415E-01 -1.1131E+02  2.5270E+00 -1.4803E+01  0.0000E+00 -2.9624E+01 -4.1894E+01
            -6.4023E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2648.34923677554        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      650
 NPARAMETR:  9.3473E-01  1.7096E-01  1.1466E-01  8.8060E-01  1.5853E-01  8.7323E-01  1.4388E+00  1.0000E-02  1.4284E+00  6.3961E-01
             2.5437E+00
 PARAMETER:  3.2502E-02 -1.6663E+00 -2.0658E+00 -2.7152E-02 -1.7418E+00 -3.5551E-02  4.6382E-01 -5.8398E+00  4.5656E-01 -3.4689E-01
             1.0336E+00
 GRADIENT:   4.0840E+00 -1.5086E+01  3.6240E+01 -7.3598E+00 -5.8727E+01 -1.0039E-01 -8.0302E+00  0.0000E+00 -1.6869E+01 -2.0835E+01
            -1.6360E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2648.72274741044        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  9.3388E-01  1.6887E-01  1.0816E-01  8.6450E-01  1.5522E-01  8.7289E-01  1.4834E+00  1.0000E-02  1.5151E+00  6.8024E-01
             2.5536E+00
 PARAMETER:  3.1593E-02 -1.6786E+00 -2.1242E+00 -4.5607E-02 -1.7629E+00 -3.5940E-02  4.9437E-01 -5.7415E+00  5.1550E-01 -2.8531E-01
             1.0375E+00
 GRADIENT:   1.3520E+00 -5.2694E+00  9.2391E+00 -5.1739E+00 -1.5660E+01 -4.1996E-01 -2.0812E+00  0.0000E+00 -5.0007E+00 -6.3091E+00
            -2.4114E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2649.57042438515        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  9.3328E-01  1.7095E-01  1.0898E-01  8.7362E-01  1.5659E-01  8.7479E-01  1.4917E+00  1.0000E-02  1.5323E+00  6.9349E-01
             2.5493E+00
 PARAMETER:  3.0953E-02 -1.6664E+00 -2.1166E+00 -3.5110E-02 -1.7541E+00 -3.3770E-02  4.9994E-01 -5.6044E+00  5.2674E-01 -2.6602E-01
             1.0358E+00
 GRADIENT:   2.2409E-01 -1.3891E+00  1.3767E+00 -1.5268E+00 -2.2892E+00  1.1124E-02 -9.4887E-01  0.0000E+00 -8.7780E-01 -5.2533E-01
            -6.6917E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2660.58789393227        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      957
 NPARAMETR:  9.3413E-01  2.2979E-01  1.7367E-01  1.0509E+00  2.1174E-01  8.7990E-01  1.4581E+00  1.5319E-02  1.3009E+00  6.3364E-01
             2.5631E+00
 PARAMETER:  3.1860E-02 -1.3706E+00 -1.6506E+00  1.4964E-01 -1.4524E+00 -2.7942E-02  4.7713E-01 -4.0787E+00  3.6309E-01 -3.5628E-01
             1.0412E+00
 GRADIENT:  -1.7555E-01 -5.6963E-01 -6.8568E-01  2.5897E+00 -3.5439E+00  3.3967E-01 -2.1845E-01 -3.9547E-04 -1.0697E+00  3.0568E-01
             1.5956E-01

0ITERATION NO.:   64    OBJECTIVE VALUE:  -2660.62922099811        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1084
 NPARAMETR:  9.3417E-01  2.3455E-01  1.7894E-01  1.0582E+00  2.1625E-01  8.7896E-01  1.4585E+00  1.6902E-02  1.2926E+00  6.2905E-01
             2.5661E+00
 PARAMETER:  3.1901E-02 -1.3501E+00 -1.6207E+00  1.5653E-01 -1.4313E+00 -2.9014E-02  4.7743E-01 -3.9803E+00  3.5664E-01 -3.6354E-01
             1.0424E+00
 GRADIENT:  -3.7982E-04 -7.2230E-02 -3.9226E-02  9.5049E-02  6.3454E-02 -2.4598E-02  1.8964E-02 -6.4264E-04 -3.7607E-02  5.3254E-02
             5.7084E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1084
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0708E-03  8.7132E-03  1.5708E-04 -7.4302E-03  6.5329E-03
 SE:             2.9215E-02  2.4524E-02  4.6135E-04  2.7581E-02  2.4801E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4349E-01  7.2237E-01  7.3349E-01  7.8763E-01  7.9223E-01

 ETASHRINKSD(%)  2.1246E+00  1.7841E+01  9.8454E+01  7.5986E+00  1.6915E+01
 ETASHRINKVR(%)  4.2040E+00  3.2499E+01  9.9976E+01  1.4620E+01  3.0969E+01
 EBVSHRINKSD(%)  2.2041E+00  1.6307E+01  9.8478E+01  5.8429E+00  1.8023E+01
 EBVSHRINKVR(%)  4.3596E+00  2.9955E+01  9.9977E+01  1.1344E+01  3.2798E+01
 RELATIVEINF(%)  9.5613E+01  1.5087E+01  1.2980E-03  4.1628E+01  3.4177E+00
 EPSSHRINKSD(%)  1.8716E+01
 EPSSHRINKVR(%)  3.3929E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2660.6292209981061     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1006.5398612296954     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.95
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2660.629       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.34E-01  2.35E-01  1.79E-01  1.06E+00  2.16E-01  8.79E-01  1.46E+00  1.69E-02  1.29E+00  6.29E-01  2.57E+00
 


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
+        1.58E+03
 
 TH 2
+        2.12E+00  6.73E+03
 
 TH 3
+        4.91E+01 -1.11E+03  2.58E+04
 
 TH 4
+       -1.84E+00 -8.26E+01 -4.70E+02  5.24E+02
 
 TH 5
+       -2.15E+01 -5.45E+03 -2.74E+04 -4.35E+02  4.03E+04
 
 TH 6
+        9.77E+00 -9.68E+00  4.08E+01 -7.62E+00 -2.81E+01  2.36E+02
 
 TH 7
+       -1.80E+00  4.82E+01  3.07E+01 -1.53E+00 -6.84E+01  2.67E-01  4.49E+01
 
 TH 8
+       -1.73E+00 -2.18E+00  3.75E+00 -5.29E+00 -1.12E+00  1.96E+00 -3.57E-01 -5.08E+00
 
 TH 9
+        1.16E+01 -2.67E+00  1.90E+02 -7.02E+00  4.00E+01 -1.82E+00 -1.19E-01 -3.68E-01  8.57E+01
 
 TH10
+       -9.00E+00 -2.93E+00  9.60E+01  1.05E+01  6.38E+01 -2.69E+00  7.82E+00  1.36E+00  1.05E+01  2.31E+02
 
 TH11
+       -2.20E+01 -1.56E+01 -1.56E+02 -5.12E+00  1.36E+02  3.89E+00  7.62E+00  5.68E-01  4.39E+00  1.63E+01  1.65E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.782
Stop Time:
Sat Sep 18 01:49:01 CDT 2021
