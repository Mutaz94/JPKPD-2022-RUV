Sat Sep 18 09:32:44 CDT 2021
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
$DATA ../../../../data/spa/A1/dat95.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1500.26775570113        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0726E+01 -8.2300E+01 -3.4773E+01 -1.1722E+02  7.3531E+01  1.3952E+01 -2.3923E+01  1.0087E+01 -8.1364E+01  8.6631E+00
            -2.3253E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1553.07487779627        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0006E+00  1.0294E+00  1.1404E+00  1.1250E+00  9.7296E-01  9.3775E-01  1.0597E+00  8.9727E-01  1.3509E+00  8.0938E-01
             1.5827E+00
 PARAMETER:  1.0056E-01  1.2896E-01  2.3140E-01  2.1780E-01  7.2585E-02  3.5723E-02  1.5802E-01 -8.3976E-03  4.0076E-01 -1.1148E-01
             5.5914E-01
 GRADIENT:  -2.5884E+00  4.0625E+01  2.2353E+01  4.8876E+01 -3.0905E+01 -9.3679E+00  4.9833E+00  1.0147E+00  1.8387E+01  1.9425E+00
             3.1270E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1556.04294636654        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0065E+00  9.5911E-01  6.8179E-01  1.1505E+00  7.4310E-01  9.4632E-01  1.3402E+00  4.3528E-01  1.1893E+00  5.4413E-01
             1.5567E+00
 PARAMETER:  1.0653E-01  5.8250E-02 -2.8303E-01  2.4020E-01 -1.9693E-01  4.4827E-02  3.9279E-01 -7.3177E-01  2.7334E-01 -5.0857E-01
             5.4254E-01
 GRADIENT:   4.9455E+00  2.4896E+01 -3.8930E+01  8.7038E+01  3.5731E+01 -7.3367E+00  8.4587E+00  2.0528E+00  1.3622E+01  5.0835E+00
             3.1429E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1561.45219693521        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0016E+00  8.0655E-01  6.6610E-01  1.1500E+00  6.7159E-01  9.6458E-01  1.4525E+00  3.7876E-01  1.0542E+00  4.7233E-01
             1.4365E+00
 PARAMETER:  1.0164E-01 -1.1500E-01 -3.0632E-01  2.3973E-01 -2.9811E-01  6.3934E-02  4.7332E-01 -8.7086E-01  1.5279E-01 -6.5008E-01
             4.6218E-01
 GRADIENT:   1.2524E+00  3.2360E+00  6.2321E+00 -3.0791E+00 -8.2205E+00  5.6179E-01 -9.3468E-01 -2.9231E-01 -7.7167E-01 -9.0348E-01
            -4.2494E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1562.25708410667        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0003E+00  6.6170E-01  7.1313E-01  1.2446E+00  6.4931E-01  9.6006E-01  1.6737E+00  4.2937E-01  1.0065E+00  5.0726E-01
             1.4513E+00
 PARAMETER:  1.0030E-01 -3.1294E-01 -2.3809E-01  3.1879E-01 -3.3185E-01  5.9237E-02  6.1501E-01 -7.4542E-01  1.0643E-01 -5.7873E-01
             4.7245E-01
 GRADIENT:   2.3789E-01  4.7457E+00  4.7428E+00  7.6857E+00 -8.1361E+00 -5.0468E-01 -5.5616E-01 -2.0637E-01 -2.3248E-01 -3.9610E-01
             1.5010E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1562.43908353797        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.9955E-01  6.0349E-01  7.1775E-01  1.2791E+00  6.3420E-01  9.5926E-01  1.7852E+00  4.3480E-01  9.8729E-01  5.0742E-01
             1.4472E+00
 PARAMETER:  9.9547E-02 -4.0502E-01 -2.3164E-01  3.4618E-01 -3.5540E-01  5.8405E-02  6.7953E-01 -7.3287E-01  8.7209E-02 -5.7842E-01
             4.6960E-01
 GRADIENT:  -4.7661E-02  4.7976E+00  4.8684E+00  1.2640E+01 -8.2052E+00 -5.4022E-01 -5.5290E-01 -5.2010E-01 -3.6730E-01 -1.2300E+00
             9.1134E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1562.45529036294        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      470
 NPARAMETR:  9.9941E-01  6.0767E-01  7.2245E-01  1.2718E+00  6.3952E-01  9.6012E-01  1.7767E+00  4.4432E-01  9.9115E-01  5.1442E-01
             1.4444E+00
 PARAMETER:  9.9414E-02 -3.9812E-01 -2.2510E-01  3.4046E-01 -3.4704E-01  5.9305E-02  6.7478E-01 -7.1122E-01  9.1115E-02 -5.6471E-01
             4.6771E-01
 GRADIENT:  -1.7769E+01 -3.0243E-01  1.5271E+00 -1.0621E+01 -5.6210E+00 -1.6533E+00 -1.3655E+00 -2.3250E-01 -5.8875E-01 -5.5718E-01
            -3.8424E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1563.55605655095        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  9.9935E-01  3.4601E-01  8.9695E-01  1.4657E+00  6.5666E-01  9.5731E-01  2.4916E+00  6.2980E-01  9.4585E-01  6.1304E-01
             1.4456E+00
 PARAMETER:  9.9351E-02 -9.6130E-01 -8.7564E-03  4.8234E-01 -3.2059E-01  5.6374E-02  1.0129E+00 -3.6236E-01  4.4327E-02 -3.8933E-01
             4.6854E-01
 GRADIENT:  -3.2877E+00  5.4986E+00  1.8194E+01  1.8161E+01 -2.4248E+01 -3.6616E-02 -2.3500E+00  1.7344E-01  2.0030E+00  5.5649E-01
             9.1608E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1564.43463373320        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      826
 NPARAMETR:  9.9858E-01  2.1566E-01  8.5085E-01  1.5124E+00  6.1407E-01  9.5447E-01  3.2015E+00  5.8777E-01  9.1971E-01  6.2009E-01
             1.4415E+00
 PARAMETER:  9.8577E-02 -1.4340E+00 -6.1514E-02  5.1370E-01 -3.8765E-01  5.3404E-02  1.2636E+00 -4.3141E-01  1.6301E-02 -3.7790E-01
             4.6569E-01
 GRADIENT:   1.5999E+00 -3.2102E-01 -1.0513E+00 -3.1270E+00  1.7009E+00  1.2169E-02 -1.0580E+00  6.0814E-03  2.0413E+00  7.8951E-01
             5.7723E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1564.72172506911        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  9.9552E-01  1.1518E-01  8.5694E-01  1.5687E+00  5.9733E-01  9.5154E-01  4.4058E+00  6.3727E-01  8.8992E-01  6.0747E-01
             1.4371E+00
 PARAMETER:  9.5508E-02 -2.0613E+00 -5.4383E-02  5.5023E-01 -4.1528E-01  5.0330E-02  1.5829E+00 -3.5056E-01 -1.6621E-02 -3.9846E-01
             4.6263E-01
 GRADIENT:   1.5334E+00  4.8944E+00 -3.5942E+00 -3.7690E+00  3.9742E+00 -4.2911E-01  9.6935E+00 -4.6409E-01 -8.1991E+00 -7.8288E-01
            -1.0634E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1564.88886415622        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1179
 NPARAMETR:  9.9275E-01  5.8832E-02  8.1982E-01  1.5868E+00  5.6944E-01  9.5074E-01  6.0067E+00  5.9957E-01  8.7771E-01  5.9667E-01
             1.4384E+00
 PARAMETER:  9.2727E-02 -2.7331E+00 -9.8669E-02  5.6170E-01 -4.6310E-01  4.9481E-02  1.8929E+00 -4.1155E-01 -3.0444E-02 -4.1638E-01
             4.6352E-01
 GRADIENT:  -2.3882E+00  2.3240E+00  2.5932E+00  2.0211E+00 -5.0197E+00 -3.8789E-01  4.2864E+00 -1.0650E+00 -1.6986E+00 -1.0577E+00
            -5.0678E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1565.04002220419        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1355
 NPARAMETR:  9.9294E-01  3.9685E-02  8.4155E-01  1.5938E+00  5.7840E-01  9.5080E-01  7.1134E+00  6.5224E-01  8.6916E-01  5.7717E-01
             1.4378E+00
 PARAMETER:  9.2911E-02 -3.1268E+00 -7.2512E-02  5.6610E-01 -4.4750E-01  4.9553E-02  2.0620E+00 -3.2735E-01 -4.0224E-02 -4.4961E-01
             4.6309E-01
 GRADIENT:  -1.2680E+00 -3.9973E+00  9.9596E+00  1.0928E+01 -1.1908E+01  4.7146E-01 -1.0127E+01  1.1160E+00  4.6096E+00  1.8487E+00
             1.0505E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1565.20869903015        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1534
 NPARAMETR:  9.9209E-01  2.0061E-02  7.9602E-01  1.5907E+00  5.5320E-01  9.5009E-01  9.3574E+00  6.0483E-01  8.6799E-01  5.5230E-01
             1.4349E+00
 PARAMETER:  9.2059E-02 -3.8090E+00 -1.2813E-01  5.6420E-01 -4.9203E-01  4.8806E-02  2.3362E+00 -4.0281E-01 -4.1574E-02 -4.9366E-01
             4.6109E-01
 GRADIENT:  -2.5079E+00 -4.4988E+00  1.6812E+01  1.5827E+01 -2.2999E+01  3.0047E-01 -1.0798E+01  5.9918E-01  5.2111E+00  1.4337E+00
             3.5321E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1565.57550702085        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1711
 NPARAMETR:  9.9311E-01  1.6769E-02  7.7947E-01  1.5865E+00  5.4777E-01  9.4978E-01  1.0382E+01  5.9614E-01  8.5772E-01  5.6214E-01
             1.4346E+00
 PARAMETER:  9.3084E-02 -3.9882E+00 -1.4914E-01  5.6155E-01 -5.0189E-01  4.8470E-02  2.4400E+00 -4.1727E-01 -5.3483E-02 -4.7600E-01
             4.6090E-01
 GRADIENT:   2.7900E-01 -2.8322E+00  3.2923E+00  1.0333E+01 -6.6176E+00  1.8178E-01 -6.1764E+00  1.3774E+00  7.6427E-01  2.4815E+00
             1.0955E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1565.60508822093        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1879
 NPARAMETR:  9.9283E-01  1.5769E-02  7.7813E-01  1.5843E+00  5.4718E-01  9.4957E-01  1.0685E+01  5.8071E-01  8.5733E-01  5.6547E-01
             1.4365E+00
 PARAMETER:  9.2801E-02 -4.0494E+00 -1.5081E-01  5.6012E-01 -5.0294E-01  4.8277E-02  2.4691E+00 -4.4308E-01 -5.3836E-02 -4.6995E-01
             4.6236E-01
 GRADIENT:  -7.2230E-02  2.3150E+02  7.2783E-01 -1.6885E+03  1.8812E+03  2.6595E-02  3.7923E+02  2.8107E-01  1.1625E-01  5.2962E-01
             3.2491E-01
 NUMSIGDIG:         3.5         3.3         2.7         3.3         3.3         2.8         3.3         2.2         2.2         2.7
                    2.8

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1879
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.7921E-04  2.7424E-02 -9.8578E-03 -1.1441E-02 -7.6920E-03
 SE:             2.9670E-02  1.0080E-02  1.4030E-02  2.8536E-02  1.8418E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8711E-01  6.5143E-03  4.8230E-01  6.8846E-01  6.7621E-01

 ETASHRINKSD(%)  6.0321E-01  6.6232E+01  5.2997E+01  4.4015E+00  3.8297E+01
 ETASHRINKVR(%)  1.2028E+00  8.8597E+01  7.7907E+01  8.6092E+00  6.1928E+01
 EBVSHRINKSD(%)  8.9555E-01  7.3724E+01  5.2228E+01  4.8550E+00  3.6049E+01
 EBVSHRINKVR(%)  1.7831E+00  9.3096E+01  7.7178E+01  9.4742E+00  5.9103E+01
 RELATIVEINF(%)  9.8118E+01  5.9117E+00  1.4562E+00  6.5013E+01  2.5974E+00
 EPSSHRINKSD(%)  4.0267E+01
 EPSSHRINKVR(%)  6.4319E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1565.6050882209345     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -830.45426165719630     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.08
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1565.605       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.93E-01  1.58E-02  7.78E-01  1.58E+00  5.47E-01  9.50E-01  1.07E+01  5.81E-01  8.57E-01  5.66E-01  1.44E+00
 


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
+        1.24E+03
 
 TH 2
+       -1.59E+02  1.23E+08
 
 TH 3
+       -1.86E+02 -7.48E+03  4.26E+03
 
 TH 4
+        3.22E+01 -2.63E+03 -1.00E+02  3.01E+05
 
 TH 5
+       -1.41E+02  1.02E+04 -7.33E+06 -9.70E+05  3.13E+06
 
 TH 6
+       -2.28E+01 -3.57E+01  8.18E+01 -2.03E+01  5.60E+01  2.14E+02
 
 TH 7
+        3.61E-01  1.58E+05 -2.60E+01 -1.34E+01  4.21E+01 -5.30E-01  7.24E+02
 
 TH 8
+       -6.15E+01 -7.39E+02  8.58E+02 -1.10E+02  4.09E+02  2.41E+01 -3.70E+00  3.12E+02
 
 TH 9
+       -2.21E+01  4.02E+02  4.26E+02 -2.22E+02  7.33E+02  9.52E+00 -3.02E+00  1.27E+02  2.97E+02
 
 TH10
+       -1.62E+02 -1.13E+03  2.50E+03  1.00E+06 -3.24E+06  6.08E+01 -6.66E+00  8.00E+02  3.28E+02  2.27E+03
 
 TH11
+       -4.46E+01  6.18E+02  5.33E+02 -1.46E+02 -1.29E+06  1.62E+01  4.72E-01  1.86E+02  8.03E+01  4.90E+02  2.12E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.033
Stop Time:
Sat Sep 18 09:33:16 CDT 2021
