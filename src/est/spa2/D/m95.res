Thu Sep 30 10:07:39 CDT 2021
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
$DATA ../../../../data/spa2/D/dat95.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   35692.5217824547        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.1658E+02  7.6182E+02  4.6465E+01  6.2686E+02  1.5344E+02 -3.6716E+03 -1.6620E+03 -6.1496E+01 -2.3490E+03 -1.0615E+03
            -6.6653E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -324.090612965157        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3236E+00  1.2167E+00  9.7243E-01  1.4629E+00  1.1214E+00  2.1292E+00  1.4034E+00  9.7038E-01  1.2863E+00  9.4816E-01
             1.4325E+01
 PARAMETER:  3.8035E-01  2.9611E-01  7.2043E-02  4.8044E-01  2.1460E-01  8.5576E-01  4.3893E-01  6.9928E-02  3.5176E-01  4.6770E-02
             2.7620E+00
 GRADIENT:   3.0650E+01  2.3674E+01 -4.3126E+00  4.9759E+01 -9.9573E-01  3.2772E+01 -3.2445E+01  4.3460E+00 -4.0889E+01  8.3562E+00
            -9.1180E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -440.377010924874        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.3211E+00  2.4430E+00  2.2362E+00  1.0753E+00  2.7565E+00  2.7426E+00  5.0948E+00  3.5335E-01  2.5679E+00  3.8817E-02
             1.4894E+01
 PARAMETER:  3.7845E-01  9.9323E-01  9.0477E-01  1.7257E-01  1.1140E+00  1.1089E+00  1.7282E+00 -9.4029E-01  1.0431E+00 -3.1489E+00
             2.8010E+00
 GRADIENT:   1.4644E+01  1.0599E+01 -1.3135E+01 -6.1147E+00  1.2279E+01  6.2900E+01  6.5074E+01  3.8217E-02  2.7682E+01  4.9417E-03
             1.5153E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -463.848450510342        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      282             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2586E+00  1.3424E+00  2.4141E+00  1.5233E+00  1.8363E+00  2.6144E+00  5.7782E+00  1.7289E-01  1.0339E+00  4.3412E-02
             1.4247E+01
 PARAMETER:  3.2999E-01  3.9445E-01  9.8134E-01  5.2090E-01  7.0776E-01  1.0610E+00  1.8541E+00 -1.6551E+00  1.3330E-01 -3.0370E+00
             2.7566E+00
 GRADIENT:   1.5061E+01  7.1458E+00 -8.3425E+00  4.4557E+01 -1.8073E+00  5.9784E+01  2.3179E+01  4.1292E-02 -3.7664E+00  1.1321E-02
             7.7436E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -471.138672715592        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.2590E+00  1.3756E+00  6.9168E+00  1.2602E+00  2.4807E+00  2.6048E+00  5.7025E+00  5.2859E-02  9.1215E-01  1.5654E-02
             1.2800E+01
 PARAMETER:  3.3028E-01  4.1889E-01  2.0339E+00  3.3128E-01  1.0086E+00  1.0574E+00  1.8409E+00 -2.8401E+00  8.0462E-03 -4.0570E+00
             2.6494E+00
 GRADIENT:   3.4614E+01  2.3864E+00 -4.4226E-01 -2.3056E+00  2.0595E-01  5.7622E+01  4.3270E+01  2.1321E-04 -1.8137E+00  1.1822E-03
             8.7862E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -472.204557739388        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  1.2582E+00  1.3068E+00  8.1741E+00  1.3714E+00  2.6106E+00  2.5872E+00  5.7256E+00  5.3094E-02  1.1290E+00  1.3049E-02
             1.3073E+01
 PARAMETER:  3.2970E-01  3.6756E-01  2.2010E+00  4.1584E-01  1.0596E+00  1.0506E+00  1.8449E+00 -2.8357E+00  2.2137E-01 -4.2390E+00
             2.6705E+00
 GRADIENT:   1.8859E+01  4.3576E-01 -6.6455E-01 -3.6380E-01  1.8606E+00  4.2513E+01 -2.1716E+01  1.5051E-04  1.6835E+00  6.9451E-04
            -3.7831E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -485.792170303520        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      666
 NPARAMETR:  1.2339E+00  7.0190E-01  1.4658E+05  1.7185E+00  2.8550E+00  1.9437E+00  7.7091E+00  1.0000E-02  1.2319E+00  1.0000E-02
             1.3452E+01
 PARAMETER:  3.1021E-01 -2.5396E-01  1.1995E+01  6.4145E-01  1.1491E+00  7.6458E-01  2.1424E+00 -6.0638E+00  3.0856E-01 -1.5245E+01
             2.6992E+00
 GRADIENT:   2.3338E+01 -1.3332E+00  6.4448E-06  4.0968E+00 -1.4247E+00  8.8544E+00  4.8673E+00  0.0000E+00 -2.3122E-01  0.0000E+00
             9.8446E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -489.938553329655        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      787
 NPARAMETR:  1.0627E+00  7.1858E-01  2.4674E+05  1.7207E+00  2.8570E+00  1.8882E+00  7.5926E+00  1.0000E-02  1.1315E+00  1.0000E-02
             1.3445E+01
 PARAMETER:  1.6081E-01 -2.3048E-01  1.2516E+01  6.4272E-01  1.1498E+00  7.3564E-01  2.1272E+00 -6.2518E+00  2.2357E-01 -1.5790E+01
             2.6986E+00
 GRADIENT:  -4.1542E+02 -2.3682E+01  1.5623E-05 -2.2876E+01  5.1817E+00  2.4486E+02  9.4542E+01  0.0000E+00  5.9787E+00  0.0000E+00
             2.5942E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -491.648558786101        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:      960
 NPARAMETR:  1.0791E+00  7.0369E-01  2.5735E+05  1.6887E+00  3.0043E+00  1.7026E+00  7.8527E+00  1.0000E-02  1.1774E+00  1.0000E-02
             1.3492E+01
 PARAMETER:  1.7617E-01 -2.5141E-01  1.2558E+01  6.2394E-01  1.2001E+00  6.3216E-01  2.1609E+00 -6.2518E+00  2.6330E-01 -1.5790E+01
             2.7021E+00
 GRADIENT:  -2.6443E+01  3.2979E-02 -2.7352E-06  5.2070E-01  2.2274E-01  6.3813E-01  4.3650E+00  0.0000E+00  5.9218E-01  0.0000E+00
             9.7142E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -491.896362256225        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1156             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0867E+00  6.8703E-01  2.2167E+05  1.7022E+00  3.0258E+00  1.7368E+00  8.2112E+00  1.0000E-02  1.1533E+00  1.0000E-02
             1.3546E+01
 PARAMETER:  1.8317E-01 -2.7538E-01  1.2409E+01  6.3195E-01  1.2072E+00  6.5207E-01  2.2055E+00 -6.2518E+00  2.4260E-01 -1.5790E+01
             2.7061E+00
 GRADIENT:  -2.4178E+01  2.3785E+00 -5.4873E-06  1.4802E+01  6.6283E-01  8.2071E+00  1.1731E+02  0.0000E+00 -1.0005E+00  0.0000E+00
             4.6588E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -491.965182989924        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     1359             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0881E+00  6.4556E-01  1.9889E+05  1.7158E+00  3.0211E+00  1.7354E+00  8.2507E+00  1.0000E-02  1.1615E+00  1.0000E-02
             1.3547E+01
 PARAMETER:  1.8444E-01 -3.3763E-01  1.2300E+01  6.3987E-01  1.2056E+00  6.5122E-01  2.2103E+00 -6.2518E+00  2.4975E-01 -1.5790E+01
             2.7061E+00
 GRADIENT:  -5.3402E+01  1.5332E+00 -1.3514E-06  2.1724E+01  3.9259E-01  2.0519E+01  1.1789E+02  0.0000E+00 -3.4357E+00  0.0000E+00
             5.7561E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -491.997730119524        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1558             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0879E+00  6.3680E-01  2.3900E+05  1.7316E+00  3.0032E+00  1.7335E+00  8.3263E+00  1.0000E-02  1.1594E+00  1.0000E-02
             1.3535E+01
 PARAMETER:  1.8423E-01 -3.5131E-01  1.2484E+01  6.4905E-01  1.1997E+00  6.5017E-01  2.2194E+00 -6.2518E+00  2.4787E-01 -1.5790E+01
             2.7053E+00
 GRADIENT:  -4.4806E+01  1.6729E+00 -5.4306E-08  2.3313E+01  1.6073E-02  1.8970E+01  1.1966E+02  0.0000E+00 -4.1297E+00  0.0000E+00
             4.8019E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -492.011656833872        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1757             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0879E+00  6.2269E-01  2.1471E+05  1.7348E+00  3.0096E+00  1.7327E+00  8.3854E+00  1.0000E-02  1.1768E+00  1.0000E-02
             1.3544E+01
 PARAMETER:  1.8424E-01 -3.7370E-01  1.2377E+01  6.5086E-01  1.2018E+00  6.4968E-01  2.2265E+00 -6.2518E+00  2.6277E-01 -1.5790E+01
             2.7059E+00
 GRADIENT:  -1.1410E+02 -3.5450E-01  5.3471E-06  3.2236E+01  1.6559E+00  4.2922E+01  1.2139E+02  0.0000E+00 -7.6632E-01  0.0000E+00
             8.0510E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -492.020499004285        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1945             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0838E+00  6.0864E-01  2.4630E+05  1.7441E+00  3.0114E+00  1.7320E+00  8.4320E+00  1.0000E-02  1.1751E+00  1.0000E-02
             1.3558E+01
 PARAMETER:  1.8045E-01 -3.9653E-01  1.2514E+01  6.5623E-01  1.2024E+00  6.4926E-01  2.2320E+00 -6.2518E+00  2.6135E-01 -1.5790E+01
             2.7070E+00
 GRADIENT:  -6.3001E+01  9.7879E-01 -5.4121E-07  2.2717E+01  2.8007E-01  2.4351E+01  1.2330E+02  0.0000E+00 -4.0377E+00  0.0000E+00
             6.2563E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -492.031348194957        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2136             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0882E+00  5.9977E-01  2.9729E+05  1.7491E+00  3.0102E+00  1.7318E+00  8.4679E+00  1.0000E-02  1.1838E+00  1.0000E-02
             1.3547E+01
 PARAMETER:  1.8455E-01 -4.1121E-01  1.2702E+01  6.5910E-01  1.2020E+00  6.4918E-01  2.2363E+00 -6.2518E+00  2.6870E-01 -1.5790E+01
             2.7062E+00
 GRADIENT:  -4.1074E+01  7.4129E-01 -3.7779E-06  2.1908E+01  3.1771E-01  1.3419E+01  1.2580E+02  0.0000E+00  2.6061E-01  0.0000E+00
             5.1758E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -492.032998132829        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2327
 NPARAMETR:  1.0883E+00  5.9606E-01  4.4260E+05  1.7517E+00  3.0168E+00  1.7319E+00  8.4924E+00  1.0000E-02  1.1840E+00  1.0000E-02
             1.3555E+01
 PARAMETER:  1.8462E-01 -4.1741E-01  1.3100E+01  6.6057E-01  1.2042E+00  6.4919E-01  2.2392E+00 -6.2518E+00  2.6888E-01 -1.5790E+01
             2.7068E+00
 GRADIENT:  -3.0892E+01 -3.4942E+00 -3.2049E-06  1.1980E-02  7.4223E-01  2.1047E+00  1.4740E+01  0.0000E+00  2.1876E+00  0.0000E+00
             2.1071E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -492.035046101117        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:     2528             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0885E+00  5.9270E-01  1.4027E+05  1.7556E+00  3.0134E+00  1.7314E+00  8.5063E+00  1.0000E-02  1.1803E+00  1.0000E-02
             1.3553E+01
 PARAMETER:  1.8484E-01 -4.2307E-01  1.1951E+01  6.6282E-01  1.2031E+00  6.4896E-01  2.2408E+00 -6.2518E+00  2.6581E-01 -1.5790E+01
             2.7066E+00
 GRADIENT:  -2.0856E+01  1.4004E+00 -4.7648E-06  1.8061E+01  1.9261E-01  1.2489E+01  1.2628E+02  0.0000E+00 -1.5240E+00  0.0000E+00
             4.1632E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -492.035413641103        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2718
 NPARAMETR:  1.0886E+00  5.9071E-01  1.6893E+05  1.7567E+00  3.0137E+00  1.7314E+00  8.5116E+00  1.0000E-02  1.1816E+00  1.0000E-02
             1.3553E+01
 PARAMETER:  1.8487E-01 -4.2643E-01  1.2137E+01  6.6344E-01  1.2032E+00  6.4892E-01  2.2414E+00 -6.2518E+00  2.6686E-01 -1.5790E+01
             2.7066E+00
 GRADIENT:  -7.0563E+01 -2.7227E+00 -1.0212E-07  1.6215E+01 -5.6346E-01  2.3587E+01  1.2873E+01  0.0000E+00 -3.4270E+00  0.0000E+00
             2.7658E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -492.036105814162        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     2918             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0886E+00  5.8852E-01  1.5674E+05  1.7578E+00  3.0143E+00  1.7312E+00  8.5216E+00  1.0000E-02  1.1845E+00  1.0000E-02
             1.3554E+01
 PARAMETER:  1.8488E-01 -4.3014E-01  1.2062E+01  6.6406E-01  1.2034E+00  6.4881E-01  2.2426E+00 -6.2518E+00  2.6934E-01 -1.5790E+01
             2.7067E+00
 GRADIENT:  -4.8234E+01  1.7919E+00 -5.7821E-07  2.5295E+01 -2.7697E-02  3.0679E+01  1.2598E+02  0.0000E+00 -6.9204E+00  0.0000E+00
             6.4385E+01

0ITERATION NO.:   92    OBJECTIVE VALUE:  -492.036105814162        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     2984
 NPARAMETR:  1.0889E+00  5.8755E-01  1.7684E+05  1.7607E+00  3.0114E+00  1.7308E+00  8.5391E+00  1.0000E-02  1.1855E+00  1.0000E-02
             1.3554E+01
 PARAMETER:  1.8488E-01 -4.3014E-01  1.2062E+01  6.6406E-01  1.2034E+00  6.4881E-01  2.2426E+00 -6.2518E+00  2.6934E-01 -1.5790E+01
             2.7067E+00
 GRADIENT:  -4.0356E-02  1.1031E-02 -1.9154E-04 -1.8418E-01  1.1034E-02  1.5397E-02 -8.1480E-02  0.0000E+00 -1.5757E-02  0.0000E+00
            -2.0086E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2984
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3116E-02  5.4863E-02  5.7051E-10 -8.1305E-02  9.8018E-07
 SE:             2.7141E-02  2.1653E-02  1.7137E-10  1.4299E-02  2.5517E-05
 N:                     100         100         100         100         100

 P VAL.:         3.9437E-01  1.1285E-02  8.7106E-04  1.3032E-08  9.6936E-01

 ETASHRINKSD(%)  9.0752E+00  2.7460E+01  1.0000E+02  5.2097E+01  9.9915E+01
 ETASHRINKVR(%)  1.7327E+01  4.7380E+01  1.0000E+02  7.7053E+01  1.0000E+02
 EBVSHRINKSD(%)  1.2689E+01  2.3239E+01  1.0000E+02  5.1349E+01  9.9845E+01
 EBVSHRINKVR(%)  2.3767E+01  4.1078E+01  1.0000E+02  7.6331E+01  1.0000E+02
 RELATIVEINF(%)  7.5380E+01  3.4683E+01  0.0000E+00  1.3393E+01  2.0219E-04
 EPSSHRINKSD(%)  5.0208E+00
 EPSSHRINKVR(%)  9.7895E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -492.03610581416234     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       610.69013403144481     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    80.18
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    12.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -492.036       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.09E+00  5.89E-01  1.57E+05  1.76E+00  3.01E+00  1.73E+00  8.52E+00  1.00E-02  1.18E+00  1.00E-02  1.36E+01
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        4.45E+02
 
 TH 2
+       -5.04E+01  1.24E+01
 
 TH 3
+        1.43E-06 -2.33E-07  7.13E-15
 
 TH 4
+       -1.93E+02  4.10E+01 -5.96E-07  1.71E+02
 
 TH 5
+        9.27E+00 -2.16E+00  3.24E-08 -8.64E+00  4.46E-01
 
 TH 6
+       -3.82E+01  1.56E+01 -3.67E-07  3.87E+01 -2.52E+00  4.43E+01
 
 TH 7
+        5.89E+00 -9.61E-01  7.54E-09 -5.41E+00  2.60E-01 -3.22E-01  2.17E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.64E+01 -1.32E+01  1.54E-07 -6.16E+01  3.04E+00 -8.88E+00  2.16E+00  0.00E+00  2.33E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        6.96E-01 -1.81E+00  8.58E-09 -6.44E+00  3.30E-01 -1.09E+00  1.61E-01  0.00E+00  2.28E+00  0.00E+00  6.14E-01
 
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
+        2.39E+02
 
 TH 2
+       -4.49E+00  3.45E+01
 
 TH 3
+        7.05E-07 -4.50E-07  2.44E-14
 
 TH 4
+       -1.16E+01  1.55E+01 -1.12E-08  8.80E+01
 
 TH 5
+        1.16E-01 -9.08E-01 -2.42E-08 -4.12E+00  3.00E+00
 
 TH 6
+        5.06E+00 -3.02E+00 -2.68E-08  2.61E+00 -9.18E-01  5.11E+01
 
 TH 7
+        5.13E-01  3.44E+00 -4.93E-09 -3.61E+00  9.06E-02 -2.37E-01  1.33E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.40E+00 -1.93E+00 -3.97E-07 -2.79E+01  1.54E+00 -2.74E+00  9.52E-01  0.00E+00  3.04E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.21E+00 -1.34E+00  2.50E-09 -6.72E+00  2.96E-01  6.71E-01  1.87E-01  0.00E+00  2.42E+00  0.00E+00  3.35E+00
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.44E+02
 
 TH 2
+        4.28E+01  2.87E+01
 
 TH 3
+        1.57E-09  9.22E-10  1.78E-18
 
 TH 4
+        8.77E+01  1.80E+01  3.43E-10  8.81E+01
 
 TH 5
+       -4.24E+00  3.61E-01 -2.40E-10 -3.81E+00  9.21E-01
 
 TH 6
+       -8.75E+00  2.88E+00 -1.70E-12 -1.82E+01  2.37E-01  4.84E+01
 
 TH 7
+       -3.42E+00  3.27E+00  2.41E-11 -3.91E+00  4.65E-01  6.30E-01  1.37E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.16E+01 -1.47E+00  1.68E-11 -3.07E+01  1.30E+00  1.90E+01  9.37E-01  0.00E+00  2.72E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.02E+01 -6.69E+00  1.91E-09 -2.05E+01  1.08E+00  1.35E+01  7.21E-01  0.00E+00  8.66E+00  0.00E+00  3.87E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       92.552
Stop Time:
Thu Sep 30 10:09:13 CDT 2021
