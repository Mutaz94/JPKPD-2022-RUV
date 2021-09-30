Thu Sep 30 03:29:24 CDT 2021
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
$DATA ../../../../data/spa1/D/dat69.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17847.9979952951        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6408E+02  3.6804E+02 -9.7538E+01  1.9826E+02  3.6581E+02 -2.0161E+03 -9.5822E+02 -3.4484E+01 -1.7415E+03 -6.7581E+02
            -3.4059E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -644.380564794216        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3832E+00  8.8079E-01  8.2551E-01  2.4913E+00  1.3982E+00  3.6668E+00  1.5003E+00  9.1471E-01  3.3153E+00  1.1156E+00
             1.1820E+01
 PARAMETER:  4.2443E-01 -2.6933E-02 -9.1752E-02  1.0128E+00  4.3518E-01  1.3993E+00  5.0564E-01  1.0851E-02  1.2985E+00  2.0935E-01
             2.5698E+00
 GRADIENT:   1.3744E+01  1.4438E+01 -7.2836E+01  9.6328E+01  9.9879E+00  1.4313E+02  5.5244E+00  7.4508E+00  3.8365E+01  3.0362E+00
             1.1625E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -695.483937261029        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3871E+00  9.9604E-01  2.7514E+00  2.3181E+00  8.3396E+00  2.6208E+00  9.2436E+00  1.5922E-01  2.9674E+00  7.5046E+00
             1.0680E+01
 PARAMETER:  4.2723E-01  9.6029E-02  1.1121E+00  9.4076E-01  2.2210E+00  1.0635E+00  2.3239E+00 -1.7374E+00  1.1877E+00  2.1155E+00
             2.4684E+00
 GRADIENT:   2.5859E+01  1.7327E+01 -2.9315E-01  6.2539E+01 -1.0952E-01  3.0083E+01  5.5621E+01  1.5916E-03  5.6008E+01 -1.7264E-01
             1.0874E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -710.971613756635        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.2408E+00  1.0225E+00  1.9971E+00  1.5196E+00  4.9116E+00  2.4027E+00  5.3005E+00  3.5108E-02  1.7582E+00  6.8849E+00
             1.1210E+01
 PARAMETER:  3.1577E-01  1.2225E-01  7.9167E-01  5.1844E-01  1.6916E+00  9.7660E-01  1.7678E+00 -3.2493E+00  6.6431E-01  2.0293E+00
             2.5168E+00
 GRADIENT:  -2.1097E+01  6.2176E+00  5.1439E+00  8.8208E+00 -4.7082E+00  9.5550E+00  2.6622E+01 -7.9912E-05  1.3614E+01  7.1759E+00
             1.5003E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -715.637453309252        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.2482E+00  1.2752E+00  1.6885E+00  1.1546E+00  5.0793E+00  2.3735E+00  3.9273E+00  2.5476E-02  1.6384E+00  6.1204E+00
             1.0004E+01
 PARAMETER:  3.2172E-01  3.4307E-01  6.2384E-01  2.4373E-01  1.7252E+00  9.6438E-01  1.4679E+00 -3.5700E+00  5.9369E-01  1.9116E+00
             2.4030E+00
 GRADIENT:  -1.2642E+00 -9.1222E+00  2.7574E+00 -6.5321E+00 -5.7730E+00  1.1445E+01 -1.0623E+01 -1.9002E-05  5.5177E+00  8.3920E+00
             6.7592E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -716.437904630793        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.2412E+00  1.4296E+00  1.6021E+00  1.0564E+00  5.9346E+00  2.3192E+00  3.7876E+00  1.9728E-02  1.6184E+00  5.9698E+00
             9.5894E+00
 PARAMETER:  3.1607E-01  4.5737E-01  5.7133E-01  1.5485E-01  1.8808E+00  9.4122E-01  1.4317E+00 -3.8257E+00  5.8143E-01  1.8867E+00
             2.3607E+00
 GRADIENT:  -2.2463E-01 -5.0780E+00  8.8709E-01 -4.7657E-01 -1.9150E+00  1.7396E+00 -7.5516E+00  2.2767E-05  7.0481E+00 -4.8147E-01
             3.6407E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -730.983770989597        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  1.3273E+00  1.0625E+00  1.6287E+00  1.2730E+00  7.1769E+00  2.8014E+00  6.0064E+00  1.0000E-02  1.2302E+00  6.5635E+00
             9.3176E+00
 PARAMETER:  3.8318E-01  1.6061E-01  5.8781E-01  3.4139E-01  2.0709E+00  1.1301E+00  1.8928E+00 -5.2907E+00  3.0719E-01  1.9815E+00
             2.3319E+00
 GRADIENT:  -1.2151E+00  5.5920E+00  2.3404E-01  4.5871E+00 -9.6113E-01 -1.3126E+00  7.1332E+00  0.0000E+00  2.6572E+00 -2.7304E-01
            -6.1606E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -732.560616535035        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  1.3382E+00  1.1225E+00  1.1650E+00  1.1131E+00  1.9972E+01  2.8184E+00  5.4664E+00  1.0000E-02  9.0306E-01  9.3069E+00
             9.4583E+00
 PARAMETER:  3.9136E-01  2.1551E-01  2.5273E-01  2.0719E-01  3.0944E+00  1.1362E+00  1.7986E+00 -8.9000E+00 -1.9633E-03  2.3308E+00
             2.3469E+00
 GRADIENT:   3.1463E+00 -7.9388E-01  1.7162E-01  1.4848E+00 -2.4066E-01 -2.5258E-01 -3.0657E+00  0.0000E+00 -8.4888E-01 -1.1144E-01
             5.5309E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -732.793020842314        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.3099E+00  1.1828E+00  1.0585E+00  1.0716E+00  4.4038E+01  2.8308E+00  5.4496E+00  1.0000E-02  8.5048E-01  1.2951E+01
             9.3985E+00
 PARAMETER:  3.6993E-01  2.6785E-01  1.5686E-01  1.6918E-01  3.8851E+00  1.1405E+00  1.7955E+00 -1.0554E+01 -6.1959E-02  2.6612E+00
             2.3405E+00
 GRADIENT:  -1.1595E+00  5.4745E-01  9.1956E-02  2.6004E-01 -1.2224E-01  3.4271E-01  2.3847E-01  0.0000E+00 -6.7874E-01 -4.3954E-02
             1.4027E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -733.007897967811        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  1.3156E+00  1.1739E+00  1.0343E+00  1.0747E+00  1.4829E+03  2.8765E+00  5.3884E+00  1.0000E-02  8.9017E-01  1.4024E+01
             9.3745E+00
 PARAMETER:  3.7429E-01  2.6030E-01  1.3376E-01  1.7204E-01  7.4017E+00  1.1566E+00  1.7843E+00 -1.0554E+01 -1.6341E-02  2.7408E+00
             2.3380E+00
 GRADIENT:   7.8804E-01  2.4475E-01 -1.9803E-01  1.8288E+00 -4.6315E-03  6.0157E+00 -2.5505E+00  0.0000E+00 -3.7808E-01 -5.2512E-05
            -7.8324E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -733.017634647496        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1261             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3135E+00  1.1707E+00  1.0412E+00  1.0727E+00  1.6754E+04  2.8684E+00  5.3815E+00  1.0000E-02  9.0161E-01  1.6015E+01
             9.3720E+00
 PARAMETER:  3.7272E-01  2.5763E-01  1.4040E-01  1.7018E-01  9.8264E+00  1.1538E+00  1.7830E+00 -1.0554E+01 -3.5690E-03  2.8735E+00
             2.3377E+00
 GRADIENT:   2.9992E+01  2.8237E+00  4.0673E-02  3.1526E+00 -4.0225E-04  7.8028E+01  5.3989E+01  0.0000E+00  1.2526E-01 -4.6822E-07
             2.3463E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -733.018997002484        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1451
 NPARAMETR:  1.3135E+00  1.1748E+00  1.0461E+00  1.0711E+00  2.1355E+05  2.8718E+00  5.3706E+00  1.0000E-02  9.0369E-01  1.8109E+01
             9.3708E+00
 PARAMETER:  3.7271E-01  2.6109E-01  1.4503E-01  1.6872E-01  1.2372E+01  1.1549E+00  1.7809E+00 -1.0554E+01 -1.2690E-03  2.9964E+00
             2.3376E+00
 GRADIENT:   3.3272E-01 -1.9954E-01  7.3804E-02 -3.7702E-01 -3.3111E-05  5.5418E+00 -2.2514E+00  0.0000E+00  1.2577E-01 -4.2153E-09
             2.8399E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -733.020021935916        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     1651             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3139E+00  1.1802E+00  1.0401E+00  1.0706E+00  2.6567E+06  2.8719E+00  5.3593E+00  1.0000E-02  8.9770E-01  1.8539E+01
             9.3683E+00
 PARAMETER:  3.7300E-01  2.6567E-01  1.3930E-01  1.6823E-01  1.4893E+01  1.1550E+00  1.7788E+00 -1.0554E+01 -7.9139E-03  3.0199E+00
             2.3373E+00
 GRADIENT:   3.0036E+01  3.1148E+00 -4.0402E-02  3.9673E+00 -2.4441E-06  7.8559E+01  5.3156E+01  0.0000E+00  9.7650E-03 -2.2588E-11
             2.2747E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -733.020558131389        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1839
 NPARAMETR:  1.3140E+00  1.1825E+00  1.0393E+00  1.0695E+00  3.3607E+06  2.8720E+00  5.3556E+00  1.0000E-02  8.9623E-01  1.8818E+01
             9.3682E+00
 PARAMETER:  3.7307E-01  2.6765E-01  1.3854E-01  1.6715E-01  1.5128E+01  1.1550E+00  1.7781E+00 -1.0554E+01 -9.5621E-03  3.0348E+00
             2.3373E+00
 GRADIENT:   4.3247E-01  3.9687E-02 -6.7883E-02  7.5667E-01 -2.0863E-06  5.5315E+00 -2.6272E+00  0.0000E+00 -8.2393E-02 -1.8731E-11
            -6.2652E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -733.021225525564        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2035             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3140E+00  1.1840E+00  1.0440E+00  1.0669E+00  6.9842E+06  2.8722E+00  5.3519E+00  1.0000E-02  8.9760E-01  1.7213E+01
             9.3700E+00
 PARAMETER:  3.7309E-01  2.6890E-01  1.4304E-01  1.6477E-01  1.5859E+01  1.1551E+00  1.7775E+00 -1.0554E+01 -8.0314E-03  2.9457E+00
             2.3375E+00
 GRADIENT:   2.9981E+01  2.9779E+00  9.5284E-02  2.7867E+00 -9.2898E-07  7.8632E+01  5.3490E+01  0.0000E+00  2.0741E-01  0.0000E+00
             2.3527E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -733.021692051164        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2222
 NPARAMETR:  1.3141E+00  1.1859E+00  1.0429E+00  1.0660E+00  1.3000E+07  2.8724E+00  5.3482E+00  1.0000E-02  8.9593E-01  1.8840E+01
             9.3699E+00
 PARAMETER:  3.7315E-01  2.7049E-01  1.4196E-01  1.6396E-01  1.6480E+01  1.1551E+00  1.7768E+00 -1.0554E+01 -9.8933E-03  3.0360E+00
             2.3375E+00
 GRADIENT:   3.6422E-01 -1.4710E-01  5.6358E-02 -2.4361E-01 -5.4170E-07  5.5421E+00 -2.3440E+00  0.0000E+00  9.3833E-02  0.0000E+00
             1.0922E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -733.022309861266        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2416             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3143E+00  1.1892E+00  1.0381E+00  1.0655E+00  1.1236E+11  2.8726E+00  5.3417E+00  1.0000E-02  8.9155E-01  1.6576E+01
             9.3686E+00
 PARAMETER:  3.7328E-01  2.7328E-01  1.3743E-01  1.6346E-01  2.5545E+01  1.1552E+00  1.7756E+00 -1.0554E+01 -1.4791E-02  2.9080E+00
             2.3374E+00
 GRADIENT:   3.0051E+01  3.1851E+00 -7.9222E-03  3.5268E+00 -5.3406E-11  7.8634E+01  5.2978E+01  0.0000E+00  5.5427E-02  0.0000E+00
             2.2941E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -733.022608133717        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2606
 NPARAMETR:  1.3143E+00  1.1903E+00  1.0384E+00  1.0648E+00  1.0499E+11  2.8727E+00  5.3392E+00  1.0000E-02  8.9096E-01  1.6814E+01
             9.3688E+00
 PARAMETER:  3.7332E-01  2.7421E-01  1.3764E-01  1.6276E-01  2.5477E+01  1.1552E+00  1.7751E+00 -1.0554E+01 -1.5456E-02  2.9222E+00
             2.3374E+00
 GRADIENT:   4.1400E-01 -3.1072E-02 -2.2320E-02  3.6024E-01 -6.7381E-11  5.5487E+00 -2.5620E+00  0.0000E+00 -2.4144E-02  0.0000E+00
            -3.5850E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -733.022927014391        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2804             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3144E+00  1.1920E+00  1.0401E+00  1.0633E+00  2.1153E+09  2.8728E+00  5.3358E+00  1.0000E-02  8.9103E-01  1.6075E+01
             9.3696E+00
 PARAMETER:  3.7336E-01  2.7560E-01  1.3930E-01  1.6134E-01  2.1572E+01  1.1553E+00  1.7744E+00 -1.0554E+01 -1.5381E-02  2.8773E+00
             2.3375E+00
 GRADIENT:   3.0026E+01  3.1230E+00  6.2170E-02  2.9199E+00 -2.9712E-09  7.8681E+01  5.3099E+01  0.0000E+00  1.5448E-01  0.0000E+00
             2.3352E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -733.023166684026        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2994
 NPARAMETR:  1.3144E+00  1.1932E+00  1.0391E+00  1.0629E+00  1.7594E+09  2.8729E+00  5.3335E+00  1.0000E-02  8.8975E-01  1.5633E+01
             9.3694E+00
 PARAMETER:  3.7340E-01  2.7665E-01  1.3831E-01  1.6097E-01  2.1388E+01  1.1553E+00  1.7740E+00 -1.0554E+01 -1.6818E-02  2.8494E+00
             2.3375E+00
 GRADIENT:   3.9108E-01 -8.0806E-02  1.6909E-02  4.4826E-02 -3.9982E-09  5.5540E+00 -2.4731E+00  0.0000E+00  3.1355E-02  0.0000E+00
            -1.1790E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -733.023460359192        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3191             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3146E+00  1.1954E+00  1.0369E+00  1.0621E+00  5.7054E+08  2.8730E+00  5.3290E+00  1.0000E-02  8.8729E-01  1.6343E+01
             9.3688E+00
 PARAMETER:  3.7349E-01  2.7846E-01  1.3627E-01  1.6025E-01  2.0262E+01  1.1554E+00  1.7732E+00 -1.0554E+01 -1.9584E-02  2.8938E+00
             2.3374E+00
 GRADIENT:   3.0067E+01  3.2341E+00  1.1407E-02  3.2731E+00 -1.1056E-08  7.8689E+01  5.2812E+01  0.0000E+00  7.7434E-02  0.0000E+00
             2.3052E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -733.024233040379        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     3389             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3148E+00  1.2025E+00  1.0367E+00  1.0585E+00  1.0636E+08  2.8735E+00  5.3158E+00  1.0000E-02  8.8525E-01  1.2857E+01
             9.3694E+00
 PARAMETER:  3.7369E-01  2.8439E-01  1.3604E-01  1.5681E-01  1.8582E+01  1.1555E+00  1.7707E+00 -1.0554E+01 -2.1880E-02  2.6539E+00
             2.3375E+00
 GRADIENT:   3.0062E+01  3.2957E+00  5.2883E-02  2.7927E+00 -5.9754E-08  7.8761E+01  5.2733E+01  0.0000E+00  1.6867E-01  0.0000E+00
             2.3336E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -733.024397737047        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3578
 NPARAMETR:  1.3148E+00  1.2028E+00  1.0362E+00  1.0582E+00  8.9566E+07  2.8735E+00  5.3147E+00  1.0000E-02  8.8307E-01  1.2667E+01
             9.3693E+00
 PARAMETER:  3.7371E-01  2.8467E-01  1.3555E-01  1.5654E-01  1.8410E+01  1.1555E+00  1.7705E+00 -1.0554E+01 -2.4356E-02  2.6390E+00
             2.3374E+00
 GRADIENT:   3.9597E-01 -5.8973E-02  1.6962E-02  1.2848E-02 -7.8104E-08  5.5656E+00 -2.4982E+00  0.0000E+00  2.9135E-02  0.0000E+00
            -1.2595E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -733.024530273675        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     3776             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3149E+00  1.2041E+00  1.0342E+00  1.0577E+00  1.2009E+08  2.8736E+00  5.3123E+00  1.0000E-02  8.8107E-01  1.3548E+01
             9.3688E+00
 PARAMETER:  3.7377E-01  2.8574E-01  1.3360E-01  1.5608E-01  1.8704E+01  1.1556E+00  1.7700E+00 -1.0554E+01 -2.6620E-02  2.7062E+00
             2.3374E+00
 GRADIENT:   3.0093E+01  3.3495E+00  1.5221E-02  3.0919E+00 -5.2785E-08  7.8751E+01  5.2546E+01  0.0000E+00  8.1066E-02  0.0000E+00
             2.3080E+01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -733.024593353022        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3961
 NPARAMETR:  1.3149E+00  1.2050E+00  1.0335E+00  1.0573E+00  1.3134E+08  2.8737E+00  5.3105E+00  1.0000E-02  8.8013E-01  1.3473E+01
             9.3687E+00
 PARAMETER:  3.7379E-01  2.8652E-01  1.3294E-01  1.5569E-01  1.8793E+01  1.1556E+00  1.7697E+00 -1.0554E+01 -2.7690E-02  2.7007E+00
             2.3374E+00
 GRADIENT:   4.2187E-01 -1.2849E-02 -1.9189E-02  2.6075E-01 -5.3202E-08  5.5662E+00 -2.5835E+00  0.0000E+00 -2.4839E-02  0.0000E+00
            -3.4649E-01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -733.024682433847        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     4161             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3149E+00  1.2057E+00  1.0348E+00  1.0565E+00  3.8519E+10  2.8738E+00  5.3095E+00  1.0000E-02  8.8054E-01  1.3091E+01
             9.3694E+00
 PARAMETER:  3.7379E-01  2.8707E-01  1.3420E-01  1.5497E-01  2.4474E+01  1.1556E+00  1.7695E+00 -1.0554E+01 -2.7223E-02  2.6719E+00
             2.3374E+00
 GRADIENT:   3.0074E+01  3.3323E+00  4.6407E-02  2.8127E+00 -1.5979E-10  7.8776E+01  5.2607E+01  0.0000E+00  1.2469E-01  0.0000E+00
             2.3272E+01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -733.024753131342        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     4331
 NPARAMETR:  1.3150E+00  1.2063E+00  1.0339E+00  1.0563E+00  2.7845E+07  2.8738E+00  5.3084E+00  1.0000E-02  8.7971E-01  1.2733E+01
             9.3692E+00
 PARAMETER:  3.7381E-01  2.8754E-01  1.3332E-01  1.5478E-01  1.7242E+01  1.1556E+00  1.7693E+00 -1.0554E+01 -2.8167E-02  2.6442E+00
             2.3374E+00
 GRADIENT:   3.0084E+01  3.3538E+00  3.3055E-02  2.8977E+00 -2.2936E-07  7.8775E+01  5.2549E+01  0.0000E+00  1.0612E-01  8.5990E-12
             2.3202E+01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -733.024775846636        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     4425
 NPARAMETR:  1.3151E+00  1.2091E+00  1.0316E+00  1.0541E+00  1.7638E+06  2.8746E+00  5.3046E+00  1.0000E-02  8.7825E-01  1.2521E+01
             9.3702E+00
 PARAMETER:  3.7395E-01  2.8988E-01  1.3108E-01  1.5273E-01  1.4483E+01  1.1559E+00  1.7686E+00 -1.0554E+01 -2.9826E-02  2.6274E+00
             2.3375E+00
 GRADIENT:   4.1469E-01 -9.0668E-02  2.6394E-02 -3.0621E-01 -3.9696E-06  5.6273E+00 -2.3626E+00  0.0000E+00  7.7749E-02 -3.0288E-11
             1.0520E-01

0ITERATION NO.:  140    OBJECTIVE VALUE:  -733.025029082311        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     4623             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3151E+00  1.2102E+00  1.0310E+00  1.0543E+00  5.0302E+06  2.8742E+00  5.3012E+00  1.0000E-02  8.7581E-01  1.3682E+01
             9.3688E+00
 PARAMETER:  3.7392E-01  2.9080E-01  1.3056E-01  1.5291E-01  1.5531E+01  1.1558E+00  1.7679E+00 -1.0554E+01 -3.2604E-02  2.7161E+00
             2.3374E+00
 GRADIENT:   3.0105E+01  3.4268E+00  1.2042E-02  2.9596E+00 -1.2789E-06  7.8812E+01  5.2382E+01  0.0000E+00  7.5091E-02  0.0000E+00
             2.3091E+01

0ITERATION NO.:  145    OBJECTIVE VALUE:  -733.025067795042        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     4817
 NPARAMETR:  1.3151E+00  1.2108E+00  1.0314E+00  1.0541E+00  8.0754E+06  2.8742E+00  5.2997E+00  1.0000E-02  8.7582E-01  1.4430E+01
             9.3687E+00
 PARAMETER:  3.7394E-01  2.9132E-01  1.3088E-01  1.5271E-01  1.6004E+01  1.1558E+00  1.7677E+00 -1.0554E+01 -3.2590E-02  2.7693E+00
             2.3374E+00
 GRADIENT:   4.1559E-01 -1.8366E-02 -1.0801E-02  1.3860E-01 -8.6326E-07  5.5756E+00 -2.5538E+00  0.0000E+00 -1.2233E-02  0.0000E+00
            -2.8588E-01

0ITERATION NO.:  146    OBJECTIVE VALUE:  -733.025067795042        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     4847
 NPARAMETR:  1.3151E+00  1.2108E+00  1.0326E+00  1.0537E+00  9.4770E+06  2.8742E+00  5.2999E+00  1.0000E-02  8.7640E-01  1.4036E+01
             9.3695E+00
 PARAMETER:  3.7394E-01  2.9132E-01  1.3088E-01  1.5271E-01  1.6004E+01  1.1558E+00  1.7677E+00 -1.0554E+01 -3.2590E-02  2.7693E+00
             2.3374E+00
 GRADIENT:   4.4031E-03 -1.1694E-04 -7.1453E-03  1.4543E-01 -2.0929E-06 -4.7297E-04 -5.3801E-03  0.0000E+00 -6.9882E-03  6.0663E-06
            -7.7329E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4847
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1178E-02  2.6380E-02  9.5629E-05 -6.1348E-02 -3.6365E-11
 SE:             2.9148E-02  2.4063E-02  1.1173E-05  1.1143E-02  4.1773E-09
 N:                     100         100         100         100         100

 P VAL.:         7.0135E-01  2.7296E-01  1.1570E-17  3.6941E-08  9.9305E-01

 ETASHRINKSD(%)  2.3522E+00  1.9385E+01  9.9963E+01  6.2668E+01  1.0000E+02
 ETASHRINKVR(%)  4.6490E+00  3.5012E+01  1.0000E+02  8.6064E+01  1.0000E+02
 EBVSHRINKSD(%)  2.3221E+00  1.4591E+01  9.9949E+01  6.8103E+01  1.0000E+02
 EBVSHRINKVR(%)  4.5903E+00  2.7054E+01  1.0000E+02  8.9826E+01  1.0000E+02
 RELATIVEINF(%)  3.3741E+01  2.8186E+01  7.8033E-06  1.5454E+00  0.0000E+00
 EPSSHRINKSD(%)  1.0856E+01
 EPSSHRINKVR(%)  2.0534E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -733.02506779504154     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       185.91346540963116     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   115.70
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    10.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -733.025       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.32E+00  1.21E+00  1.03E+00  1.05E+00  8.08E+06  2.87E+00  5.30E+00  1.00E-02  8.76E-01  1.44E+01  9.37E+00
 


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
+        4.27E+01
 
 TH 2
+       -9.55E-01  5.53E+00
 
 TH 3
+        3.50E-02 -2.96E+00  1.59E+00
 
 TH 4
+       -1.83E+01  3.07E+01 -1.63E+01  1.75E+02
 
 TH 5
+        5.21E-12 -5.53E-12  2.91E-12 -3.20E-11  5.96E-24
 
 TH 6
+        3.69E+00  8.31E-02 -7.53E-02 -7.69E-01  2.58E-13  4.21E-01
 
 TH 7
+        2.75E+00 -2.12E+00  1.11E+00 -1.25E+01  2.35E-12  1.82E-01  9.43E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.86E+00 -6.01E+00  3.17E+00 -3.49E+01  6.47E-12  3.82E-01  2.57E+00  0.00E+00  7.10E+00
 
 TH10
+       -2.46E-05  7.01E-06 -3.52E-06  4.63E-05 -9.28E-18 -2.16E-06 -4.01E-06  0.00E+00 -1.04E-05  2.23E-11
 
 TH11
+       -4.71E-01 -1.53E+00  8.45E-01 -8.43E+00  1.40E-12  8.27E-02  5.55E-01  0.00E+00  1.67E+00 -1.92E-06  7.29E-01
 
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
+        7.20E+01
 
 TH 2
+       -2.54E+00  1.99E+01
 
 TH 3
+       -3.23E+00 -3.16E+00  5.59E+00
 
 TH 4
+        2.25E+00  3.28E+01 -1.75E+01  1.79E+02
 
 TH 5
+       -3.36E-13  3.20E-12 -6.64E-13 -3.41E-11 -4.62E-21
 
 TH 6
+       -1.26E+00 -8.66E-01  1.90E-01  3.19E+00  7.87E-13  1.96E+01
 
 TH 7
+        3.72E-04  2.61E+00  1.20E+00 -1.21E+01 -4.39E-13 -7.85E-01  3.57E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.77E-01 -1.93E+00  1.66E+00 -3.36E+01 -3.19E-12  2.57E-01  2.30E+00  0.00E+00  1.38E+01
 
 TH10
+       -1.36E-05  1.81E-05 -1.45E-05 -9.07E-06  5.25E-15 -2.41E-06 -9.01E-08  0.00E+00 -2.24E-04 -4.39E-08
 
 TH11
+       -3.12E+00 -2.06E+00  1.30E+00 -1.08E+01  1.07E-13  8.35E-01  5.18E-01  0.00E+00  2.44E+00  2.15E-07  5.74E+00
 
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
+        7.44E+01
 
 TH 2
+        3.36E+01  2.00E+01
 
 TH 3
+       -4.05E+00 -2.41E+00  2.37E+00
 
 TH 4
+        5.34E+01  3.37E+01 -1.95E+01  1.99E+02
 
 TH 5
+        3.26E-13  1.73E-13 -4.42E-14  4.80E-13  2.78E-27
 
 TH 6
+        1.53E+01  7.77E+00  1.08E+00 -1.27E+01  5.56E-14  1.84E+01
 
 TH 7
+        3.46E+00  1.69E+00  1.44E+00 -1.33E+01 -1.26E-14  4.44E+00  2.57E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.88E+00 -2.65E+00  3.26E+00 -3.72E+01 -5.24E-14  3.82E+00  2.40E+00  0.00E+00  1.43E+01
 
 TH10
+       -1.14E-08  1.88E-08 -1.19E-08  1.20E-07 -1.63E-22 -1.24E-08 -7.11E-09  0.00E+00  5.18E-09  1.23E-14
 
 TH11
+       -1.98E+01 -8.14E+00  2.44E+00 -4.15E+01 -1.44E-13  7.07E+00  1.52E+00  0.00E+00  9.89E+00 -2.43E-07  9.56E+01
 
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
 #CPUT: Total CPU Time in Seconds,      126.463
Stop Time:
Thu Sep 30 03:31:32 CDT 2021
