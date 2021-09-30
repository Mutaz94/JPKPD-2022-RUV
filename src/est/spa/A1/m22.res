Wed Sep 29 11:59:20 CDT 2021
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
$DATA ../../../../data/spa/A1/dat22.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m22.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1129.68316942528        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9701E+02  5.3565E+01 -1.1576E+01  1.0343E+02  1.2745E+02  7.9801E+01 -3.5066E+01 -7.2948E+00 -6.0458E+01 -3.1601E+01
            -8.9683E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1400.82521449816        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1114E+00  1.0420E+00  1.0943E+00  1.0313E+00  9.4759E-01  9.9734E-01  9.8628E-01  9.3014E-01  9.9541E-01  8.6729E-01
             2.6256E+00
 PARAMETER:  2.0561E-01  1.4111E-01  1.9009E-01  1.3082E-01  4.6169E-02  9.7335E-02  8.6184E-02  2.7575E-02  9.5402E-02 -4.2381E-02
             1.0653E+00
 GRADIENT:   3.9670E+02  5.1910E+01  1.9145E+01  4.6681E+01 -4.5954E+01  2.6827E+01  2.9010E+00  3.2775E+00  2.6041E+00  2.8078E+00
             4.2420E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1415.87073040544        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0402E+00  1.0194E+00  6.2996E-01  9.9702E-01  7.8071E-01  9.0678E-01  6.7763E-01  2.0946E-01  1.2447E+00  6.2879E-01
             2.4522E+00
 PARAMETER:  1.3941E-01  1.1917E-01 -3.6210E-01  9.7019E-02 -1.4756E-01  2.1472E-03 -2.8916E-01 -1.4632E+00  3.1887E-01 -3.6395E-01
             9.9699E-01
 GRADIENT:   2.5483E+02 -2.1678E+01 -4.6322E+01  4.6565E+01  7.6316E+01  2.3470E+01 -7.2393E+00  2.5516E-01  2.8207E+01 -7.7383E+00
            -7.5788E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1429.78476415370        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.4740E-01  7.7011E-01  5.8279E-01  1.0878E+00  6.1661E-01  8.0853E-01  1.2168E+00  1.5536E-01  8.9825E-01  5.5405E-01
             2.3403E+00
 PARAMETER:  4.5963E-02 -1.6122E-01 -4.3993E-01  1.8415E-01 -3.8352E-01 -1.1254E-01  2.9624E-01 -1.7620E+00 -7.3078E-03 -4.9049E-01
             9.5030E-01
 GRADIENT:  -2.1342E+01  1.3398E+01  5.4470E+00  1.8394E+01  2.6832E+00 -5.4991E+00 -1.2029E+00  3.0804E-01 -2.3032E+00 -3.0300E-01
            -1.1455E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1432.61867703912        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  9.7915E-01  6.9924E-01  7.9499E-01  1.1752E+00  7.1547E-01  8.2760E-01  1.3123E+00  1.1205E-01  8.9100E-01  7.1850E-01
             2.4017E+00
 PARAMETER:  7.8930E-02 -2.5776E-01 -1.2942E-01  2.6147E-01 -2.3481E-01 -8.9223E-02  3.7178E-01 -2.0888E+00 -1.5416E-02 -2.3059E-01
             9.7617E-01
 GRADIENT:   2.0012E+01  9.2824E+00  1.0555E+00  7.4428E+00 -4.7429E+00  7.4338E-01  5.6093E-01  1.0044E-01 -1.7911E+00  1.6088E+00
             2.6148E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1433.97781213035        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  9.6372E-01  4.6211E-01  9.8262E-01  1.3309E+00  7.3368E-01  8.2197E-01  1.4800E+00  2.0178E-02  8.5436E-01  7.9728E-01
             2.4002E+00
 PARAMETER:  6.3049E-02 -6.7196E-01  8.2468E-02  3.8586E-01 -2.0968E-01 -9.6057E-02  4.9206E-01 -3.8032E+00 -5.7404E-02 -1.2655E-01
             9.7554E-01
 GRADIENT:  -1.3810E+01  7.0576E+00  6.5854E+00  1.1156E+01 -1.0235E+01 -2.3072E-01  1.1979E-02  2.3258E-03 -1.3120E+00  2.2140E+00
            -6.2031E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1435.34738774843        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      759
 NPARAMETR:  9.6453E-01  2.0894E-01  8.6010E-01  1.4542E+00  6.1433E-01  8.2107E-01  1.3976E+00  1.0000E-02  8.2573E-01  7.1427E-01
             2.4149E+00
 PARAMETER:  6.3888E-02 -1.4657E+00 -5.0703E-02  4.7444E-01 -3.8722E-01 -9.7152E-02  4.3476E-01 -8.3645E+00 -9.1486E-02 -2.3649E-01
             9.8166E-01
 GRADIENT:   8.1543E-01  2.6023E+00  3.2894E+00  9.2463E+00 -8.3555E+00  7.7884E-03 -5.7710E-01  0.0000E+00 -9.5108E-01  3.1154E-01
             1.7154E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1436.15244819535        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  9.6057E-01  5.8075E-02  9.2534E-01  1.5426E+00  6.1457E-01  8.1928E-01  1.3210E+00  1.0000E-02  7.8482E-01  7.1534E-01
             2.4192E+00
 PARAMETER:  5.9772E-02 -2.7460E+00  2.2405E-02  5.3345E-01 -3.8684E-01 -9.9327E-02  3.7838E-01 -1.6797E+01 -1.4230E-01 -2.3499E-01
             9.8344E-01
 GRADIENT:   1.3098E+00  5.2106E-01  3.4197E+00  3.1615E+00 -3.8295E+00  3.0329E-01 -3.1289E-02  0.0000E+00 -1.1950E+00 -6.0838E-01
            -4.0352E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1436.44619101418        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1110
 NPARAMETR:  9.5923E-01  1.2796E-02  8.5602E-01  1.5523E+00  5.7510E-01  8.1816E-01  1.1266E+00  1.0000E-02  7.8191E-01  6.9404E-01
             2.4153E+00
 PARAMETER:  5.8379E-02 -4.2587E+00 -5.5467E-02  5.3972E-01 -4.5321E-01 -1.0070E-01  2.1922E-01 -2.7475E+01 -1.4601E-01 -2.6523E-01
             9.8183E-01
 GRADIENT:   9.0947E-01  8.1261E-02  9.6274E-01 -1.6720E+00 -1.7353E+00  5.5331E-02 -1.0436E-03  0.0000E+00 -1.4479E-01  4.9062E-03
            -7.0739E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1436.46142058542        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1302             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5890E-01  1.0000E-02  8.6609E-01  1.5554E+00  5.7972E-01  8.1794E-01  1.1242E+00  1.0000E-02  7.8067E-01  6.9759E-01
             2.4155E+00
 PARAMETER:  5.8034E-02 -4.6157E+00 -4.3772E-02  5.4175E-01 -4.4520E-01 -1.0097E-01  2.1708E-01 -3.0009E+01 -1.4760E-01 -2.6013E-01
             9.8191E-01
 GRADIENT:   6.4126E+01  0.0000E+00  9.2748E-01  1.4380E+02  1.1139E+01  4.2197E+00 -3.9116E-04  0.0000E+00  3.8088E+00  5.8512E-01
             7.9020E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1436.46147411214        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     1497
 NPARAMETR:  9.5890E-01  1.0000E-02  8.6591E-01  1.5554E+00  5.7965E-01  8.1794E-01  1.1425E+00  1.0000E-02  7.8071E-01  6.9768E-01
             2.4152E+00
 PARAMETER:  5.8029E-02 -4.6157E+00 -4.3975E-02  5.4173E-01 -4.4534E-01 -1.0097E-01  2.3320E-01 -3.0009E+01 -1.4756E-01 -2.5999E-01
             9.8177E-01
 GRADIENT:   2.0197E-01  0.0000E+00  1.7248E-01 -2.4454E+00 -6.0927E-02  8.7920E-03 -6.1637E-04  0.0000E+00  1.8190E-02  5.7011E-03
            -2.0557E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1436.46153224188        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     1697             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5890E-01  1.0000E-02  8.6557E-01  1.5553E+00  5.7958E-01  8.1794E-01  1.1580E+00  1.0000E-02  7.8072E-01  6.9760E-01
             2.4151E+00
 PARAMETER:  5.8030E-02 -4.6157E+00 -4.4372E-02  5.4169E-01 -4.4546E-01 -1.0097E-01  2.4669E-01 -3.0009E+01 -1.4753E-01 -2.6011E-01
             9.8173E-01
 GRADIENT:   6.4150E+01  0.0000E+00  7.3930E-01  1.4380E+02  1.1416E+01  4.2172E+00 -3.8762E-04  0.0000E+00  3.8060E+00  5.8586E-01
             7.8107E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1436.46157761729        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1889
 NPARAMETR:  9.5890E-01  1.0000E-02  8.6541E-01  1.5553E+00  5.7950E-01  8.1794E-01  1.1858E+00  1.0000E-02  7.8074E-01  6.9755E-01
             2.4151E+00
 PARAMETER:  5.8030E-02 -4.6157E+00 -4.4552E-02  5.4167E-01 -4.4559E-01 -1.0097E-01  2.7038E-01 -3.0009E+01 -1.4751E-01 -2.6019E-01
             9.8173E-01
 GRADIENT:   2.1497E-01  0.0000E+00 -2.8753E-03 -2.4638E+00  1.8275E-01  6.8362E-03 -6.6312E-04  0.0000E+00  1.4112E-02  2.7000E-03
            -4.2931E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1436.46163035879        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     2088             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5890E-01  1.0000E-02  8.6526E-01  1.5553E+00  5.7936E-01  8.1794E-01  1.2017E+00  1.0000E-02  7.8077E-01  6.9747E-01
             2.4151E+00
 PARAMETER:  5.8029E-02 -4.6157E+00 -4.4730E-02  5.4165E-01 -4.4584E-01 -1.0096E-01  2.8377E-01 -3.0009E+01 -1.4748E-01 -2.6029E-01
             9.8173E-01
 GRADIENT:   6.4138E+01  0.0000E+00  8.8608E-01  1.4379E+02  1.1203E+01  4.2188E+00 -3.7848E-04  0.0000E+00  3.8097E+00  5.8934E-01
             7.8173E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1436.46168946260        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2281
 NPARAMETR:  9.5890E-01  1.0000E-02  8.6510E-01  1.5552E+00  5.7922E-01  8.1795E-01  1.2691E+00  1.0000E-02  7.8079E-01  6.9741E-01
             2.4151E+00
 PARAMETER:  5.8028E-02 -4.6157E+00 -4.4914E-02  5.4162E-01 -4.4607E-01 -1.0096E-01  3.3827E-01 -3.0009E+01 -1.4745E-01 -2.6039E-01
             9.8174E-01
 GRADIENT:   1.9651E-01  0.0000E+00  2.5761E-01 -2.4366E+00 -2.1646E-01  1.0078E-02 -7.6226E-04  0.0000E+00  2.0977E-02  8.2656E-03
            -3.0973E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1436.46181974378        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2474
 NPARAMETR:  9.5890E-01  1.0000E-02  8.6493E-01  1.5552E+00  5.7917E-01  8.1795E-01  1.4217E+00  1.0000E-02  7.8080E-01  6.9735E-01
             2.4151E+00
 PARAMETER:  5.8034E-02 -4.6157E+00 -4.5104E-02  5.4160E-01 -4.4615E-01 -1.0095E-01  4.5187E-01 -3.0009E+01 -1.4743E-01 -2.6047E-01
             9.8174E-01
 GRADIENT:   2.1762E-01  0.0000E+00  2.0394E-01 -2.4516E+00 -1.4030E-01  1.2615E-02 -9.5350E-04  0.0000E+00  2.3024E-02  7.7547E-03
            -2.8222E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1436.46216977265        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2668
 NPARAMETR:  9.5891E-01  1.0000E-02  8.6478E-01  1.5552E+00  5.7906E-01  8.1796E-01  1.8564E+00  1.0000E-02  7.8082E-01  6.9730E-01
             2.4151E+00
 PARAMETER:  5.8039E-02 -4.6157E+00 -4.5275E-02  5.4158E-01 -4.4635E-01 -1.0094E-01  7.1865E-01 -3.0009E+01 -1.4742E-01 -2.6055E-01
             9.8175E-01
 GRADIENT:   2.2868E-01  0.0000E+00  2.9060E-01 -2.4551E+00 -2.7190E-01  1.6640E-02 -1.6179E-03  0.0000E+00  3.0572E-02  1.1303E-02
            -1.6353E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1436.46326226706        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2846
 NPARAMETR:  9.5887E-01  1.0000E-02  8.6493E-01  1.5577E+00  5.7942E-01  8.1800E-01  4.9113E+00  1.0000E-02  7.8064E-01  6.9726E-01
             2.4157E+00
 PARAMETER:  5.8000E-02 -4.6157E+00 -4.5106E-02  5.4318E-01 -4.4573E-01 -1.0089E-01  1.6915E+00 -3.0009E+01 -1.4765E-01 -2.6060E-01
             9.8200E-01
 GRADIENT:  -1.9357E-01  0.0000E+00 -4.7109E-01  2.1089E+00  1.7160E-01 -5.2687E-03 -1.0548E-02  0.0000E+00  4.6900E-02 -1.0692E-02
             3.8120E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1436.47354786147        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     3010
 NPARAMETR:  9.5885E-01  1.0000E-02  8.6516E-01  1.5588E+00  5.7958E-01  8.1798E-01  1.0589E+01  1.0000E-02  7.8011E-01  6.9709E-01
             2.4158E+00
 PARAMETER:  5.7983E-02 -4.6157E+00 -4.4841E-02  5.4392E-01 -4.4545E-01 -1.0091E-01  2.4598E+00 -3.0009E+01 -1.4831E-01 -2.6084E-01
             9.8202E-01
 GRADIENT:  -3.7136E-01  0.0000E+00 -6.4287E-01  4.1002E+00  2.2329E-01 -9.3758E-03  5.6748E-03  0.0000E+00  1.6304E-01 -6.9101E-03
             8.0591E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3010
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.1254E-04  2.0410E-04  8.0722E-05 -1.3105E-02 -1.5939E-02
 SE:             2.8960E-02  1.8089E-03  1.7501E-04  2.7122E-02  1.8702E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9139E-01  9.1016E-01  6.4462E-01  6.2897E-01  3.9408E-01

 ETASHRINKSD(%)  2.9797E+00  9.3940E+01  9.9414E+01  9.1378E+00  3.7346E+01
 ETASHRINKVR(%)  5.8707E+00  9.9633E+01  9.9997E+01  1.7441E+01  6.0744E+01
 EBVSHRINKSD(%)  2.9894E+00  9.4349E+01  9.9387E+01  9.0090E+00  3.7269E+01
 EBVSHRINKVR(%)  5.8894E+00  9.9681E+01  9.9996E+01  1.7206E+01  6.0649E+01
 RELATIVEINF(%)  8.1164E+01  7.9461E-03  1.7262E-04  3.2709E+00  1.3047E+00
 EPSSHRINKSD(%)  3.0782E+01
 EPSSHRINKVR(%)  5.2089E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1436.4735478614707     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -701.32272129773253     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.46
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1436.474       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.59E-01  1.00E-02  8.65E-01  1.56E+00  5.80E-01  8.18E-01  1.06E+01  1.00E-02  7.80E-01  6.97E-01  2.42E+00
 


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
+        1.70E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.28E-01  0.00E+00  5.48E+02
 
 TH 4
+       -6.54E+01  0.00E+00 -4.36E+01  6.22E+02
 
 TH 5
+        4.39E+01  0.00E+00 -1.14E+03 -1.61E+02  2.58E+03
 
 TH 6
+       -2.04E+00  0.00E+00  5.97E+00 -1.28E+01 -2.40E+00  2.63E+02
 
 TH 7
+        3.08E-03  0.00E+00 -4.79E-03 -2.61E-02  2.53E-02  6.34E-04  1.98E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.37E+01  0.00E+00  1.76E+01 -1.32E+01  9.40E+00  7.60E+00  1.78E-02  0.00E+00  2.20E+02
 
 TH10
+       -1.07E+01  0.00E+00 -7.27E+00 -7.23E+00 -3.77E+01 -1.67E+00 -3.80E-03  0.00E+00  2.13E+00  7.06E+01
 
 TH11
+       -1.85E+01  0.00E+00 -2.23E+00 -8.90E+00 -1.08E+01  3.53E+00 -2.92E-04  0.00E+00  1.07E+01  2.50E+01  5.32E+01
 
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
 #CPUT: Total CPU Time in Seconds,       46.716
Stop Time:
Wed Sep 29 12:00:08 CDT 2021
