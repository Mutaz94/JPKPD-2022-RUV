Sat Sep 25 09:30:00 CDT 2021
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
$DATA ../../../../data/spa/A3/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   10.1272602380245        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6001E+02  4.1588E+01  4.4831E+01 -1.8396E+01  2.4856E+02  4.7515E+00 -9.1829E+01 -1.8988E+01 -1.7393E+02 -2.1206E+02
            -2.8078E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1210.00662429001        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0074E+00  9.4658E-01  9.8659E-01  1.0993E+00  9.2456E-01  9.0174E-01  1.0338E+00  9.5795E-01  1.1890E+00  9.6451E-01
             4.3901E+00
 PARAMETER:  1.0740E-01  4.5099E-02  8.6497E-02  1.9471E-01  2.1564E-02 -3.4300E-03  1.3326E-01  5.7045E-02  2.7308E-01  6.3860E-02
             1.5794E+00
 GRADIENT:   3.0820E+01 -2.8539E+01 -2.2991E+01 -2.2837E+01  2.1277E+01 -7.7369E+00  8.2799E+00  6.0312E+00  1.3517E+01  2.1331E+01
             8.6241E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1220.16121065781        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  1.0028E+00  6.7539E-01  5.8752E-01  1.2729E+00  5.4833E-01  9.5008E-01  9.3334E-01  3.5545E-01  1.1976E+00  4.9096E-01
             4.3064E+00
 PARAMETER:  1.0285E-01 -2.9246E-01 -4.3185E-01  3.4131E-01 -5.0087E-01  4.8792E-02  3.1012E-02 -9.3437E-01  2.8036E-01 -6.1140E-01
             1.5601E+00
 GRADIENT:   3.6027E+00  2.1590E+01 -3.0138E+00  4.6745E+01 -1.6354E+01  3.6842E+00  2.7998E+00  1.9666E+00  1.3433E+01  8.3113E+00
             8.6028E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1230.55496388959        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  9.8806E-01  4.8443E-01  5.6268E-01  1.3213E+00  4.8240E-01  9.4060E-01  8.0416E-01  7.9997E-02  1.1297E+00  3.5182E-01
             3.7938E+00
 PARAMETER:  8.7991E-02 -6.2479E-01 -4.7505E-01  3.7861E-01 -6.2899E-01  3.8766E-02 -1.1795E-01 -2.4258E+00  2.2193E-01 -9.4464E-01
             1.4334E+00
 GRADIENT:   3.8162E+00  2.0268E+01  2.0332E+01  3.4978E+01 -3.4481E+01  3.5710E-02 -9.8943E-01  5.6177E-02  4.8370E-01 -1.1971E-01
            -6.7584E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1233.79707003055        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  9.7250E-01  2.3122E-01  3.7665E-01  1.3412E+00  3.2800E-01  9.2789E-01  1.0393E+00  1.3480E-02  1.1641E+00  1.0185E-01
             3.7460E+00
 PARAMETER:  7.2112E-02 -1.3644E+00 -8.7643E-01  3.9353E-01 -1.0147E+00  2.5159E-02  1.3854E-01 -4.2066E+00  2.5198E-01 -2.1842E+00
             1.4207E+00
 GRADIENT:  -2.1321E+01  1.1072E+01  2.3666E+01  4.2170E+01 -4.3310E+01 -5.7800E+00 -1.3902E+00 -2.9785E-03 -2.9109E+00 -1.0076E+00
            -8.1649E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1236.39703660247        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  9.7561E-01  1.2608E-01  3.4188E-01  1.3114E+00  3.0021E-01  9.3825E-01  2.4894E+00  1.0000E-02  1.1864E+00  2.9946E-02
             3.7554E+00
 PARAMETER:  7.5310E-02 -1.9708E+00 -9.7330E-01  3.7113E-01 -1.1033E+00  3.6257E-02  1.0120E+00 -5.9843E+00  2.7095E-01 -3.4084E+00
             1.4232E+00
 GRADIENT:   4.5382E+00  4.1680E-01 -1.0226E+01 -2.5311E+00  1.7507E+01  9.5326E-03 -9.0270E-01  0.0000E+00  6.3937E+00 -7.1040E-02
             2.8119E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1237.38435871455        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  9.6911E-01  5.8247E-02  3.4237E-01  1.3402E+00  2.9281E-01  9.3583E-01  4.6143E+00  1.0000E-02  1.1363E+00  1.0000E-02
             3.7343E+00
 PARAMETER:  6.8620E-02 -2.7431E+00 -9.7186E-01  3.9278E-01 -1.1282E+00  3.3677E-02  1.6292E+00 -8.4300E+00  2.2782E-01 -4.6468E+00
             1.4176E+00
 GRADIENT:  -1.9686E+00  1.6486E+00 -7.0191E-01  4.5014E+00 -1.8258E+00 -8.2687E-01  1.5255E+00  0.0000E+00 -1.4798E+00  0.0000E+00
            -3.1771E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1237.44981659625        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.6855E-01  4.7451E-02  3.3990E-01  1.3389E+00  2.9090E-01  9.3595E-01  5.1845E+00  1.0000E-02  1.1383E+00  1.0000E-02
             3.7332E+00
 PARAMETER:  6.8041E-02 -2.9480E+00 -9.7909E-01  3.9184E-01 -1.1348E+00  3.3802E-02  1.7457E+00 -9.1301E+00  2.2957E-01 -5.0007E+00
             1.4173E+00
 GRADIENT:   1.4758E+00  1.0845E+01 -4.1601E+00 -7.4391E+00  8.6609E-01 -4.4505E+00  1.8942E+01  0.0000E+00 -1.6293E+00  0.0000E+00
            -1.4444E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1237.46405355580        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      598
 NPARAMETR:  9.6869E-01  4.3121E-02  3.3816E-01  1.3363E+00  2.8939E-01  9.3658E-01  5.4555E+00  1.0000E-02  1.1434E+00  1.0000E-02
             3.7365E+00
 PARAMETER:  6.8193E-02 -3.0437E+00 -9.8424E-01  3.8990E-01 -1.1400E+00  3.4478E-02  1.7966E+00 -9.4565E+00  2.3397E-01 -5.1732E+00
             1.4182E+00
 GRADIENT:   4.2397E+00  1.6491E+01 -3.0866E+00 -1.5474E+01 -1.3045E+00 -6.6228E+00  3.0140E+01  0.0000E+00 -1.2034E+00  0.0000E+00
            -2.0700E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1237.47268121622        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      669
 NPARAMETR:  9.6878E-01  4.1878E-02  3.3832E-01  1.3362E+00  2.8911E-01  9.3678E-01  5.5429E+00  1.0000E-02  1.1443E+00  1.0000E-02
             3.7388E+00
 PARAMETER:  6.8285E-02 -3.0730E+00 -9.8378E-01  3.8984E-01 -1.1409E+00  3.4690E-02  1.8125E+00 -9.5534E+00  2.3478E-01 -5.2236E+00
             1.4188E+00
 GRADIENT:   4.4471E+00  1.6045E+01 -1.2252E+00 -1.5759E+01 -3.4868E+00 -6.3667E+00  2.9680E+01  0.0000E+00 -9.5126E-01  0.0000E+00
            -1.9713E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1237.64015476594        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      847
 NPARAMETR:  9.6945E-01  4.1588E-02  3.6306E-01  1.3663E+00  3.0340E-01  9.3401E-01  5.6723E+00  1.0000E-02  1.1188E+00  1.0000E-02
             3.7505E+00
 PARAMETER:  6.8971E-02 -3.0800E+00 -9.1318E-01  4.1211E-01 -1.0927E+00  3.1729E-02  1.8356E+00 -9.5683E+00  2.1228E-01 -5.1420E+00
             1.4219E+00
 GRADIENT:  -6.5766E-02  4.2270E+00  1.2006E+00 -1.9415E+00 -3.5145E+00 -1.9280E+00  7.3734E+00  0.0000E+00 -8.4745E-01  0.0000E+00
            -5.4907E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1237.65359367868        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1025            RESET HESSIAN, TYPE II
 NPARAMETR:  9.6957E-01  4.0962E-02  3.6237E-01  1.3653E+00  3.0330E-01  9.3490E-01  5.7074E+00  1.0000E-02  1.1222E+00  1.0000E-02
             3.7553E+00
 PARAMETER:  6.9099E-02 -3.0951E+00 -9.1508E-01  4.1139E-01 -1.0930E+00  3.2686E-02  1.8418E+00 -9.6176E+00  2.1532E-01 -5.1707E+00
             1.4232E+00
 GRADIENT:   1.9846E+00  1.0238E+00 -1.6155E-01  3.9188E+00  4.1602E+00 -7.7605E-02  1.2938E+00  0.0000E+00  5.0251E-01  0.0000E+00
             1.1755E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1237.65497418176        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1181
 NPARAMETR:  9.7010E-01  4.0927E-02  3.6263E-01  1.3653E+00  3.0329E-01  9.3505E-01  5.7052E+00  1.0000E-02  1.1196E+00  1.0000E-02
             3.7552E+00
 PARAMETER:  6.9639E-02 -3.0960E+00 -9.1438E-01  4.1138E-01 -1.0931E+00  3.2846E-02  1.8414E+00 -9.6176E+00  2.1299E-01 -5.1707E+00
             1.4231E+00
 GRADIENT:   6.5819E-01  8.4828E-01 -7.1454E-02  1.9725E-01 -2.3388E-01 -1.5758E-01  1.0034E+00  0.0000E+00 -3.3214E-01  0.0000E+00
            -3.9986E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1237.66514911279        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1273
 NPARAMETR:  9.6869E-01  3.9288E-02  3.6218E-01  1.3623E+00  3.0289E-01  9.3500E-01  5.6396E+00  1.0000E-02  1.1221E+00  1.0000E-02
             3.7494E+00
 PARAMETER:  6.8185E-02 -3.1368E+00 -9.1560E-01  4.0918E-01 -1.0944E+00  3.2797E-02  1.8298E+00 -9.6176E+00  2.1523E-01 -5.1707E+00
             1.4216E+00
 GRADIENT:   4.6945E-01  1.9272E-01  1.2538E+00  1.6161E+00  3.3565E+00  3.5696E-01 -2.5374E-01  0.0000E+00  4.4085E-01  0.0000E+00
             1.0436E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1237.66797240379        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:     1389
 NPARAMETR:  9.6911E-01  3.8687E-02  3.6144E-01  1.3625E+00  3.0232E-01  9.3395E-01  5.6320E+00  1.0000E-02  1.1211E+00  1.0000E-02
             3.7460E+00
 PARAMETER:  6.8625E-02 -3.1522E+00 -9.1766E-01  4.0932E-01 -1.0963E+00  3.1667E-02  1.8285E+00 -9.6176E+00  2.1427E-01 -5.1707E+00
             1.4207E+00
 GRADIENT:  -8.6552E-01  1.4575E-01  7.6679E-01 -1.2656E+00 -8.1406E-01 -2.5353E-01 -3.0966E-01  0.0000E+00 -3.1851E-01  0.0000E+00
            -1.0440E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1237.71976428456        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1568
 NPARAMETR:  9.6739E-01  2.5166E-02  3.5790E-01  1.3614E+00  2.9884E-01  9.3236E-01  6.2095E+00  1.0000E-02  1.1194E+00  1.0000E-02
             3.7321E+00
 PARAMETER:  6.6844E-02 -3.5822E+00 -9.2750E-01  4.0850E-01 -1.1079E+00  2.9966E-02  1.9261E+00 -9.6176E+00  2.1283E-01 -5.1707E+00
             1.4170E+00
 GRADIENT:  -2.2893E+00  1.1491E-01  2.4355E+00 -3.0374E+00 -3.9664E+00 -8.3307E-01 -2.3195E-01  0.0000E+00 -1.2656E+00  0.0000E+00
            -3.9566E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1237.87551786562        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1745
 NPARAMETR:  9.6654E-01  1.0000E-02  3.6025E-01  1.3709E+00  2.9903E-01  9.3482E-01  1.1587E+01  1.0000E-02  1.1184E+00  1.0000E-02
             3.7484E+00
 PARAMETER:  6.5972E-02 -5.5405E+00 -9.2095E-01  4.1544E-01 -1.1072E+00  3.2594E-02  2.5498E+00 -9.6176E+00  2.1189E-01 -5.1707E+00
             1.4213E+00
 GRADIENT:  -2.5607E+00  0.0000E+00  1.3170E+00 -1.8296E+00 -1.7478E+00  2.6672E-01  1.0618E-01  0.0000E+00 -1.2309E-01  0.0000E+00
             8.9766E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1237.87990777979        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1909
 NPARAMETR:  9.6766E-01  1.0000E-02  3.6176E-01  1.3744E+00  3.0009E-01  9.3404E-01  1.1506E+01  1.0000E-02  1.1173E+00  1.0000E-02
             3.7491E+00
 PARAMETER:  6.7129E-02 -5.4744E+00 -9.1678E-01  4.1799E-01 -1.1037E+00  3.1762E-02  2.5429E+00 -9.6176E+00  2.1095E-01 -5.1707E+00
             1.4215E+00
 GRADIENT:  -2.1391E-03  0.0000E+00  5.4161E-04  7.3650E-03  2.0932E-03  2.7031E-03 -1.4050E-02  0.0000E+00  5.5494E-04  0.0000E+00
             8.9087E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1909
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.9863E-04 -7.9512E-04  1.3934E-04 -1.4596E-02  8.7914E-05
 SE:             2.8133E-02  1.7188E-03  2.4596E-04  2.6361E-02  3.8470E-04
 N:                     100         100         100         100         100

 P VAL.:         9.8302E-01  6.4365E-01  5.7105E-01  5.7978E-01  8.1924E-01

 ETASHRINKSD(%)  5.7512E+00  9.4242E+01  9.9176E+01  1.1688E+01  9.8711E+01
 ETASHRINKVR(%)  1.1172E+01  9.9668E+01  9.9993E+01  2.2009E+01  9.9983E+01
 EBVSHRINKSD(%)  5.2168E+00  9.4808E+01  9.9197E+01  1.0635E+01  9.8861E+01
 EBVSHRINKVR(%)  1.0161E+01  9.9730E+01  9.9994E+01  2.0139E+01  9.9987E+01
 RELATIVEINF(%)  7.3889E+01  2.5082E-02  2.7830E-04  1.9960E+01  3.8312E-04
 EPSSHRINKSD(%)  2.1739E+01
 EPSSHRINKVR(%)  3.8752E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1237.8799077797883     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -502.72908121605008     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.91
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.66
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1237.880       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  1.00E-02  3.62E-01  1.37E+00  3.00E-01  9.34E-01  1.15E+01  1.00E-02  1.12E+00  1.00E-02  3.75E+00
 


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
+        6.10E+00
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.46E+02  0.00E+00  3.51E+03
 
 TH 4
+       -5.16E+00  0.00E+00  1.24E+02  4.37E+00
 
 TH 5
+        2.67E+02  0.00E+00 -6.40E+03 -2.26E+02  1.17E+04
 
 TH 6
+       -5.94E-01  0.00E+00  1.43E+01  5.03E-01 -2.60E+01  5.79E-02
 
 TH 7
+       -1.37E-04  0.00E+00  3.28E-03  1.16E-04 -6.00E-03  1.34E-05  3.08E-09
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.51E+00  0.00E+00 -3.61E+01 -1.27E+00  6.59E+01 -1.47E-01 -3.38E-05  0.00E+00  3.72E-01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        7.16E-01  0.00E+00 -1.72E+01 -6.06E-01  3.13E+01 -6.98E-02 -1.61E-05  0.00E+00  1.77E-01  0.00E+00  8.41E-02
 
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
+        1.21E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.16E+02  0.00E+00  5.03E+03
 
 TH 4
+       -5.02E+01  0.00E+00 -1.15E+02  3.76E+02
 
 TH 5
+        3.26E+02  0.00E+00 -8.24E+03 -4.34E+02  1.55E+04
 
 TH 6
+       -2.72E+00  0.00E+00  4.71E+01 -1.65E+01 -1.87E+01  1.81E+02
 
 TH 7
+        3.86E-02  0.00E+00  3.87E-03 -1.38E-01 -1.36E-02 -3.66E-02  2.73E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.67E+01  0.00E+00  2.62E+01 -1.09E+01  1.27E+02  1.56E-01 -2.47E-02  0.00E+00  9.95E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.38E+01  0.00E+00 -2.41E+01 -5.23E+00  4.17E+01  2.66E+00 -4.72E-02  0.00E+00  5.14E+00  0.00E+00  3.16E+01
 
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
+        1.22E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.07E+02  0.00E+00  5.93E+03
 
 TH 4
+       -2.99E+02  0.00E+00  1.13E+02  3.78E+02
 
 TH 5
+        1.12E+03  0.00E+00 -1.07E+04 -5.19E+02  2.07E+04
 
 TH 6
+       -1.42E+02  0.00E+00  1.07E+02  9.69E+00 -8.81E+01  1.41E+02
 
 TH 7
+        1.06E-02  0.00E+00 -1.22E-01  1.13E-03  2.96E-01 -8.65E-03  1.96E-05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.11E+01  0.00E+00  2.44E+02 -7.55E+01 -3.67E+02 -3.75E+00  8.36E-03  0.00E+00  1.08E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        6.05E+01  0.00E+00 -1.01E+02 -5.90E+01  2.21E+02  2.30E+01  5.72E-03  0.00E+00  3.02E+01  0.00E+00  9.89E+01
 
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
 #CPUT: Total CPU Time in Seconds,       29.653
Stop Time:
Sat Sep 25 09:30:31 CDT 2021
