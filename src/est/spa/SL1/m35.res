Sat Sep 25 10:28:45 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1687.64029345745        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.5193E+01 -3.8417E+01 -5.4176E+01  9.2663E+00  9.3187E+01  1.6058E+01 -1.0883E+00  1.1017E+01 -3.3408E+00 -1.5586E+01
            -1.6473E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1696.24222636416        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      176
 NPARAMETR:  1.0590E+00  1.0165E+00  1.0792E+00  1.0015E+00  9.6902E-01  9.3500E-01  9.9221E-01  9.3609E-01  1.0251E+00  1.0382E+00
             1.0640E+00
 PARAMETER:  1.5734E-01  1.1632E-01  1.7625E-01  1.0147E-01  6.8530E-02  3.2795E-02  9.2181E-02  3.3961E-02  1.2476E-01  1.3745E-01
             1.6206E-01
 GRADIENT:   5.7561E+00  1.6508E+00 -2.4377E+00 -3.4071E+00 -7.4828E+00 -8.4217E+00  2.6646E+00  6.8329E+00  1.6346E-01 -3.8639E+00
             9.1798E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1697.84398214632        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  1.0560E+00  9.2380E-01  1.1673E+00  1.0660E+00  9.7439E-01  9.4931E-01  8.5908E-01  7.6295E-01  1.0150E+00  1.1224E+00
             1.0575E+00
 PARAMETER:  1.5449E-01  2.0744E-02  2.5474E-01  1.6389E-01  7.4052E-02  4.7982E-02 -5.1896E-02 -1.7056E-01  1.1488E-01  2.1547E-01
             1.5588E-01
 GRADIENT:   1.0795E+00  5.6365E+00  8.7484E-01  4.9037E+00 -7.2543E+00 -1.7383E+00 -1.5836E+00  1.4092E+00 -1.7805E+00 -1.7617E+00
             7.4075E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1698.15591965091        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  1.0549E+00  8.6900E-01  1.1592E+00  1.0983E+00  9.5795E-01  9.5398E-01  9.8881E-01  6.5914E-01  9.7941E-01  1.1283E+00
             1.0353E+00
 PARAMETER:  1.5349E-01 -4.0411E-02  2.4776E-01  1.9373E-01  5.7040E-02  5.2886E-02  8.8751E-02 -3.1682E-01  7.9190E-02  2.2070E-01
             1.3465E-01
 GRADIENT:  -1.9883E-01  2.0021E+00  1.0519E+00  1.5860E+00 -1.3374E+00  1.1080E-01  5.4248E-02 -1.4375E-01 -1.2260E-01 -5.1134E-02
            -6.9657E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1698.26225254811        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  1.0529E+00  6.7724E-01  1.2547E+00  1.2222E+00  9.2596E-01  9.5164E-01  1.1128E+00  7.3354E-01  9.0563E-01  1.1292E+00
             1.0369E+00
 PARAMETER:  1.5151E-01 -2.8974E-01  3.2686E-01  3.0068E-01  2.3077E-02  5.0430E-02  2.0688E-01 -2.0987E-01  8.7652E-04  2.2154E-01
             1.3626E-01
 GRADIENT:   4.6281E-02  2.5328E+00  7.7632E-01  4.9802E+00 -1.9688E+00 -6.4768E-03  4.3866E-02  6.3219E-02  1.1607E-01  1.0498E-01
            -1.8912E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1698.31980119663        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  1.0513E+00  5.6310E-01  1.3112E+00  1.2922E+00  9.1153E-01  9.5014E-01  1.1768E+00  7.7507E-01  8.6915E-01  1.1346E+00
             1.0371E+00
 PARAMETER:  1.4999E-01 -4.7430E-01  3.7091E-01  3.5631E-01  7.3731E-03  4.8855E-02  2.6283E-01 -1.5480E-01 -4.0236E-02  2.2632E-01
             1.3641E-01
 GRADIENT:   2.1150E-02  7.0699E-01  3.9497E-01  1.5023E+00 -5.7902E-01 -2.1083E-02 -1.4948E-01 -9.4230E-02 -1.3486E-01 -3.0411E-02
            -2.3847E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1698.36424074460        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0490E+00  4.1647E-01  1.3647E+00  1.3833E+00  8.8605E-01  9.4821E-01  1.4673E+00  8.3853E-01  8.1910E-01  1.1264E+00
             1.0377E+00
 PARAMETER:  1.4784E-01 -7.7594E-01  4.1093E-01  4.2448E-01 -2.0984E-02  4.6821E-02  4.8340E-01 -7.6107E-02 -9.9543E-02  2.1905E-01
             1.3703E-01
 GRADIENT:  -1.1759E-01  4.7120E-01  2.5565E-01  1.0540E+00 -9.1292E-01 -8.9210E-02  1.2082E-01  1.8938E-01  1.6393E-01  2.7675E-01
             2.4455E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1698.37590078191        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1236
 NPARAMETR:  1.0476E+00  3.2982E-01  1.4025E+00  1.4375E+00  8.7487E-01  9.4743E-01  1.6442E+00  8.7898E-01  7.9516E-01  1.1260E+00
             1.0372E+00
 PARAMETER:  1.4650E-01 -1.0092E+00  4.3825E-01  4.6294E-01 -3.3680E-02  4.5999E-02  5.9724E-01 -2.8995E-02 -1.2921E-01  2.1870E-01
             1.3654E-01
 GRADIENT:  -8.3750E-02  3.6447E-01 -2.6842E-03  1.9581E+00 -5.3386E-01  6.4225E-02  2.0996E-02  7.9285E-02  4.1800E-02  2.0487E-01
            -1.2469E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1698.37876065391        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1414
 NPARAMETR:  1.0468E+00  2.8032E-01  1.4374E+00  1.4688E+00  8.7299E-01  9.4654E-01  1.7073E+00  9.1386E-01  7.8349E-01  1.1286E+00
             1.0378E+00
 PARAMETER:  1.4570E-01 -1.1718E+00  4.6284E-01  4.8446E-01 -3.5826E-02  4.5058E-02  6.3492E-01  9.9237E-03 -1.4400E-01  2.2098E-01
             1.3707E-01
 GRADIENT:  -2.0476E-02  2.4324E-01  1.7772E-01  1.6664E+00 -3.1747E-01 -8.7211E-05 -5.6688E-02 -6.6061E-02 -1.7336E-01 -9.0142E-02
            -8.8063E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1698.38969517633        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1595
 NPARAMETR:  1.0451E+00  1.8419E-01  1.5115E+00  1.5315E+00  8.7054E-01  9.4545E-01  1.9212E+00  9.9135E-01  7.6074E-01  1.1321E+00
             1.0376E+00
 PARAMETER:  1.4411E-01 -1.5918E+00  5.1309E-01  5.2628E-01 -3.8646E-02  4.3910E-02  7.5295E-01  9.1315E-02 -1.7347E-01  2.2405E-01
             1.3693E-01
 GRADIENT:   4.0896E-02  4.5206E-01  5.1690E-01  4.7301E+00 -5.6072E-01  1.1529E-01 -5.6478E-02 -2.1718E-01 -3.6770E-01 -4.6490E-01
            -4.4149E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1698.42324588491        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1773
 NPARAMETR:  1.0435E+00  9.9977E-02  1.5503E+00  1.5844E+00  8.6009E-01  9.4415E-01  2.4117E+00  1.0387E+00  7.4055E-01  1.1315E+00
             1.0378E+00
 PARAMETER:  1.4260E-01 -2.2028E+00  5.3843E-01  5.6024E-01 -5.0719E-02  4.2534E-02  9.8032E-01  1.3799E-01 -2.0036E-01  2.2351E-01
             1.3713E-01
 GRADIENT:  -1.8900E-01  3.4817E-01  2.4608E-01  6.5224E+00 -1.7207E+00  4.1312E-02 -3.1375E-02 -6.9464E-02 -5.0235E-01  1.2850E-01
            -3.7320E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1698.46135186823        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1949
 NPARAMETR:  1.0427E+00  4.7303E-02  1.6210E+00  1.6189E+00  8.6729E-01  9.4277E-01  3.3586E+00  1.1088E+00  7.2669E-01  1.1364E+00
             1.0394E+00
 PARAMETER:  1.4178E-01 -2.9512E+00  5.8306E-01  5.8176E-01 -4.2386E-02  4.1071E-02  1.3115E+00  2.0329E-01 -2.1926E-01  2.2790E-01
             1.3866E-01
 GRADIENT:  -4.7171E-02  1.3461E-01  6.0693E-01  3.7379E+00 -1.2164E+00 -2.1383E-01 -4.4982E-03  1.5651E-02 -6.6900E-01  2.2289E-02
             3.0201E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1698.48961751067        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2126
 NPARAMETR:  1.0420E+00  1.5052E-02  1.6145E+00  1.6366E+00  8.5879E-01  9.4287E-01  6.1258E+00  1.1055E+00  7.1953E-01  1.1328E+00
             1.0387E+00
 PARAMETER:  1.4110E-01 -4.0962E+00  5.7901E-01  5.9262E-01 -5.2231E-02  4.1168E-02  1.9125E+00  2.0034E-01 -2.2916E-01  2.2471E-01
             1.3798E-01
 GRADIENT:  -2.8276E-01  1.9838E-02 -3.5897E-01  1.0684E+00  5.9129E-02 -8.2229E-03 -1.5134E-03 -4.1402E-03 -5.3443E-01  1.5976E-01
             7.7247E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1698.49517836976        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2303
 NPARAMETR:  1.0420E+00  1.0000E-02  1.6312E+00  1.6402E+00  8.6200E-01  9.4280E-01  7.6355E+00  1.1214E+00  7.1915E-01  1.1337E+00
             1.0385E+00
 PARAMETER:  1.4114E-01 -4.5258E+00  5.8934E-01  5.9481E-01 -4.8504E-02  4.1104E-02  2.1328E+00  2.1455E-01 -2.2969E-01  2.2548E-01
             1.3776E-01
 GRADIENT:   1.0257E-02  0.0000E+00 -1.1672E-02 -4.3191E-02  4.3403E-02 -7.8447E-03  3.3892E-04 -6.3607E-03 -9.6692E-04 -1.2232E-02
            -2.6179E-03

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1698.49517836976        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     2325
 NPARAMETR:  1.0420E+00  1.0000E-02  1.6312E+00  1.6402E+00  8.6200E-01  9.4280E-01  7.6355E+00  1.1214E+00  7.1915E-01  1.1337E+00
             1.0385E+00
 PARAMETER:  1.4114E-01 -4.5258E+00  5.8934E-01  5.9481E-01 -4.8504E-02  4.1104E-02  2.1328E+00  2.1455E-01 -2.2969E-01  2.2548E-01
             1.3776E-01
 GRADIENT:   1.0257E-02  0.0000E+00 -1.1672E-02 -4.3191E-02  4.3403E-02 -7.8447E-03  3.3892E-04 -6.3607E-03 -9.6692E-04 -1.2232E-02
            -2.6179E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2325
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.6054E-04 -1.4715E-03 -2.4541E-02 -6.9401E-03 -3.3334E-02
 SE:             2.9846E-02  1.2498E-03  1.4734E-02  2.9139E-02  2.2575E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9571E-01  2.3905E-01  9.5806E-02  8.1175E-01  1.3980E-01

 ETASHRINKSD(%)  1.3060E-02  9.5813E+01  5.0638E+01  2.3790E+00  2.4369E+01
 ETASHRINKVR(%)  2.6118E-02  9.9825E+01  7.5634E+01  4.7015E+00  4.2800E+01
 EBVSHRINKSD(%)  4.7761E-01  9.5896E+01  5.4122E+01  2.6667E+00  2.0406E+01
 EBVSHRINKVR(%)  9.5293E-01  9.9832E+01  7.8952E+01  5.2624E+00  3.6647E+01
 RELATIVEINF(%)  9.4318E+01  4.3205E-03  4.3756E+00  2.8557E+00  7.3910E+00
 EPSSHRINKSD(%)  4.3993E+01
 EPSSHRINKVR(%)  6.8632E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1698.4951783697618     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -963.34435180602361     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1698.495       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.00E-02  1.63E+00  1.64E+00  8.62E-01  9.43E-01  7.64E+00  1.12E+00  7.19E-01  1.13E+00  1.04E+00
 


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
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.70E-01  0.00E+00  7.55E+01
 
 TH 4
+       -1.24E+01  0.00E+00 -3.07E+01  7.57E+02
 
 TH 5
+       -1.54E+00  0.00E+00 -1.86E+02 -5.59E+01  7.50E+02
 
 TH 6
+       -8.32E+00  0.00E+00 -3.46E+00  1.11E+00 -3.74E+01  2.32E+02
 
 TH 7
+        1.13E-02  0.00E+00  1.60E-02 -3.10E-02 -1.71E-01  2.88E-01  1.32E-03
 
 TH 8
+        1.30E+01  0.00E+00 -1.91E+01 -1.82E+00  6.18E+00 -2.29E+01 -4.93E-02  2.07E+01
 
 TH 9
+        4.44E-01  0.00E+00  6.54E+00 -1.60E+00 -2.41E+01 -2.37E+01 -7.00E-02 -4.58E+00  3.50E+02
 
 TH10
+        1.39E+00  0.00E+00 -3.93E+00 -7.89E-01 -9.40E+01  1.77E+01  1.41E-01  2.04E+01 -3.82E-01  6.70E+01
 
 TH11
+        9.36E+00  0.00E+00 -7.49E+00 -1.01E+01  2.28E+01  1.86E+01  4.66E-02  1.43E+01  1.54E+01  1.47E+01  2.10E+02
 
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
 #CPUT: Total CPU Time in Seconds,       34.483
Stop Time:
Sat Sep 25 10:29:21 CDT 2021
