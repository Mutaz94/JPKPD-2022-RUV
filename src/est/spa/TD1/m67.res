Sat Sep 25 13:02:23 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1652.56499016208        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1779E+01 -1.0239E+02  3.0254E+00 -1.6678E+02 -1.2403E+01 -3.7995E+00 -1.7532E+01  4.8145E+00 -1.9091E+01  5.1853E+00
             5.7528E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1663.93721909812        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.8287E-01  1.1031E+00  1.0480E+00  1.0561E+00  1.0785E+00  1.0193E+00  1.1387E+00  9.6219E-01  1.0001E+00  9.8749E-01
             9.7517E-01
 PARAMETER:  8.2724E-02  1.9811E-01  1.4692E-01  1.5460E-01  1.7560E-01  1.1907E-01  2.2987E-01  6.1458E-02  1.0008E-01  8.7411E-02
             7.4852E-02
 GRADIENT:   1.5563E+01  7.3634E+00  1.7694E+00  1.3233E+01  1.2452E+01  4.4677E+00 -4.1095E+00 -1.1749E+00  1.6289E+00 -9.1997E+00
            -8.3818E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1664.63775449533        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  9.8788E-01  1.0519E+00  1.0549E+00  1.0882E+00  1.0503E+00  1.0261E+00  1.2623E+00  8.7445E-01  9.6418E-01  1.0319E+00
             9.8753E-01
 PARAMETER:  8.7810E-02  1.5059E-01  1.5347E-01  1.8454E-01  1.4903E-01  1.2578E-01  3.3290E-01 -3.4157E-02  6.3519E-02  1.3140E-01
             8.7447E-02
 GRADIENT:   2.6210E+01  1.1099E+01  7.5276E+00  1.2007E+01 -5.3500E+00  7.5951E+00  2.8404E+00 -1.1964E+00  2.6128E+00 -1.0491E+00
            -2.4717E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.86227380324        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      321
 NPARAMETR:  9.9660E-01  9.8259E-01  1.0635E+00  1.1312E+00  1.0250E+00  1.0228E+00  1.3079E+00  9.0169E-01  9.2800E-01  1.0147E+00
             9.9033E-01
 PARAMETER:  9.6598E-02  8.2435E-02  1.6154E-01  2.2328E-01  1.2473E-01  1.2257E-01  3.6841E-01 -3.4804E-03  2.5272E-02  1.1457E-01
             9.0280E-02
 GRADIENT:   4.2220E+00  2.2248E+00  4.4412E-01  5.0606E-01 -1.3297E+00 -6.1364E-02 -7.2894E-02  4.3327E-01 -8.0658E-01 -4.8467E-01
            -2.7289E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1665.04189483607        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      498
 NPARAMETR:  9.9140E-01  7.4054E-01  1.2096E+00  1.2892E+00  9.9761E-01  1.0199E+00  1.5424E+00  9.6187E-01  8.7392E-01  1.0349E+00
             9.9402E-01
 PARAMETER:  9.1364E-02 -2.0037E-01  2.9032E-01  3.5402E-01  9.7607E-02  1.1971E-01  5.3332E-01  6.1123E-02 -3.4763E-02  1.3430E-01
             9.4006E-02
 GRADIENT:  -6.0623E-01  2.1692E+00  1.2608E-01  4.9280E+00 -3.5279E-01  2.8715E-02  6.9690E-02 -1.9676E-01 -2.5255E-01  3.7889E-02
             8.1632E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1665.14358669092        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      675
 NPARAMETR:  9.8719E-01  4.9711E-01  1.4129E+00  1.4553E+00  1.0005E+00  1.0161E+00  1.8179E+00  1.1608E+00  8.4130E-01  1.0748E+00
             9.9808E-01
 PARAMETER:  8.7110E-02 -5.9895E-01  4.4563E-01  4.7521E-01  1.0053E-01  1.1594E-01  6.9769E-01  2.4910E-01 -7.2810E-02  1.7213E-01
             9.8081E-02
 GRADIENT:  -1.4256E+00  4.8523E+00 -8.6532E-01  1.6508E+01 -2.6420E+00  1.3732E-02  6.3107E-01  1.1883E+00  3.2952E-01  1.2055E+00
             2.7160E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1665.54195888083        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      852
 NPARAMETR:  9.8571E-01  3.2404E-01  1.5196E+00  1.5566E+00  9.8739E-01  1.0132E+00  2.0717E+00  1.2457E+00  8.1507E-01  1.0718E+00
             9.8749E-01
 PARAMETER:  8.5608E-02 -1.0269E+00  5.1846E-01  5.4252E-01  8.7311E-02  1.1316E-01  8.2838E-01  3.1974E-01 -1.0449E-01  1.6932E-01
             8.7413E-02
 GRADIENT:   1.6748E+00  1.3053E+00  2.8346E+00  3.5609E+00 -2.9171E+00 -4.5431E-02 -4.0555E-01 -7.4207E-01 -1.3508E-01 -8.0463E-01
            -2.0956E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1665.90088271218        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1027
 NPARAMETR:  9.8172E-01  1.4744E-01  1.6061E+00  1.6642E+00  9.7335E-01  1.0103E+00  2.7370E+00  1.3385E+00  7.8400E-01  1.0792E+00
             9.8755E-01
 PARAMETER:  8.1556E-02 -1.8143E+00  5.7383E-01  6.0933E-01  7.2985E-02  1.1026E-01  1.1069E+00  3.9151E-01 -1.4335E-01  1.7625E-01
             8.7469E-02
 GRADIENT:  -3.8677E-01  4.9347E-02 -2.7787E-01 -9.7418E-01  1.0086E+00 -8.4477E-02 -3.8706E-01 -6.7781E-01 -7.1689E-01 -4.2156E-01
            -1.6436E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1666.08217801215        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  9.8014E-01  5.5560E-02  1.7117E+00  1.7237E+00  9.8510E-01  1.0092E+00  4.2249E+00  1.4465E+00  7.6558E-01  1.0853E+00
             9.9230E-01
 PARAMETER:  7.9942E-02 -2.7903E+00  6.3749E-01  6.4447E-01  8.4984E-02  1.0915E-01  1.5410E+00  4.6916E-01 -1.6712E-01  1.8189E-01
             9.2269E-02
 GRADIENT:  -4.8547E-01 -2.6830E-02 -1.7564E-01 -2.9939E+00  2.2053E+00  3.9957E-02 -1.4034E-01 -3.2800E-01 -5.8684E-01 -8.1161E-01
             2.4508E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1666.15528853773        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1377
 NPARAMETR:  9.7941E-01  1.6519E-02  1.7075E+00  1.7499E+00  9.7097E-01  1.0084E+00  7.7164E+00  1.4498E+00  7.5868E-01  1.0859E+00
             9.9145E-01
 PARAMETER:  7.9197E-02 -4.0033E+00  6.3505E-01  6.5958E-01  7.0544E-02  1.0841E-01  2.1434E+00  4.7139E-01 -1.7618E-01  1.8241E-01
             9.1418E-02
 GRADIENT:  -5.4514E-01  3.9916E-02 -1.4382E-02  4.1448E+00 -1.4968E+00 -8.4228E-03 -4.1790E-02  2.3272E-01 -4.6210E-01  4.1253E-01
             2.2434E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1666.17463523363        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1552
 NPARAMETR:  9.7973E-01  1.0000E-02  1.7269E+00  1.7531E+00  9.7607E-01  1.0084E+00  1.2275E+01  1.4642E+00  7.5768E-01  1.0868E+00
             9.9111E-01
 PARAMETER:  7.9517E-02 -4.8619E+00  6.4630E-01  6.6136E-01  7.5781E-02  1.0840E-01  2.6075E+00  4.8128E-01 -1.7750E-01  1.8327E-01
             9.1067E-02
 GRADIENT:   4.3678E-01  0.0000E+00  3.4311E-01  2.8082E-01 -4.6436E-01  1.5145E-02  9.9321E-03 -4.1560E-02 -2.2807E-03  2.1882E-02
            -5.7939E-02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1666.17498030372        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1644
 NPARAMETR:  9.7951E-01  1.0000E-02  1.7238E+00  1.7528E+00  9.7564E-01  1.0083E+00  1.2079E+01  1.4622E+00  7.5776E-01  1.0865E+00
             9.9120E-01
 PARAMETER:  7.9301E-02 -4.8315E+00  6.4454E-01  6.6119E-01  7.5335E-02  1.0831E-01  2.5915E+00  4.7993E-01 -1.7739E-01  1.8295E-01
             9.1157E-02
 GRADIENT:  -3.1361E-02  0.0000E+00 -1.1713E-02 -7.9265E-02  2.5774E-02 -8.9426E-03 -8.5980E-06  2.7673E-03 -8.9203E-04 -1.0359E-02
             6.7068E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1644
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4877E-04  4.5402E-04 -3.5228E-02 -7.2536E-03 -4.8557E-02
 SE:             2.9845E-02  1.8212E-03  1.8954E-02  2.9234E-02  1.9789E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9335E-01  8.0313E-01  6.3090E-02  8.0404E-01  1.4135E-02

 ETASHRINKSD(%)  1.6197E-02  9.3899E+01  3.6501E+01  2.0628E+00  3.3706E+01
 ETASHRINKVR(%)  3.2391E-02  9.9628E+01  5.9679E+01  4.0830E+00  5.6051E+01
 EBVSHRINKSD(%)  4.0218E-01  9.4112E+01  3.9805E+01  2.3631E+00  3.0340E+01
 EBVSHRINKVR(%)  8.0274E-01  9.9653E+01  6.3766E+01  4.6704E+00  5.1474E+01
 RELATIVEINF(%)  9.3562E+01  9.2580E-03  9.8326E+00  3.0048E+00  7.0105E+00
 EPSSHRINKSD(%)  4.5275E+01
 EPSSHRINKVR(%)  7.0052E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1666.1749803037237     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -931.02415373998554     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1666.175       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.00E-02  1.72E+00  1.75E+00  9.76E-01  1.01E+00  1.21E+01  1.46E+00  7.58E-01  1.09E+00  9.91E-01
 


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
+       -5.82E-01  0.00E+00  6.70E+01
 
 TH 4
+       -5.51E+00  0.00E+00 -1.51E+01  6.01E+02
 
 TH 5
+        2.03E+01  0.00E+00 -1.37E+02 -4.42E+01  5.13E+02
 
 TH 6
+       -5.56E+00  0.00E+00 -2.12E-01 -4.41E+00  1.75E+01  2.08E+02
 
 TH 7
+       -1.72E-02  0.00E+00 -8.48E-04 -4.70E-03  3.48E-02  3.07E-02  1.91E-03
 
 TH 8
+       -1.43E+00  0.00E+00 -1.86E+01 -2.96E+00 -6.32E+00  1.62E+00 -1.89E-03  2.25E+01
 
 TH 9
+        1.02E+01  0.00E+00  5.82E+00 -2.17E+00  2.81E+00 -7.48E+00  9.21E-02 -2.14E-01  3.26E+02
 
 TH10
+        1.03E+00  0.00E+00  7.75E-01  3.17E-01 -7.11E+01 -2.75E-01 -2.95E-02  1.36E+01  6.02E+00  5.89E+01
 
 TH11
+       -2.68E+01  0.00E+00 -6.37E+00 -8.50E+00 -2.02E+01  3.16E+00  5.75E-02  1.26E+01  2.96E+01  2.77E+01  2.39E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.621
Stop Time:
Sat Sep 25 13:02:52 CDT 2021
