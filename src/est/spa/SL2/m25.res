Sat Sep 18 12:09:36 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1679.62430576556        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.4861E+01 -5.8724E+01  3.5293E+00 -1.0646E+02 -3.0705E+01  3.5700E+01  2.7737E+00  8.5221E+00 -2.1616E+00  7.9299E+00
             8.0894E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1685.91476427070        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9236E-01  1.0474E+00  1.0975E+00  1.0484E+00  1.1086E+00  8.7266E-01  9.5564E-01  9.1619E-01  9.7931E-01  1.0036E+00
             1.0165E+00
 PARAMETER:  9.2328E-02  1.4633E-01  1.9302E-01  1.4729E-01  2.0311E-01 -3.6210E-02  5.4621E-02  1.2465E-02  7.9095E-02  1.0359E-01
             1.1638E-01
 GRADIENT:   4.8670E+01 -1.1618E-01 -5.9331E+00  1.0208E+01  2.7422E+01 -1.5090E+01 -1.1601E-01 -1.8029E+00 -4.7352E+00 -1.1763E+01
             1.4781E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1686.67115047978        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.8453E-01  1.1367E+00  1.0836E+00  9.9344E-01  1.1418E+00  8.8746E-01  8.2234E-01  9.6427E-01  1.0703E+00  1.0907E+00
             1.0054E+00
 PARAMETER:  8.4405E-02  2.2817E-01  1.8026E-01  9.3418E-02  2.3263E-01 -1.9395E-02 -9.5607E-02  6.3614E-02  1.6791E-01  1.8684E-01
             1.0535E-01
 GRADIENT:   2.6215E+01  8.5536E+00 -4.6119E+00  1.0508E+01  1.1644E+01 -8.1932E+00 -1.6053E+00 -1.0982E+00  4.2728E-01 -8.9410E-01
            -1.1760E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1687.05664148458        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      317
 NPARAMETR:  9.8928E-01  1.0309E+00  1.2135E+00  1.0636E+00  1.1396E+00  9.1285E-01  8.6858E-01  1.0949E+00  1.0203E+00  1.1079E+00
             1.0041E+00
 PARAMETER:  8.9220E-02  1.3040E-01  2.9349E-01  1.6168E-01  2.3068E-01  8.8163E-03 -4.0896E-02  1.9068E-01  1.2005E-01  2.0248E-01
             1.0405E-01
 GRADIENT:   3.1724E+00  2.0302E+00  7.4961E-01  1.1322E+00 -1.3542E+00  3.9923E-01  5.7003E-02  3.6851E-02  1.3798E-02  1.1155E-01
            -8.1514E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1687.18453307441        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      493
 NPARAMETR:  9.8521E-01  8.1555E-01  1.3796E+00  1.2062E+00  1.1169E+00  9.0898E-01  9.0755E-01  1.1685E+00  9.3819E-01  1.1181E+00
             1.0057E+00
 PARAMETER:  8.5100E-02 -1.0389E-01  4.2176E-01  2.8750E-01  2.1060E-01  4.5709E-03  2.9956E-03  2.5575E-01  3.6198E-02  2.1161E-01
             1.0568E-01
 GRADIENT:  -1.7526E+00  3.2885E+00  1.0548E+00  4.9196E+00 -1.8447E+00 -3.4706E-01 -5.3430E-02 -2.4759E-01 -2.6558E-01  1.3779E-02
             5.3207E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1687.32435656584        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      669
 NPARAMETR:  9.8369E-01  6.0356E-01  1.5605E+00  1.3443E+00  1.1077E+00  9.0829E-01  9.1267E-01  1.3033E+00  8.7228E-01  1.1363E+00
             1.0038E+00
 PARAMETER:  8.3559E-02 -4.0491E-01  5.4501E-01  3.9584E-01  2.0232E-01  3.8030E-03  8.6181E-03  3.6492E-01 -3.6640E-02  2.2782E-01
             1.0376E-01
 GRADIENT:   1.3881E+00  2.5610E+00  3.9151E-01  6.4826E+00 -2.0314E+00  4.7596E-01 -1.7283E-02  4.3041E-01  2.9898E-01  2.4979E-01
            -3.5199E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1687.51337155592        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  9.8143E-01  3.8739E-01  1.7374E+00  1.4799E+00  1.0974E+00  9.0543E-01  9.1313E-01  1.4534E+00  8.0298E-01  1.1431E+00
             1.0032E+00
 PARAMETER:  8.1252E-02 -8.4831E-01  6.5238E-01  4.9200E-01  1.9297E-01  6.5278E-04  9.1183E-03  4.7388E-01 -1.1942E-01  2.3378E-01
             1.0319E-01
 GRADIENT:   3.4081E+00  4.0814E-01 -8.2628E-02 -1.2656E+00 -4.7628E-01  5.6321E-01 -7.1153E-02  3.3995E-01 -9.6396E-01 -7.0025E-02
            -4.3631E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1687.79560984476        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1027
 NPARAMETR:  9.7431E-01  1.8376E-01  2.0226E+00  1.6252E+00  1.1164E+00  8.9689E-01  8.7066E-01  1.7272E+00  7.3800E-01  1.1686E+00
             1.0062E+00
 PARAMETER:  7.3974E-02 -1.5941E+00  8.0440E-01  5.8561E-01  2.1013E-01 -8.8241E-03 -3.8501E-02  6.4648E-01 -2.0382E-01  2.5578E-01
             1.0617E-01
 GRADIENT:  -8.6855E+00  2.3669E+00 -2.2814E-01  2.1439E+01 -3.8291E+00 -1.8405E+00 -5.7526E-02  9.3786E-01 -2.1553E+00  5.4736E-01
             9.0686E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1688.18910952907        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  9.7652E-01  7.7984E-02  2.2841E+00  1.6908E+00  1.1539E+00  9.0209E-01  8.0281E-01  1.9527E+00  7.0797E-01  1.1906E+00
             1.0017E+00
 PARAMETER:  7.6241E-02 -2.4512E+00  9.2599E-01  6.2520E-01  2.4316E-01 -3.0434E-03 -1.1963E-01  7.6920E-01 -2.4536E-01  2.7447E-01
             1.0169E-01
 GRADIENT:   1.4889E+00  3.3235E-01  3.3189E-02  1.0108E+00 -3.7969E-02  1.2538E+00 -4.8082E-03  5.7436E-01 -4.4100E-02 -5.3845E-01
            -8.3109E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1688.41244217228        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1377
 NPARAMETR:  9.7510E-01  1.6187E-02  2.3369E+00  1.7305E+00  1.1480E+00  8.9881E-01  6.9934E-01  1.9842E+00  6.8916E-01  1.1919E+00
             1.0028E+00
 PARAMETER:  7.4781E-02 -4.0235E+00  9.4883E-01  6.4841E-01  2.3803E-01 -6.6816E-03 -2.5762E-01  7.8520E-01 -2.7228E-01  2.7559E-01
             1.0283E-01
 GRADIENT:   7.6533E-03  6.4099E-02  5.5646E-01 -7.5289E-02 -7.9975E-01  2.3466E-01 -1.7410E-04 -3.1940E-01 -5.1237E-01  6.0623E-02
            -2.8632E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1688.43600753288        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1552
 NPARAMETR:  9.7491E-01  1.0000E-02  2.3498E+00  1.7347E+00  1.1507E+00  8.9817E-01  6.6222E-01  2.0008E+00  6.8805E-01  1.1940E+00
             1.0033E+00
 PARAMETER:  7.4588E-02 -4.5467E+00  9.5431E-01  6.5081E-01  2.4033E-01 -7.4010E-03 -3.1215E-01  7.9352E-01 -2.7390E-01  2.7729E-01
             1.0329E-01
 GRADIENT:  -3.1798E-01  0.0000E+00 -3.9741E-02 -5.9930E-01 -6.6489E-02  2.5823E-02 -5.0474E-05  1.3461E-03 -5.8404E-02  6.5599E-02
            -1.6224E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1688.43619892379        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1722
 NPARAMETR:  9.7504E-01  1.0000E-02  2.3624E+00  1.7352E+00  1.1533E+00  8.9811E-01  6.6279E-01  2.0110E+00  6.8794E-01  1.1955E+00
             1.0034E+00
 PARAMETER:  7.4724E-02 -4.5487E+00  9.5962E-01  6.5110E-01  2.4276E-01 -7.3643E-03 -3.1140E-01  7.9867E-01 -2.7407E-01  2.7849E-01
             1.0334E-01
 GRADIENT:  -3.3223E-03  0.0000E+00 -3.4236E-03 -3.0204E-03  2.7197E-02  2.2473E-02 -1.2270E-03  1.9260E-03 -1.6136E-03 -5.0675E-03
            -1.6094E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1722
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.7657E-05 -2.2814E-04 -4.0206E-02 -9.0915E-03 -5.8575E-02
 SE:             2.9823E-02  1.2038E-04  1.8790E-02  2.9203E-02  1.9349E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9819E-01  5.8058E-02  3.2371E-02  7.5555E-01  2.4681E-03

 ETASHRINKSD(%)  8.9145E-02  9.9597E+01  3.7052E+01  2.1678E+00  3.5178E+01
 ETASHRINKVR(%)  1.7821E-01  9.9998E+01  6.0375E+01  4.2886E+00  5.7980E+01
 EBVSHRINKSD(%)  5.2429E-01  9.9634E+01  4.0435E+01  2.7006E+00  3.0827E+01
 EBVSHRINKVR(%)  1.0458E+00  9.9999E+01  6.4520E+01  5.3283E+00  5.2151E+01
 RELATIVEINF(%)  9.6632E+01  7.4521E-05  1.4396E+01  5.8050E+00  1.3897E+01
 EPSSHRINKSD(%)  4.4890E+01
 EPSSHRINKVR(%)  6.9629E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1688.4361989237921     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -953.28537236005388     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.86
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1688.436       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.00E-02  2.36E+00  1.74E+00  1.15E+00  8.98E-01  6.63E-01  2.01E+00  6.88E-01  1.20E+00  1.00E+00
 


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
+        1.46E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.99E-01  0.00E+00  1.98E+01
 
 TH 4
+       -1.19E+01  0.00E+00 -1.03E+01  7.40E+02
 
 TH 5
+       -1.22E+01  0.00E+00 -4.54E+01 -3.86E+01  3.08E+02
 
 TH 6
+        8.18E-01  0.00E+00 -8.71E-01 -4.47E+00  3.17E+00  2.45E+02
 
 TH 7
+       -1.85E+01  0.00E+00 -3.04E-02 -3.58E-01 -2.65E+00 -1.78E+00 -4.00E+00
 
 TH 8
+        2.17E-01  0.00E+00 -1.03E+01 -3.10E+00 -1.36E+01  7.44E-01  2.37E-01  1.60E+01
 
 TH 9
+        2.12E+00  0.00E+00  4.83E+00 -8.09E-01  8.99E-01 -2.69E+01  1.16E+00  1.66E+00  3.87E+02
 
 TH10
+        1.90E+00  0.00E+00  9.13E-01 -5.94E-01 -5.32E+01  1.89E+00  2.85E+00  3.50E+00  1.39E+00  5.53E+01
 
 TH11
+       -2.14E+01  0.00E+00 -1.42E+00 -1.11E+01 -1.03E+01  1.63E+01 -7.95E+00  3.67E+00  1.20E+01 -6.46E+00  2.18E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.004
Stop Time:
Sat Sep 18 12:10:04 CDT 2021
