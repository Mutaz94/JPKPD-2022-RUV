Sat Sep 18 07:26:17 CDT 2021
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
$DATA ../../../../data/int/D/dat72.csv ignore=@
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m72.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26170.8462928341        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.1698E+02  3.9529E+02 -1.0685E+01  4.4898E+00  2.2578E+02 -3.1189E+03 -1.1473E+03 -7.8427E+01 -2.0433E+03 -9.2383E+02
            -5.2872E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1004.89173783694        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  2.5679E+00  1.3657E+00  7.9185E-01  2.4307E+00  9.5052E-01  6.5944E+00  5.0411E+00  1.0365E+00  4.9296E+00  2.5776E+00
             1.0034E+01
 PARAMETER:  1.0431E+00  4.1167E-01 -1.3339E-01  9.8820E-01  4.9251E-02  1.9862E+00  1.7176E+00  1.3584E-01  1.6953E+00  1.0468E+00
             2.4060E+00
 GRADIENT:   3.4920E+01 -1.7222E+01 -4.7007E+01  6.2380E+01 -4.1060E+01  1.4495E+02  3.3593E+01  4.7845E+00  9.9650E+01  6.4249E+01
             3.4176E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1101.90024687771        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  2.2438E+00  3.2665E+00  1.3821E+01  1.6575E+00  3.2379E+00  6.2751E+00  7.2715E+00  7.7183E-01  2.3079E+01  8.8622E-01
             9.5709E+00
 PARAMETER:  9.0819E-01  1.2837E+00  2.7262E+00  6.0534E-01  1.2749E+00  1.9366E+00  2.0840E+00 -1.5899E-01  3.2389E+00 -2.0795E-02
             2.3587E+00
 GRADIENT:   3.0841E+01 -1.6872E+01 -3.5512E+00  2.0963E+01  1.8077E+01  1.4211E+02  1.2240E+02  1.7616E-01  6.0001E+01  1.1049E+01
             3.6082E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1221.99756592262        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.6320E+00  4.3264E+00  8.6627E+00  1.3023E+00  2.6604E+00  3.3380E+00  3.5397E+00  6.0007E-01  1.9212E+01  1.0652E+00
             9.7717E+00
 PARAMETER:  5.8980E-01  1.5647E+00  2.2590E+00  3.6412E-01  1.0785E+00  1.3054E+00  1.3640E+00 -4.1070E-01  3.0555E+00  1.6313E-01
             2.3795E+00
 GRADIENT:   3.1947E+01 -2.0377E+00  1.0646E+00  2.0149E+01 -3.0784E+01  3.5944E+01 -1.0039E+01  9.8592E-01  3.6413E+01  2.0080E+01
             3.7269E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1282.24612039621        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.5251E+00  2.0789E+00  6.9472E+00  3.9432E-01  2.5630E+00  3.1257E+00  2.8518E+00  6.8305E+00  6.7242E+00  1.4061E+00
             7.7282E+00
 PARAMETER:  5.2203E-01  8.3183E-01  2.0383E+00 -8.3060E-01  1.0412E+00  1.2397E+00  1.1480E+00  2.0214E+00  2.0057E+00  4.4081E-01
             2.1449E+00
 GRADIENT:   2.9229E+01 -4.4663E+01 -1.8356E+01 -6.2987E+00  8.9848E+00 -2.7556E+01  2.8382E+01 -1.1552E+01  1.9645E+01  2.9468E+01
             4.1425E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1323.86729051049        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.3166E+00  2.3209E+00  2.0392E+01  3.9165E-01  2.4728E+00  2.8427E+00  2.4338E+00  7.3328E+00  6.2566E+00  4.7340E-01
             7.4075E+00
 PARAMETER:  3.7506E-01  9.4194E-01  3.1151E+00 -8.3740E-01  1.0054E+00  1.1448E+00  9.8945E-01  2.0924E+00  1.9336E+00 -6.4781E-01
             2.1025E+00
 GRADIENT:   8.7779E+00  2.8480E+00 -1.0707E+01 -1.0295E+01 -2.4666E+01 -2.7700E+01  1.6654E+01  4.9672E-01 -1.5317E+01  9.4009E-01
            -8.3456E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1353.63992990790        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  1.2536E+00  1.9062E+00  1.7787E+03  7.4362E-01  2.6991E+00  3.1128E+00  1.9201E+00  1.1941E+00  5.5996E+00  2.3041E-02
             7.6816E+00
 PARAMETER:  3.2600E-01  7.4513E-01  7.5836E+00 -1.9623E-01  1.0929E+00  1.2355E+00  7.5238E-01  2.7740E-01  1.8227E+00 -3.6705E+00
             2.1388E+00
 GRADIENT:  -5.7620E+00  7.8681E+00  6.9189E-01 -2.0576E+00  4.5669E-01  1.0635E+01  3.3848E+00 -3.8316E-03 -1.8277E+00  2.7191E-03
            -1.5126E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1355.33770313417        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  1.3175E+00  2.0313E+00  6.6694E+02  6.0592E-01  2.6709E+00  3.0898E+00  1.8072E+00  1.5121E+00  6.5687E+00  3.9583E-02
             7.7398E+00
 PARAMETER:  3.7576E-01  8.0867E-01  6.6027E+00 -4.0102E-01  1.0824E+00  1.2281E+00  6.9180E-01  5.1347E-01  1.9823E+00 -3.1294E+00
             2.1464E+00
 GRADIENT:   1.5877E+00 -5.9787E-01  8.8955E-02 -4.1667E-01 -2.1368E+00 -2.6997E+00  1.0580E+00  3.0544E-01  1.5530E+00  7.3217E-03
            -6.3834E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1355.46421015630        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  1.3060E+00  1.9425E+00  6.3693E+02  6.6926E-01  2.6799E+00  3.1147E+00  1.7184E+00  1.3391E+00  6.3097E+00  4.1693E-02
             7.7886E+00
 PARAMETER:  3.6698E-01  7.6399E-01  6.5567E+00 -3.0159E-01  1.0858E+00  1.2361E+00  6.4141E-01  3.9203E-01  1.9421E+00 -3.0774E+00
             2.1527E+00
 GRADIENT:  -7.2577E-01  1.3322E-01  2.1755E-02 -3.9272E-01 -3.8804E-01  5.3688E-01 -3.1015E-01  2.5737E-01 -3.6969E-01  8.8400E-03
             5.5129E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1355.57020032871        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      960
 NPARAMETR:  1.3127E+00  1.9329E+00  5.7408E+02  6.8898E-01  2.6759E+00  3.1201E+00  1.7292E+00  7.3719E-01  6.2549E+00  4.7548E-02
             7.7728E+00
 PARAMETER:  3.7212E-01  7.5905E-01  6.4528E+00 -2.7254E-01  1.0843E+00  1.2379E+00  6.4767E-01 -2.0491E-01  1.9334E+00 -2.9460E+00
             2.1506E+00
 GRADIENT:   5.2040E-01  9.9914E-01 -7.0682E-03  2.2442E-01 -8.0906E-01  1.0883E+00 -1.0405E-01  9.5519E-02 -5.3259E-01  1.1217E-02
             1.1537E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1355.62310131224        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1135
 NPARAMETR:  1.3100E+00  1.9263E+00  5.5160E+02  6.8537E-01  2.6772E+00  3.1112E+00  1.7218E+00  9.3952E-02  6.2717E+00  5.3327E-02
             7.7660E+00
 PARAMETER:  3.7003E-01  7.5560E-01  6.4128E+00 -2.7780E-01  1.0848E+00  1.2350E+00  6.4338E-01 -2.2650E+00  1.9360E+00 -2.8313E+00
             2.1498E+00
 GRADIENT:   1.1772E-01  1.1549E-03 -1.8450E-03 -2.5796E-02 -8.0213E-02 -5.6481E-03 -1.8783E-02  1.7016E-03 -5.9889E-03  1.4230E-02
            -1.9784E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1355.62347333581        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1311
 NPARAMETR:  1.3094E+00  1.9253E+00  5.5492E+02  6.8666E-01  2.6778E+00  3.1111E+00  1.7213E+00  5.0356E-02  6.2680E+00  5.4286E-02
             7.7671E+00
 PARAMETER:  3.6956E-01  7.5507E-01  6.4188E+00 -2.7591E-01  1.0850E+00  1.2350E+00  6.4311E-01 -2.8886E+00  1.9355E+00 -2.8135E+00
             2.1499E+00
 GRADIENT:  -4.7523E-03  2.1871E-02  1.1473E-02 -2.6909E-03 -5.5898E-04 -1.0110E-02 -1.1356E-02  4.8290E-04 -5.8040E-03  1.4813E-02
             7.6025E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1355.62738209697        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1492
 NPARAMETR:  1.3124E+00  1.9255E+00  5.4711E+02  6.8598E-01  2.6760E+00  3.1096E+00  1.7265E+00  4.4596E-02  6.2674E+00  3.1083E-02
             7.7664E+00
 PARAMETER:  3.7188E-01  7.5516E-01  6.4047E+00 -2.7691E-01  1.0843E+00  1.2345E+00  6.4608E-01 -3.0101E+00  1.9354E+00 -3.3711E+00
             2.1498E+00
 GRADIENT:   5.2290E-01 -2.0070E-01 -1.7233E-02 -5.1645E-02 -2.9934E-01 -1.9054E-01  1.6433E-01  3.8991E-04  6.6384E-03  4.7851E-03
            -1.0530E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1355.63085728832        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1661
 NPARAMETR:  1.3094E+00  1.9248E+00  5.5156E+02  6.8741E-01  2.6781E+00  3.1132E+00  1.7212E+00  1.0035E-02  6.2688E+00  1.0000E-02
             7.7667E+00
 PARAMETER:  3.6958E-01  7.5459E-01  6.4140E+00 -2.7510E-01  1.0851E+00  1.2351E+00  6.4311E-01 -4.5042E+00  1.9352E+00 -4.7990E+00
             2.1499E+00
 GRADIENT:   1.3116E-03 -1.6984E-02  8.0938E-04 -1.4263E-02 -4.8792E-03 -5.8252E-02  1.8268E-03 -1.8880E-04 -2.5429E-02  0.0000E+00
             1.7668E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1661
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8045E-03 -8.0400E-02  1.5275E-06  4.1082E-02 -3.2729E-05
 SE:             2.9340E-02  1.5175E-02  2.5025E-06  2.3824E-02  1.4548E-04
 N:                     100         100         100         100         100

 P VAL.:         8.4317E-01  1.1721E-07  5.4161E-01  8.4639E-02  8.2200E-01

 ETASHRINKSD(%)  1.7083E+00  4.9162E+01  9.9992E+01  2.0186E+01  9.9513E+01
 ETASHRINKVR(%)  3.3874E+00  7.4154E+01  1.0000E+02  3.6298E+01  9.9998E+01
 EBVSHRINKSD(%)  1.5528E+00  4.8779E+01  9.9968E+01  1.6326E+01  9.9487E+01
 EBVSHRINKVR(%)  3.0814E+00  7.3764E+01  1.0000E+02  2.9986E+01  9.9997E+01
 RELATIVEINF(%)  9.6835E+01  1.2036E+01  1.0160E-05  3.2288E+01  2.5194E-03
 EPSSHRINKSD(%)  7.8878E+00
 EPSSHRINKVR(%)  1.5153E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1355.6308572883231     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       298.45850248008765     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    49.97
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1355.631       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.31E+00  1.92E+00  5.52E+02  6.87E-01  2.68E+00  3.11E+00  1.72E+00  1.00E-02  6.27E+00  1.00E-02  7.77E+00
 


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
+        6.80E+01
 
 TH 2
+       -7.53E-01  3.47E+01
 
 TH 3
+        6.77E-04  1.51E-05  3.83E-06
 
 TH 4
+        7.96E+00  2.06E+01 -4.63E-04  4.00E+01
 
 TH 5
+       -3.93E-01 -5.60E+00 -1.82E-03 -4.45E+00  3.74E+01
 
 TH 6
+        6.58E-01  2.43E-01 -4.38E-06 -7.14E-01 -3.87E-02  1.90E+01
 
 TH 7
+        8.86E-01 -1.19E+01 -2.25E-04  2.50E+00  2.08E+00 -8.39E-02  1.14E+01
 
 TH 8
+       -9.17E+01 -4.72E-01 -8.13E-04 -8.82E+01 -6.74E+00  2.38E+00 -2.26E+01  2.55E+02
 
 TH 9
+        1.52E-01 -3.99E+00  9.16E-05  4.63E+00  4.09E-01 -8.66E-02  1.64E+00 -1.31E+00  2.76E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.91E+00 -3.81E+00 -3.64E-05 -2.94E+00 -7.39E-02  8.58E-01  2.37E+00  1.60E+00  3.08E-01  0.00E+00  1.69E+01
 
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
 #CPUT: Total CPU Time in Seconds,       67.310
Stop Time:
Sat Sep 18 07:27:26 CDT 2021
