Wed Sep 29 12:42:07 CDT 2021
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
$DATA ../../../../data/spa/A2/dat27.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m27.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -887.858676421124        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7088E+02  5.0840E+01  1.5172E+01  5.0015E+01  1.3236E+02  5.0803E+01 -1.8313E+01 -6.4477E+00 -6.7080E+01 -5.8143E+01
            -1.3582E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1354.87570130231        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1104E+00  9.5001E-01  9.9050E-01  1.0576E+00  8.7855E-01  1.1100E+00  9.8170E-01  9.5834E-01  1.2052E+00  8.9507E-01
             2.3008E+00
 PARAMETER:  2.0472E-01  4.8721E-02  9.0454E-02  1.5601E-01 -2.9485E-02  2.0433E-01  8.1528E-02  5.7446E-02  2.8664E-01 -1.0858E-02
             9.3326E-01
 GRADIENT:   3.6619E+02  2.1717E+01  2.2235E+01  1.9329E+01 -2.0455E+01  4.5654E+01  6.5327E+00  3.6503E+00  1.2489E+01  9.6017E+00
            -8.2011E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1367.82761350093        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0827E+00  7.1448E-01  4.4547E-01  1.1759E+00  5.1059E-01  1.0212E+00  1.0776E+00  4.1069E-01  1.1137E+00  4.3677E-01
             2.2626E+00
 PARAMETER:  1.7947E-01 -2.3620E-01 -7.0862E-01  2.6204E-01 -5.7219E-01  1.2101E-01  1.7476E-01 -7.8992E-01  2.0766E-01 -7.2834E-01
             9.1650E-01
 GRADIENT:   3.0998E+02  1.5324E+01 -3.4487E+01  1.2934E+02  7.3579E+01  2.2072E+01 -2.5101E+00  2.0258E+00  1.5159E+01  2.3159E-01
            -8.2661E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1378.94249453315        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0196E+00  7.1442E-01  3.2685E-01  1.1019E+00  4.2711E-01  9.4374E-01  9.5952E-01  2.3028E-01  1.0720E+00  2.7144E-01
             2.3630E+00
 PARAMETER:  1.1946E-01 -2.3628E-01 -1.0183E+00  1.9706E-01 -7.5071E-01  4.2093E-02  5.8679E-02 -1.3685E+00  1.6955E-01 -1.2040E+00
             9.5992E-01
 GRADIENT:   6.0373E+01  2.1961E+01 -5.2901E+00  7.2986E+01 -1.5425E+01  1.5909E+00 -1.6534E+01  1.0173E-02  2.8643E+00 -3.2610E+00
            -6.2816E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1386.30965791457        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      477
 NPARAMETR:  1.0152E+00  6.0325E-01  4.4437E-01  1.1595E+00  4.6589E-01  9.3003E-01  1.2723E+00  4.7292E-01  9.8815E-01  1.6915E-01
             2.6087E+00
 PARAMETER:  1.1507E-01 -4.0542E-01 -7.1109E-01  2.4801E-01 -6.6381E-01  2.7458E-02  3.4083E-01 -6.4882E-01  8.8080E-02 -1.6770E+00
             1.0589E+00
 GRADIENT:   4.4315E+01  8.0274E+00  1.1958E+01 -3.2423E+00 -3.1453E+01  1.9006E+00 -2.4733E+00  7.5397E-01  9.3971E-01 -5.2331E-01
            -1.0700E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1389.95993206881        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  9.8957E-01  4.9563E-01  7.3468E-01  1.2925E+00  6.1754E-01  9.0794E-01  1.4277E+00  4.8060E-01  9.3861E-01  1.1738E-01
             2.8031E+00
 PARAMETER:  8.9519E-02 -6.0192E-01 -2.0831E-01  3.5657E-01 -3.8202E-01  3.4233E-03  4.5605E-01 -6.3271E-01  3.6646E-02 -2.0424E+00
             1.1307E+00
 GRADIENT:  -3.4504E+00  2.6663E+00  3.0592E+00  4.8737E-01 -4.3242E+00  3.7262E-01  1.4544E+00 -3.1673E-01 -2.5966E-01 -8.3892E-02
            -2.6377E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1390.02740412037        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      829
 NPARAMETR:  9.9064E-01  4.8203E-01  7.4317E-01  1.3020E+00  6.2031E-01  9.0605E-01  1.3522E+00  4.5702E-01  9.3639E-01  1.2577E-01
             2.8311E+00
 PARAMETER:  9.0601E-02 -6.2975E-01 -1.9683E-01  3.6390E-01 -3.7753E-01  1.3361E-03  4.0177E-01 -6.8302E-01  3.4281E-02 -1.9733E+00
             1.1407E+00
 GRADIENT:  -2.9664E-01  1.1808E+00  1.0284E+00  1.4932E+00 -1.6672E+00 -1.9480E-02 -2.1611E-01 -1.0158E-01 -4.1214E-01 -9.1349E-02
             6.2582E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1390.35296741456        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1009
 NPARAMETR:  9.8922E-01  2.3842E-01  6.8985E-01  1.4296E+00  5.3067E-01  9.1000E-01  2.1534E+00  5.4846E-01  8.7893E-01  2.8129E-01
             2.7539E+00
 PARAMETER:  8.9160E-02 -1.3337E+00 -2.7128E-01  4.5742E-01 -5.3362E-01  5.6890E-03  8.6703E-01 -5.0064E-01 -2.9045E-02 -1.1684E+00
             1.1130E+00
 GRADIENT:   6.1547E+00  2.5317E+00  4.8331E+00  1.6033E+01 -1.0043E+01  8.0132E-01 -1.2010E-01  7.0996E-01 -6.1300E-02 -4.9489E-01
             1.2554E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1391.27423630447        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1185
 NPARAMETR:  9.8028E-01  1.5288E-01  6.2004E-01  1.4457E+00  4.7830E-01  9.0720E-01  2.7939E+00  3.1401E-01  8.6884E-01  5.2556E-01
             2.6807E+00
 PARAMETER:  8.0084E-02 -1.7781E+00 -3.7796E-01  4.6862E-01 -6.3753E-01  2.6079E-03  1.1274E+00 -1.0583E+00 -4.0596E-02 -5.4329E-01
             1.0861E+00
 GRADIENT:  -1.1809E+01  9.3119E-01 -2.4646E+00  4.7005E+00  2.2933E+00 -1.2459E+00  2.1362E-01  7.2129E-01 -1.0367E+00  1.4076E+00
             6.2894E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1391.72416142312        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1360
 NPARAMETR:  9.8293E-01  1.1340E-01  5.7608E-01  1.4489E+00  4.4638E-01  9.1245E-01  3.2846E+00  4.9195E-02  8.6889E-01  5.6956E-01
             2.6309E+00
 PARAMETER:  8.2785E-02 -2.0768E+00 -4.5152E-01  4.7080E-01 -7.0660E-01  8.3759E-03  1.2892E+00 -2.9120E+00 -4.0533E-02 -4.6289E-01
             1.0673E+00
 GRADIENT:  -1.9007E+00  8.2663E-01  2.0160E+00  5.0201E+00 -4.3368E+00  1.7597E-01  7.5188E-02  1.3142E-02 -1.3593E+00  3.7109E-01
            -1.1006E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1391.88777966285        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1535
 NPARAMETR:  9.8032E-01  3.8717E-02  5.7882E-01  1.4830E+00  4.3776E-01  9.1153E-01  4.9493E+00  1.0000E-02  8.5847E-01  5.7684E-01
             2.6366E+00
 PARAMETER:  8.0123E-02 -3.1515E+00 -4.4677E-01  4.9404E-01 -7.2607E-01  7.3728E-03  1.6993E+00 -1.0013E+01 -5.2598E-02 -4.5019E-01
             1.0695E+00
 GRADIENT:  -6.6753E-01  9.8677E-02  1.0315E+00  3.1911E+00 -1.9878E+00  4.3010E-01 -1.6649E-01  0.0000E+00 -3.0444E-01  1.1210E-01
             2.8338E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1391.94838896702        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1710
 NPARAMETR:  9.7922E-01  1.0000E-02  5.7365E-01  1.4927E+00  4.3104E-01  9.1019E-01  9.5160E+00  1.0000E-02  8.5479E-01  5.8011E-01
             2.6313E+00
 PARAMETER:  7.9005E-02 -4.6962E+00 -4.5574E-01  5.0058E-01 -7.4156E-01  5.8968E-03  2.3530E+00 -2.0817E+01 -5.6895E-02 -4.4453E-01
             1.0675E+00
 GRADIENT:  -2.6406E-01  0.0000E+00  8.9981E-01  1.2710E+00 -1.5541E+00  5.8849E-02 -4.5557E-02  0.0000E+00 -3.5942E-02 -4.4800E-02
            -9.5768E-02

0ITERATION NO.:   59    OBJECTIVE VALUE:  -1391.95505521672        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1838
 NPARAMETR:  9.7928E-01  1.0000E-02  5.7341E-01  1.4918E+00  4.3117E-01  9.1001E-01  1.1182E+01  1.0000E-02  8.5470E-01  5.7959E-01
             2.6316E+00
 PARAMETER:  7.9064E-02 -5.0616E+00 -4.5615E-01  4.9999E-01 -7.4124E-01  5.7009E-03  2.5143E+00 -2.3429E+01 -5.7002E-02 -4.4543E-01
             1.0676E+00
 GRADIENT:   6.0685E-03  0.0000E+00 -1.4461E-03 -2.1906E-01  3.4907E-02  1.4788E-02  7.2883E-03  0.0000E+00  2.2984E-02 -3.4233E-02
            -5.3603E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1838
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.4621E-04  6.9358E-04  7.9663E-05 -1.1353E-02 -8.6515E-03
 SE:             2.8958E-02  1.8239E-03  2.2272E-04  2.7028E-02  1.8029E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9322E-01  7.0375E-01  7.2058E-01  6.7445E-01  6.3131E-01

 ETASHRINKSD(%)  2.9864E+00  9.3890E+01  9.9254E+01  9.4532E+00  3.9602E+01
 ETASHRINKVR(%)  5.8836E+00  9.9627E+01  9.9994E+01  1.8013E+01  6.3521E+01
 EBVSHRINKSD(%)  2.9170E+00  9.4285E+01  9.9238E+01  9.0830E+00  3.9546E+01
 EBVSHRINKVR(%)  5.7490E+00  9.9673E+01  9.9994E+01  1.7341E+01  6.3453E+01
 RELATIVEINF(%)  8.0906E+01  1.6193E-02  2.2316E-04  8.0692E+00  1.0148E+00
 EPSSHRINKSD(%)  2.9952E+01
 EPSSHRINKVR(%)  5.0932E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1391.9550552167238     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -656.80422865298567     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.35
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1391.955       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.00E-02  5.73E-01  1.49E+00  4.31E-01  9.10E-01  1.12E+01  1.00E-02  8.55E-01  5.80E-01  2.63E+00
 


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
+        1.32E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.00E+01  0.00E+00  1.92E+03
 
 TH 4
+       -5.01E+01  0.00E+00 -1.17E+02  5.64E+02
 
 TH 5
+        1.07E+02  0.00E+00 -3.33E+03 -2.07E+02  6.28E+03
 
 TH 6
+       -8.06E-01  0.00E+00  1.24E+01 -1.19E+01 -1.21E+00  2.12E+02
 
 TH 7
+        5.69E-04  0.00E+00  3.61E-03 -1.49E-02  3.14E-02  2.55E-03  3.88E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.09E+01  0.00E+00  3.28E+01 -1.46E+01  8.81E+00  8.72E-01  1.92E-02  0.00E+00  1.88E+02
 
 TH10
+       -1.12E+01  0.00E+00 -3.85E+01 -5.93E+00  5.32E+01 -3.72E-01  4.76E-03  0.00E+00  1.49E-01  7.10E+01
 
 TH11
+       -1.62E+01  0.00E+00 -7.41E+00 -6.18E+00 -1.30E+01  4.06E+00  1.17E-03  0.00E+00  9.70E+00  3.05E+01  4.46E+01
 
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
 #CPUT: Total CPU Time in Seconds,       29.792
Stop Time:
Wed Sep 29 12:42:40 CDT 2021
