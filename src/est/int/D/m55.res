Wed Sep 29 09:09:09 CDT 2021
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
$DATA ../../../../data/int/D/dat55.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m55.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   37222.0362110942        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.0087E+02  6.0992E+02 -2.9276E+01  4.3714E+02  2.2009E+02 -2.3372E+03 -1.0928E+03 -9.4519E+01 -1.7753E+03 -7.2503E+02
            -7.5287E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -807.136261706297        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1685E+00  1.6428E+00  8.6424E-01  1.9797E+00  9.1870E-01  3.9435E+00  3.3516E+00  9.7881E-01  2.1855E+00  1.4275E+00
             1.3251E+01
 PARAMETER:  2.5570E-01  5.9642E-01 -4.5899E-02  7.8292E-01  1.5202E-02  1.4721E+00  1.3094E+00  7.8579E-02  8.8183E-01  4.5590E-01
             2.6841E+00
 GRADIENT:  -1.9379E+01  1.3904E+01 -3.7183E+01  1.2825E+02 -3.2407E+00  1.6047E+02 -8.4534E+00  4.7225E+00  1.7661E+01  3.4992E+01
             3.7306E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -879.633291099375        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.0977E+00  1.1959E+00  1.9781E+01  2.4658E+00  2.1227E+00  4.3887E+00  7.5827E+00  5.8045E-01  2.3939E+00  1.1009E+00
             1.2814E+01
 PARAMETER:  1.9321E-01  2.7886E-01  3.0847E+00  1.0025E+00  8.5271E-01  1.5790E+00  2.1259E+00 -4.4396E-01  9.7292E-01  1.9615E-01
             2.6506E+00
 GRADIENT:  -2.0584E+01  2.4510E+01 -2.7120E+00  7.6297E+01 -7.2411E+00  1.7374E+02  7.7456E+01  4.5452E-02  3.9096E+01  1.7510E+01
             3.9722E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -979.364003533965        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.2289E+00  1.1878E+00  1.4421E+01  1.4066E+00  2.2662E+00  1.8507E+00  4.4143E+00  3.0734E+00  1.1686E+00  1.2079E+00
             1.2126E+01
 PARAMETER:  3.0612E-01  2.7207E-01  2.7687E+00  4.4118E-01  9.1811E-01  7.1556E-01  1.5849E+00  1.2228E+00  2.5578E-01  2.8890E-01
             2.5954E+00
 GRADIENT:   2.5569E+01 -1.3568E+01 -2.7186E+00 -5.4892E+00  3.9959E+00  1.6302E+01  1.3005E+01  1.4222E+00  8.1089E+00  2.0602E+01
             3.3573E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1011.47126371688        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.0560E+00  1.7033E+00  6.6668E+00  8.9599E-01  2.1985E+00  1.7710E+00  3.5225E+00  8.3840E-01  9.5346E-01  3.8675E-01
             9.6143E+00
 PARAMETER:  1.5449E-01  6.3255E-01  1.9971E+00 -9.8296E-03  8.8779E-01  6.7155E-01  1.3592E+00 -7.6258E-02  5.2337E-02 -8.4998E-01
             2.3633E+00
 GRADIENT:  -1.3622E+01 -5.1366E+00 -1.6189E+00 -1.3749E+01  3.4128E+00 -8.1081E+00  3.1572E+00  1.4298E-01  5.3678E+00  1.6773E+00
             7.8086E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1021.28116753037        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      413
 NPARAMETR:  1.1299E+00  1.3628E+00  1.1911E+01  1.1554E+00  2.2880E+00  1.9404E+00  4.3583E+00  4.1240E-01  7.4127E-01  2.2171E-01
             9.8287E+00
 PARAMETER:  2.2210E-01  4.0956E-01  2.5775E+00  2.4446E-01  9.2769E-01  7.6292E-01  1.5721E+00 -7.8576E-01 -1.9940E-01 -1.4064E+00
             2.3853E+00
 GRADIENT:   2.3273E+00 -2.4403E+00 -1.2754E+00  6.0131E+00  2.4883E+00 -2.0729E+00 -2.4033E+01  2.4938E-02 -1.6038E+00  4.8929E-01
            -1.2358E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1026.90551076420        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  1.0925E+00  1.0455E+00  2.6179E+01  1.3782E+00  2.4036E+00  1.8421E+00  5.3283E+00  3.1765E-01  9.8288E-01  2.4121E-01
             9.9226E+00
 PARAMETER:  1.8847E-01  1.4450E-01  3.3650E+00  4.2076E-01  9.7696E-01  7.1093E-01  1.7730E+00 -1.0468E+00  8.2732E-02 -1.3221E+00
             2.3948E+00
 GRADIENT:  -5.9250E+00  1.1335E+00 -6.5876E-01 -2.9951E+00  2.1976E+00 -2.0530E+00 -1.0671E+00  3.0776E-03 -6.5627E-01  5.5062E-01
             2.9631E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1028.67822221327        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      766
 NPARAMETR:  1.1081E+00  7.1261E-01  5.5373E+01  1.5984E+00  2.4552E+00  1.8649E+00  6.0552E+00  1.8877E-01  1.1941E+00  2.1548E-01
             9.8921E+00
 PARAMETER:  2.0261E-01 -2.3882E-01  4.1141E+00  5.6898E-01  9.9820E-01  7.2321E-01  1.9009E+00 -1.5672E+00  2.7743E-01 -1.4349E+00
             2.3917E+00
 GRADIENT:   1.1264E+00 -7.8183E-03 -1.3973E-01  7.7164E-01 -5.5036E-04  3.0406E-01  1.6122E-01  1.3110E-04 -6.2842E-02  4.1599E-01
            -7.8744E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1028.88234180940        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  1.0991E+00  6.9668E-01  6.6243E+01  1.5979E+00  2.4722E+00  1.8550E+00  6.0795E+00  1.0211E-01  1.1857E+00  8.9948E-02
             9.8888E+00
 PARAMETER:  1.9452E-01 -2.6144E-01  4.2933E+00  5.6869E-01  1.0051E+00  7.1788E-01  1.9049E+00 -2.1817E+00  2.7036E-01 -2.3085E+00
             2.3914E+00
 GRADIENT:  -2.2859E+00 -7.4442E-01 -9.4980E-02 -7.4058E-01 -3.1230E-04 -1.8539E-01 -1.3197E-01  2.4843E-05 -2.9152E-01  7.0965E-02
            -4.9176E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1028.88615128923        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1119
 NPARAMETR:  1.1026E+00  7.0678E-01  6.9190E+01  1.5970E+00  2.4790E+00  1.8499E+00  6.0579E+00  7.9987E-02  1.1877E+00  6.3401E-02
             9.9018E+00
 PARAMETER:  1.9763E-01 -2.4704E-01  4.3369E+00  5.6812E-01  1.0078E+00  7.1512E-01  1.9014E+00 -2.4259E+00  2.7201E-01 -2.6583E+00
             2.3927E+00
 GRADIENT:  -7.3332E-01 -4.9898E-01 -9.5232E-02 -2.2870E-01  2.3034E-01 -7.9605E-01 -2.6067E-01  1.3786E-05 -1.1377E-01  3.5362E-02
             6.6289E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1028.88690190434        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1294
 NPARAMETR:  1.1036E+00  7.1445E-01  7.2304E+01  1.5941E+00  2.4827E+00  1.8529E+00  6.0443E+00  6.3077E-02  1.1869E+00  4.5669E-02
             9.9017E+00
 PARAMETER:  1.9859E-01 -2.3624E-01  4.3809E+00  5.6632E-01  1.0094E+00  7.1678E-01  1.8991E+00 -2.6634E+00  2.7136E-01 -2.9863E+00
             2.3927E+00
 GRADIENT:  -3.6837E-01 -3.1782E-01 -8.5277E-02 -1.0020E-01  2.1655E-01 -5.1250E-01 -1.3766E-01  7.9615E-06 -4.8244E-02  1.8326E-02
             4.7340E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1028.88695424878        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1469
 NPARAMETR:  1.1041E+00  7.1913E-01  7.5815E+01  1.5922E+00  2.4862E+00  1.8550E+00  6.0362E+00  5.0611E-02  1.1861E+00  3.3767E-02
             9.9011E+00
 PARAMETER:  1.9899E-01 -2.2972E-01  4.4283E+00  5.6510E-01  1.0108E+00  7.1788E-01  1.8978E+00 -2.8836E+00  2.7066E-01 -3.2883E+00
             2.3926E+00
 GRADIENT:  -2.2172E-01 -2.1567E-01 -7.5363E-02 -5.6599E-02  1.8739E-01 -3.0156E-01 -5.5731E-02  4.5477E-06 -1.6183E-02  1.0019E-02
             2.8744E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1029.03167147084        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1632
 NPARAMETR:  1.1042E+00  7.2346E-01  1.0314E+02  1.5909E+00  2.4863E+00  1.8583E+00  6.2728E+00  4.7752E-02  1.1851E+00  1.0500E-02
             9.8942E+00
 PARAMETER:  1.9913E-01 -2.2371E-01  4.7361E+00  5.6430E-01  1.0108E+00  7.1964E-01  1.9362E+00 -2.9417E+00  2.6985E-01 -4.4564E+00
             2.3919E+00
 GRADIENT:   5.4912E-01  1.5546E+00  1.2934E-02 -4.1751E+00 -2.4905E+00  4.9927E-01  8.3560E+00  1.1447E-06  9.2245E-01  9.4897E-04
            -9.5103E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1029.13058319199        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1810
 NPARAMETR:  1.1046E+00  6.7910E-01  1.3261E+02  1.6247E+00  2.5069E+00  1.8562E+00  6.2540E+00  4.7862E-02  1.2060E+00  1.0000E-02
             9.9099E+00
 PARAMETER:  1.9945E-01 -2.8699E-01  4.9874E+00  5.8532E-01  1.0191E+00  7.1854E-01  1.9332E+00 -2.9394E+00  2.8729E-01 -4.5537E+00
             2.3935E+00
 GRADIENT:   2.4717E-01  4.0285E-01  3.1038E-03  1.7548E-01 -1.3061E+00  4.1650E-01  3.5454E+00  7.6930E-07  1.3663E-01  0.0000E+00
             1.5883E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1029.21579709277        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1969
 NPARAMETR:  1.1037E+00  6.6754E-01  1.6543E+02  1.6279E+00  2.5188E+00  1.8536E+00  6.4207E+00  4.7816E-02  1.2083E+00  1.0000E-02
             9.8967E+00
 PARAMETER:  1.9864E-01 -3.0415E-01  5.2086E+00  5.8726E-01  1.0238E+00  7.1712E-01  1.9595E+00 -2.9404E+00  2.8921E-01 -4.5568E+00
             2.3922E+00
 GRADIENT:   1.0523E+01  3.1050E+00 -1.5572E-03  1.5474E+01  2.9515E+00  1.2733E+01  1.2235E+02  7.6854E-07  1.8093E+00  0.0000E+00
             4.2677E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1029.26279699836        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2152
 NPARAMETR:  1.1037E+00  6.3879E-01  1.8227E+02  1.6473E+00  2.5234E+00  1.8528E+00  6.4255E+00  3.2742E-02  1.2096E+00  1.0000E-02
             9.9057E+00
 PARAMETER:  1.9862E-01 -3.4817E-01  5.3055E+00  5.9912E-01  1.0256E+00  7.1671E-01  1.9603E+00 -3.3191E+00  2.9029E-01 -4.5568E+00
             2.3931E+00
 GRADIENT:   3.3992E-01  2.2564E-01 -1.8075E-03 -2.5067E-01 -3.6538E-01  4.8864E-01  5.6172E+00  2.6614E-07 -2.9292E-01  0.0000E+00
             8.5341E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1029.33162242504        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2332
 NPARAMETR:  1.1036E+00  6.0779E-01  2.1483E+02  1.6652E+00  2.5289E+00  1.8533E+00  6.5685E+00  3.0694E-02  1.2132E+00  1.0000E-02
             9.9168E+00
 PARAMETER:  1.9862E-01 -3.9792E-01  5.4699E+00  6.0993E-01  1.0278E+00  7.1698E-01  1.9823E+00 -3.3837E+00  2.9323E-01 -4.5568E+00
             2.3942E+00
 GRADIENT:   3.8673E-01  9.3649E-02 -6.8332E-04 -7.1536E-01 -1.4625E-01  9.9660E-01  7.3566E+00  6.2365E-08 -5.6761E-01  0.0000E+00
             3.1324E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1029.36588678240        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2510
 NPARAMETR:  1.1030E+00  5.9284E-01  2.3141E+02  1.6790E+00  2.5317E+00  1.8513E+00  6.6455E+00  3.1573E-02  1.2249E+00  1.0000E-02
             9.9141E+00
 PARAMETER:  1.9806E-01 -4.2283E-01  5.5442E+00  6.1820E-01  1.0289E+00  7.1587E-01  1.9939E+00 -3.3555E+00  3.0287E-01 -4.5568E+00
             2.3940E+00
 GRADIENT:   3.4789E-01  2.7598E-01  1.4153E-05 -9.8574E-02 -1.4662E-01  9.3668E-01  8.2315E+00 -1.0783E-07 -4.5750E-01  0.0000E+00
             2.4410E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1029.39888990932        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2696             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1003E+00  5.6298E-01  2.4454E+02  1.6892E+00  2.5325E+00  1.8429E+00  6.7496E+00  2.6666E-02  1.2351E+00  1.0000E-02
             9.8939E+00
 PARAMETER:  1.9558E-01 -4.7451E-01  5.5994E+00  6.2428E-01  1.0292E+00  7.1136E-01  2.0095E+00 -3.5244E+00  3.1115E-01 -4.5568E+00
             2.3919E+00
 GRADIENT:   9.7442E+00  2.9397E+00  1.2144E-03  1.8801E+01  3.4711E+00  1.2378E+01  1.3638E+02 -1.2260E-07  1.0806E+00  0.0000E+00
             4.2879E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1029.40964709259        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2875            RESET HESSIAN, TYPE II
 NPARAMETR:  1.1004E+00  5.6202E-01  2.4740E+02  1.6981E+00  2.5323E+00  1.8432E+00  6.7815E+00  2.6840E-02  1.2403E+00  1.0000E-02
             9.8975E+00
 PARAMETER:  1.9572E-01 -4.7622E-01  5.6110E+00  6.2953E-01  1.0291E+00  7.1151E-01  2.0142E+00 -3.5179E+00  3.1532E-01 -4.5568E+00
             2.3923E+00
 GRADIENT:   9.7512E+00  3.3783E+00  3.9958E-03  2.0952E+01  3.0041E+00  1.2506E+01  1.3694E+02  3.8076E-08  9.2305E-01  0.0000E+00
             4.2668E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1029.41918877403        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3058             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1006E+00  5.4768E-01  2.4896E+02  1.7052E+00  2.5328E+00  1.8399E+00  6.8307E+00  2.2727E-02  1.2463E+00  1.0000E-02
             9.8932E+00
 PARAMETER:  1.9582E-01 -5.0206E-01  5.6173E+00  6.3369E-01  1.0293E+00  7.0972E-01  2.0214E+00 -3.6842E+00  3.2015E-01 -4.5568E+00
             2.3918E+00
 GRADIENT:   1.0112E+01  3.3105E+00  3.0783E-03  2.0721E+01  3.1802E+00  1.2152E+01  1.3940E+02  1.2899E-07  1.0726E+00  0.0000E+00
             4.2125E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1029.42229644145        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3242             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1002E+00  5.4600E-01  2.5592E+02  1.7109E+00  2.5324E+00  1.8416E+00  6.8532E+00  2.3643E-02  1.2484E+00  1.0000E-02
             9.8958E+00
 PARAMETER:  1.9547E-01 -5.0513E-01  5.6449E+00  6.3700E-01  1.0292E+00  7.1064E-01  2.0247E+00 -3.6447E+00  3.2187E-01 -4.5568E+00
             2.3921E+00
 GRADIENT:   9.8024E+00  3.5518E+00  5.8448E-03  2.2072E+01  2.7987E+00  1.2450E+01  1.3986E+02  2.4762E-07  8.9342E-01  0.0000E+00
             4.2275E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1029.42562801457        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3425             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1000E+00  5.4122E-01  2.5401E+02  1.7146E+00  2.5325E+00  1.8418E+00  6.8788E+00  2.2813E-02  1.2528E+00  1.0000E-02
             9.8958E+00
 PARAMETER:  1.9533E-01 -5.1393E-01  5.6374E+00  6.3915E-01  1.0292E+00  7.1074E-01  2.0284E+00 -3.6804E+00  3.2541E-01 -4.5568E+00
             2.3921E+00
 GRADIENT:   9.7469E+00  3.6112E+00  5.2881E-03  2.2068E+01  2.8723E+00  1.2517E+01  1.4096E+02 -3.1575E-08  1.0615E+00  0.0000E+00
             4.2448E+01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1029.42738109263        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     3606
 NPARAMETR:  1.1000E+00  5.3751E-01  2.5211E+02  1.7163E+00  2.5323E+00  1.8410E+00  6.8919E+00  2.2678E-02  1.2551E+00  1.0000E-02
             9.8962E+00
 PARAMETER:  1.9531E-01 -5.2081E-01  5.6299E+00  6.4017E-01  1.0291E+00  7.1033E-01  2.0303E+00 -3.6864E+00  3.2719E-01 -4.5568E+00
             2.3922E+00
 GRADIENT:   1.0179E-01  1.5852E-01  3.2946E-03 -7.1289E-02 -3.7986E-01  3.2723E-01  1.0605E+01  4.4542E-08 -1.5602E-01  0.0000E+00
            -4.0706E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1029.42920972274        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3796             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0998E+00  5.3165E-01  2.3272E+02  1.7189E+00  2.5332E+00  1.8406E+00  6.9103E+00  2.1359E-02  1.2585E+00  1.0000E-02
             9.8973E+00
 PARAMETER:  1.9517E-01 -5.3176E-01  5.5498E+00  6.4169E-01  1.0295E+00  7.1012E-01  2.0330E+00 -3.7463E+00  3.2988E-01 -4.5568E+00
             2.3923E+00
 GRADIENT:   9.7501E+00  3.5059E+00 -3.7096E-04  2.1439E+01  3.4476E+00  1.2475E+01  1.4261E+02  2.8316E-07  1.3465E+00  0.0000E+00
             4.3187E+01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1029.42981828289        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3980
 NPARAMETR:  1.0998E+00  5.3020E-01  2.3447E+02  1.7203E+00  2.5331E+00  1.8405E+00  6.9183E+00  1.9416E-02  1.2593E+00  1.0000E-02
             9.8971E+00
 PARAMETER:  1.9514E-01 -5.3450E-01  5.5573E+00  6.4252E-01  1.0295E+00  7.1006E-01  2.0342E+00 -3.8417E+00  3.3059E-01 -4.5568E+00
             2.3922E+00
 GRADIENT:   6.5352E-02  3.9036E-02 -1.9749E-03 -4.6010E-01  1.0021E-01  3.5031E-01  1.0782E+01  2.1422E-08 -1.3701E-02  0.0000E+00
             5.4623E-02

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1029.43059325036        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     4168             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0997E+00  5.2731E-01  2.4767E+02  1.7234E+00  2.5320E+00  1.8401E+00  6.9330E+00  1.6053E-02  1.2610E+00  1.0000E-02
             9.8964E+00
 PARAMETER:  1.9507E-01 -5.3997E-01  5.6121E+00  6.4429E-01  1.0290E+00  7.0982E-01  2.0363E+00 -4.0318E+00  3.3188E-01 -4.5568E+00
             2.3922E+00
 GRADIENT:   9.7578E+00  3.6126E+00  4.9044E-03  2.2153E+01  2.9274E+00  1.2439E+01  1.4342E+02 -1.1814E-07  1.2737E+00  0.0000E+00
             4.2857E+01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1029.43090986782        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     4347
 NPARAMETR:  1.0997E+00  5.2602E-01  2.4350E+02  1.7240E+00  2.5324E+00  1.8401E+00  6.9373E+00  1.4667E-02  1.2618E+00  1.0000E-02
             9.8965E+00
 PARAMETER:  1.9504E-01 -5.4242E-01  5.5951E+00  6.4467E-01  1.0292E+00  7.0980E-01  2.0369E+00 -4.1222E+00  3.3252E-01 -4.5568E+00
             2.3922E+00
 GRADIENT:   6.5380E-02  5.8501E-02  1.3770E-03 -1.4035E-01 -2.0713E-01  3.2227E-01  1.0899E+01 -1.5134E-07 -5.5116E-02  0.0000E+00
            -9.6976E-02

0ITERATION NO.:  136    OBJECTIVE VALUE:  -1029.43090986782        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     4375
 NPARAMETR:  1.0997E+00  5.2416E-01  2.3025E+02  1.7243E+00  2.5348E+00  1.8403E+00  6.9436E+00  1.5235E-02  1.2635E+00  1.0000E-02
             9.8974E+00
 PARAMETER:  1.9504E-01 -5.4242E-01  5.5951E+00  6.4467E-01  1.0292E+00  7.0980E-01  2.0369E+00 -4.1222E+00  3.3252E-01 -4.5568E+00
             2.3922E+00
 GRADIENT:   9.9632E-03  4.8755E-02  1.4375E-03 -4.6830E-02 -1.8516E-01 -8.2509E-03 -8.3317E-02 -2.9145E-05 -6.0688E-02  0.0000E+00
            -7.7182E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4375
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8318E-02  5.5081E-02  3.0712E-07 -8.5096E-02 -2.1455E-06
 SE:             2.8765E-02  2.1884E-02  1.3572E-06  1.5986E-02  1.2223E-04
 N:                     100         100         100         100         100

 P VAL.:         5.2423E-01  1.1838E-02  8.2098E-01  1.0226E-07  9.8600E-01

 ETASHRINKSD(%)  3.6342E+00  2.6685E+01  9.9995E+01  4.6444E+01  9.9591E+01
 ETASHRINKVR(%)  7.1363E+00  4.6250E+01  1.0000E+02  7.1317E+01  9.9998E+01
 EBVSHRINKSD(%)  6.7870E+00  2.4824E+01  9.9994E+01  4.1672E+01  9.9549E+01
 EBVSHRINKVR(%)  1.3113E+01  4.3486E+01  1.0000E+02  6.5978E+01  9.9998E+01
 RELATIVEINF(%)  8.6654E+01  2.8431E+01  5.6669E-08  1.6821E+01  3.5012E-04
 EPSSHRINKSD(%)  5.6921E+00
 EPSSHRINKVR(%)  1.1060E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1029.4309098678152     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       624.65844990059554     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   161.25
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    18.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1029.431       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.10E+00  5.26E-01  2.43E+02  1.72E+00  2.53E+00  1.84E+00  6.94E+00  1.47E-02  1.26E+00  1.00E-02  9.90E+00
 


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
+        1.93E+02
 
 TH 2
+        1.23E+00  5.01E+01
 
 TH 3
+        1.17E-04  2.47E-04  4.33E-07
 
 TH 4
+       -5.28E+00  3.40E+01  6.13E-04  1.01E+02
 
 TH 5
+        2.47E+00 -6.10E+00 -3.16E-03 -1.18E+01  3.17E+01
 
 TH 6
+       -2.84E+01  1.39E+00 -1.06E-05  1.02E-01  2.84E+00  2.30E+01
 
 TH 7
+        1.45E+00  5.10E+00 -1.91E-05 -4.32E+00  1.46E-01  7.80E-01  1.92E+00
 
 TH 8
+       -9.81E-01 -7.19E-02 -4.04E-05  1.23E-02 -1.35E-02  8.67E-02  1.10E-03 -2.14E-01
 
 TH 9
+        9.69E-01 -2.08E+00 -1.87E-04 -3.16E+01  3.44E+00 -1.28E+00  1.49E+00 -1.79E-01  2.82E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  2.50E+01
 
 TH11
+       -4.31E+00 -4.69E+00 -1.49E-05 -6.90E+00  1.10E-01  4.88E+00 -1.10E-01  2.33E-03  3.03E+00  0.00E+00  9.39E+00
 
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
 #CPUT: Total CPU Time in Seconds,      179.849
Stop Time:
Wed Sep 29 09:12:11 CDT 2021
