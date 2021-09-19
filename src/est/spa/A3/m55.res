Sat Sep 18 10:33:19 CDT 2021
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
$DATA ../../../../data/spa/A3/dat55.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -103.622808301138        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.0540E+01  1.3788E+01  7.4327E+01 -5.2046E+01  2.3270E+02  5.6584E+01 -6.1510E+01 -4.3174E+01 -9.4337E+01 -2.0546E+02
            -2.7147E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1253.64661937288        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0532E+00  1.0000E+00  9.1237E-01  1.1254E+00  8.7676E-01  7.5801E-01  1.0034E+00  1.0126E+00  9.9269E-01  1.1212E+00
             4.8944E+00
 PARAMETER:  1.5187E-01  1.0000E-01  8.2937E-03  2.1817E-01 -3.1518E-02 -1.7706E-01  1.0335E-01  1.1251E-01  9.2666E-02  2.1442E-01
             1.6881E+00
 GRADIENT:   4.8642E+01 -1.3600E+01 -1.1649E+01 -1.0609E+01 -1.2307E+01 -9.8301E-01  1.3699E+01  5.9751E+00  2.7105E+01  2.9070E+01
             1.5823E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1269.55711030165        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0454E+00  6.7360E-01  4.9307E-01  1.3241E+00  5.0567E-01  7.9338E-01  8.1333E-01  3.4463E-01  7.2916E-01  6.0222E-01
             4.6429E+00
 PARAMETER:  1.4443E-01 -2.9513E-01 -6.0710E-01  3.8072E-01 -5.8186E-01 -1.3145E-01 -1.0661E-01 -9.6528E-01 -2.1586E-01 -4.0714E-01
             1.6353E+00
 GRADIENT:  -1.1790E+01  3.1835E+01 -2.4511E+01  1.1811E+02 -8.0569E+00 -3.4180E+00  1.9179E+00  2.1485E+00  4.5506E+00  1.4964E+01
             1.1004E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1293.18707175019        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0207E+00  7.7771E-01  1.0325E+00  1.2450E+00  8.9456E-01  7.8778E-01  8.8111E-01  1.9571E-01  7.1142E-01  4.4273E-01
             4.0138E+00
 PARAMETER:  1.2048E-01 -1.5141E-01  1.3195E-01  3.1916E-01 -1.1428E-02 -1.3854E-01 -2.6570E-02 -1.5311E+00 -2.4049E-01 -7.1480E-01
             1.4897E+00
 GRADIENT:   1.6974E+01 -5.9436E+00 -1.2092E+01  2.1036E-01  1.6651E+01  2.1745E+00  2.1097E+00  2.2196E-01  2.7170E+00  3.3944E+00
             2.8677E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1295.00024972987        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0134E+00  7.4968E-01  1.3430E+00  1.2652E+00  9.9161E-01  7.8229E-01  6.8823E-01  1.2681E-01  6.7391E-01  1.8263E-01
             4.0397E+00
 PARAMETER:  1.1333E-01 -1.8810E-01  3.9493E-01  3.3526E-01  9.1570E-02 -1.4553E-01 -2.7363E-01 -1.9651E+00 -2.9466E-01 -1.6003E+00
             1.4962E+00
 GRADIENT:   6.2990E+00 -1.3271E+00  1.6695E+00 -3.0755E+00 -2.4948E+00  1.2378E+00 -9.5577E-02  5.1568E-02 -1.5075E-01  3.2019E-01
            -1.6093E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1295.27588947702        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0141E+00  1.0007E+00  1.2358E+00  1.1107E+00  1.0525E+00  7.7986E-01  5.6571E-01  1.0196E-01  7.4807E-01  4.3974E-02
             4.0506E+00
 PARAMETER:  1.1405E-01  1.0071E-01  3.1173E-01  2.0499E-01  1.5114E-01 -1.4864E-01 -4.6967E-01 -2.1832E+00 -1.9025E-01 -3.0242E+00
             1.4989E+00
 GRADIENT:  -4.7181E-01  2.5192E+00  8.6310E-01  3.7466E+00 -2.1113E+00 -1.6079E-01 -3.3681E-01  2.9267E-02 -4.4084E-01  1.6566E-02
            -1.8567E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1295.31477531927        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.0150E+00  1.0953E+00  1.1542E+00  1.0474E+00  1.0637E+00  7.8013E-01  5.7587E-01  9.0418E-02  7.7599E-01  2.7100E-02
             4.0542E+00
 PARAMETER:  1.1487E-01  1.9105E-01  2.4340E-01  1.4630E-01  1.6178E-01 -1.4829E-01 -4.5188E-01 -2.3033E+00 -1.5361E-01 -3.5082E+00
             1.4998E+00
 GRADIENT:  -1.6545E-02 -6.8281E-02 -2.7797E-01  1.9848E-01  3.7516E-01 -1.6006E-02 -2.7291E-03  2.2551E-02 -4.7002E-02  6.7466E-03
            -1.9240E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1295.32764143548        NO. OF FUNC. EVALS.: 113
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  1.0157E+00  1.1360E+00  1.1772E+00  1.0235E+00  1.0899E+00  7.8037E-01  5.6320E-01  7.0602E-02  7.8492E-01  1.8959E-02
             4.0618E+00
 PARAMETER:  1.1563E-01  2.2752E-01  2.6316E-01  1.2327E-01  1.8613E-01 -1.4799E-01 -4.7411E-01 -2.5507E+00 -1.4217E-01 -3.8655E+00
             1.5016E+00
 GRADIENT:  -1.0266E+00  4.1664E-01  1.8536E-01  4.1919E-01 -3.7591E-01 -9.1673E-02 -2.5977E-03  1.2191E-02  3.1699E-02  3.1559E-03
            -4.9702E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1295.33346619824        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      729
 NPARAMETR:  1.0165E+00  1.2067E+00  1.1493E+00  9.7812E-01  1.1120E+00  7.8067E-01  5.4588E-01  4.6193E-02  8.0794E-01  1.1126E-02
             4.0663E+00
 PARAMETER:  1.1635E-01  2.8792E-01  2.3911E-01  7.7874E-02  2.0618E-01 -1.4760E-01 -5.0535E-01 -2.9749E+00 -1.1327E-01 -4.3985E+00
             1.5027E+00
 GRADIENT:  -7.3134E-02 -9.0231E-03  2.3761E-02 -6.1879E-02 -1.2971E-02 -8.3264E-04 -9.8188E-04  4.7845E-03  1.2300E-02  1.0781E-03
            -5.0187E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1295.33475555912        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  1.0164E+00  1.1855E+00  1.1553E+00  9.9171E-01  1.1045E+00  7.8061E-01  5.5119E-01  3.5500E-02  8.0086E-01  1.1916E-02
             4.0654E+00
 PARAMETER:  1.1624E-01  2.7014E-01  2.4438E-01  9.1677E-02  1.9941E-01 -1.4768E-01 -4.9568E-01 -3.2382E+00 -1.2207E-01 -4.3299E+00
             1.5025E+00
 GRADIENT:  -1.0345E-01  1.0861E-02  3.0258E-02 -2.1459E-02 -3.4616E-02 -1.6479E-02  3.0193E-03  2.9188E-03  1.3560E-02  1.2428E-03
            -7.5294E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1295.33635298190        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  1.0164E+00  1.1733E+00  1.1564E+00  9.9941E-01  1.0995E+00  7.8067E-01  5.5410E-01  1.0000E-02  7.9683E-01  1.0000E-02
             4.0654E+00
 PARAMETER:  1.1623E-01  2.5982E-01  2.4528E-01  9.9409E-02  1.9481E-01 -1.4761E-01 -4.9041E-01 -4.5129E+00 -1.2711E-01 -4.5749E+00
             1.5025E+00
 GRADIENT:   1.6243E-02 -2.9363E-02  3.2858E-03 -3.8679E-02 -3.8461E-03  4.6364E-03 -1.3903E-03  0.0000E+00 -3.4811E-03  0.0000E+00
             4.3700E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1295.33635856681        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1245
 NPARAMETR:  1.0164E+00  1.1763E+00  1.1553E+00  9.9754E-01  1.1004E+00  7.8066E-01  5.5346E-01  1.0000E-02  7.9787E-01  1.0000E-02
             4.0655E+00
 PARAMETER:  1.1624E-01  2.6236E-01  2.4432E-01  9.7532E-02  1.9568E-01 -1.4761E-01 -4.9157E-01 -4.5396E+00 -1.2581E-01 -4.5994E+00
             1.5025E+00
 GRADIENT:   3.5714E-04 -1.3092E-03 -4.0054E-04 -1.0055E-03  1.0391E-03 -7.6620E-05  8.5918E-05  0.0000E+00  1.2437E-04  0.0000E+00
             2.6046E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1245
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8520E-03 -1.0201E-02  5.5757E-05 -9.9847E-03 -2.7536E-05
 SE:             2.7894E-02  1.1125E-02  5.3036E-05  1.7745E-02  1.3269E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4706E-01  3.5917E-01  2.9312E-01  5.7365E-01  8.3560E-01

 ETASHRINKSD(%)  6.5517E+00  6.2730E+01  9.9822E+01  4.0553E+01  9.9555E+01
 ETASHRINKVR(%)  1.2674E+01  8.6109E+01  1.0000E+02  6.4661E+01  9.9998E+01
 EBVSHRINKSD(%)  6.4933E+00  6.3018E+01  9.9791E+01  4.0540E+01  9.9508E+01
 EBVSHRINKVR(%)  1.2565E+01  8.6324E+01  1.0000E+02  6.4646E+01  9.9998E+01
 RELATIVEINF(%)  7.8882E+01  6.4963E-02  1.8598E-05  2.3255E-01  6.8629E-05
 EPSSHRINKSD(%)  1.6790E+01
 EPSSHRINKVR(%)  3.0760E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1295.3363585668130     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -560.18553200307485     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.93
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1295.336       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.18E+00  1.16E+00  9.98E-01  1.10E+00  7.81E-01  5.53E-01  1.00E-02  7.98E-01  1.00E-02  4.07E+00
 


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
+        1.55E+03
 
 TH 2
+       -1.42E+02  3.48E+02
 
 TH 3
+        1.00E+01  5.02E+01  3.24E+01
 
 TH 4
+       -2.12E+02  4.28E+02  2.25E+01  5.99E+02
 
 TH 5
+        1.24E+01 -1.69E+02 -7.87E+01 -1.15E+02  2.14E+02
 
 TH 6
+        3.23E+00 -2.39E+01  2.92E+00 -3.28E+01 -1.33E+00  2.47E+02
 
 TH 7
+        1.83E+00 -1.44E+01 -1.62E+00 -1.23E+01  1.23E+01  3.18E+00  1.57E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.25E-01 -1.23E+01  2.68E+00 -1.48E+00  9.88E+00  5.82E+00  1.64E+01  0.00E+00  3.50E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.25E+01 -1.30E+01 -8.33E-01 -1.43E+01  3.50E+00  4.92E+00  5.78E+00  0.00E+00  9.16E+00  0.00E+00  2.90E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.548
Stop Time:
Sat Sep 18 10:33:41 CDT 2021
