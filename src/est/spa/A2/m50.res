Sat Sep 18 09:54:15 CDT 2021
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
$DATA ../../../../data/spa/A2/dat50.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1047.26345071865        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.7567E+01  9.9614E+00  2.8650E+01 -5.1619E+01  7.9226E+01  9.1303E+00 -1.3956E+01  9.1631E-01 -6.0933E+01 -4.6224E+01
            -1.0982E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1420.36496169225        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9946E-01  8.7722E-01  1.0832E+00  1.1211E+00  9.4234E-01  8.9075E-01  8.6208E-01  8.4661E-01  1.1534E+00  7.4930E-01
             2.5217E+00
 PARAMETER:  9.9462E-02 -3.0994E-02  1.7990E-01  2.1430E-01  4.0613E-02 -1.5686E-02 -4.8408E-02 -6.6510E-02  2.4275E-01 -1.8861E-01
             1.0249E+00
 GRADIENT:   1.2289E+01  7.8446E-02  2.9851E+00  1.1233E+00  1.1491E+01 -2.7931E+01  8.2100E+00  5.9659E+00  1.3242E+01  4.4719E+00
             8.2360E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1428.91558072044        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0084E+00  7.9433E-01  7.2877E-01  1.1646E+00  7.2433E-01  9.9181E-01  4.9089E-01  4.2954E-01  1.1426E+00  4.7539E-01
             2.5456E+00
 PARAMETER:  1.0841E-01 -1.3026E-01 -2.1640E-01  2.5235E-01 -2.2251E-01  9.1777E-02 -6.1153E-01 -7.4503E-01  2.3327E-01 -6.4362E-01
             1.0344E+00
 GRADIENT:   2.2615E+01  8.6872E+00 -1.3133E+01  3.0120E+01  1.6654E+01  9.6404E+00  3.8160E-02  1.5257E+00  6.2141E+00 -7.1338E-01
             6.6347E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1433.29676204168        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.8712E-01  4.8734E-01  5.7629E-01  1.3007E+00  5.2197E-01  9.6270E-01  7.2275E-01  6.6363E-02  9.4805E-01  5.7349E-01
             2.3268E+00
 PARAMETER:  8.7039E-02 -6.1880E-01 -4.5115E-01  3.6293E-01 -5.5015E-01  6.1991E-02 -2.2469E-01 -2.6126E+00  4.6648E-02 -4.5602E-01
             9.4448E-01
 GRADIENT:  -2.3313E+01  1.7830E+01 -4.1534E-02  3.8811E+01 -3.2456E+00 -3.3957E+00 -1.5898E+00  4.3729E-02 -1.1706E+01 -2.1966E+00
            -1.2553E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1435.04073858329        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0040E+00  3.5269E-01  3.6079E-01  1.2731E+00  3.5162E-01  1.0091E+00  1.4838E+00  1.0000E-02  9.2607E-01  5.6248E-01
             2.2197E+00
 PARAMETER:  1.0397E-01 -9.4217E-01 -9.1947E-01  3.4144E-01 -9.4520E-01  1.0902E-01  4.9459E-01 -4.7343E+00  2.3199E-02 -4.7540E-01
             8.9739E-01
 GRADIENT:   7.3378E+00  1.0682E+01 -2.8684E-01  4.3895E+01 -1.1601E+01  1.0767E+01 -1.4207E+00  0.0000E+00 -5.7718E+00  7.4585E+00
             3.9837E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1435.45432596589        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0041E+00  3.2211E-01  3.2557E-01  1.2387E+00  3.2483E-01  9.9453E-01  1.7347E+00  1.0000E-02  9.3762E-01  5.0947E-01
             2.2150E+00
 PARAMETER:  1.0406E-01 -1.0329E+00 -1.0222E+00  3.1409E-01 -1.0244E+00  9.4515E-02  6.5081E-01 -5.3235E+00  3.5590E-02 -5.7439E-01
             8.9525E-01
 GRADIENT:   8.5120E+00  4.4327E+00  2.4292E+00  1.2189E+01 -9.9244E+00  5.3615E+00  3.8382E-01  0.0000E+00 -5.9655E-01  3.0593E+00
             3.5123E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1435.82663021398        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  9.9703E-01  3.0664E-01  3.3766E-01  1.2420E+00  3.3244E-01  9.7189E-01  1.8812E+00  1.0000E-02  9.2695E-01  4.7809E-01
             2.2082E+00
 PARAMETER:  9.7025E-02 -1.0821E+00 -9.8572E-01  3.1673E-01 -1.0013E+00  7.1486E-02  7.3189E-01 -5.4091E+00  2.4145E-02 -6.3796E-01
             8.9220E-01
 GRADIENT:  -3.8603E+00 -1.1851E+00 -2.8384E-01 -5.4729E+00  3.0194E+00 -2.2921E+00 -1.4690E-01  0.0000E+00  3.3307E-03 -1.4887E+00
            -1.9817E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1435.83298190506        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  9.9813E-01  3.0810E-01  3.3590E-01  1.2420E+00  3.3117E-01  9.7530E-01  1.8625E+00  1.0000E-02  9.2854E-01  4.8387E-01
             2.2094E+00
 PARAMETER:  9.8133E-02 -1.0773E+00 -9.9093E-01  3.1669E-01 -1.0051E+00  7.4992E-02  7.2192E-01 -5.4081E+00  2.5860E-02 -6.2595E-01
             8.9270E-01
 GRADIENT:  -1.8443E+00 -4.0573E-01  2.5795E-02 -2.7683E+00  1.1716E+00 -1.0891E+00 -6.4216E-02  0.0000E+00 -3.8098E-02 -7.3411E-01
            -9.8285E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1436.86044329081        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      676
 NPARAMETR:  1.0005E+00  3.2513E-01  4.2208E-01  1.2988E+00  3.8957E-01  9.7277E-01  1.8865E+00  1.0000E-02  9.0303E-01  5.1275E-01
             2.2629E+00
 PARAMETER:  1.0046E-01 -1.0235E+00 -7.6256E-01  3.6144E-01 -8.4271E-01  7.2391E-02  7.3474E-01 -4.8106E+00 -2.0030E-03 -5.6796E-01
             9.1664E-01
 GRADIENT:   4.4192E-01  2.1812E+00  3.5022E+00 -3.7024E-01 -5.3333E+00 -6.1787E-01 -2.4152E-02  0.0000E+00 -9.2900E-01 -3.2116E-01
            -2.5287E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1436.96700684485        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      852
 NPARAMETR:  9.9860E-01  2.8249E-01  4.3595E-01  1.3247E+00  3.9153E-01  9.7256E-01  2.0763E+00  1.0000E-02  8.9752E-01  5.3609E-01
             2.2637E+00
 PARAMETER:  9.8599E-02 -1.1641E+00 -7.3023E-01  3.8115E-01 -8.3770E-01  7.2180E-02  8.3059E-01 -5.3057E+00 -8.1231E-03 -5.2345E-01
             9.1702E-01
 GRADIENT:   2.5497E-02  1.6110E-02 -2.3926E-02 -2.1306E-02 -1.8065E-02 -3.4621E-02  2.3225E-02  0.0000E+00 -2.1611E-02  1.1471E-02
             4.0050E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1436.96703406597        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  9.9856E-01  2.8227E-01  4.3644E-01  1.3250E+00  3.9179E-01  9.7257E-01  2.0758E+00  1.0000E-02  8.9760E-01  5.3643E-01
             2.2639E+00
 PARAMETER:  9.8561E-02 -1.1649E+00 -7.2911E-01  3.8144E-01 -8.3702E-01  7.2191E-02  8.3034E-01 -5.3065E+00 -8.0322E-03 -5.2282E-01
             9.1708E-01
 GRADIENT:  -1.5685E-02  5.3525E-03 -9.7037E-03  2.0852E-02 -1.0448E-03 -2.1326E-02  1.0215E-02  0.0000E+00  5.0701E-03  7.6816E-04
             1.7452E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      909
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.7709E-04  2.7254E-02 -2.2041E-04 -1.7245E-02  5.5101E-03
 SE:             2.9232E-02  1.5441E-02  2.3363E-04  2.6287E-02  1.8057E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8152E-01  7.7546E-02  3.4546E-01  5.1180E-01  7.6025E-01

 ETASHRINKSD(%)  2.0707E+00  4.8272E+01  9.9217E+01  1.1934E+01  3.9508E+01
 ETASHRINKVR(%)  4.0984E+00  7.3242E+01  9.9994E+01  2.2444E+01  6.3407E+01
 EBVSHRINKSD(%)  2.0608E+00  5.7929E+01  9.9125E+01  1.0443E+01  3.6198E+01
 EBVSHRINKVR(%)  4.0791E+00  8.2300E+01  9.9992E+01  1.9795E+01  5.9293E+01
 RELATIVEINF(%)  9.4481E+01  3.9082E+00  2.6381E-04  2.7660E+01  1.2858E+00
 EPSSHRINKSD(%)  3.4615E+01
 EPSSHRINKVR(%)  5.7248E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1436.9670340659691     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -701.81620750223090     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.17
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1436.967       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  2.82E-01  4.36E-01  1.33E+00  3.92E-01  9.73E-01  2.08E+00  1.00E-02  8.98E-01  5.36E-01  2.26E+00
 


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
+       -5.94E+01  6.15E+02
 
 TH 3
+        1.06E+01  1.01E+03  4.49E+03
 
 TH 4
+       -3.08E+01  2.61E+02 -3.68E+02  6.34E+02
 
 TH 5
+        7.92E+01 -1.82E+03 -6.54E+03  4.40E+01  1.02E+04
 
 TH 6
+       -1.19E+00 -7.44E+00  1.03E+01 -1.01E+01  5.37E+00  1.93E+02
 
 TH 7
+        9.63E-01  2.56E+01  2.24E+01  4.73E-01 -3.81E+01  4.19E-01  3.79E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.43E+00 -5.57E+00 -6.47E+00 -1.44E+01  6.15E+01  3.40E+00  7.89E+00  0.00E+00  1.61E+02
 
 TH10
+       -5.78E+00  4.52E+01 -1.61E+02 -1.41E+01  1.62E+02  4.15E-01  7.40E+00  0.00E+00 -8.80E+00  1.23E+02
 
 TH11
+       -1.36E+01  6.43E+00 -5.13E+01 -9.61E+00  3.20E+01  2.37E+00  2.74E+00  0.00E+00  6.81E+00  2.66E+01  5.37E+01
 
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
 #CPUT: Total CPU Time in Seconds,       15.865
Stop Time:
Sat Sep 18 09:54:33 CDT 2021
