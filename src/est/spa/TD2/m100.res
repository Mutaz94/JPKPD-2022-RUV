Wed Sep 29 19:33:33 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat100.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1656.37730619524        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0978E+02 -1.8357E+01 -2.7489E+01  4.5620E+01  3.8831E+01  3.3646E+01  2.0933E-01  5.0439E+00  2.9050E+01 -6.3787E+00
             9.7717E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1661.89355696640        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.6023E-01  1.1235E+00  1.0833E+00  9.5280E-01  1.0651E+00  1.0747E+00  1.0168E+00  9.7260E-01  8.5058E-01  1.0494E+00
             9.5454E-01
 PARAMETER:  5.9419E-02  2.1640E-01  1.7998E-01  5.1652E-02  1.6309E-01  1.7205E-01  1.1662E-01  7.2216E-02 -6.1835E-02  1.4819E-01
             5.3476E-02
 GRADIENT:   1.1887E-01  1.5988E+00  9.1927E+00 -2.5993E+00 -4.9457E+00  1.6846E+01 -6.4169E+00 -3.7706E+00 -9.5948E+00 -1.1845E+01
            -1.4855E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1663.04905704592        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.6019E-01  1.2566E+00  1.2194E+00  8.8435E-01  1.1847E+00  1.0509E+00  9.6522E-01  1.1587E+00  9.2319E-01  1.2194E+00
             9.6624E-01
 PARAMETER:  5.9373E-02  3.2840E-01  2.9836E-01 -2.2904E-02  2.6947E-01  1.4969E-01  6.4602E-02  2.4732E-01  2.0079E-02  2.9840E-01
             6.5659E-02
 GRADIENT:   1.6353E-01  2.1206E+01  1.2300E+01  9.2345E+00 -7.4093E+00  8.3563E+00  6.0999E-01 -5.8349E+00 -5.8001E+00 -1.8617E+00
            -8.6278E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1663.87629817017        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.6065E-01  1.2813E+00  1.2206E+00  8.6437E-01  1.1973E+00  1.0288E+00  9.0604E-01  1.3512E+00  9.8913E-01  1.2166E+00
             9.8205E-01
 PARAMETER:  5.9858E-02  3.4790E-01  2.9930E-01 -4.5759E-02  2.8010E-01  1.2840E-01  1.3282E-03  4.0100E-01  8.9070E-02  2.9606E-01
             8.1882E-02
 GRADIENT:   4.8974E-01  8.4632E+00  2.7865E+00  7.9186E+00 -4.4631E+00 -8.9530E-02  4.7917E-01 -3.5154E-01 -3.9286E-01 -6.0839E-01
             1.2170E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.45890050793        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  9.6268E-01  1.6050E+00  7.4515E-01  6.4587E-01  1.1877E+00  1.0302E+00  8.0570E-01  1.0065E+00  1.1590E+00  1.1596E+00
             9.8107E-01
 PARAMETER:  6.1961E-02  5.7310E-01 -1.9416E-01 -3.3715E-01  2.7204E-01  1.2971E-01 -1.1604E-01  1.0648E-01  2.4755E-01  2.4811E-01
             8.0890E-02
 GRADIENT:  -2.9182E-01  7.4691E+00 -2.2307E+00  1.1788E+01  2.3895E+00 -5.3926E-01 -6.7174E-01  2.7765E-02 -4.5147E-02 -9.3676E-01
             4.8528E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1664.50859826015        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.6315E-01  1.7537E+00  6.2520E-01  5.5390E-01  1.2202E+00  1.0329E+00  7.6638E-01  9.4573E-01  1.2604E+00  1.1780E+00
             9.7872E-01
 PARAMETER:  6.2450E-02  6.6175E-01 -3.6968E-01 -4.9078E-01  2.9901E-01  1.3239E-01 -1.6608E-01  4.4198E-02  3.3143E-01  2.6380E-01
             7.8489E-02
 GRADIENT:  -3.3287E-01  2.2693E+01 -6.2069E-01  1.5651E+01  4.3937E-01  2.6164E-01 -1.3512E+00 -3.0907E-01 -1.4603E+00 -7.2323E-01
            -1.5317E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1664.56923924562        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.6336E-01  1.8716E+00  5.4213E-01  4.7629E-01  1.2594E+00  1.0347E+00  7.3832E-01  9.2202E-01  1.3764E+00  1.2082E+00
             9.7811E-01
 PARAMETER:  6.2671E-02  7.2682E-01 -5.1224E-01 -6.4174E-01  3.3067E-01  1.3412E-01 -2.0338E-01  1.8813E-02  4.1944E-01  2.8915E-01
             7.7864E-02
 GRADIENT:  -2.3705E-01  2.5693E+01  8.9951E-01  1.2743E+01 -1.0988E+00  8.4555E-01 -1.4032E+00 -4.7731E-01 -2.0325E+00 -2.1129E-01
            -2.2320E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1664.57432574320        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  9.6347E-01  1.9509E+00  4.8303E-01  4.2280E-01  1.2868E+00  1.0354E+00  7.2261E-01  8.9983E-01  1.4743E+00  1.2287E+00
             9.7856E-01
 PARAMETER:  6.2790E-02  7.6832E-01 -6.2767E-01 -7.6085E-01  3.5213E-01  1.3476E-01 -2.2489E-01 -5.5484E-03  4.8818E-01  3.0599E-01
             7.8328E-02
 GRADIENT:  -1.5702E-01  2.4690E+01  1.6654E+00  9.9906E+00 -2.1954E+00  1.0558E+00 -1.0471E+00 -4.7560E-01 -2.0569E+00  1.1272E-01
            -2.3345E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1664.83937872894        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1397
 NPARAMETR:  9.6024E-01  1.9500E+00  4.6715E-01  4.0711E-01  1.2922E+00  1.0314E+00  7.1955E-01  9.3293E-01  1.5022E+00  1.2316E+00
             9.8039E-01
 PARAMETER:  5.9428E-02  7.6784E-01 -6.6111E-01 -7.9868E-01  3.5635E-01  1.3090E-01 -2.2912E-01  3.0577E-02  5.0691E-01  3.0835E-01
             8.0193E-02
 GRADIENT:  -6.9081E+00 -1.2992E+01  1.9044E-01  1.9363E-01  6.4528E-01 -4.2519E-01 -1.7246E-01 -1.1645E-01 -1.2367E+00  9.0388E-01
            -7.8777E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1664.86317506182        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1560
 NPARAMETR:  9.6282E-01  1.9464E+00  4.6765E-01  4.0527E-01  1.2923E+00  1.0327E+00  7.1756E-01  9.6330E-01  1.5240E+00  1.2289E+00
             9.8202E-01
 PARAMETER:  6.2110E-02  7.6599E-01 -6.6003E-01 -8.0319E-01  3.5642E-01  1.3221E-01 -2.3190E-01  6.2605E-02  5.2131E-01  3.0611E-01
             8.1858E-02
 GRADIENT:   4.4769E+02  1.0132E+03  1.9283E+00  1.0605E+02  2.1511E+01  7.8962E+01  1.6431E+01  1.0591E-01  1.3951E+01  3.9518E+00
             8.8616E-01

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1664.86317506182        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1619
 NPARAMETR:  9.6282E-01  1.9464E+00  4.6765E-01  4.0527E-01  1.2923E+00  1.0327E+00  7.1756E-01  9.6330E-01  1.5240E+00  1.2289E+00
             9.8202E-01
 PARAMETER:  6.2110E-02  7.6599E-01 -6.6003E-01 -8.0319E-01  3.5642E-01  1.3221E-01 -2.3190E-01  6.2605E-02  5.2131E-01  3.0611E-01
             8.1858E-02
 GRADIENT:  -1.2449E+00 -2.3526E+01 -4.5340E-01 -1.6192E+00  9.5676E-01  1.5363E-01  7.7374E-02  7.5789E-02  2.0960E-01  9.8162E-01
             8.6245E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1619
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5500E-03 -3.0926E-02 -2.0794E-02  3.5170E-02 -4.4755E-02
 SE:             2.9861E-02  2.4463E-02  6.7583E-03  2.0622E-02  2.2529E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5860E-01  2.0616E-01  2.0922E-03  8.8104E-02  4.6976E-02

 ETASHRINKSD(%)  1.0000E-10  1.8045E+01  7.7359E+01  3.0914E+01  2.4524E+01
 ETASHRINKVR(%)  1.0000E-10  3.2833E+01  9.4874E+01  5.2272E+01  4.3034E+01
 EBVSHRINKSD(%)  3.9107E-01  1.6904E+01  8.1041E+01  3.5619E+01  1.9780E+01
 EBVSHRINKVR(%)  7.8061E-01  3.0951E+01  9.6405E+01  5.8550E+01  3.5648E+01
 RELATIVEINF(%)  9.9159E+01  3.9822E+00  3.2737E-01  2.0401E+00  2.1036E+01
 EPSSHRINKSD(%)  4.4861E+01
 EPSSHRINKVR(%)  6.9597E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1664.8631750618194     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -929.71234849808127     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.78
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1664.863       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.63E-01  1.95E+00  4.68E-01  4.05E-01  1.29E+00  1.03E+00  7.18E-01  9.63E-01  1.52E+00  1.23E+00  9.82E-01
 


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
+        1.12E+03
 
 TH 2
+       -1.87E+06  1.21E+05
 
 TH 3
+        2.54E+02  2.89E+03  2.33E+02
 
 TH 4
+       -2.49E+02  5.54E+05 -5.35E+06  2.55E+06
 
 TH 5
+        1.21E+07 -3.28E+01  3.78E+06 -5.78E+03  1.27E+06
 
 TH 6
+       -6.50E+02 -1.64E+00  2.05E+02 -1.97E+02  1.35E+02  1.84E+02
 
 TH 7
+       -1.49E+03  5.69E+00 -5.22E+06 -4.96E+06  3.00E+02  2.88E-01  1.99E+02
 
 TH 8
+        7.87E+02  1.87E+06 -2.11E+04  1.72E+07 -1.21E+07  6.50E+02  3.17E+00  2.39E+00
 
 TH 9
+       -3.50E+06 -1.21E+01 -8.22E+03  7.83E+03 -5.47E+03  3.56E-02  2.00E+02 -2.61E+04  2.48E+01
 
 TH10
+        7.43E+06 -3.30E+01  4.63E+06 -7.31E+03  1.80E+02  1.67E+02  3.84E+02 -1.72E+04 -6.70E+03  1.90E+06
 
 TH11
+       -7.81E+02 -1.84E+06 -8.85E+06 -8.40E+06  1.19E+07 -6.38E+02  1.65E+07  1.52E+00  3.46E+06  7.28E+06  2.79E+07
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       28.461
Stop Time:
Wed Sep 29 19:34:03 CDT 2021
