Fri Sep 24 22:38:58 CDT 2021
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
$DATA ../../../../data/int/A3/dat72.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 (2E4.0,E20.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -950.888298357345        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0046E+01  2.3499E+02  3.7846E+02 -1.2432E+02  1.6561E+02  3.0780E+01 -2.4598E+02 -4.1226E+02 -4.5735E+01 -1.6046E+02
            -4.9612E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2797.81243591701        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0431E+00  9.1473E-01  8.3217E-01  1.0935E+00  8.9899E-01  8.5900E-01  1.0501E+00  9.7334E-01  8.6464E-01  1.0029E+00
             2.4918E+00
 PARAMETER:  1.4215E-01  1.0873E-02 -8.3719E-02  1.8938E-01 -6.4861E-03 -5.1983E-02  1.4885E-01  7.2983E-02 -4.5444E-02  1.0289E-01
             1.0130E+00
 GRADIENT:   4.5921E+01  1.8822E+01  7.0404E+00 -1.7477E+01  1.0258E+01 -2.3787E+01 -8.1229E+00  1.1713E+01 -7.5091E+00  2.5556E+00
            -7.2521E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2812.30700047388        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0491E+00  6.1242E-01  5.7216E-01  1.3095E+00  5.7625E-01  8.8116E-01  1.3350E+00  2.5637E-01  7.3842E-01  6.4594E-01
             2.4595E+00
 PARAMETER:  1.4790E-01 -3.9034E-01 -4.5833E-01  3.6968E-01 -4.5121E-01 -2.6513E-02  3.8891E-01 -1.2612E+00 -2.0324E-01 -3.3704E-01
             9.9996E-01
 GRADIENT:   5.3908E+01  7.4519E+01  5.5135E+01  2.3257E+02 -3.2949E+00 -1.7083E+01 -1.2179E+01 -1.1019E+00 -7.9768E+01 -3.1758E+01
            -5.9462E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2888.44142623336        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0011E+00  2.8483E-01  2.0236E-01  1.3510E+00  2.4895E-01  1.0974E+00  1.4512E+00  1.0000E-02  1.3454E+00  7.5460E-01
             2.2473E+00
 PARAMETER:  1.0107E-01 -1.1559E+00 -1.4977E+00  4.0084E-01 -1.2905E+00  1.9297E-01  4.7237E-01 -6.1627E+00  3.9670E-01 -1.8157E-01
             9.0972E-01
 GRADIENT:  -4.3797E+01  7.5270E+01 -1.8115E+01  2.2645E+02  1.5765E+01  5.1902E+01  1.1831E+00  0.0000E+00 -7.8019E+00  2.5970E+01
             3.7574E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2923.80991814973        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0222E+00  1.4730E-01  7.8865E-02  8.3629E-01  1.3685E-01  9.0423E-01  1.6542E+00  1.0000E-02  1.8633E+00  8.2312E-01
             2.0252E+00
 PARAMETER:  1.2196E-01 -1.8153E+00 -2.4400E+00 -7.8782E-02 -1.8889E+00 -6.6858E-04  6.0331E-01 -1.5781E+01  7.2236E-01 -9.4648E-02
             8.0567E-01
 GRADIENT:   1.0265E+01 -8.9863E+00 -3.7241E+00  2.6846E+01 -5.5836E+01 -9.4134E+00  1.0764E+00  0.0000E+00 -4.5114E+01 -1.6205E+01
            -3.0582E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2931.28537476976        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  1.0238E+00  1.5474E-01  8.2974E-02  8.3985E-01  1.4337E-01  9.0802E-01  1.6683E+00  1.0000E-02  1.8811E+00  7.9126E-01
             2.0370E+00
 PARAMETER:  1.2352E-01 -1.7660E+00 -2.3892E+00 -7.4535E-02 -1.8423E+00  3.5157E-03  6.1179E-01 -1.5460E+01  7.3184E-01 -1.3412E-01
             8.1148E-01
 GRADIENT:   6.5592E+00 -3.9551E+00 -2.0450E+01  9.8935E+00 -1.0220E+02 -8.1785E+00  5.7094E+00  0.0000E+00 -3.7478E+01 -1.9803E+01
            -1.2298E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2934.05658656489        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      597
 NPARAMETR:  1.0194E+00  1.5596E-01  8.3533E-02  8.2540E-01  1.4818E-01  9.3474E-01  1.6222E+00  1.0000E-02  1.8907E+00  8.4752E-01
             2.0368E+00
 PARAMETER:  1.1918E-01 -1.7582E+00 -2.3825E+00 -9.1881E-02 -1.8094E+00  3.2517E-02  5.8377E-01 -1.5415E+01  7.3694E-01 -6.5436E-02
             8.1137E-01
 GRADIENT:   3.7982E+00  5.2138E+00 -3.4769E+01 -2.7505E+00  9.5072E+01  3.4069E+00 -1.3391E+00  0.0000E+00 -3.0494E+01  6.4917E+00
            -6.3120E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2934.21908233633        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      775
 NPARAMETR:  1.0210E+00  1.5747E-01  8.3574E-02  8.3423E-01  1.4861E-01  9.2856E-01  1.6319E+00  1.0000E-02  1.8908E+00  8.3164E-01
             2.0368E+00
 PARAMETER:  1.2074E-01 -1.7485E+00 -2.3820E+00 -8.1251E-02 -1.8064E+00  2.5885E-02  5.8975E-01 -1.5415E+01  7.3700E-01 -8.4354E-02
             8.1138E-01
 GRADIENT:  -4.4975E-01  7.0758E-02 -5.6431E+01  1.3430E+00 -4.5551E-01  3.1351E-01 -4.2447E-01  0.0000E+00 -3.3003E+01  4.1312E-01
            -8.3097E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2934.30414332420        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:      943
 NPARAMETR:  1.0209E+00  1.5712E-01  8.3699E-02  8.3395E-01  1.4876E-01  9.3065E-01  1.6264E+00  1.0000E-02  1.8913E+00  8.2914E-01
             2.0367E+00
 PARAMETER:  1.2072E-01 -1.7508E+00 -2.3805E+00 -8.1580E-02 -1.8054E+00  2.8129E-02  5.8637E-01 -1.5415E+01  7.3727E-01 -8.7362E-02
             8.1133E-01
 GRADIENT:  -4.2169E-01 -2.1006E+00 -5.5628E+01  7.4037E-01  1.1346E+00  1.1623E+00 -1.3399E+00  0.0000E+00 -3.2870E+01 -3.8271E-01
            -8.4084E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2934.35952107133        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1132             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0212E+00  1.5766E-01  8.3754E-02  8.3379E-01  1.4898E-01  9.2836E-01  1.6327E+00  1.0000E-02  1.8922E+00  8.3025E-01
             2.0361E+00
 PARAMETER:  1.2094E-01 -1.7473E+00 -2.3799E+00 -8.1769E-02 -1.8039E+00  2.5669E-02  5.9026E-01 -1.5415E+01  7.3775E-01 -8.6029E-02
             8.1104E-01
 GRADIENT:   8.4976E+00  1.1240E+01 -4.4141E+01  1.0745E+00  1.0413E+02  7.8153E-01  1.0097E+00  0.0000E+00 -2.9974E+01  3.9265E-01
            -7.2603E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2934.38165364844        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  1.0212E+00  1.5767E-01  8.3787E-02  8.3190E-01  1.4899E-01  9.2539E-01  1.6388E+00  1.0000E-02  1.8924E+00  8.3007E-01
             2.0361E+00
 PARAMETER:  1.2096E-01 -1.7473E+00 -2.3795E+00 -8.4048E-02 -1.8039E+00  2.2364E-02  5.9398E-01 -1.5415E+01  7.3782E-01 -8.6249E-02
             8.1103E-01
 GRADIENT:   1.7271E+04 -1.2011E+03  8.5517E+02 -1.0467E+04 -1.1625E+03 -9.6403E-01 -3.5294E+03  0.0000E+00 -1.4557E+03 -1.0442E+04
             2.5814E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1315
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5811E-03  1.1627E-02  6.9447E-04 -3.0281E-03  7.2412E-03
 SE:             2.9392E-02  2.6693E-02  2.8688E-04  2.9153E-02  2.8318E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5710E-01  6.6313E-01  1.5487E-02  9.1727E-01  7.9817E-01

 ETASHRINKSD(%)  1.5343E+00  1.0576E+01  9.9039E+01  2.3329E+00  5.1325E+00
 ETASHRINKVR(%)  3.0451E+00  2.0034E+01  9.9991E+01  4.6115E+00  1.0001E+01
 EBVSHRINKSD(%)  1.6371E+00  9.3721E+00  9.9245E+01  7.9062E+00  6.2401E+00
 EBVSHRINKVR(%)  3.2473E+00  1.7866E+01  9.9994E+01  1.5187E+01  1.2091E+01
 RELATIVEINF(%)  9.6691E+01  5.7153E+01  1.8230E-03  5.8938E+01  3.2784E+01
 EPSSHRINKSD(%)  2.0598E+01
 EPSSHRINKVR(%)  3.6954E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2934.3816536484392     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1280.2922938800284     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    41.35
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2934.382       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.58E-01  8.38E-02  8.32E-01  1.49E-01  9.25E-01  1.64E+00  1.00E-02  1.89E+00  8.30E-01  2.04E+00
 


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
+        3.43E+07
 
 TH 2
+       -2.84E+03  6.94E+06
 
 TH 3
+       -4.55E+03 -1.35E+04  1.39E+07
 
 TH 4
+       -4.34E+04  2.38E+04  2.53E+03  7.57E+07
 
 TH 5
+       -9.88E+02  2.42E+04 -6.87E+04  1.56E+04  7.37E+06
 
 TH 6
+        4.57E+07 -2.37E+03 -8.94E+02 -2.36E+04 -1.46E+03  2.05E+02
 
 TH 7
+       -1.39E+03  6.99E+03 -8.40E+02  9.46E+03  4.82E+03 -9.51E+02  5.54E+05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.65E+01  4.05E+03  7.20E+03  2.13E+03  9.38E+03 -1.79E+02  7.27E+02  0.00E+00  2.70E+05
 
 TH10
+       -5.10E+07  1.08E+04 -2.41E+03  8.26E+04  8.43E+03 -6.81E+07  3.86E+03  0.00E+00  1.31E+03  7.59E+07
 
 TH11
+       -1.44E+02 -2.90E+03 -8.60E+03 -1.38E+03 -6.63E+03  1.03E+02 -5.12E+02  0.00E+00 -5.15E+02 -9.61E+02  1.93E+05
 
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
 #CPUT: Total CPU Time in Seconds,       59.375
Stop Time:
Fri Sep 24 22:39:59 CDT 2021
