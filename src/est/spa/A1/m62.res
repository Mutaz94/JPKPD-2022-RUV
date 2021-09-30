Wed Sep 29 12:15:22 CDT 2021
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
$DATA ../../../../data/spa/A1/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1367.35844100540        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9274E+02  5.2269E+01  2.7656E+01  8.1719E+01  1.0359E+02  7.2528E+01 -2.6229E+01 -1.8174E+01 -2.4916E+01 -5.3285E+01
            -4.9500E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1481.64201149504        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1413E+00  9.6736E-01  8.9146E-01  1.0488E+00  8.7766E-01  1.1364E+00  1.0762E+00  1.0269E+00  1.0769E+00  1.0768E+00
             2.1486E+00
 PARAMETER:  2.3217E-01  6.6818E-02 -1.4898E-02  1.4767E-01 -3.0498E-02  2.2789E-01  1.7345E-01  1.2653E-01  1.7406E-01  1.7398E-01
             8.6481E-01
 GRADIENT:   3.8614E+02  2.1039E+01  2.6406E-01  4.4078E+01 -6.7553E+00  6.9893E+01  5.7478E-01  6.3155E+00  1.2767E+01  1.3539E+01
             7.2261E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1492.85259631837        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      199
 NPARAMETR:  1.1283E+00  7.0509E-01  5.9149E-01  1.1909E+00  5.9028E-01  1.1301E+00  1.6919E+00  7.2568E-01  9.1328E-01  6.4206E-01
             1.9833E+00
 PARAMETER:  2.2074E-01 -2.4943E-01 -4.2511E-01  2.7472E-01 -4.2715E-01  2.2231E-01  6.2584E-01 -2.2065E-01  9.2838E-03 -3.4308E-01
             7.8476E-01
 GRADIENT:   1.3732E+02  4.2742E+01  2.9915E+00  6.9185E+01 -1.4678E+01  4.9334E+01  1.1222E+01  7.2918E+00  8.5437E+00 -3.0481E+00
             4.0469E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1508.02130294054        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0774E+00  6.5667E-01  6.5033E-01  1.1895E+00  6.2489E-01  9.9222E-01  1.8536E+00  4.3618E-01  8.4683E-01  8.0838E-01
             1.7083E+00
 PARAMETER:  1.7459E-01 -3.2058E-01 -3.3027E-01  2.7350E-01 -3.7018E-01  9.2193E-02  7.1713E-01 -7.2970E-01 -6.6255E-02 -1.1272E-01
             6.3551E-01
 GRADIENT:   8.7096E+01  2.7417E+01  1.2746E+01  2.1031E+01 -4.4846E+00  1.7815E+01  1.1067E+01  4.0144E-01 -8.2276E+00 -3.6582E+00
            -2.1266E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1516.09599786245        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0335E+00  4.2290E-01  4.9935E-01  1.2399E+00  4.6861E-01  9.4797E-01  1.9717E+00  2.8016E-01  8.5189E-01  7.2310E-01
             1.7844E+00
 PARAMETER:  1.3292E-01 -7.6062E-01 -5.9444E-01  3.1501E-01 -6.5799E-01  4.6565E-02  7.7891E-01 -1.1724E+00 -6.0294E-02 -2.2421E-01
             6.7910E-01
 GRADIENT:  -8.7117E+00 -6.5025E-01  4.0355E+00 -1.2789E+01 -3.9192E+00  5.4982E+00 -2.3831E+00  6.3111E-01  1.8960E+00 -7.4990E-01
             9.2408E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1516.98521506507        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0322E+00  2.6530E-01  5.4099E-01  1.3462E+00  4.5861E-01  9.2975E-01  2.5726E+00  8.8473E-02  8.2931E-01  8.4819E-01
             1.7420E+00
 PARAMETER:  1.3170E-01 -1.2269E+00 -5.1435E-01  3.9731E-01 -6.7956E-01  2.7164E-02  1.0449E+00 -2.3251E+00 -8.7157E-02 -6.4656E-02
             6.5501E-01
 GRADIENT:  -4.9908E-02  5.0191E+00  8.5554E+00  1.6263E+01 -1.5188E+01  2.2380E-01  1.3433E+00  3.1022E-02 -2.4306E+00 -1.0110E+00
            -1.1982E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1518.55210820541        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  1.0265E+00  1.7133E-01  5.4058E-01  1.3880E+00  4.4565E-01  9.3023E-01  3.0779E+00  2.4064E-02  8.3214E-01  8.7427E-01
             1.7363E+00
 PARAMETER:  1.2617E-01 -1.6641E+00 -5.1511E-01  4.2788E-01 -7.0822E-01  2.7674E-02  1.2243E+00 -3.6270E+00 -8.3757E-02 -3.4362E-02
             6.5177E-01
 GRADIENT:  -4.8245E+00  2.9654E+00  1.0301E+00  1.7626E+01 -5.6588E+00  1.4769E+00 -1.7726E+00  3.3074E-03  7.2701E-01 -2.4626E-01
            -4.5989E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1519.92950974832        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1078
 NPARAMETR:  1.0232E+00  6.8424E-02  5.9002E-01  1.4437E+00  4.6111E-01  9.2614E-01  5.1156E+00  1.0000E-02  7.9257E-01  9.2922E-01
             1.7534E+00
 PARAMETER:  1.2292E-01 -2.5820E+00 -4.2759E-01  4.6722E-01 -6.7412E-01  2.3267E-02  1.7323E+00 -6.8613E+00 -1.3247E-01  2.6591E-02
             6.6155E-01
 GRADIENT:   1.6458E+00  1.3214E+01 -9.8166E+00 -3.2485E+01  1.8822E+01 -7.6260E-02  2.9420E+01  0.0000E+00 -1.1396E+01 -1.0627E+01
            -2.3093E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1520.98403433231        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1255
 NPARAMETR:  1.0216E+00  3.2840E-02  5.5728E-01  1.4408E+00  4.3896E-01  9.2100E-01  6.8379E+00  1.0000E-02  7.9528E-01  8.9637E-01
             1.7399E+00
 PARAMETER:  1.2134E-01 -3.3161E+00 -4.8469E-01  4.6520E-01 -7.2335E-01  1.7702E-02  2.0225E+00 -9.3968E+00 -1.2906E-01 -9.3985E-03
             6.5381E-01
 GRADIENT:  -9.0363E-01 -2.3330E+00  5.0022E+00  3.1748E+00 -6.5329E+00 -2.9371E-01 -5.6007E+00  0.0000E+00  8.5965E-01  3.5771E+00
             1.1363E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1520.99872170828        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1432
 NPARAMETR:  1.0224E+00  3.3777E-02  5.4933E-01  1.4408E+00  4.3537E-01  9.2389E-01  6.7662E+00  1.0000E-02  8.0000E-01  8.8391E-01
             1.7385E+00
 PARAMETER:  1.2214E-01 -3.2880E+00 -4.9906E-01  4.6522E-01 -7.3157E-01  2.0839E-02  2.0119E+00 -9.3006E+00 -1.2314E-01 -2.3398E-02
             6.5303E-01
 GRADIENT:   7.8285E-01 -2.1198E+00  1.8335E-01  1.0520E+01 -1.8404E+00  7.5713E-01 -5.2237E+00  0.0000E+00  2.1750E+00  1.4089E+00
             4.5204E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1521.23310227542        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1596
 NPARAMETR:  1.0216E+00  2.9677E-02  5.4918E-01  1.4393E+00  4.3432E-01  9.2204E-01  7.2479E+00  1.0000E-02  7.9589E-01  8.8435E-01
             1.7386E+00
 PARAMETER:  1.2135E-01 -3.4174E+00 -4.9932E-01  4.6414E-01 -7.3398E-01  1.8829E-02  2.0807E+00 -9.8844E+00 -1.2830E-01 -2.2897E-02
             6.5310E-01
 GRADIENT:  -4.6705E-01 -6.3666E-01  1.5041E+00  1.7355E+00 -2.2174E+00 -5.3756E-02 -1.3086E+00  0.0000E+00  2.3156E-01  2.3930E-01
            -9.1365E-03

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1521.23556299054        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1689
 NPARAMETR:  1.0217E+00  2.9761E-02  5.4894E-01  1.4390E+00  4.3446E-01  9.2220E-01  7.2477E+00  1.0000E-02  7.9629E-01  8.8612E-01
             1.7391E+00
 PARAMETER:  1.2143E-01 -3.4146E+00 -4.9976E-01  4.6394E-01 -7.3364E-01  1.9006E-02  2.0807E+00 -9.8844E+00 -1.2779E-01 -2.0897E-02
             6.5339E-01
 GRADIENT:  -1.3501E-01  3.1538E-01 -3.7397E-01 -6.5065E-01  5.0333E-01 -8.2233E-02  7.6070E-01  0.0000E+00 -1.9260E-01 -2.9068E-01
            -4.9916E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1689
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3510E-04  1.4072E-02 -1.2412E-04 -1.1109E-02 -7.7102E-03
 SE:             2.9519E-02  8.3214E-03  2.2318E-04  2.7796E-02  2.4773E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8824E-01  9.0820E-02  5.7810E-01  6.8940E-01  7.5562E-01

 ETASHRINKSD(%)  1.1077E+00  7.2122E+01  9.9252E+01  6.8783E+00  1.7008E+01
 ETASHRINKVR(%)  2.2031E+00  9.2228E+01  9.9994E+01  1.3283E+01  3.1122E+01
 EBVSHRINKSD(%)  1.3464E+00  8.1931E+01  9.9258E+01  6.0122E+00  1.4875E+01
 EBVSHRINKVR(%)  2.6747E+00  9.6735E+01  9.9994E+01  1.1663E+01  2.7537E+01
 RELATIVEINF(%)  9.6590E+01  2.0385E+00  2.9962E-04  3.3520E+01  4.0887E+00
 EPSSHRINKSD(%)  3.9449E+01
 EPSSHRINKVR(%)  6.3336E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1521.2355629905362     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -786.08473642679803     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1521.236       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.98E-02  5.49E-01  1.44E+00  4.34E-01  9.22E-01  7.25E+00  1.00E-02  7.96E-01  8.86E-01  1.74E+00
 


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
+        1.23E+03
 
 TH 2
+       -1.27E+02  1.95E+05
 
 TH 3
+        3.84E+01  7.15E+02  2.94E+04
 
 TH 4
+       -1.62E+01  8.24E+02 -3.98E+02  5.30E+03
 
 TH 5
+        1.24E+01 -1.96E+03 -2.88E+03  1.76E+02  2.53E+04
 
 TH 6
+       -8.58E+00  3.74E+01 -2.40E+01 -1.47E+01  2.08E+01  2.28E+02
 
 TH 7
+       -3.26E-01 -2.37E+01  3.26E+00  3.52E+00 -9.12E+00  3.00E-01  9.21E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.27E+01  1.27E+02 -1.03E+02 -5.86E+01  8.95E+01  5.68E+01  1.17E+00  0.00E+00  2.01E+05
 
 TH10
+       -1.24E+02  2.28E+02 -2.28E+02 -5.75E+01  6.74E+01  1.02E+02  1.29E+00  0.00E+00  2.30E+05  2.64E+05
 
 TH11
+       -2.88E+00  5.86E+01 -7.19E+01 -2.51E+01  4.67E+01 -4.83E+00  4.28E-01  0.00E+00 -2.42E+01 -3.03E+01  1.68E+03
 
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
 #CPUT: Total CPU Time in Seconds,       27.573
Stop Time:
Wed Sep 29 12:15:51 CDT 2021
