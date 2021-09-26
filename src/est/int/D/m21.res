Sat Sep 25 05:29:58 CDT 2021
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
$DATA ../../../../data/int/D/dat21.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1063.29633414775        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.5847E+02  4.6577E+01 -5.3152E+01 -2.9358E+02  2.1797E+02 -1.0494E+03 -3.0056E+02 -7.9863E+01 -8.8686E+02 -4.7347E+02
            -7.5993E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2382.91567620735        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.4643E+00  1.2049E+00  1.1742E+00  1.5294E+00  8.1966E-01  3.8122E+00  1.7649E+00  9.4744E-01  4.5433E+00  2.9683E+00
             3.2777E+00
 PARAMETER:  4.8138E-01  2.8643E-01  2.6055E-01  5.2487E-01 -9.8867E-02  1.4382E+00  6.6811E-01  4.6003E-02  1.6137E+00  1.1880E+00
             1.2871E+00
 GRADIENT:   4.8553E+01  8.8981E+00 -1.7824E+01  4.2387E+01 -6.1558E+01  1.3941E+02  4.5301E+01  8.5404E+00  1.0189E+02  5.6726E+01
             2.9540E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2471.95565926885        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0789E+00  3.5583E+00  1.1540E+01  8.2559E-01  2.3721E+00  2.5581E+00  1.8349E+00  2.1464E-01  1.8439E+01  2.3096E+00
             2.8994E+00
 PARAMETER:  1.7590E-01  1.3693E+00  2.5458E+00 -9.1663E-02  9.6380E-01  1.0393E+00  7.0698E-01 -1.4388E+00  3.0144E+00  9.3708E-01
             1.1645E+00
 GRADIENT:  -2.4126E+01  1.5349E+02  3.7057E+00  1.9358E+01 -1.9108E+00  7.2807E+01 -3.7802E+01 -2.4981E-01  6.5981E+01  8.1836E+01
            -2.7152E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2558.35755497604        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.1981E+00  1.3681E+00  8.7253E+00  6.1764E-01  1.9047E+00  2.7906E+00  1.3545E+00  1.4532E-01  4.5584E+00  1.6178E+00
             2.9213E+00
 PARAMETER:  2.8075E-01  4.1342E-01  2.2662E+00 -3.8185E-01  7.4434E-01  1.1262E+00  4.0345E-01 -1.8288E+00  1.6170E+00  5.8107E-01
             1.1720E+00
 GRADIENT:   1.0333E+01 -6.0836E+01 -5.3029E+01 -3.1732E+01 -9.1825E+00  7.5874E+01  2.0811E+01 -1.7782E-01  4.4530E+01 -2.3025E+01
             6.8974E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2618.97830530487        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.1744E+00  9.5503E-01  6.3058E+01  1.2383E+00  2.3375E+00  2.4115E+00  1.1183E+00  1.0000E-02  3.4831E+00  1.3688E+00
             2.8780E+00
 PARAMETER:  2.6080E-01  5.3987E-02  4.2441E+00  3.1375E-01  9.4910E-01  9.8024E-01  2.1184E-01 -5.6381E+00  1.3479E+00  4.1393E-01
             1.1571E+00
 GRADIENT:  -1.1756E+00 -2.2229E+01 -4.4407E+00 -3.5853E+00  1.4790E+01  2.3662E+01  8.9917E+00  0.0000E+00  3.5476E+01 -8.0689E-01
             3.3664E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2619.40735638639        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.1793E+00  9.2069E-01  7.2093E+01  1.2984E+00  2.3193E+00  2.3040E+00  1.0651E+00  1.0000E-02  3.2016E+00  1.3524E+00
             2.8516E+00
 PARAMETER:  2.6489E-01  1.7364E-02  4.3780E+00  3.6113E-01  9.4128E-01  9.3465E-01  1.6308E-01 -5.5959E+00  1.2637E+00  4.0186E-01
             1.1479E+00
 GRADIENT:  -5.5089E-01 -1.1058E+01 -2.3176E+00 -3.4889E-01  2.7943E+00  4.6963E+00  7.0421E+00  0.0000E+00  1.0867E+01 -2.6172E+00
             9.0546E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2620.21263944929        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  1.1818E+00  9.9204E-01  1.1440E+02  1.2506E+00  2.3408E+00  2.2411E+00  8.4744E-01  1.0000E-02  3.1431E+00  1.3598E+00
             2.8441E+00
 PARAMETER:  2.6700E-01  9.2009E-02  4.8397E+00  3.2361E-01  9.5051E-01  9.0696E-01 -6.5534E-02 -6.0654E+00  1.2452E+00  4.0736E-01
             1.1452E+00
 GRADIENT:  -5.4849E-01  3.5326E-01  1.3168E+00  1.0659E+00 -6.7351E+00 -7.5659E+00  4.6545E+00  0.0000E+00 -1.1892E+01 -2.4851E-01
            -2.7494E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2624.64642895811        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  1.1826E+00  1.0693E+00  1.2891E+02  1.1525E+00  2.3642E+00  2.2745E+00  5.1494E-01  1.0000E-02  3.4311E+00  1.3766E+00
             2.8501E+00
 PARAMETER:  2.6775E-01  1.6697E-01  4.9591E+00  2.4196E-01  9.6046E-01  9.2175E-01 -5.6370E-01 -6.2521E+00  1.3329E+00  4.1964E-01
             1.1473E+00
 GRADIENT:   1.5705E-01  2.2111E+00  1.5717E+00 -1.0724E-01 -2.2796E+00 -7.6675E-01  2.1219E+00  0.0000E+00 -6.9672E-01  2.1561E+00
            -5.4183E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2627.18266540015        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      637
 NPARAMETR:  1.2137E+00  1.0402E+00  1.0151E+02  1.1654E+00  2.3580E+00  2.3388E+00  2.3153E-01  1.0000E-02  3.5950E+00  1.3640E+00
             2.8433E+00
 PARAMETER:  2.9371E-01  1.3944E-01  4.7202E+00  2.5311E-01  9.5781E-01  9.4966E-01 -1.3631E+00 -6.1356E+00  1.3796E+00  4.1041E-01
             1.1450E+00
 GRADIENT:  -5.0286E+00 -3.3618E+00  1.7667E-01 -2.7584E+00 -4.3434E-01 -9.5429E+00  2.2104E-01  0.0000E+00 -1.1344E+01 -6.8884E-01
            -8.4332E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2627.91574177988        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  1.2255E+00  1.1637E+00  9.8395E+01  1.0736E+00  2.3625E+00  2.3894E+00  4.2660E-01  1.0000E-02  4.0238E+00  1.3672E+00
             2.8566E+00
 PARAMETER:  3.0333E-01  2.5160E-01  4.6890E+00  1.7101E-01  9.5970E-01  9.7105E-01 -7.5191E-01 -5.9230E+00  1.4922E+00  4.1277E-01
             1.1496E+00
 GRADIENT:  -1.1412E+00 -5.6156E+00 -2.7621E-01  1.4034E+00  2.7336E-01 -9.0166E-01  8.6707E-01  0.0000E+00  4.0272E+00  1.0586E-01
             1.7292E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2628.00834156433        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  1.2292E+00  1.1931E+00  1.0292E+02  1.0323E+00  2.3665E+00  2.3946E+00  4.2256E-01  1.0000E-02  4.0642E+00  1.3668E+00
             2.8562E+00
 PARAMETER:  3.0636E-01  2.7653E-01  4.7339E+00  1.3181E-01  9.6141E-01  9.7323E-01 -7.6144E-01 -5.9218E+00  1.5022E+00  4.1244E-01
             1.1495E+00
 GRADIENT:   1.4860E-02  6.3027E-02 -2.1017E-03  2.7935E-02  4.6251E-04  2.3244E-03 -2.4074E-03  0.0000E+00 -9.1464E-03 -1.3847E-02
            -8.6824E-03

0ITERATION NO.:   52    OBJECTIVE VALUE:  -2628.00834444452        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1045
 NPARAMETR:  1.2291E+00  1.1928E+00  1.0292E+02  1.0324E+00  2.3665E+00  2.3946E+00  4.2220E-01  1.0000E-02  4.0637E+00  1.3668E+00
             2.8562E+00
 PARAMETER:  3.0632E-01  2.7628E-01  4.7340E+00  1.3189E-01  9.6141E-01  9.7321E-01 -7.6228E-01 -5.9230E+00  1.5021E+00  4.1250E-01
             1.1495E+00
 GRADIENT:  -2.0747E-02  1.2723E-02 -9.5794E-04  3.9410E-03 -4.6726E-03  5.9952E-03 -4.6596E-04  0.0000E+00 -1.8253E-03  3.9133E-04
            -2.5288E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1045
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5078E-03 -3.4465E-02 -1.1716E-06  7.2985E-03 -1.2495E-02
 SE:             2.9775E-02  7.8916E-03  9.4544E-06  2.8261E-02  2.5631E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5961E-01  1.2587E-05  9.0138E-01  7.9621E-01  6.2591E-01

 ETASHRINKSD(%)  2.5025E-01  7.3562E+01  9.9968E+01  5.3222E+00  1.4134E+01
 ETASHRINKVR(%)  4.9988E-01  9.3010E+01  1.0000E+02  1.0361E+01  2.6271E+01
 EBVSHRINKSD(%)  3.3322E-01  8.2180E+01  9.9939E+01  2.5960E+00  1.3259E+01
 EBVSHRINKVR(%)  6.6533E-01  9.6825E+01  1.0000E+02  5.1247E+00  2.4760E+01
 RELATIVEINF(%)  9.9327E+01  1.6672E+00  3.3398E-05  4.9926E+01  6.7323E+01
 EPSSHRINKSD(%)  1.5260E+01
 EPSSHRINKVR(%)  2.8192E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2628.0083444445249     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -973.91898467611418     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.36
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2628.008       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.23E+00  1.19E+00  1.03E+02  1.03E+00  2.37E+00  2.39E+00  4.22E-01  1.00E-02  4.06E+00  1.37E+00  2.86E+00
 


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
+        1.30E+02
 
 TH 2
+       -1.53E+00  2.03E+02
 
 TH 3
+       -1.52E-03 -7.62E-04  3.51E-04
 
 TH 4
+        9.05E+00  5.86E+01 -1.36E-04  5.52E+01
 
 TH 5
+       -3.93E-01 -7.11E+00 -5.69E-02 -3.57E+00  8.03E+01
 
 TH 6
+        3.17E-01  8.97E-01  2.03E-04 -1.38E-01 -6.94E-02  3.38E+01
 
 TH 7
+       -1.47E+00 -5.06E+01  2.91E-03 -6.58E+00  1.54E+00 -3.96E-01  2.80E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.01E-01 -2.73E+01 -3.97E-04  2.27E+00  5.08E-01 -1.16E-01  5.18E+00  0.00E+00  9.77E+00
 
 TH10
+       -6.01E-01 -4.92E-01  6.09E-03 -1.26E+00 -7.59E+00 -1.79E-02 -2.08E+00  0.00E+00  3.38E-02  6.06E+01
 
 TH11
+       -2.70E+00 -1.36E+01 -4.32E-03 -3.16E+00  3.83E-01  1.24E+00  3.83E+00  0.00E+00  1.27E+00  7.76E+00  1.44E+02
 
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
 #CPUT: Total CPU Time in Seconds,       35.526
Stop Time:
Sat Sep 25 05:30:35 CDT 2021
