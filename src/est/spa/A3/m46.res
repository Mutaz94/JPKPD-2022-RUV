Sat Sep 18 10:29:55 CDT 2021
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
$DATA ../../../../data/spa/A3/dat46.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   89.4676758156713        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6736E+02  2.6132E+00  1.3634E+02 -1.8173E+02  1.4614E+02  3.7565E+00 -6.6979E+01 -5.0484E+01 -1.5483E+02 -2.1085E+02
            -2.9380E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1181.30057913643        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9552E-01  1.0141E+00  8.5378E-01  1.3127E+00  9.3540E-01  8.2919E-01  9.8775E-01  1.0032E+00  1.0120E+00  1.1506E+00
             5.2163E+00
 PARAMETER:  9.5513E-02  1.1400E-01 -5.8085E-02  3.7207E-01  3.3215E-02 -8.7312E-02  8.7678E-02  1.0324E-01  1.1195E-01  2.4029E-01
             1.7518E+00
 GRADIENT:  -3.8992E+01  1.4518E+01 -2.6301E+01  6.0298E+01 -7.6468E+00 -2.1086E+01  1.3576E+01  8.7328E+00  3.1714E+01  2.6713E+01
             1.5955E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1218.49025778737        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  9.8550E-01  8.5499E-01  3.5944E-01  1.2190E+00  4.5734E-01  9.7785E-01  4.8111E-01  3.8745E-02  9.1452E-01  4.7041E-01
             4.4370E+00
 PARAMETER:  8.5390E-02 -5.6664E-02 -9.2322E-01  2.9801E-01 -6.8233E-01  7.7601E-02 -6.3165E-01 -3.1508E+00  1.0644E-02 -6.5415E-01
             1.5900E+00
 GRADIENT:  -2.6575E+01  8.1907E+01  3.8987E+01  5.8093E+01 -9.6491E+01  1.3156E+01  1.7935E+00  1.0950E-02  8.1092E+00  7.5237E+00
             7.1883E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1234.01622568487        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.8033E-01  6.2097E-01  2.3616E-01  1.1986E+00  3.2677E-01  9.3221E-01  2.0435E-01  4.3153E-02  1.1031E+00  1.8562E-01
             3.7450E+00
 PARAMETER:  8.0129E-02 -3.7647E-01 -1.3432E+00  2.8119E-01 -1.0185E+00  2.9798E-02 -1.4879E+00 -3.0430E+00  1.9811E-01 -1.5841E+00
             1.4204E+00
 GRADIENT:   1.0930E+01  7.1641E+01  5.1033E+01  5.0853E+01 -9.9128E+01 -1.2438E+01 -1.2233E+00 -5.1819E-02 -2.4608E+00 -3.4954E+00
            -2.6776E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1240.47806180358        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.5077E-01  5.7118E-01  1.5349E-01  1.0699E+00  2.7042E-01  9.8990E-01  9.0421E-02  1.2371E-02  1.3321E+00  1.1480E-01
             3.6965E+00
 PARAMETER:  4.9515E-02 -4.6005E-01 -1.7741E+00  1.6757E-01 -1.2078E+00  8.9847E-02 -2.3033E+00 -4.2924E+00  3.8676E-01 -2.0645E+00
             1.4074E+00
 GRADIENT:  -7.2630E+00  1.5771E+01  7.8277E+00  6.5676E+00 -1.9199E+01 -5.8896E-01 -3.3872E-01 -5.5062E-03 -2.5639E+00 -2.0026E+00
             1.5985E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1245.31257266392        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  9.6623E-01  4.1531E-01  1.5363E-01  1.0688E+00  2.2540E-01  1.0094E+00  3.9375E-01  1.0000E-02  1.4573E+00  5.5616E-01
             3.3365E+00
 PARAMETER:  6.5649E-02 -7.7874E-01 -1.7732E+00  1.6657E-01 -1.3899E+00  1.0936E-01 -8.3203E-01 -7.1442E+00  4.7661E-01 -4.8670E-01
             1.3049E+00
 GRADIENT:   8.2101E+00 -3.0068E+01  8.8154E+00 -2.4440E+01  3.9689E+01  8.9207E+00  1.3682E+00  0.0000E+00  1.2520E+01  4.4112E+00
             2.4402E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1249.29602956475        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      449
 NPARAMETR:  9.4479E-01  4.5378E-01  1.1754E-01  1.0220E+00  2.0880E-01  1.0041E+00  2.0812E-01  1.0000E-02  1.5928E+00  6.0533E-01
             3.1238E+00
 PARAMETER:  4.3209E-02 -6.9015E-01 -2.0410E+00  1.2172E-01 -1.4664E+00  1.0405E-01 -1.4696E+00 -7.9184E+00  5.6552E-01 -4.0198E-01
             1.2390E+00
 GRADIENT:  -1.1897E+00  5.5057E+00  1.9166E-01 -3.2498E+00 -5.8740E+00 -5.4349E-02  6.6698E-01  0.0000E+00  1.3868E+00  4.3792E-02
             3.0759E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1249.67253841269        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      520
 NPARAMETR:  9.4705E-01  4.5482E-01  1.2079E-01  1.0330E+00  2.1128E-01  1.0045E+00  1.3552E-01  1.0000E-02  1.5666E+00  5.9337E-01
             3.1273E+00
 PARAMETER:  4.5592E-02 -6.8786E-01 -2.0137E+00  1.3248E-01 -1.4546E+00  1.0445E-01 -1.8987E+00 -7.4841E+00  5.4893E-01 -4.2193E-01
             1.2402E+00
 GRADIENT:   6.3315E-01  3.7448E+00  2.2828E+00  5.7878E-01 -5.9921E+00  4.7594E-01  2.3757E-01  0.0000E+00  3.3469E-01 -1.8613E+00
             2.0829E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1250.10775982186        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  9.4915E-01  4.8682E-01  1.2942E-01  1.0413E+00  2.2542E-01  9.9942E-01  8.4362E-02  1.0000E-02  1.5221E+00  5.8562E-01
             3.1504E+00
 PARAMETER:  4.7807E-02 -6.1987E-01 -1.9447E+00  1.4051E-01 -1.3898E+00  9.9416E-02 -2.3726E+00 -6.8910E+00  5.2009E-01 -4.3508E-01
             1.2475E+00
 GRADIENT:   9.0569E-01 -1.9086E-01  3.5875E-01 -1.0176E+00  5.7321E-02  1.1169E-01  7.9331E-02  0.0000E+00  3.3108E-02  2.0066E-01
            -4.7164E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1250.14711258273        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  9.4867E-01  4.8769E-01  1.2924E-01  1.0418E+00  2.2559E-01  9.9894E-01  1.4872E-02  1.0000E-02  1.5231E+00  5.8521E-01
             3.1531E+00
 PARAMETER:  4.7305E-02 -6.1808E-01 -1.9461E+00  1.4098E-01 -1.3890E+00  9.8942E-02 -4.1083E+00 -5.8666E+00  5.2076E-01 -4.3578E-01
             1.2484E+00
 GRADIENT:   1.5282E-01 -1.8138E-01 -1.3737E-01 -3.1417E-01  5.2826E-01 -1.0840E-01  2.4320E-03  0.0000E+00  1.6385E-02  1.6184E-02
            -1.9587E-01

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1250.14800762098        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      960
 NPARAMETR:  9.4866E-01  4.8714E-01  1.2921E-01  1.0424E+00  2.2540E-01  9.9931E-01  1.0000E-02  1.0000E-02  1.5231E+00  5.8515E-01
             3.1542E+00
 PARAMETER:  4.7300E-02 -6.1920E-01 -1.9463E+00  1.4149E-01 -1.3899E+00  9.9307E-02 -4.9235E+00 -5.3808E+00  5.2074E-01 -4.3588E-01
             1.2487E+00
 GRADIENT:  -1.0187E-02  1.2847E-02  7.8288E-03  2.2134E-02 -2.9800E-02  1.3856E-03  0.0000E+00  0.0000E+00 -6.8127E-04  3.8182E-03
             3.9524E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      960
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7394E-03 -9.8514E-05  2.5637E-04 -1.7633E-02  8.3042E-03
 SE:             2.8554E-02  1.9516E-04  1.8603E-04  2.5630E-02  2.1324E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5142E-01  6.1371E-01  1.6816E-01  4.9147E-01  6.9696E-01

 ETASHRINKSD(%)  4.3407E+00  9.9346E+01  9.9377E+01  1.4135E+01  2.8561E+01
 ETASHRINKVR(%)  8.4929E+00  9.9996E+01  9.9996E+01  2.6271E+01  4.8965E+01
 EBVSHRINKSD(%)  3.6372E+00  9.9309E+01  9.9477E+01  1.1094E+01  2.8731E+01
 EBVSHRINKVR(%)  7.1421E+00  9.9995E+01  9.9997E+01  2.0957E+01  4.9208E+01
 RELATIVEINF(%)  8.7068E+01  1.9770E-04  4.1500E-04  6.3837E+01  1.2520E+00
 EPSSHRINKSD(%)  3.0159E+01
 EPSSHRINKVR(%)  5.1222E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1250.1480076209828     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -514.99718105724457     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.07
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
 





 #OBJV:********************************************    -1250.148       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.49E-01  4.87E-01  1.29E-01  1.04E+00  2.25E-01  9.99E-01  1.00E-02  1.00E-02  1.52E+00  5.85E-01  3.15E+00
 


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
+        1.15E+03
 
 TH 2
+       -2.18E+01  2.26E+03
 
 TH 3
+       -7.09E+02  4.37E+03  1.85E+04
 
 TH 4
+       -1.12E+01  8.39E+01 -5.75E+02  3.52E+02
 
 TH 5
+        4.77E+02 -8.04E+03 -1.99E+04 -7.45E+01  3.26E+04
 
 TH 6
+        4.96E+00 -7.21E+00  1.02E+02 -1.27E+01  1.53E+01  1.74E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.29E+01 -2.53E+01  2.25E+02 -8.96E+00  6.35E+01 -1.69E+00  0.00E+00  0.00E+00  5.06E+01
 
 TH10
+       -8.03E+00 -7.43E+01 -1.45E+02  5.83E+00  4.12E+02  9.93E+00  0.00E+00  0.00E+00 -3.28E-01  1.35E+02
 
 TH11
+       -1.94E+01 -1.23E+01 -2.39E+01 -4.77E+00  2.69E+01  2.27E+00  0.00E+00  0.00E+00  5.60E+00  2.54E+01  3.06E+01
 
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
 #CPUT: Total CPU Time in Seconds,       16.683
Stop Time:
Sat Sep 18 10:30:13 CDT 2021
