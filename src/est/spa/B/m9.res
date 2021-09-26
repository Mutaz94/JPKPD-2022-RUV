Sat Sep 25 07:06:18 CDT 2021
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
$DATA ../../../../data/spa/B/dat9.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1632.87432452960        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8291E+02 -5.5453E+01 -2.6558E+01 -3.9799E+01  6.5381E+01 -2.1687E+00  9.4769E-01  1.8271E+00  5.9364E+00 -5.7505E+00
             5.1659E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1638.04155822640        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9880E-01  9.6948E-01  9.6781E-01  1.0132E+00  9.4806E-01  9.9775E-01  9.9495E-01  1.0003E+00  9.7823E-01  1.0192E+00
             8.6535E-01
 PARAMETER:  9.8795E-02  6.9001E-02  6.7281E-02  1.1313E-01  4.6666E-02  9.7750E-02  9.4932E-02  1.0032E-01  7.7990E-02  1.1898E-01
            -4.4623E-02
 GRADIENT:   1.9466E+02 -4.7066E+01 -1.6795E+01 -4.8976E+01  3.2257E+01 -2.7180E+00 -4.3962E-01  4.0508E+00  1.4393E+00  3.4659E+00
             2.5627E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1640.01389313233        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  9.9162E-01  7.3909E-01  9.8496E-01  1.2141E+00  8.2727E-01  1.0470E+00  1.0412E+00  8.8663E-01  9.0040E-01  9.5415E-01
             8.6500E-01
 PARAMETER:  9.1581E-02 -2.0234E-01  8.4850E-02  2.9403E-01 -8.9627E-02  1.4591E-01  1.4042E-01 -2.0327E-02 -4.9179E-03  5.3061E-02
            -4.5031E-02
 GRADIENT:   1.2068E+02  1.3432E+01 -2.2004E+00  3.6969E+01 -1.0448E+01  1.4404E+01 -2.6196E+00  1.0361E+00  1.9581E+00  1.0279E+00
             1.8126E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1644.99349502883        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  9.3984E-01  7.8280E-01  9.5857E-01  1.1700E+00  8.4048E-01  9.8976E-01  1.3450E+00  8.3917E-01  8.5896E-01  9.4409E-01
             8.5213E-01
 PARAMETER:  3.7954E-02 -1.4487E-01  5.7689E-02  2.5703E-01 -7.3779E-02  8.9707E-02  3.9642E-01 -7.5344E-02 -5.2033E-02  4.2466E-02
            -6.0014E-02
 GRADIENT:   1.3501E+01  4.3203E+00 -3.2512E+00  8.8761E+00  2.4242E+00  2.9238E-01  3.1800E+00  8.4817E-01  3.0419E+00  2.1104E+00
            -3.1009E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1645.21885760769        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      622
 NPARAMETR:  9.3370E-01  6.7993E-01  1.0360E+00  1.2365E+00  8.3546E-01  9.9045E-01  1.3493E+00  8.9223E-01  8.3579E-01  9.6106E-01
             8.5836E-01
 PARAMETER:  3.1402E-02 -2.8577E-01  1.3539E-01  3.1226E-01 -7.9771E-02  9.0400E-02  3.9957E-01 -1.4034E-02 -7.9375E-02  6.0279E-02
            -5.2733E-02
 GRADIENT:   1.8541E+00  3.6763E+00  1.4843E+00  5.9330E+00 -5.3557E+00  1.2768E+00 -3.6483E-01  9.9705E-02 -5.5828E-01  1.4859E+00
            -6.1455E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1645.47827288219        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      797
 NPARAMETR:  9.2834E-01  4.7550E-01  1.2579E+00  1.3761E+00  8.6501E-01  9.7964E-01  1.3850E+00  1.0838E+00  8.1076E-01  1.0055E+00
             8.6234E-01
 PARAMETER:  2.5647E-02 -6.4339E-01  3.2942E-01  4.1925E-01 -4.5014E-02  7.9434E-02  4.2570E-01  1.8047E-01 -1.0978E-01  1.0547E-01
            -4.8101E-02
 GRADIENT:  -2.9672E+00  4.3142E+00  3.1497E+00  1.2902E+01 -2.7018E+00 -1.7115E+00  2.8464E-01 -6.8870E-01  3.5972E-01 -1.5228E+00
             2.9334E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1645.84969946036        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      972
 NPARAMETR:  9.2584E-01  2.8702E-01  1.3249E+00  1.4856E+00  8.4073E-01  9.7902E-01  1.2835E+00  1.1600E+00  7.7852E-01  1.0261E+00
             8.6029E-01
 PARAMETER:  2.2943E-02 -1.1482E+00  3.8132E-01  4.9580E-01 -7.3482E-02  7.8793E-02  3.4956E-01  2.4842E-01 -1.5037E-01  1.2574E-01
            -5.0488E-02
 GRADIENT:  -2.2081E+00  2.3899E-01 -1.0605E+00  1.9995E-01 -9.3870E-01 -9.2304E-01  2.5543E-01  6.8621E-01  1.2469E+00  9.3668E-01
             6.9001E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1645.95463111451        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1149
 NPARAMETR:  9.2444E-01  1.7194E-01  1.4259E+00  1.5616E+00  8.4477E-01  9.7904E-01  1.0671E+00  1.2684E+00  7.4811E-01  1.0259E+00
             8.6011E-01
 PARAMETER:  2.1434E-02 -1.6606E+00  4.5483E-01  5.4573E-01 -6.8690E-02  7.8812E-02  1.6493E-01  3.3775E-01 -1.9021E-01  1.2560E-01
            -5.0691E-02
 GRADIENT:  -9.3821E-01  3.9251E-01  3.0165E-01  4.5271E+00 -7.8746E-01 -3.0385E-01  1.2999E-01  1.6239E-01  1.3982E-01 -9.8099E-02
            -2.7802E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1645.97259029685        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1324
 NPARAMETR:  9.2390E-01  1.2309E-01  1.4437E+00  1.5899E+00  8.3873E-01  9.7890E-01  8.8009E-01  1.2880E+00  7.3484E-01  1.0239E+00
             8.6059E-01
 PARAMETER:  2.0851E-02 -1.9948E+00  4.6724E-01  5.6367E-01 -7.5863E-02  7.8676E-02 -2.7729E-02  3.5309E-01 -2.0810E-01  1.2361E-01
            -5.0139E-02
 GRADIENT:  -2.4712E-01 -4.1926E-02  5.3799E-02 -2.0621E-01  8.5693E-02 -8.9566E-02  5.9629E-02 -6.3054E-02 -2.2730E-02 -1.2439E-02
             2.6470E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1645.97268436988        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1501
 NPARAMETR:  9.2394E-01  1.1741E-01  1.4456E+00  1.5933E+00  8.3787E-01  9.7909E-01  8.4684E-01  1.2920E+00  7.3359E-01  1.0230E+00
             8.6050E-01
 PARAMETER:  2.0896E-02 -2.0421E+00  4.6849E-01  5.6582E-01 -7.6895E-02  7.8873E-02 -6.6238E-02  3.5617E-01 -2.0981E-01  1.2277E-01
            -5.0246E-02
 GRADIENT:   9.1372E-02 -6.4232E-02 -6.9400E-02 -4.0929E-01  1.0314E-01  1.6824E-02  5.2065E-02  3.2981E-02  1.1300E-01  4.8720E-03
             9.4824E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1645.98435868833        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1684
 NPARAMETR:  9.2404E-01  1.1198E-01  1.4440E+00  1.5943E+00  8.3693E-01  9.7936E-01  4.1028E-01  1.2920E+00  7.3395E-01  1.0232E+00
             8.6038E-01
 PARAMETER:  2.1005E-02 -2.0895E+00  4.6739E-01  5.6641E-01 -7.8017E-02  7.9141E-02 -7.9092E-01  3.5619E-01 -2.0931E-01  1.2293E-01
            -5.0381E-02
 GRADIENT:   6.0628E-01 -4.0456E-01 -7.2242E-01 -5.7657E+00  1.4542E+00  1.6138E-01  1.3643E-02  1.2194E-01  6.7875E-01  1.6632E-02
             4.6463E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1646.05563720016        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1863
 NPARAMETR:  9.2511E-01  1.9087E-01  1.4218E+00  1.5482E+00  8.4922E-01  9.7995E-01  3.3650E-02  1.2586E+00  7.5750E-01  1.0329E+00
             8.6066E-01
 PARAMETER:  2.2156E-02 -1.5561E+00  4.5191E-01  5.3708E-01 -6.3436E-02  7.9748E-02 -3.2917E+00  3.3002E-01 -1.7773E-01  1.3232E-01
            -5.0051E-02
 GRADIENT:  -5.5822E-02  4.6178E-02  1.9265E-01  1.7251E+00 -3.2613E-01 -1.8066E-02  2.6438E-04 -2.9170E-04  9.6588E-03 -6.0745E-02
            -3.7335E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1646.10164220463        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2039
 NPARAMETR:  9.2675E-01  3.0052E-01  1.3802E+00  1.4799E+00  8.6506E-01  9.8131E-01  1.0000E-02  1.2040E+00  7.9424E-01  1.0469E+00
             8.6098E-01
 PARAMETER:  2.3923E-02 -1.1022E+00  4.2224E-01  4.9196E-01 -4.4958E-02  8.1130E-02 -9.9884E+00  2.8564E-01 -1.3037E-01  1.4586E-01
            -4.9684E-02
 GRADIENT:  -2.9535E-01  5.3346E-01  4.0278E-01  3.6073E+00 -8.9270E-01 -3.0628E-02  0.0000E+00 -7.4931E-02 -3.6965E-01  2.4705E-02
            -5.2132E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1646.11034044098        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2214
 NPARAMETR:  9.2759E-01  3.5141E-01  1.3637E+00  1.4464E+00  8.7429E-01  9.8198E-01  1.0000E-02  1.1854E+00  8.1364E-01  1.0533E+00
             8.6113E-01
 PARAMETER:  2.4834E-02 -9.4581E-01  4.1018E-01  4.6908E-01 -3.4342E-02  8.1818E-02 -1.4204E+01  2.7009E-01 -1.0624E-01  1.5191E-01
            -4.9511E-02
 GRADIENT:  -5.0191E-02  8.0983E-02  3.9329E-02  3.9708E-01 -8.7328E-02 -5.9952E-03  0.0000E+00 -1.7096E-02 -4.9817E-02 -1.4163E-02
            -8.1865E-04

0ITERATION NO.:   68    OBJECTIVE VALUE:  -1646.11037384163        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:     2320
 NPARAMETR:  9.2761E-01  3.5135E-01  1.3638E+00  1.4463E+00  8.7425E-01  9.8200E-01  1.0000E-02  1.1855E+00  8.1366E-01  1.0534E+00
             8.6110E-01
 PARAMETER:  2.4849E-02 -9.4636E-01  4.1015E-01  4.6901E-01 -3.4294E-02  8.1828E-02 -1.4204E+01  2.7022E-01 -1.0612E-01  1.5197E-01
            -4.9508E-02
 GRADIENT:  -2.1437E-02 -2.9466E-02 -4.9734E-02 -4.0888E-02  4.8763E-02 -2.7484E-03  0.0000E+00  3.7742E-03  1.3757E-02 -8.2246E-03
             1.1490E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2320
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3519E-04 -1.4229E-04 -3.5439E-02 -4.9760E-03 -3.5980E-02
 SE:             2.9870E-02  6.5816E-05  1.7581E-02  2.9394E-02  2.2045E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8838E-01  3.0619E-02  4.3831E-02  8.6557E-01  1.0266E-01

 ETASHRINKSD(%)  1.0000E-10  9.9780E+01  4.1100E+01  1.5267E+00  2.6146E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.5308E+01  3.0301E+00  4.5456E+01
 EBVSHRINKSD(%)  3.2522E-01  9.9792E+01  4.5363E+01  1.9737E+00  2.2340E+01
 EBVSHRINKVR(%)  6.4938E-01  1.0000E+02  7.0148E+01  3.9084E+00  3.9690E+01
 RELATIVEINF(%)  9.8121E+01  2.6868E-05  7.1003E+00  7.3199E+00  9.4264E+00
 EPSSHRINKSD(%)  4.6403E+01
 EPSSHRINKVR(%)  7.1274E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1646.1103738416293     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -910.95954727789115     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.92
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.47
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1646.110       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.28E-01  3.51E-01  1.36E+00  1.45E+00  8.74E-01  9.82E-01  1.00E-02  1.19E+00  8.14E-01  1.05E+00  8.61E-01
 


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
+        1.33E+03
 
 TH 2
+       -2.74E+01  4.10E+02
 
 TH 3
+        1.79E+00  6.24E+01  1.36E+02
 
 TH 4
+       -8.50E+00  4.87E+02 -2.22E+01  7.71E+02
 
 TH 5
+        1.58E+00 -2.47E+02 -2.61E+02 -4.65E+01  8.15E+02
 
 TH 6
+       -5.77E+00 -1.72E+00  5.79E-01 -1.44E+00 -5.10E+00  2.07E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.43E+00 -4.14E+00 -3.31E+01 -3.87E+00  8.19E-02 -9.07E-01  0.00E+00  3.09E+01
 
 TH 9
+        1.23E+00 -1.02E+02  3.53E+00  5.24E-01  4.56E+00  4.48E+00  0.00E+00  2.70E-01  2.81E+02
 
 TH10
+        1.42E+00  9.28E+00 -6.09E+00 -4.89E-01 -7.93E+01  6.18E-02  0.00E+00  1.68E+01  3.79E+00  7.20E+01
 
 TH11
+       -8.05E+00 -1.34E+01 -9.70E+00 -7.55E+00 -9.52E+00  9.33E+00  0.00E+00  9.36E+00  6.77E+00  1.14E+01  2.83E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.442
Stop Time:
Sat Sep 25 07:06:53 CDT 2021
