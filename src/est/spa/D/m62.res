Sat Sep 18 15:30:27 CDT 2021
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
$DATA ../../../../data/spa/D/dat62.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:   14468.2409008517        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.5419E+01  3.0711E+02 -9.8265E+01  4.4900E+01  3.5865E+02 -1.9570E+03 -7.1379E+02 -6.7228E+01 -1.6441E+03 -6.4679E+02
            -2.7253E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -616.177716386005        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.6890E+00  1.0944E+00  9.4969E-01  1.7410E+00  1.2086E+00  2.1419E+00  1.1779E+00  9.7457E-01  1.4059E+00  1.0580E+00
             1.4150E+01
 PARAMETER:  6.2411E-01  1.9017E-01  4.8380E-02  6.5448E-01  2.8946E-01  8.6171E-01  2.6370E-01  7.4245E-02  4.4064E-01  1.5642E-01
             2.7497E+00
 GRADIENT:   5.9526E+01  3.2173E+01  2.9179E-01  4.1166E+01 -1.1094E+01  3.4350E+01 -3.0519E+00  5.3196E+00 -1.0802E+01  2.6294E+00
             2.7513E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -632.575962140314        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.5846E+00  8.5671E-01  1.4096E+00  1.9408E+00  2.0082E+00  2.1218E+00  1.9904E+00  3.4417E-01  2.0832E+00  2.1069E+00
             1.2554E+01
 PARAMETER:  5.6034E-01 -5.4659E-02  4.4333E-01  7.6309E-01  7.9724E-01  8.5225E-01  7.8831E-01 -9.6663E-01  8.3389E-01  8.4524E-01
             2.6301E+00
 GRADIENT:   4.0608E+01  1.6945E+01 -2.5782E+00  2.8571E+01 -7.0741E+00 -2.0812E+00  3.8412E+00  3.0184E-01  5.9175E+00  3.8256E+00
             4.0857E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -652.684240002619        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1694E+00  3.4243E-01  1.1347E+00  1.9318E+00  3.2018E+00  2.2004E+00  1.1761E+00  1.3732E-01  1.4154E+00  5.2625E+00
             1.1202E+01
 PARAMETER:  2.5650E-01 -9.7169E-01  2.2639E-01  7.5846E-01  1.2637E+00  8.8864E-01  2.6223E-01 -1.8854E+00  4.4742E-01  1.7606E+00
             2.5160E+00
 GRADIENT:  -4.6829E+01  1.1977E+01  1.8518E+01  4.2582E+01 -6.7383E+00  2.1877E+01  2.4476E-01 -4.9081E-02 -2.2430E+01 -2.1409E+00
             6.8030E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -733.157731199446        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  8.6623E-01  9.5222E-02  1.1131E-01  9.6187E-01  4.2146E+01  1.8947E+00  1.0000E-02  1.0000E-02  1.3048E+00  8.4686E+00
             9.2159E+00
 PARAMETER: -4.3602E-02 -2.2515E+00 -2.0954E+00  6.1125E-02  3.8411E+00  7.3909E-01 -4.6138E+00 -7.4984E+00  3.6604E-01  2.2364E+00
             2.3209E+00
 GRADIENT:   2.3027E+01  1.2258E+01 -7.7642E+01  1.4905E+02 -6.4874E-01 -1.7866E+01  0.0000E+00  0.0000E+00  7.6232E+00  4.9967E-01
            -9.0333E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -781.233831826889        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  3.3660E-01  1.0000E-02  1.0242E-02  1.3975E-01  1.6490E+04  1.7105E+00  1.0000E-02  1.0000E-02  5.1780E-01  6.3369E+01
             1.0172E+01
 PARAMETER: -9.8887E-01 -5.0832E+00 -4.4813E+00 -1.8679E+00  9.8105E+00  6.3679E-01 -1.9035E+01 -1.6856E+01 -5.5816E-01  4.2490E+00
             2.4197E+00
 GRADIENT:   6.6348E+00  0.0000E+00 -3.1413E+01  4.5071E+01 -3.9300E-04  7.8306E+00  0.0000E+00  0.0000E+00  5.1336E+00  1.3409E-05
            -1.7190E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -782.941746849930        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      453
 NPARAMETR:  3.2308E-01  1.0000E-02  1.0000E-02  1.3554E-01  1.6906E+04  1.5856E+00  1.0000E-02  1.0000E-02  3.1545E-01  7.8115E+01
             1.0196E+01
 PARAMETER: -1.0299E+00 -5.1390E+00 -4.5289E+00 -1.8985E+00  9.8354E+00  5.6098E-01 -1.8136E+01 -1.6841E+01 -1.0538E+00  4.4582E+00
             2.4220E+00
 GRADIENT:  -6.6632E+00  0.0000E+00  0.0000E+00  9.6948E+00 -5.0058E-04 -6.4553E+00  0.0000E+00  0.0000E+00  2.3859E+00  2.1400E-05
            -1.3923E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -784.605298768682        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  3.2491E-01  1.0000E-02  1.0000E-02  1.3600E-01  9.1420E+03  1.6095E+00  1.0000E-02  1.0000E-02  6.5402E-02  1.0536E+02
             1.0525E+01
 PARAMETER: -1.0242E+00 -5.1627E+00 -4.5536E+00 -1.8951E+00  9.2206E+00  5.7591E-01 -1.4386E+01 -1.7384E+01 -2.6272E+00  4.7574E+00
             2.4538E+00
 GRADIENT:  -8.2134E-01  0.0000E+00  0.0000E+00  3.3715E+00 -9.3364E-04 -7.0370E-01  0.0000E+00  0.0000E+00  1.2791E-01  1.2781E-04
             2.1478E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -784.721385879478        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  3.2799E-01  1.0000E-02  1.0000E-02  1.3613E-01  5.1131E+03  1.6182E+00  1.0000E-02  1.0000E-02  1.3561E-02  1.4821E+02
             1.0558E+01
 PARAMETER: -1.0148E+00 -5.1856E+00 -4.5767E+00 -1.8941E+00  8.6396E+00  5.8131E-01 -1.0344E+01 -1.7777E+01 -4.2005E+00  5.0986E+00
             2.4569E+00
 GRADIENT:  -1.5170E+00  0.0000E+00  0.0000E+00 -2.0579E-01 -2.5678E-03 -3.7034E-01  0.0000E+00  0.0000E+00  5.2309E-03  1.4720E-03
            -1.2549E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -784.724991524504        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      818
 NPARAMETR:  3.2863E-01  1.0000E-02  1.0000E-02  1.3616E-01  4.9408E+03  1.6204E+00  1.0000E-02  1.0000E-02  1.2511E-02  1.5012E+02
             1.0575E+01
 PARAMETER: -1.0128E+00 -5.1861E+00 -4.5776E+00 -1.8939E+00  8.6053E+00  5.8265E-01 -1.0146E+01 -1.7803E+01 -4.2811E+00  5.1114E+00
             2.4585E+00
 GRADIENT:  -3.1571E-01  0.0000E+00  0.0000E+00 -2.6325E-02 -5.8643E-03 -8.7208E-02  0.0000E+00  0.0000E+00  4.4576E-03  7.7998E-03
            -2.6895E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -784.725499818608        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      995
 NPARAMETR:  3.2881E-01  1.0000E-02  1.0000E-02  1.3610E-01  4.9415E+03  1.6210E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.4968E+02
             1.0586E+01
 PARAMETER: -1.0123E+00 -5.1859E+00 -4.5776E+00 -1.8944E+00  8.6054E+00  5.8306E-01 -1.0144E+01 -1.7806E+01 -4.5061E+00  5.1085E+00
             2.4595E+00
 GRADIENT:  -4.6961E-02  0.0000E+00  0.0000E+00 -1.2757E+00 -3.4893E-03  4.4298E-02  0.0000E+00  0.0000E+00  1.1730E-03  3.1044E-03
             5.4838E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -784.727108657355        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1174
 NPARAMETR:  3.2884E-01  1.0000E-02  1.0000E-02  1.3618E-01  1.0761E+04  1.6217E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.3748E+02
             1.0576E+01
 PARAMETER: -1.0122E+00 -5.1859E+00 -4.5776E+00 -1.8938E+00  9.3837E+00  5.8347E-01 -1.0144E+01 -1.7806E+01 -4.6239E+00  5.0234E+00
             2.4586E+00
 GRADIENT:   1.2077E-01  0.0000E+00  0.0000E+00  2.4803E-01 -7.9785E-04  8.0224E-02  0.0000E+00  0.0000E+00  0.0000E+00  1.5770E-04
            -2.3892E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -784.727811790869        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1349
 NPARAMETR:  3.2892E-01  1.0000E-02  1.0000E-02  1.3617E-01  2.1389E+05  1.6213E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0125E+02
             1.0581E+01
 PARAMETER: -1.0119E+00 -5.1859E+00 -4.5776E+00 -1.8938E+00  1.2373E+01  5.8321E-01 -1.0144E+01 -1.7806E+01 -4.6239E+00  4.7176E+00
             2.4591E+00
 GRADIENT:   2.4557E-01  0.0000E+00  0.0000E+00  9.3319E-02 -3.2500E-05  2.4682E-02  0.0000E+00  0.0000E+00  0.0000E+00  2.0657E-07
             7.4821E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -784.727895972021        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1528
 NPARAMETR:  3.2880E-01  1.0000E-02  1.0000E-02  1.3617E-01  3.3657E+07  1.6210E+00  1.0000E-02  1.0000E-02  1.0000E-02  6.0209E+01
             1.0579E+01
 PARAMETER: -1.0123E+00 -5.1859E+00 -4.5776E+00 -1.8939E+00  1.7432E+01  5.8306E-01 -1.0144E+01 -1.7806E+01 -4.6239E+00  4.1978E+00
             2.4589E+00
 GRADIENT:   5.0033E-04  0.0000E+00  0.0000E+00  2.9126E-03 -2.0574E-07  4.6536E-05  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -6.1089E-04

0ITERATION NO.:   67    OBJECTIVE VALUE:  -784.727896007319        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1599
 NPARAMETR:  3.2882E-01  1.0000E-02  1.0000E-02  1.3616E-01  3.4524E+07  1.6210E+00  1.0000E-02  1.0000E-02  1.0000E-02  5.9695E+01
             1.0579E+01
 PARAMETER: -1.0123E+00 -5.1859E+00 -4.5776E+00 -1.8939E+00  1.7475E+01  5.8305E-01 -1.0144E+01 -1.7806E+01 -4.6239E+00  4.1934E+00
             2.4589E+00
 GRADIENT:  -1.4903E-02  0.0000E+00  0.0000E+00  1.1013E-02  4.3704E-07 -8.2522E-04  0.0000E+00  0.0000E+00  0.0000E+00  1.1328E-06
             1.1052E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1599
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4634E-04 -5.4952E-07  7.6116E-05 -1.3779E-04  1.4494E-09
 SE:             2.9152E-02  1.9926E-05  3.3534E-04  3.8052E-04  3.4616E-09
 N:                     100         100         100         100         100

 P VAL.:         9.7684E-01  9.7800E-01  8.2044E-01  7.1726E-01  6.7544E-01

 ETASHRINKSD(%)  2.3385E+00  9.9933E+01  9.8877E+01  9.8725E+01  1.0000E+02
 ETASHRINKVR(%)  4.6223E+00  1.0000E+02  9.9987E+01  9.9984E+01  1.0000E+02
 EBVSHRINKSD(%)  2.7304E+00  9.9897E+01  9.8866E+01  9.8685E+01  1.0000E+02
 EBVSHRINKVR(%)  5.3863E+00  1.0000E+02  9.9987E+01  9.9983E+01  1.0000E+02
 RELATIVEINF(%)  5.6374E+00  1.0000E-10  3.3594E-05  4.4374E-05  1.0000E-10
 EPSSHRINKSD(%)  7.6106E+00
 EPSSHRINKVR(%)  1.4642E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -784.72789600731880     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -49.577069443580626     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.21
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -784.728       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.29E-01  1.00E-02  1.00E-02  1.36E-01  3.51E+07  1.62E+00  1.00E-02  1.00E-02  1.00E-02  5.99E+01  1.06E+01
 


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
+        3.63E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        0.00E+00  0.00E+00  0.00E+00
 
 TH 4
+       -1.08E+03  0.00E+00  0.00E+00  6.88E+04
 
 TH 5
+       -1.46E-11  0.00E+00  0.00E+00 -2.81E-11  8.82E-21
 
 TH 6
+       -2.95E+01  0.00E+00  0.00E+00 -1.28E+02  1.21E-12  4.44E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        5.22E-05  0.00E+00  0.00E+00 -4.72E-05  1.67E-14  8.88E-07  0.00E+00  0.00E+00  0.00E+00  1.84E-08
 
 TH11
+       -2.52E+01  0.00E+00  0.00E+00 -2.33E+01  2.52E-13  7.48E-01  0.00E+00  0.00E+00  0.00E+00 -1.29E-07  3.69E+00
 
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
 #CPUT: Total CPU Time in Seconds,       26.819
Stop Time:
Sat Sep 18 15:30:56 CDT 2021
