Wed Sep 29 12:45:41 CDT 2021
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
$DATA ../../../../data/spa/A2/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1125.23262084870        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2108E+02  3.3245E+01  5.1130E+01 -7.2370E+00  8.5679E+01  5.1443E+01 -1.5362E+01 -7.8011E+00 -5.1690E+01 -7.2975E+01
            -9.5372E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1444.59884716410        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1170E+00  9.6167E-01  9.3738E-01  1.0814E+00  8.9350E-01  9.3491E-01  9.5887E-01  9.3534E-01  1.0685E+00  9.5941E-01
             2.3705E+00
 PARAMETER:  2.1060E-01  6.0918E-02  3.5337E-02  1.7824E-01 -1.2608E-02  3.2695E-02  5.8005E-02  3.3155E-02  1.6628E-01  5.8560E-02
             9.6309E-01
 GRADIENT:   2.8233E+02  1.5190E+01  1.6039E+01  1.1650E+01 -1.6303E+01  4.6452E+00  4.1531E+00  4.5951E+00  3.2962E+00  7.1959E+00
             2.8283E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1449.92987063865        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.1032E+00  7.6996E-01  3.4451E-01  1.1347E+00  4.5278E-01  9.6393E-01  7.6109E-01  4.3269E-01  1.0902E+00  4.1484E-01
             2.2502E+00
 PARAMETER:  1.9818E-01 -1.6142E-01 -9.6562E-01  2.2636E-01 -6.9235E-01  6.3260E-02 -1.7300E-01 -7.3772E-01  1.8639E-01 -7.7987E-01
             9.1101E-01
 GRADIENT:   2.2372E+02  6.1326E+01  3.9067E+01  1.3264E+02 -1.6052E+01  1.2902E+01 -2.2950E+01 -1.0341E+00  2.5944E+01 -8.9750E+00
            -1.5900E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1456.92875314033        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      336
 NPARAMETR:  1.0822E+00  7.9794E-01  4.0033E-01  1.1186E+00  5.1805E-01  9.4446E-01  1.1698E+00  3.3217E-01  1.0954E+00  4.2014E-01
             2.2750E+00
 PARAMETER:  1.7904E-01 -1.2572E-01 -8.1547E-01  2.1206E-01 -5.5769E-01  4.2857E-02  2.5686E-01 -1.0021E+00  1.9110E-01 -7.6718E-01
             9.2198E-01
 GRADIENT:   3.4622E+01  1.6718E+00 -2.6115E+01  6.5834E+01  4.8108E+01  2.5912E+00  1.9171E+00  1.0495E+00  3.0255E+01  1.5939E+00
            -1.0065E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1462.83081030494        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      511
 NPARAMETR:  1.0696E+00  5.1758E-01  4.2918E-01  1.2538E+00  4.4326E-01  9.3688E-01  1.6362E+00  3.9712E-01  8.4605E-01  3.5759E-01
             2.2881E+00
 PARAMETER:  1.6728E-01 -5.5859E-01 -7.4588E-01  3.2615E-01 -7.1359E-01  3.4795E-02  5.9236E-01 -8.2351E-01 -6.7172E-02 -9.2836E-01
             9.2771E-01
 GRADIENT:   3.3878E+00  1.1870E+01  9.5698E+00  6.5337E+01 -2.4527E+00  1.6514E+00 -4.7065E+00 -1.7458E+00 -1.8075E+00 -4.5688E+00
            -7.2958E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1468.05558451462        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  1.0685E+00  3.0032E-01  3.3683E-01  1.2649E+00  3.3562E-01  9.3654E-01  2.0949E+00  2.1164E-01  8.5773E-01  5.9318E-01
             2.1868E+00
 PARAMETER:  1.6629E-01 -1.1029E+00 -9.8819E-01  3.3501E-01 -9.9177E-01  3.4439E-02  8.3950E-01 -1.4529E+00 -5.3461E-02 -4.2225E-01
             8.8243E-01
 GRADIENT:   5.5324E+00 -5.2030E-01 -5.5657E+00  2.4239E+01  9.2808E+00  8.3918E-01 -3.0273E-01  1.9955E-01  7.4722E+00  1.8291E-01
             2.4950E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1468.59297975075        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      862
 NPARAMETR:  1.0646E+00  2.3477E-01  3.0442E-01  1.2447E+00  3.0274E-01  9.3273E-01  2.2927E+00  1.1585E-01  8.3826E-01  6.3962E-01
             2.1662E+00
 PARAMETER:  1.6257E-01 -1.3491E+00 -1.0894E+00  3.1893E-01 -1.0949E+00  3.0360E-02  9.2973E-01 -2.0554E+00 -7.6425E-02 -3.4687E-01
             8.7298E-01
 GRADIENT:  -9.3343E-01 -8.1982E-01 -1.2814E+00  8.1621E-01  3.0402E+00  1.0559E-01 -3.6628E-01  5.2993E-02  6.8993E-01  4.0396E-01
             8.0795E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1468.62869630651        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1040
 NPARAMETR:  1.0647E+00  2.3515E-01  3.0120E-01  1.2415E+00  3.0036E-01  9.3237E-01  2.2848E+00  4.0025E-02  8.3742E-01  6.4408E-01
             2.1650E+00
 PARAMETER:  1.6269E-01 -1.3475E+00 -1.1000E+00  3.1632E-01 -1.1028E+00  2.9970E-02  9.2630E-01 -3.1183E+00 -7.7431E-02 -3.3994E-01
             8.7240E-01
 GRADIENT:  -1.1139E+00 -3.1003E-01 -4.5634E-02  1.3869E+00  3.3611E-01 -9.9132E-02 -2.5246E-02  5.9896E-03 -6.3598E-02  5.7477E-01
             5.8736E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1471.77650471499        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1225
 NPARAMETR:  1.0573E+00  5.1890E-01  2.4758E-01  1.0960E+00  3.2527E-01  9.4946E-01  1.4160E+00  1.0000E-02  8.7406E-01  4.1018E-01
             2.1941E+00
 PARAMETER:  1.5568E-01 -5.5605E-01 -1.2960E+00  1.9166E-01 -1.0231E+00  4.8141E-02  4.4782E-01 -4.8218E+01 -3.4606E-02 -7.9115E-01
             8.8576E-01
 GRADIENT:  -1.8585E+01  3.9235E+00  3.0969E+00  1.7262E+01  2.2106E+00  2.3397E+00  3.8603E+00  0.0000E+00 -1.4480E+00  1.7806E+00
             1.3742E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1473.58200006828        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1401
 NPARAMETR:  1.0524E+00  6.6922E-01  2.0425E-01  9.8267E-01  3.4181E-01  9.4511E-01  1.1367E+00  1.0000E-02  9.3095E-01  2.5200E-01
             2.1397E+00
 PARAMETER:  1.5111E-01 -3.0164E-01 -1.4884E+00  8.2517E-02 -9.7350E-01  4.3549E-02  2.2812E-01 -7.7538E+01  2.8446E-02 -1.2783E+00
             8.6065E-01
 GRADIENT:  -5.3274E+00 -2.2267E+00 -1.4261E+00  1.0265E+00  1.2025E+01  3.0680E-01 -8.3407E-01  0.0000E+00 -2.6650E-02  5.5605E-01
             2.4039E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1473.95981578360        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1578            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0559E+00  6.1610E-01  1.9833E-01  9.9429E-01  3.2173E-01  9.4758E-01  1.2062E+00  1.0000E-02  9.2554E-01  2.6230E-01
             2.1218E+00
 PARAMETER:  1.5441E-01 -3.8435E-01 -1.5178E+00  9.4274E-02 -1.0340E+00  4.6157E-02  2.8745E-01 -7.4023E+01  2.2620E-02 -1.2383E+00
             8.5228E-01
 GRADIENT:   1.1635E+02  1.4192E+01  2.5429E+01  1.5683E+01  9.0034E+01  6.0392E+00  3.1402E+00  0.0000E+00  1.3400E+00  7.4488E-01
             8.1640E+00

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1473.97621396143        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1705
 NPARAMETR:  1.0564E+00  6.1579E-01  1.9867E-01  9.9504E-01  3.2145E-01  9.4759E-01  1.1898E+00  1.0000E-02  9.2589E-01  2.7385E-01
             2.1160E+00
 PARAMETER:  1.5489E-01 -3.8485E-01 -1.5161E+00  9.5025E-02 -1.0349E+00  4.6169E-02  2.7376E-01 -7.4023E+01  2.2997E-02 -1.1952E+00
             8.4951E-01
 GRADIENT:   1.5887E+00  1.7874E-01 -2.7616E-01 -1.1228E-01  2.5163E+00  4.2846E-02 -8.4083E-02  0.0000E+00 -1.3376E-01  5.8841E-03
            -1.0186E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1705
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.8248E-04  1.1959E-02 -2.9833E-04 -1.0344E-02  1.0480E-02
 SE:             2.9418E-02  2.3465E-02  2.1084E-04  2.6634E-02  1.0262E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9505E-01  6.1029E-01  1.5708E-01  6.9772E-01  3.0717E-01

 ETASHRINKSD(%)  1.4463E+00  2.1389E+01  9.9294E+01  1.0774E+01  6.5619E+01
 ETASHRINKVR(%)  2.8718E+00  3.8204E+01  9.9995E+01  2.0386E+01  8.8180E+01
 EBVSHRINKSD(%)  1.6599E+00  2.0192E+01  9.9358E+01  1.0480E+01  6.6722E+01
 EBVSHRINKVR(%)  3.2923E+00  3.6306E+01  9.9996E+01  1.9862E+01  8.8926E+01
 RELATIVEINF(%)  9.4171E+01  4.5606E+00  4.5690E-04  3.2230E+01  4.1385E-01
 EPSSHRINKSD(%)  3.7150E+01
 EPSSHRINKVR(%)  6.0499E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1473.9762139614304     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -738.82538739769223     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1473.976       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  6.16E-01  1.99E-01  9.95E-01  3.21E-01  9.48E-01  1.19E+00  1.00E-02  9.26E-01  2.74E-01  2.12E+00
 


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
+        1.07E+03
 
 TH 2
+       -3.12E+01  1.33E+03
 
 TH 3
+       -2.55E+02  2.47E+03  1.24E+04
 
 TH 4
+       -3.54E+01  1.37E+02 -1.81E+03  1.05E+03
 
 TH 5
+        2.08E+02 -4.15E+03 -1.23E+04  7.65E+02  1.64E+04
 
 TH 6
+       -1.27E+00 -3.06E+00  9.09E+00 -8.64E+00  3.13E+01  2.07E+02
 
 TH 7
+       -5.78E-01  4.73E+01 -1.91E+02 -1.82E+00  9.52E+01  1.51E-01  5.72E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.81E+00 -1.73E+01  1.02E+02 -8.47E+00 -7.61E+00 -3.73E+00  6.42E+00  0.00E+00  1.47E+02
 
 TH10
+        7.43E-01 -1.77E+01 -3.14E+02 -1.02E+01  3.32E+02  1.78E+00  3.17E+01  0.00E+00  8.15E-01  5.21E+01
 
 TH11
+       -1.25E+01 -1.29E+01 -7.51E+01 -9.38E+00  5.81E+01  4.23E+00  9.29E+00  0.00E+00  1.37E+01  1.53E+01  5.51E+01
 
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
 #CPUT: Total CPU Time in Seconds,       26.905
Stop Time:
Wed Sep 29 12:46:09 CDT 2021
