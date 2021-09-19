Sat Sep 18 14:27:59 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat9.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1633.58713144690        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8272E+02 -5.6289E+01 -2.7956E+01 -3.9797E+01  6.7085E+01 -1.7952E+00  1.0863E+00  2.5987E+00  5.5371E+00 -5.1481E+00
             5.1640E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1638.80522089295        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9853E-01  9.6449E-01  9.6658E-01  1.0158E+00  9.4449E-01  9.9623E-01  9.9215E-01  9.9212E-01  9.8093E-01  1.0159E+00
             8.6380E-01
 PARAMETER:  9.8526E-02  6.3843E-02  6.6006E-02  1.1564E-01  4.2895E-02  9.6223E-02  9.2120E-02  9.2093E-02  8.0748E-02  1.1580E-01
            -4.6416E-02
 GRADIENT:   1.9452E+02 -4.7916E+01 -1.7901E+01 -4.9310E+01  3.3498E+01 -2.8805E+00 -3.3386E-01  4.5986E+00  1.6812E+00  3.8692E+00
             1.7849E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1640.85759608104        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  9.9399E-01  7.6484E-01  9.7997E-01  1.1869E+00  8.3976E-01  1.0459E+00  1.0150E+00  8.7826E-01  9.2028E-01  9.6632E-01
             8.6396E-01
 PARAMETER:  9.3971E-02 -1.6809E-01  7.9772E-02  2.7132E-01 -7.4638E-02  1.4488E-01  1.1488E-01 -2.9808E-02  1.6918E-02  6.5742E-02
            -4.6227E-02
 GRADIENT:   1.2528E+02  2.7278E+00 -7.0271E+00  1.6431E+01 -1.0671E+00  1.3754E+01 -1.9399E+00  1.7549E+00  4.2288E+00  2.3823E+00
             1.5585E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1645.98777761202        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.3753E-01  7.6431E-01  9.3561E-01  1.1781E+00  8.2198E-01  9.8246E-01  1.3085E+00  7.4941E-01  8.5241E-01  9.3596E-01
             8.5572E-01
 PARAMETER:  3.5490E-02 -1.6879E-01  3.3446E-02  2.6391E-01 -9.6039E-02  8.2301E-02  3.6886E-01 -1.8847E-01 -5.9692E-02  3.3813E-02
            -5.5818E-02
 GRADIENT:   7.5350E+00  2.1239E+00 -1.0115E+00  7.2530E+00  2.7445E+00 -2.1968E+00 -1.7104E-01 -1.1494E+00  1.8816E-01  1.3105E+00
            -1.9360E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1646.07629339428        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      621
 NPARAMETR:  9.3456E-01  7.5018E-01  9.2228E-01  1.1815E+00  8.0702E-01  9.8798E-01  1.3430E+00  7.6907E-01  8.4288E-01  9.0502E-01
             8.5945E-01
 PARAMETER:  3.2326E-02 -1.8744E-01  1.9094E-02  2.6676E-01 -1.1440E-01  8.7904E-02  3.9489E-01 -1.6257E-01 -7.0932E-02  2.0208E-04
            -5.1464E-02
 GRADIENT:   2.1361E-01  3.7515E-01  2.3327E-01 -4.7634E-02 -7.9354E-01  9.2770E-02 -1.1277E-02 -2.2814E-02 -8.3237E-02  9.3662E-02
            -1.1324E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1646.14954920635        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  9.3169E-01  6.0090E-01  1.1038E+00  1.2876E+00  8.4007E-01  9.8487E-01  1.4505E+00  9.3019E-01  8.2095E-01  9.7087E-01
             8.6000E-01
 PARAMETER:  2.9242E-02 -4.0933E-01  1.9876E-01  3.5274E-01 -7.4270E-02  8.4752E-02  4.7189E-01  2.7638E-02 -9.7289E-02  7.0441E-02
            -5.0824E-02
 GRADIENT:  -2.4833E-01  2.7475E+00  2.0676E+00  3.9664E+00 -2.5443E+00 -9.3228E-02  5.5820E-01 -1.8727E-01  3.2366E-01  4.3287E-01
            -2.6305E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1646.40026573535        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      978
 NPARAMETR:  9.2877E-01  4.0362E-01  1.2076E+00  1.4136E+00  8.2497E-01  9.8067E-01  1.6045E+00  1.0058E+00  7.8649E-01  9.9068E-01
             8.6138E-01
 PARAMETER:  2.6110E-02 -8.0728E-01  2.8866E-01  4.4612E-01 -9.2407E-02  8.0481E-02  5.7284E-01  1.0582E-01 -1.4018E-01  9.0638E-02
            -4.9224E-02
 GRADIENT:  -2.4600E-01  2.2268E+00  1.7541E+00  7.9964E+00 -2.2789E+00 -6.2054E-01  1.7619E-02 -1.0353E+00 -2.0398E-01 -3.9428E-01
             5.5314E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1646.61464569049        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1153
 NPARAMETR:  9.2573E-01  2.3725E-01  1.3220E+00  1.5179E+00  8.2431E-01  9.7888E-01  1.5997E+00  1.1382E+00  7.6014E-01  1.0055E+00
             8.6234E-01
 PARAMETER:  2.2827E-02 -1.3387E+00  3.7914E-01  5.1730E-01 -9.3206E-02  7.8658E-02  5.6979E-01  2.2948E-01 -1.7425E-01  1.0552E-01
            -4.8101E-02
 GRADIENT:  -8.8620E-01  7.4744E-01  8.0987E-01  5.9864E+00 -1.1322E+00 -3.3823E-01  1.0474E-01 -6.0270E-01  4.3287E-01 -5.2920E-01
             5.3194E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1646.68361302352        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1328
 NPARAMETR:  9.2446E-01  1.4600E-01  1.3992E+00  1.5762E+00  8.2765E-01  9.7881E-01  1.4471E+00  1.2329E+00  7.3960E-01  1.0149E+00
             8.6072E-01
 PARAMETER:  2.1457E-02 -1.8241E+00  4.3593E-01  5.5503E-01 -8.9168E-02  7.8579E-02  4.6958E-01  3.0935E-01 -2.0164E-01  1.1474E-01
            -4.9983E-02
 GRADIENT:  -1.2401E-01  4.5033E-01  5.6259E-01  5.5232E+00 -1.6282E+00  1.3243E-01  1.2190E-01  7.3732E-02 -1.0853E-01  2.5746E-01
            -1.2598E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1646.71219393756        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1503
 NPARAMETR:  9.2353E-01  8.7740E-02  1.4316E+00  1.6114E+00  8.2416E-01  9.7764E-01  1.1658E+00  1.2720E+00  7.2533E-01  1.0116E+00
             8.6079E-01
 PARAMETER:  2.0450E-02 -2.3334E+00  4.5881E-01  5.7711E-01 -9.3390E-02  7.7386E-02  2.5339E-01  3.4063E-01 -2.2112E-01  1.1153E-01
            -4.9899E-02
 GRADIENT:   4.3592E-02  8.0856E-02  2.6464E-01  2.2182E+00 -6.7346E-01 -1.7510E-02  4.8526E-02  3.3295E-02 -8.1747E-02  6.2905E-02
            -4.8781E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1646.71398857929        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1664
 NPARAMETR:  9.2327E-01  8.3765E-02  1.4303E+00  1.6131E+00  8.2272E-01  9.7744E-01  1.1077E+00  1.2701E+00  7.2412E-01  1.0103E+00
             8.6080E-01
 PARAMETER:  2.0163E-02 -2.3797E+00  4.5792E-01  5.7813E-01 -9.5137E-02  7.7179E-02  2.0226E-01  3.3911E-01 -2.2279E-01  1.1025E-01
            -4.9894E-02
 GRADIENT:   4.7390E+01  7.8872E-01  1.5638E+00  1.1201E+02  3.6141E-01  3.8541E+00  4.7394E-02  6.9367E-02  2.0097E+00  7.5061E-02
             2.8933E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1646.75456484632        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1848
 NPARAMETR:  9.2402E-01  1.2020E-01  1.4204E+00  1.5903E+00  8.3023E-01  9.7751E-01  4.3762E-01  1.2641E+00  7.3502E-01  1.0187E+00
             8.6164E-01
 PARAMETER:  2.0984E-02 -2.0186E+00  4.5094E-01  5.6393E-01 -8.6053E-02  7.7250E-02 -7.2641E-01  3.3439E-01 -2.0786E-01  1.1857E-01
            -4.8916E-02
 GRADIENT:  -1.2243E-01 -7.6716E-02 -1.3932E+00 -8.1714E-02  1.0302E+00 -2.3833E-01  1.6400E-02  7.3062E-01 -4.9628E-01  4.6188E-01
             4.3109E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1646.79656498192        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2026
 NPARAMETR:  9.2482E-01  1.7177E-01  1.4089E+00  1.5592E+00  8.3808E-01  9.7891E-01  1.0134E-01  1.2314E+00  7.5218E-01  1.0233E+00
             8.6031E-01
 PARAMETER:  2.1847E-02 -1.6616E+00  4.4280E-01  5.4418E-01 -7.6641E-02  7.8681E-02 -2.1893E+00  3.0812E-01 -1.8478E-01  1.2305E-01
            -5.0467E-02
 GRADIENT:  -2.7250E-01  1.3959E-01  1.3668E+00  2.2087E+00 -1.1618E+00  4.2648E-02  1.8826E-03 -4.7798E-01 -6.8582E-02 -2.7222E-01
            -4.9183E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1646.83281145982        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2204
 NPARAMETR:  9.2616E-01  2.4131E-01  1.3695E+00  1.5132E+00  8.4509E-01  9.7957E-01  1.1008E-02  1.1930E+00  7.7577E-01  1.0312E+00
             8.6187E-01
 PARAMETER:  2.3287E-02 -1.3217E+00  4.1442E-01  5.1424E-01 -6.8312E-02  7.9360E-02 -4.4091E+00  2.7644E-01 -1.5390E-01  1.3069E-01
            -4.8656E-02
 GRADIENT:   1.9454E-01 -3.3663E-01 -9.5779E-01 -1.7370E+00  1.0550E+00 -4.2340E-02  4.4982E-05  3.2093E-01  1.2428E-01  1.0173E-01
             3.6479E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1646.84669836005        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2382
 NPARAMETR:  9.2703E-01  2.9989E-01  1.3435E+00  1.4771E+00  8.5110E-01  9.8024E-01  1.0000E-02  1.1565E+00  7.9595E-01  1.0370E+00
             8.6131E-01
 PARAMETER:  2.4231E-02 -1.1043E+00  3.9529E-01  4.9007E-01 -6.1221E-02  8.0047E-02 -6.2752E+00  2.4539E-01 -1.2822E-01  1.3633E-01
            -4.9300E-02
 GRADIENT:   4.4214E-02  1.4376E-01  2.5928E-02  9.1721E-01 -7.0765E-01 -8.1537E-02  0.0000E+00 -3.5273E-02 -6.9776E-02  5.5449E-02
            -2.1273E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1646.85031356461        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2564             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2738E-01  3.2584E-01  1.3463E+00  1.4611E+00  8.5981E-01  9.8085E-01  1.0000E-02  1.1570E+00  8.0531E-01  1.0426E+00
             8.6149E-01
 PARAMETER:  2.4613E-02 -1.0214E+00  3.9739E-01  4.7917E-01 -5.1045E-02  8.0661E-02 -7.1733E+00  2.4582E-01 -1.1653E-01  1.4170E-01
            -4.9097E-02
 GRADIENT:   4.8063E+01  4.6061E+00  1.0639E+00  7.3294E+01  9.4605E-01  4.0022E+00  0.0000E+00  8.6706E-02  1.0950E+00  1.3135E-01
             7.1302E-02

0ITERATION NO.:   78    OBJECTIVE VALUE:  -1646.85031875415        NO. OF FUNC. EVALS.:  90
 CUMULATIVE NO. OF FUNC. EVALS.:     2654
 NPARAMETR:  9.2739E-01  3.2582E-01  1.3463E+00  1.4611E+00  8.5982E-01  9.8079E-01  1.0000E-02  1.1571E+00  8.0536E-01  1.0427E+00
             8.6149E-01
 PARAMETER:  2.4619E-02 -1.0214E+00  3.9739E-01  4.7917E-01 -5.1033E-02  8.0600E-02 -7.1733E+00  2.4591E-01 -1.1647E-01  1.4184E-01
            -4.9088E-02
 GRADIENT:   2.4939E-02  9.0753E-03 -1.9070E-02  3.0270E-02 -2.8540E-02  6.2448E-03  0.0000E+00 -2.9974E-03  3.5143E-03  1.2331E-03
             2.2754E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2654
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.9599E-04 -1.2849E-04 -3.4521E-02 -5.0539E-03 -3.5081E-02
 SE:             2.9869E-02  6.0881E-05  1.7508E-02  2.9393E-02  2.2096E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8942E-01  3.4815E-02  4.8640E-02  8.6348E-01  1.1236E-01

 ETASHRINKSD(%)  1.0000E-10  9.9796E+01  4.1347E+01  1.5292E+00  2.5975E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.5598E+01  3.0351E+00  4.5203E+01
 EBVSHRINKSD(%)  3.2593E-01  9.9807E+01  4.5548E+01  1.9742E+00  2.2241E+01
 EBVSHRINKVR(%)  6.5079E-01  1.0000E+02  7.0350E+01  3.9094E+00  3.9536E+01
 RELATIVEINF(%)  9.8033E+01  2.2731E-05  6.7849E+00  7.2839E+00  8.9444E+00
 EPSSHRINKSD(%)  4.6357E+01
 EPSSHRINKVR(%)  7.1224E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1646.8503187541485     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -911.69949219041030     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.70
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1646.850       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.27E-01  3.26E-01  1.35E+00  1.46E+00  8.60E-01  9.81E-01  1.00E-02  1.16E+00  8.05E-01  1.04E+00  8.61E-01
 


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
+        1.34E+03
 
 TH 2
+       -2.66E+01  4.08E+02
 
 TH 3
+        3.11E+00  6.46E+01  1.45E+02
 
 TH 4
+       -8.19E+00  4.85E+02 -2.37E+01  7.72E+02
 
 TH 5
+        4.14E+00 -2.57E+02 -2.81E+02 -4.94E+01  8.58E+02
 
 TH 6
+        1.54E+00 -2.71E+00  2.53E-01 -2.76E+00 -2.22E+00  2.08E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -3.36E-02 -3.52E+00 -3.40E+01 -4.00E+00  1.05E+00  6.76E-01  0.00E+00  3.15E+01
 
 TH 9
+        5.39E+00 -1.04E+02  4.40E+00 -9.06E-01 -4.70E-01  4.74E-01  0.00E+00 -5.68E-01  2.85E+02
 
 TH10
+       -2.43E-01  7.90E+00 -6.31E+00 -2.12E-01 -7.73E+01  1.02E+00  0.00E+00  2.02E+01  3.85E-02  7.04E+01
 
 TH11
+       -6.94E+00 -1.16E+01 -1.00E+01 -7.01E+00 -5.21E+00  4.71E+00  0.00E+00  7.75E+00  1.19E+01  1.26E+01  2.77E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.119
Stop Time:
Sat Sep 18 14:28:38 CDT 2021
