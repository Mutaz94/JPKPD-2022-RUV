Sat Sep 18 10:18:11 CDT 2021
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
$DATA ../../../../data/spa/A3/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   33.1982556083512        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5624E+01  1.8881E+01  8.9932E+01 -1.1278E+02  1.2443E+02  1.7562E+01 -6.5640E+01 -2.9584E+01 -1.7026E+02 -1.2825E+02
            -3.0215E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1224.32998384058        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1055E+00  9.8625E-01  9.3095E-01  1.2219E+00  9.8384E-01  7.9898E-01  1.0390E+00  9.7236E-01  1.2131E+00  1.0025E+00
             5.2653E+00
 PARAMETER:  2.0026E-01  8.6154E-02  2.8447E-02  3.0040E-01  8.3712E-02 -1.2441E-01  1.3827E-01  7.1968E-02  2.9322E-01  1.0253E-01
             1.7611E+00
 GRADIENT:   6.7815E+01 -4.8201E+00 -2.2990E+01  2.3351E+01  2.1951E+00 -2.2968E+01  8.6050E+00  7.2206E+00  2.6118E+01  1.8745E+01
             1.7493E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1256.88883703856        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0647E+00  7.8184E-01  4.9422E-01  1.2467E+00  5.6556E-01  9.1004E-01  1.2188E+00  6.6546E-02  1.0858E+00  4.3770E-01
             4.4459E+00
 PARAMETER:  1.6272E-01 -1.4611E-01 -6.0477E-01  3.2049E-01 -4.6994E-01  5.7329E-03  2.9790E-01 -2.6099E+00  1.8234E-01 -7.2622E-01
             1.5920E+00
 GRADIENT:  -5.8561E+00  1.7395E+01 -2.7286E+01  4.8359E+01  1.2861E+01  5.0912E+00  9.6334E+00  7.2472E-02  1.0755E+01  7.8459E+00
             1.0638E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1272.16155900428        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0446E+00  6.4413E-01  7.2152E-01  1.2910E+00  6.6910E-01  8.9963E-01  1.0201E+00  1.5921E-01  1.0629E+00  3.0726E-01
             3.8105E+00
 PARAMETER:  1.4365E-01 -3.3985E-01 -2.2639E-01  3.5539E-01 -3.0182E-01 -5.7688E-03  1.1988E-01 -1.7375E+00  1.6099E-01 -1.0801E+00
             1.4378E+00
 GRADIENT:   2.9550E+00  6.3863E+00  7.0718E+00 -2.4244E+00 -1.3069E+01  1.5311E+00 -1.7200E-01  2.9799E-01  4.0622E-01  1.8227E+00
            -1.3247E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1275.03377023118        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0386E+00  3.2971E-01  6.8129E-01  1.4465E+00  5.7273E-01  8.8282E-01  1.6232E+00  1.7573E-02  9.3319E-01  3.4489E-02
             3.9218E+00
 PARAMETER:  1.3789E-01 -1.0096E+00 -2.8377E-01  4.6911E-01 -4.5734E-01 -2.4630E-02  5.8439E-01 -3.9414E+00  3.0850E-02 -3.2671E+00
             1.4665E+00
 GRADIENT:   1.1524E+00  4.8349E-02 -2.9946E+00 -5.1466E+00  5.8019E+00 -2.5129E+00 -5.3723E-01  5.0010E-03 -1.9165E+00  2.0921E-02
             4.0262E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1276.39639233286        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0337E+00  1.5649E-01  3.9883E-01  1.4302E+00  3.5701E-01  9.2534E-01  2.8120E+00  1.0000E-02  1.0339E+00  1.0000E-02
             3.7933E+00
 PARAMETER:  1.3311E-01 -1.7548E+00 -8.1923E-01  4.5784E-01 -9.3000E-01  2.2404E-02  1.1339E+00 -7.9059E+00  1.3333E-01 -5.8401E+00
             1.4332E+00
 GRADIENT:  -7.9477E+00  4.8486E+00  1.0168E+01  3.2510E+01 -1.8477E+01  4.4985E+00  4.0099E-01  0.0000E+00  3.4459E+00  0.0000E+00
            -5.0930E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1277.44947489538        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.0332E+00  8.8497E-02  3.2232E-01  1.3441E+00  2.9895E-01  9.0912E-01  4.2050E+00  1.0000E-02  1.0465E+00  1.0000E-02
             3.7948E+00
 PARAMETER:  1.3265E-01 -2.3248E+00 -1.0322E+00  3.9570E-01 -1.1075E+00  4.7212E-03  1.5363E+00 -1.1446E+01  1.4547E-01 -8.2963E+00
             1.4336E+00
 GRADIENT:  -1.3252E+00 -1.2654E+00  3.0418E+01 -6.7439E+00 -3.9094E+01 -1.1919E+00 -5.8339E+00  0.0000E+00 -4.8876E-01  0.0000E+00
             7.6391E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1278.05316426597        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      516
 NPARAMETR:  1.0324E+00  6.2616E-02  3.0124E-01  1.3310E+00  2.8581E-01  9.0637E-01  5.3566E+00  1.0000E-02  1.0497E+00  1.0000E-02
             3.7892E+00
 PARAMETER:  1.3185E-01 -2.6707E+00 -1.0999E+00  3.8596E-01 -1.1524E+00  1.6888E-03  1.7783E+00 -1.3594E+01  1.4852E-01 -9.8104E+00
             1.4322E+00
 GRADIENT:   2.3425E+00 -4.0274E-01  1.2764E+01 -6.6716E+00 -1.9038E+01 -3.3140E+00 -2.7040E+00  0.0000E+00 -2.7796E+00  0.0000E+00
             6.9561E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1278.28727952728        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      592
 NPARAMETR:  1.0305E+00  5.2643E-02  2.9704E-01  1.3330E+00  2.8410E-01  9.1143E-01  5.9692E+00  1.0000E-02  1.0557E+00  1.0000E-02
             3.7715E+00
 PARAMETER:  1.3007E-01 -2.8442E+00 -1.1139E+00  3.8742E-01 -1.1584E+00  7.2614E-03  1.8866E+00 -1.4653E+01  1.5422E-01 -1.0555E+01
             1.4275E+00
 GRADIENT:   8.7223E+00  1.2051E+01  3.1981E+00 -1.9100E+01 -4.7917E+00 -5.3236E+00  1.9075E+01  0.0000E+00  1.0317E+00  0.0000E+00
            -1.4365E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1278.46975020435        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  1.0290E+00  4.9066E-02  3.0172E-01  1.3447E+00  2.8783E-01  9.1827E-01  6.2418E+00  1.0000E-02  1.0579E+00  1.0000E-02
             3.7566E+00
 PARAMETER:  1.2855E-01 -2.9146E+00 -1.0983E+00  3.9620E-01 -1.1454E+00  1.4732E-02  1.9313E+00 -1.5040E+01  1.5625E-01 -1.0835E+01
             1.4235E+00
 GRADIENT:   5.4045E+00  9.1741E+00 -2.2155E+00 -1.1478E+01  3.6154E+00 -2.2337E+00  1.5120E+01  0.0000E+00  2.2162E+00  0.0000E+00
            -1.3356E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1278.77899594919        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      803
 NPARAMETR:  1.0292E+00  4.7203E-02  3.3843E-01  1.3924E+00  3.1175E-01  9.1626E-01  6.5411E+00  1.0000E-02  1.0235E+00  1.0000E-02
             3.7666E+00
 PARAMETER:  1.2878E-01 -2.9533E+00 -9.8342E-01  4.3105E-01 -1.0656E+00  1.2544E-02  1.9781E+00 -1.5036E+01  1.2318E-01 -1.0937E+01
             1.4262E+00
 GRADIENT:  -1.2337E+00  1.1823E+00 -3.3704E+00  2.8114E+00  3.3291E+00  7.8612E-01  1.3276E+00  0.0000E+00  5.2511E-02  0.0000E+00
            -2.3455E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1278.79830403351        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      979
 NPARAMETR:  1.0298E+00  4.5964E-02  3.4169E-01  1.3930E+00  3.1319E-01  9.1223E-01  6.6463E+00  1.0000E-02  1.0200E+00  1.0000E-02
             3.7735E+00
 PARAMETER:  1.2939E-01 -2.9799E+00 -9.7385E-01  4.3143E-01 -1.0609E+00  8.1415E-03  1.9941E+00 -1.5166E+01  1.1979E-01 -1.1046E+01
             1.4280E+00
 GRADIENT:   3.9668E+00  6.8870E+00  1.8623E+00 -8.9358E+00 -1.6908E+00 -2.0371E+00  1.2086E+01  0.0000E+00  1.7003E+00  0.0000E+00
            -1.0885E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1278.85505919551        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1161
 NPARAMETR:  1.0290E+00  3.8675E-02  3.3606E-01  1.3878E+00  3.0877E-01  9.1431E-01  7.2476E+00  1.0000E-02  1.0290E+00  1.0000E-02
             3.7681E+00
 PARAMETER:  1.2855E-01 -3.1526E+00 -9.9046E-01  4.2774E-01 -1.0751E+00  1.0418E-02  2.0807E+00 -1.6177E+01  1.2861E-01 -1.1770E+01
             1.4266E+00
 GRADIENT:   7.0827E+00  1.1897E+01  5.5269E+00 -1.8701E+01 -4.8454E+00 -3.1676E+00  2.2435E+01  0.0000E+00  4.2730E+00  0.0000E+00
            -2.0136E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1278.89296470399        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1339
 NPARAMETR:  1.0285E+00  3.2671E-02  3.3131E-01  1.3855E+00  3.0585E-01  9.1441E-01  7.8951E+00  1.0000E-02  1.0295E+00  1.0000E-02
             3.7643E+00
 PARAMETER:  1.2810E-01 -3.3213E+00 -1.0047E+00  4.2607E-01 -1.0847E+00  1.0520E-02  2.1662E+00 -1.7173E+01  1.2910E-01 -1.2492E+01
             1.4256E+00
 GRADIENT:   1.6144E+01  2.4117E+01  6.8558E+00 -3.7511E+01 -3.6399E+00 -7.5014E+00  4.7102E+01  0.0000E+00  7.2627E+00  0.0000E+00
            -4.2161E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1278.92052916557        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1515
 NPARAMETR:  1.0286E+00  3.0290E-02  3.4183E-01  1.3980E+00  3.1246E-01  9.1117E-01  8.2480E+00  1.0000E-02  1.0177E+00  1.0000E-02
             3.7669E+00
 PARAMETER:  1.2818E-01 -3.3969E+00 -9.7344E-01  4.3502E-01 -1.0633E+00  6.9691E-03  2.2100E+00 -1.7546E+01  1.1750E-01 -1.2798E+01
             1.4263E+00
 GRADIENT:   4.0476E+00  5.6962E+00  3.2561E-01 -7.6656E+00  8.5219E-01 -2.1453E+00  1.0916E+01  0.0000E+00  1.1477E+00  0.0000E+00
            -1.0234E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1278.93561495193        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1692
 NPARAMETR:  1.0281E+00  2.6984E-02  3.4574E-01  1.4027E+00  3.1455E-01  9.1211E-01  8.7612E+00  1.0000E-02  1.0176E+00  1.0000E-02
             3.7681E+00
 PARAMETER:  1.2776E-01 -3.5125E+00 -9.6208E-01  4.3843E-01 -1.0566E+00  8.0035E-03  2.2703E+00 -1.8185E+01  1.1744E-01 -1.3277E+01
             1.4266E+00
 GRADIENT:   2.4979E+01  3.5864E+01  1.3363E+01 -5.6641E+01 -8.1861E+00 -1.1637E+01  7.1017E+01  0.0000E+00  1.0780E+01  0.0000E+00
            -6.2920E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1278.95001935138        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1871
 NPARAMETR:  1.0282E+00  2.5767E-02  3.3679E-01  1.3932E+00  3.0855E-01  9.1234E-01  8.9247E+00  1.0000E-02  1.0244E+00  1.0000E-02
             3.7649E+00
 PARAMETER:  1.2779E-01 -3.5587E+00 -9.8829E-01  4.3160E-01 -1.0759E+00  8.2586E-03  2.2888E+00 -1.8492E+01  1.2409E-01 -1.3480E+01
             1.4257E+00
 GRADIENT:   2.6550E+00  3.9905E+00  1.5755E+00 -5.7356E+00 -1.5877E+00 -1.3418E+00  7.6508E+00  0.0000E+00  1.1961E+00  0.0000E+00
            -6.9679E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1278.97438866477        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2048            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0273E+00  1.8719E-02  3.3799E-01  1.3963E+00  3.0914E-01  9.1261E-01  1.0467E+01  1.0000E-02  1.0220E+00  1.0000E-02
             3.7658E+00
 PARAMETER:  1.2695E-01 -3.8782E+00 -9.8475E-01  4.3381E-01 -1.0740E+00  8.5570E-03  2.4482E+00 -2.0309E+01  1.2179E-01 -1.4821E+01
             1.4260E+00
 GRADIENT:   2.2282E+00 -1.4163E-01 -1.0721E+00  4.1713E+00  5.1747E+00  4.5613E-01  1.5060E-01  0.0000E+00  1.1317E-01  0.0000E+00
             2.0742E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1278.97587357222        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2225
 NPARAMETR:  1.0277E+00  1.8742E-02  3.3833E-01  1.3961E+00  3.0894E-01  9.1182E-01  1.0471E+01  1.0000E-02  1.0222E+00  1.0000E-02
             3.7639E+00
 PARAMETER:  1.2728E-01 -3.8770E+00 -9.8375E-01  4.3370E-01 -1.0746E+00  7.6869E-03  2.4486E+00 -2.0309E+01  1.2193E-01 -1.4821E+01
             1.4255E+00
 GRADIENT:   1.9663E-01  2.2612E-01  7.3815E-01 -4.2306E-01 -1.1801E+00 -1.0848E-01  4.3795E-01  0.0000E+00  4.5178E-02  0.0000E+00
            -4.1098E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1278.97623618741        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     2388
 NPARAMETR:  1.0276E+00  1.8729E-02  3.3830E-01  1.3962E+00  3.0904E-01  9.1189E-01  1.0466E+01  1.0000E-02  1.0222E+00  1.0000E-02
             3.7638E+00
 PARAMETER:  1.2727E-01 -3.8777E+00 -9.8381E-01  4.3376E-01 -1.0743E+00  7.7677E-03  2.4482E+00 -2.0309E+01  1.2194E-01 -1.4821E+01
             1.4254E+00
 GRADIENT:   5.4413E-02 -3.3984E-02 -2.4931E-02 -1.2918E-02 -1.5762E-01 -4.3753E-03 -7.0913E-02  0.0000E+00  2.2363E-03  0.0000E+00
             4.3463E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2388
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.3813E-04  5.2034E-03  1.3559E-04 -1.5842E-02  9.3346E-05
 SE:             2.8192E-02  4.8193E-03  2.4539E-04  2.5930E-02  3.7235E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7345E-01  2.8027E-01  5.8057E-01  5.4124E-01  8.0205E-01

 ETASHRINKSD(%)  5.5521E+00  8.3855E+01  9.9178E+01  1.3132E+01  9.8753E+01
 ETASHRINKVR(%)  1.0796E+01  9.7393E+01  9.9993E+01  2.4539E+01  9.9984E+01
 EBVSHRINKSD(%)  5.1157E+00  8.6911E+01  9.9184E+01  1.2451E+01  9.8870E+01
 EBVSHRINKVR(%)  9.9697E+00  9.8287E+01  9.9993E+01  2.3351E+01  9.9987E+01
 RELATIVEINF(%)  8.6802E+01  1.0846E+00  2.2910E-04  2.4083E+01  4.3234E-04
 EPSSHRINKSD(%)  2.1714E+01
 EPSSHRINKVR(%)  3.8713E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1278.9762361874061     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -543.82540962366795     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.04
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1278.976       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.87E-02  3.38E-01  1.40E+00  3.09E-01  9.12E-01  1.05E+01  1.00E-02  1.02E+00  1.00E-02  3.76E+00
 


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
+        1.17E+03
 
 TH 2
+        3.58E+03  2.54E+07
 
 TH 3
+       -9.17E+01  7.48E+03  5.91E+03
 
 TH 4
+       -1.02E+02 -6.99E+03 -2.65E+02  4.85E+02
 
 TH 5
+        3.13E+02 -6.16E+03 -8.66E+03 -3.33E+02  1.44E+04
 
 TH 6
+       -1.37E+01 -2.01E+03  2.63E+01  4.65E+00  5.93E+00  2.14E+02
 
 TH 7
+       -2.26E+01  7.19E+04  1.58E+04 -8.61E+03 -1.58E+04  1.03E+01  2.04E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.13E+01  1.37E+03  4.05E+01 -3.05E+01  1.21E+02  2.19E+00 -2.47E+01  0.00E+00  1.19E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.01E+01 -3.41E+05 -5.81E+01  2.66E+01  6.67E+01  1.51E+01 -9.67E+02  0.00E+00 -1.87E-01  0.00E+00  4.62E+03
 
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
 #CPUT: Total CPU Time in Seconds,       39.154
Stop Time:
Sat Sep 18 10:18:51 CDT 2021
