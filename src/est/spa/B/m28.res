Sat Sep 25 07:13:20 CDT 2021
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
$DATA ../../../../data/spa/B/dat28.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1622.17897332628        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.3347E+01 -1.1450E+02 -2.8008E+01 -1.5705E+02  6.0225E+01  7.4450E+00  2.5169E+00  7.8722E+00 -2.6285E+01  1.0242E+01
            -3.9957E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1636.65829057828        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9207E-01  9.5140E-01  1.0414E+00  1.1409E+00  9.4276E-01  9.7461E-01  8.2785E-01  8.9974E-01  1.0863E+00  7.9872E-01
             1.0812E+00
 PARAMETER:  9.2038E-02  5.0176E-02  1.4053E-01  2.3181E-01  4.1058E-02  7.4286E-02 -8.8925E-02 -5.6448E-03  1.8277E-01 -1.2474E-01
             1.7808E-01
 GRADIENT:   4.9492E+01 -5.7014E+00 -9.4687E-01  1.2878E+01  3.9669E+01 -2.4856E+00  2.2273E+00 -2.9734E+00  1.0235E+01 -1.7540E+01
            -1.0955E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1638.13863064598        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.8204E-01  7.9420E-01  1.0426E+00  1.2510E+00  8.9096E-01  9.8726E-01  5.9804E-01  7.0018E-01  1.0297E+00  9.4173E-01
             1.0887E+00
 PARAMETER:  8.1874E-02 -1.3041E-01  1.4173E-01  3.2397E-01 -1.5460E-02  8.7178E-02 -4.1410E-01 -2.5641E-01  1.2928E-01  3.9966E-02
             1.8500E-01
 GRADIENT:   2.7121E+01  5.5125E+00 -4.9334E+00  4.2865E+01  2.4877E+01  3.0389E+00 -2.0697E-01 -4.6635E+00  6.7417E+00  7.4882E-01
            -4.2826E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1639.66549015604        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      339
 NPARAMETR:  9.8459E-01  7.0620E-01  9.2350E-01  1.2802E+00  7.7414E-01  9.8784E-01  9.2894E-01  7.0830E-01  9.3631E-01  7.9985E-01
             1.0976E+00
 PARAMETER:  8.4466E-02 -2.4786E-01  2.0420E-02  3.4699E-01 -1.5600E-01  8.7761E-02  2.6294E-02 -2.4489E-01  3.4194E-02 -1.2333E-01
             1.9315E-01
 GRADIENT:   1.1990E+00  3.0796E+00  5.0887E+00 -6.9780E+00 -1.0367E+01  5.7408E-01 -2.8999E-01  8.1636E-01 -8.7189E-01  1.2162E-01
             7.1467E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1640.36800213139        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  9.8108E-01  4.6137E-01  9.8429E-01  1.4325E+00  7.3655E-01  9.8322E-01  1.1488E+00  7.0122E-01  8.5729E-01  8.2659E-01
             1.0960E+00
 PARAMETER:  8.0904E-02 -6.7355E-01  8.4163E-02  4.5939E-01 -2.0578E-01  8.3080E-02  2.3871E-01 -2.5493E-01 -5.3980E-02 -9.0450E-02
             1.9171E-01
 GRADIENT:   1.3472E+00  1.3760E+00  3.8630E+00 -6.3119E-01 -2.8790E+00 -3.1187E-01  4.9632E-01 -9.3649E-01  1.2178E+00  2.5137E+00
             2.3347E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1640.73881189383        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      691
 NPARAMETR:  9.7820E-01  2.9527E-01  9.6076E-01  1.5260E+00  6.8325E-01  9.8268E-01  1.3243E+00  7.2781E-01  8.0418E-01  7.7776E-01
             1.0968E+00
 PARAMETER:  7.7963E-02 -1.1199E+00  5.9969E-02  5.2264E-01 -2.8089E-01  8.2531E-02  3.8090E-01 -2.1771E-01 -1.1794E-01 -1.5134E-01
             1.9235E-01
 GRADIENT:   1.0923E+00  8.0129E-01  2.9833E-01  4.5795E+00 -2.9114E-02  2.1674E-01 -1.8568E-01 -7.8049E-01 -1.1054E+00 -7.8134E-01
            -4.0263E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1640.82798832876        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      866
 NPARAMETR:  9.7533E-01  1.8244E-01  9.9851E-01  1.5903E+00  6.7580E-01  9.8074E-01  1.4774E+00  8.1033E-01  7.7690E-01  7.7591E-01
             1.0971E+00
 PARAMETER:  7.5019E-02 -1.6013E+00  9.8511E-02  5.6389E-01 -2.9186E-01  8.0549E-02  4.9028E-01 -1.1032E-01 -1.5245E-01 -1.5372E-01
             1.9265E-01
 GRADIENT:   6.7500E-01 -7.4242E-02 -1.7374E-01 -1.8917E+00  4.1239E-01  1.3906E-01 -5.6857E-02  2.3705E-01  2.7465E-01  3.7289E-01
             5.4432E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1640.86524900679        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1041
 NPARAMETR:  9.7316E-01  1.0583E-01  9.9397E-01  1.6338E+00  6.5642E-01  9.7932E-01  2.2318E+00  8.2816E-01  7.5421E-01  7.6072E-01
             1.0980E+00
 PARAMETER:  7.2791E-02 -2.1460E+00  9.3948E-02  5.9094E-01 -3.2095E-01  7.9106E-02  9.0279E-01 -8.8551E-02 -1.8208E-01 -1.7349E-01
             1.9347E-01
 GRADIENT:  -3.5281E-01  1.8934E-01  7.0310E-01  2.4803E+00 -1.8258E+00 -3.7898E-02 -1.8697E-03  7.8423E-02 -1.7483E-01  1.2644E-01
             4.2094E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1640.87125249914        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1216
 NPARAMETR:  9.7253E-01  7.5555E-02  1.0075E+00  1.6520E+00  6.5660E-01  9.7887E-01  2.9719E+00  8.4933E-01  7.4574E-01  7.5880E-01
             1.0980E+00
 PARAMETER:  7.2144E-02 -2.4829E+00  1.0746E-01  6.0198E-01 -3.2068E-01  7.8639E-02  1.1892E+00 -6.3303E-02 -1.9337E-01 -1.7601E-01
             1.9349E-01
 GRADIENT:   6.5376E-03  1.0238E-01  2.6763E-01  1.4533E+00 -2.9566E-01 -3.9623E-02  5.8083E-02 -2.7070E-03 -5.6760E-02 -1.4644E-02
            -7.4558E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1640.88275498279        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1397
 NPARAMETR:  9.7154E-01  3.8688E-02  1.0112E+00  1.6726E+00  6.5056E-01  9.7841E-01  3.6285E+00  8.6369E-01  7.3751E-01  7.5688E-01
             1.0982E+00
 PARAMETER:  7.1131E-02 -3.1522E+00  1.1114E-01  6.1437E-01 -3.2993E-01  7.8170E-02  1.3888E+00 -4.6538E-02 -2.0448E-01 -1.7855E-01
             1.9369E-01
 GRADIENT:  -1.5780E-01  4.6679E-02  4.5123E-01  1.6494E+00 -1.1994E+00 -9.2136E-03  4.0172E-03  6.8365E-02 -6.4059E-02  1.8780E-01
             4.0542E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1640.88522833603        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1572
 NPARAMETR:  9.7125E-01  2.4324E-02  1.0198E+00  1.6813E+00  6.5169E-01  9.7821E-01  3.5210E+00  8.7589E-01  7.3413E-01  7.5607E-01
             1.0982E+00
 PARAMETER:  7.0828E-02 -3.6163E+00  1.1960E-01  6.1954E-01 -3.2819E-01  7.7966E-02  1.3587E+00 -3.2515E-02 -2.0907E-01 -1.7962E-01
             1.9364E-01
 GRADIENT:   3.4224E-02  1.2017E-02  1.6552E-01  7.5312E-01 -2.0601E-01  4.4689E-03 -6.4706E-03 -1.8913E-02 -9.6180E-02 -4.5482E-02
            -3.1491E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1640.88671261592        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1752
 NPARAMETR:  9.7084E-01  1.0000E-02  1.0197E+00  1.6890E+00  6.4873E-01  9.7797E-01  3.1754E+00  8.8015E-01  7.3079E-01  7.5450E-01
             1.0984E+00
 PARAMETER:  7.0409E-02 -4.5465E+00  1.1950E-01  6.2414E-01 -3.3274E-01  7.7722E-02  1.2554E+00 -2.7668E-02 -2.1363E-01 -1.8170E-01
             1.9381E-01
 GRADIENT:  -6.7741E-02  0.0000E+00  1.3312E-01  7.0569E-01 -4.3028E-01 -8.8026E-03 -1.9689E-03  2.1187E-02 -4.4913E-02  2.3272E-02
             2.0372E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1640.88690628727        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1929
 NPARAMETR:  9.7089E-01  1.0000E-02  1.0213E+00  1.6888E+00  6.4955E-01  9.7800E-01  3.2211E+00  8.8136E-01  7.3084E-01  7.5496E-01
             1.0983E+00
 PARAMETER:  7.0463E-02 -4.5130E+00  1.2110E-01  6.2404E-01 -3.3147E-01  7.7755E-02  1.2697E+00 -2.6292E-02 -2.1356E-01 -1.8109E-01
             1.9375E-01
 GRADIENT:   6.6633E-02  0.0000E+00  2.3133E-02 -2.9431E-01 -3.6336E-03  6.4544E-03 -1.9368E-03  2.2085E-03  1.0010E-02  3.0862E-03
             2.3068E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1640.88769673467        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2109
 NPARAMETR:  9.7074E-01  1.0000E-02  1.0214E+00  1.6900E+00  6.4968E-01  9.7782E-01  6.9633E+00  8.8091E-01  7.3060E-01  7.5484E-01
             1.0985E+00
 PARAMETER:  7.0308E-02 -4.5130E+00  1.2121E-01  6.2471E-01 -3.3127E-01  7.7573E-02  2.0407E+00 -2.6800E-02 -2.1389E-01 -1.8124E-01
             1.9391E-01
 GRADIENT:  -3.0896E-01  0.0000E+00 -1.9665E-01  2.1945E+00  1.2461E-01 -7.1052E-02  1.1093E-03 -3.8730E-02  5.1599E-02 -1.2309E-02
             5.1584E-02

0ITERATION NO.:   68    OBJECTIVE VALUE:  -1640.88845599669        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     2202
 NPARAMETR:  9.7092E-01  1.0000E-02  1.0215E+00  1.6888E+00  6.4966E-01  9.7804E-01  6.6228E+00  8.8139E-01  7.3049E-01  7.5487E-01
             1.0982E+00
 PARAMETER:  7.0485E-02 -4.5130E+00  1.2129E-01  6.2402E-01 -3.3131E-01  7.7797E-02  1.9905E+00 -2.6257E-02 -2.1404E-01 -1.8121E-01
             1.9372E-01
 GRADIENT:   1.1380E-01  0.0000E+00 -6.8621E-04 -5.3573E-01  1.2216E-01  2.2256E-02 -3.5463E-05 -3.7595E-03 -1.5134E-02  1.0064E-02
            -8.6123E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2202
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5656E-04 -7.9420E-04 -1.7532E-02 -4.7509E-03 -2.3421E-02
 SE:             2.9811E-02  1.0933E-03  1.7774E-02  2.9088E-02  2.0608E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9581E-01  4.6757E-01  3.2396E-01  8.7026E-01  2.5575E-01

 ETASHRINKSD(%)  1.2960E-01  9.6337E+01  4.0454E+01  2.5506E+00  3.0962E+01
 ETASHRINKVR(%)  2.5904E-01  9.9866E+01  6.4542E+01  5.0362E+00  5.2337E+01
 EBVSHRINKSD(%)  5.0222E-01  9.6429E+01  4.1769E+01  2.9015E+00  2.9819E+01
 EBVSHRINKVR(%)  1.0019E+00  9.9872E+01  6.6091E+01  5.7189E+00  5.0746E+01
 RELATIVEINF(%)  9.3238E+01  5.0802E-03  3.5585E+00  5.8859E+00  2.5862E+00
 EPSSHRINKSD(%)  4.3948E+01
 EPSSHRINKVR(%)  6.8582E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1640.8884559966934     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -905.73762943295526     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.38
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:     5.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1640.888       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.71E-01  1.00E-02  1.02E+00  1.69E+00  6.50E-01  9.78E-01  6.62E+00  8.81E-01  7.30E-01  7.55E-01  1.10E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         2.82E-02  0.00E+00  2.08E-01  4.40E-02  9.87E-02  5.85E-02  6.39E+01  1.29E-01  5.38E-02  2.27E-01  7.75E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        7.98E-04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        1.45E-04  0.00E+00  4.34E-02
 
 TH 4
+        9.19E-05  0.00E+00  4.76E-03  1.94E-03
 
 TH 5
+       -8.72E-06  0.00E+00  2.01E-02  2.27E-03  9.74E-03
 
 TH 6
+       -3.23E-04  0.00E+00  3.36E-03  3.19E-04  1.48E-03  3.42E-03
 
 TH 7
+        1.44E-01  0.00E+00  1.30E+01  1.44E+00  6.06E+00  1.33E+00  4.08E+03
 
 TH 8
+       -1.52E-04  0.00E+00 -2.29E-02 -2.73E-03 -1.06E-02 -3.23E-03 -7.55E+00  1.66E-02
 
 TH 9
+        1.04E-04  0.00E+00 -6.02E-04 -9.95E-05 -2.12E-04 -4.28E-04  1.33E-02  8.66E-04  2.89E-03
 
 TH10
+        5.14E-04  0.00E+00  4.11E-02  4.44E-03  1.94E-02  3.67E-03  1.30E+01 -2.21E-02 -1.25E-03  5.16E-02
 
 TH11
+       -2.68E-05  0.00E+00 -2.95E-03  2.16E-04 -1.36E-03 -3.55E-04 -8.06E-01  7.81E-04  6.51E-04 -3.96E-03  6.00E-03
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.82E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        2.46E-02  0.00E+00  2.08E-01
 
 TH 4
+        7.39E-02  0.00E+00  5.19E-01  4.40E-02
 
 TH 5
+       -3.13E-03  0.00E+00  9.77E-01  5.24E-01  9.87E-02
 
 TH 6
+       -1.96E-01  0.00E+00  2.76E-01  1.24E-01  2.56E-01  5.85E-02
 
 TH 7
+        8.01E-02  0.00E+00  9.76E-01  5.13E-01  9.62E-01  3.57E-01  6.39E+01
 
 TH 8
+       -4.19E-02  0.00E+00 -8.53E-01 -4.81E-01 -8.38E-01 -4.29E-01 -9.18E-01  1.29E-01
 
 TH 9
+        6.82E-02  0.00E+00 -5.37E-02 -4.20E-02 -4.00E-02 -1.36E-01  3.86E-03  1.25E-01  5.38E-02
 
 TH10
+        8.01E-02  0.00E+00  8.69E-01  4.44E-01  8.68E-01  2.76E-01  8.98E-01 -7.57E-01 -1.03E-01  2.27E-01
 
 TH11
+       -1.22E-02  0.00E+00 -1.83E-01  6.33E-02 -1.78E-01 -7.84E-02 -1.63E-01  7.83E-02  1.56E-01 -2.25E-01  7.75E-02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.29E+03
 
 TH 2
+        6.45E-15 -2.08E-30
 
 TH 3
+       -1.11E+02  2.50E-15  4.95E+02
 
 TH 4
+       -6.89E+01  8.13E-15 -2.17E+01  6.73E+02
 
 TH 5
+        3.05E+02  1.67E-14 -1.02E+03 -1.63E+02  2.40E+03
 
 TH 6
+       -1.43E+01 -8.81E-15 -4.41E+01 -1.38E+01  2.06E+01  3.13E+01
 
 TH 7
+       -6.38E-02 -6.34E-18 -3.17E-02  8.80E-03  8.82E-03  2.24E-02  1.87E-05
 
 TH 8
+       -9.09E+00 -7.78E-15 -4.59E+01 -9.07E+00  3.09E+01  2.73E+01  1.96E-02  2.49E+01
 
 TH 9
+       -4.48E+01  3.25E-15  5.51E+01 -1.62E+01 -8.86E+01 -3.07E+01 -2.10E-02 -1.22E+01  3.01E+02
 
 TH10
+       -9.42E+00 -1.32E-14 -1.57E+01 -1.10E+01 -8.26E+01  4.11E+01  2.87E-02  3.76E+01 -1.66E-01  6.50E+01
 
 TH11
+       -2.98E+01 -8.51E-15 -3.26E+01 -1.50E+00 -8.70E+00  2.49E+01  1.87E-02  2.62E+01  5.53E+01  4.50E+01  4.39E+01
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.22E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        9.60E+01  0.00E+00  4.70E+02
 
 TH 4
+        4.91E+01  0.00E+00 -7.60E+01  6.89E+02
 
 TH 5
+       -3.00E+02  0.00E+00 -9.19E+02 -9.92E+01  2.31E+03
 
 TH 6
+       -8.75E+01  0.00E+00  1.37E+01 -7.58E+00 -5.86E+01  1.49E+02
 
 TH 7
+        8.50E-03  0.00E+00 -6.40E-03 -1.85E-02  2.21E-02 -7.10E-03  9.68E-06
 
 TH 8
+        4.87E+00  0.00E+00 -5.24E+01  4.10E+00 -9.89E+00  1.06E+00  3.07E-03  4.45E+01
 
 TH 9
+        3.90E+01  0.00E+00 -2.50E+01 -1.62E+01  8.28E+01 -3.40E+01  4.46E-02  4.16E+00  3.32E+02
 
 TH10
+        6.75E+01  0.00E+00  1.38E+01 -2.33E+01 -1.85E+02  1.73E+01  8.30E-03  4.32E+01 -1.29E+01  1.50E+02
 
 TH11
+       -1.47E+01  0.00E+00 -1.97E+01  5.80E+01  2.04E+00  2.13E+00  4.98E-03  3.51E+00  5.33E+01  6.47E+00  1.89E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.04
 #CPUT: Total CPU Time in Seconds,       29.566
Stop Time:
Sat Sep 25 07:13:51 CDT 2021
