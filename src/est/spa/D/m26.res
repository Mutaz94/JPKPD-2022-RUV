Wed Sep 29 19:53:37 CDT 2021
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
$DATA ../../../../data/spa/D/dat26.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   8018.66937173713        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.1651E+02  3.6387E+01 -1.1401E+02 -1.3676E+01  2.1675E+02 -8.8892E+02 -4.0397E+02 -2.2820E+01 -7.5227E+02 -2.8425E+02
            -1.6827E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -719.320185335430        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.4527E+00  1.2103E+00  1.0755E+00  1.5488E+00  1.0293E+00  1.6161E+00  1.2904E+00  9.8171E-01  1.4092E+00  1.1579E+00
             1.4429E+01
 PARAMETER:  4.7341E-01  2.9087E-01  1.7275E-01  5.3750E-01  1.2886E-01  5.8001E-01  3.5497E-01  8.1538E-02  4.4302E-01  2.4662E-01
             2.7693E+00
 GRADIENT:   2.4466E+01  6.1812E+00 -1.2440E+00  9.6079E+00 -1.3905E+01  4.1804E+01  2.3278E+00  3.1607E+00  2.0123E+01  6.2470E+00
             2.4548E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -738.756714115332        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.4256E+00  1.0081E+00  6.9811E+00  2.0463E+00  4.6178E+00  1.4264E+00  3.2551E+00  4.9282E-01  1.3862E+00  8.4124E+00
             1.2450E+01
 PARAMETER:  4.5457E-01  1.0808E-01  2.0432E+00  8.1604E-01  1.6299E+00  4.5517E-01  1.2802E+00 -6.0761E-01  4.2658E-01  2.2297E+00
             2.6217E+00
 GRADIENT:   4.0810E+01  2.7308E+01  1.1316E+00  6.7788E+01 -5.4939E+00 -1.8656E+01  1.6410E+01 -2.5858E-03  1.7392E+01  1.3823E+01
             1.7656E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -795.998108895898        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1694E+00  1.4203E+00  3.5983E+01  1.0952E+00  1.3485E+01  1.3502E+00  1.2521E+00  3.1958E+00  1.3731E+00  1.0791E+01
             8.5853E+00
 PARAMETER:  2.5648E-01  4.5088E-01  3.6831E+00  1.9091E-01  2.7016E+00  4.0029E-01  3.2486E-01  1.2618E+00  4.1711E-01  2.4788E+00
             2.2501E+00
 GRADIENT:   1.3182E+01  1.1253E+01  2.8421E-01  9.9928E+00  4.4454E-01 -9.8152E+00 -2.6315E+00 -2.4971E-03  3.4017E-01 -3.6389E-01
            -4.1885E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -797.390464352690        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      383
 NPARAMETR:  1.2019E+00  1.3114E+00  2.2995E+01  1.2095E+00  9.8529E+00  1.4822E+00  1.5256E+00  2.1566E+00  1.3991E+00  1.0025E+01
             8.7137E+00
 PARAMETER:  2.8392E-01  3.7108E-01  3.2353E+00  2.9017E-01  2.3878E+00  4.9355E-01  5.2241E-01  8.6855E-01  4.3582E-01  2.4051E+00
             2.2649E+00
 GRADIENT:   3.0092E+00 -6.3662E-01  5.5949E-01 -1.2985E+00 -1.5060E-01  1.4121E+00 -4.0582E-01 -4.2300E-03  7.6583E-01 -1.3106E-01
             6.9923E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -797.723189130314        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  1.2062E+00  1.2819E+00  1.5678E+01  1.2508E+00  8.2595E+00  1.5155E+00  1.6733E+00  8.8097E+00  1.4162E+00  9.1706E+00
             8.5941E+00
 PARAMETER:  2.8747E-01  3.4834E-01  2.8522E+00  3.2380E-01  2.2114E+00  5.1577E-01  6.1479E-01  2.2759E+00  4.4798E-01  2.3160E+00
             2.2511E+00
 GRADIENT:   2.7246E+01  4.1310E+00  9.2496E-01  4.9717E+00 -3.4937E-01  1.3386E+01  1.6636E+00  6.2912E-02  2.0932E+00  3.6287E-01
             2.1217E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -803.576376209647        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      603
 NPARAMETR:  1.1905E+00  9.7972E-01  1.5511E+00  1.2347E+00  8.6501E+00  1.5198E+00  1.3291E+00  2.2488E+00  1.1537E+00  8.7703E+00
             8.0517E+00
 PARAMETER:  2.7434E-01  7.9511E-02  5.3895E-01  3.1085E-01  2.2576E+00  5.1860E-01  3.8454E-01  9.1039E-01  2.4296E-01  2.2714E+00
             2.1859E+00
 GRADIENT:   6.4038E+01 -8.4826E-01  1.0355E+01 -3.2733E+01  7.2526E-01  2.3568E+01 -3.6837E+00 -6.5734E-01 -1.4588E+01 -1.9544E+00
            -3.8030E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -821.526251898819        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      675
 NPARAMETR:  1.1796E+00  5.2084E-01  6.2414E-01  1.5399E+00  8.8768E+00  1.5188E+00  1.1607E+00  2.5143E-01  1.3058E+00  8.6057E+00
             7.8782E+00
 PARAMETER:  2.6520E-01 -5.5231E-01 -3.7138E-01  5.3174E-01  2.2834E+00  5.1790E-01  2.4901E-01 -1.2806E+00  3.6682E-01  2.2524E+00
             2.1641E+00
 GRADIENT:   7.0162E+01  3.5975E+01 -2.2390E+01  5.9583E+01 -1.0890E+01 -1.3527E+01  2.4800E+00  1.4955E-01 -4.3605E+00  1.6985E+01
            -2.3511E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -844.627720688392        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      746
 NPARAMETR:  1.1235E+00  6.8970E-02  3.5670E-01  1.4186E+00  1.0097E+01  1.5041E+00  2.4062E+00  2.0692E-02  9.2969E-01  7.6842E+00
             7.4819E+00
 PARAMETER:  2.1641E-01 -2.5741E+00 -9.3087E-01  4.4967E-01  2.4122E+00  5.0819E-01  9.7807E-01 -3.7780E+00  2.7096E-02  2.1392E+00
             2.1125E+00
 GRADIENT:   1.3325E+02  4.9093E-01 -1.2318E+01  3.1302E+01  1.0347E+01 -6.7307E+00  1.5566E-01 -7.3275E-03  4.7219E+00 -3.7277E+00
            -4.9240E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -856.036760678226        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      819
 NPARAMETR:  9.9478E-01  1.0000E-02  2.4834E-01  1.2110E+00  9.4739E+00  1.4395E+00  6.6531E+00  1.0000E-02  4.2462E-01  6.9705E+00
             9.0140E+00
 PARAMETER:  9.4769E-02 -6.4110E+00 -1.2930E+00  2.9143E-01  2.3485E+00  4.6428E-01  1.9951E+00 -6.3214E+00 -7.5656E-01  2.0417E+00
             2.2988E+00
 GRADIENT:   7.4713E+01  0.0000E+00  1.6414E+01 -2.5294E+01  1.2780E+01  1.7871E+01  3.0525E-02  0.0000E+00  5.2270E+00  1.6308E+00
             5.7934E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -858.957438666445        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      942             RESET HESSIAN, TYPE I
 NPARAMETR:  8.6602E-01  1.0000E-02  2.3211E-01  1.1813E+00  9.5709E+00  1.4295E+00  2.8550E+00  1.0000E-02  3.8971E-01  6.7923E+00
             9.1107E+00
 PARAMETER: -4.3852E-02 -6.8713E+00 -1.3606E+00  2.6659E-01  2.3587E+00  4.5732E-01  1.1491E+00 -6.7189E+00 -8.4236E-01  2.0158E+00
             2.3094E+00
 GRADIENT:  -4.8847E+01  0.0000E+00  2.1023E+01  1.3128E+01  1.2591E+01  1.8565E+01  5.5327E-03  0.0000E+00  7.8913E+00 -3.1998E+00
             8.8384E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -861.247006680919        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1014
 NPARAMETR:  9.0235E-01  1.0000E-02  2.3077E-01  1.1806E+00  9.5252E+00  1.4252E+00  1.1921E+00  1.0000E-02  3.8968E-01  6.7932E+00
             8.4434E+00
 PARAMETER: -2.7507E-03 -6.8713E+00 -1.3664E+00  2.6602E-01  2.3539E+00  4.5432E-01  2.7572E-01 -6.7189E+00 -8.4244E-01  2.0159E+00
             2.2334E+00
 GRADIENT:   5.9036E+00  0.0000E+00  8.2608E+00  2.9017E+01  1.2293E+01  1.8075E+01  8.1925E-04  0.0000E+00  3.1486E+00 -9.5261E-01
             1.3749E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -861.315842209892        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1194
 NPARAMETR:  9.0479E-01  1.0000E-02  2.3088E-01  1.1803E+00  9.4816E+00  1.4240E+00  1.1593E+00  1.0000E-02  3.8950E-01  6.8083E+00
             8.5330E+00
 PARAMETER: -5.0323E-05 -6.8713E+00 -1.3659E+00  2.6579E-01  2.3494E+00  4.5346E-01  2.4786E-01 -6.7189E+00 -8.4290E-01  2.0181E+00
             2.2439E+00
 GRADIENT:   1.7420E+00  0.0000E+00 -2.4642E+00  2.4943E+01  3.6043E+00  1.4710E+01  4.2822E-04  0.0000E+00  1.1570E+00  1.3844E+00
             1.2258E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -861.385383610615        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1376
 NPARAMETR:  8.9795E-01  1.0000E-02  2.3123E-01  1.1791E+00  9.4126E+00  1.4191E+00  9.5243E-01  1.0000E-02  3.8923E-01  6.8059E+00
             8.5673E+00
 PARAMETER: -7.6402E-03 -6.8713E+00 -1.3643E+00  2.6478E-01  2.3421E+00  4.4999E-01  5.1262E-02 -6.7189E+00 -8.4359E-01  2.0178E+00
             2.2480E+00
 GRADIENT:   9.9727E+00  0.0000E+00 -4.0918E+01  8.6799E+01 -5.4491E+01  3.1434E+01 -1.3326E-03  0.0000E+00 -1.3696E+01  2.8287E+01
            -2.2144E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -861.675485563655        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1553
 NPARAMETR:  8.9222E-01  1.0000E-02  2.3372E-01  1.1729E+00  9.1787E+00  1.3939E+00  2.3535E-01  1.0000E-02  3.8834E-01  6.7321E+00
             8.5725E+00
 PARAMETER: -1.4045E-02 -6.8713E+00 -1.3536E+00  2.5950E-01  2.3169E+00  4.3210E-01 -1.3467E+00 -6.7189E+00 -8.4588E-01  2.0069E+00
             2.2486E+00
 GRADIENT:  -1.3292E+01  0.0000E+00  1.0455E+01  5.8408E+00  4.2144E+00  6.4474E+00  2.4896E-05  0.0000E+00  3.3277E+00  2.5448E+00
             1.0610E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -861.877863476145        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1727
 NPARAMETR:  9.0492E-01  1.0000E-02  2.3544E-01  1.1683E+00  9.2388E+00  1.3819E+00  8.6765E-02  1.0000E-02  3.8957E-01  6.5511E+00
             8.5507E+00
 PARAMETER:  5.9915E-04 -6.8713E+00 -1.3523E+00  2.5684E-01  2.3131E+00  4.2370E-01 -2.3213E+00 -6.7189E+00 -8.4688E-01  1.9894E+00
             2.2357E+00
 GRADIENT:   4.0704E+02  0.0000E+00 -1.2040E+02  6.9123E+02 -9.9147E+01  4.4734E+01  5.1869E-05  0.0000E+00 -2.0779E+02  1.3046E+02
            -1.1823E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1727
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0813E-03 -7.8376E-06 -8.2527E-05 -1.1171E-02 -3.7320E-02
 SE:             2.7982E-02  4.8405E-06  1.6305E-04  1.2818E-02  9.0444E-03
 N:                     100         100         100         100         100

 P VAL.:         9.6917E-01  1.0541E-01  6.1276E-01  3.8347E-01  3.6881E-05

 ETASHRINKSD(%)  6.2569E+00  9.9984E+01  9.9454E+01  5.7058E+01  6.9700E+01
 ETASHRINKVR(%)  1.2122E+01  1.0000E+02  9.9997E+01  8.1560E+01  9.0819E+01
 EBVSHRINKSD(%)  4.6566E+00  9.9976E+01  9.9465E+01  5.5568E+01  6.2746E+01
 EBVSHRINKVR(%)  9.0963E+00  1.0000E+02  9.9997E+01  8.0258E+01  8.6122E+01
 RELATIVEINF(%)  2.8277E+01  4.4426E-07  1.0177E-04  3.8587E-01  4.0451E+00
 EPSSHRINKSD(%)  1.0262E+01
 EPSSHRINKVR(%)  1.9471E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -861.87786347614497     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -126.72703691240679     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -861.878       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.05E-01  1.00E-02  2.34E-01  1.17E+00  9.14E+00  1.38E+00  8.88E-02  1.00E-02  3.88E-01  6.62E+00  8.46E+00
 


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
+        5.44E+05
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.49E+02  0.00E+00  7.05E+04
 
 TH 4
+       -1.02E+02  0.00E+00 -2.37E+04  5.07E+04
 
 TH 5
+       -2.21E+00  0.00E+00  3.59E+02 -3.42E+02  2.75E+01
 
 TH 6
+       -1.49E+00  0.00E+00  1.11E+04 -5.00E+01  1.88E-01  1.31E+04
 
 TH 7
+        4.20E-02  0.00E+00 -7.02E-03 -1.98E-02 -2.13E-04  1.97E-02  1.02E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.11E+01  0.00E+00  1.99E+04 -2.82E+02  3.03E+02 -1.75E+01 -2.36E-02  0.00E+00  4.15E+04
 
 TH10
+        1.74E+00  0.00E+00 -5.38E+02  2.32E+01 -1.66E+01 -3.97E-01  1.50E-04  0.00E+00 -1.06E+01  3.94E+01
 
 TH11
+       -1.33E+01  0.00E+00  3.92E+02 -3.81E+02  1.11E+01  2.24E+00 -8.85E-05  0.00E+00  1.49E+01 -9.07E+00  2.51E+01
 
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
 #CPUT: Total CPU Time in Seconds,       32.239
Stop Time:
Wed Sep 29 19:54:11 CDT 2021
