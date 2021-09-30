Thu Sep 30 09:35:51 CDT 2021
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
$DATA ../../../../data/spa2/D/dat71.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   18720.5095769540        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2183E+02  3.2122E+02  8.1235E+00 -8.6252E+01  2.7784E+02 -2.4889E+03 -1.0841E+03 -1.0006E+02 -2.0388E+03 -7.6776E+02
            -3.5526E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -582.974644747040        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4229E+00  2.2945E+00  9.2575E-01  3.8611E+00  1.0882E+00  6.2287E+00  5.0704E+00  9.8192E-01  4.5842E+00  1.3217E+00
             9.2465E+00
 PARAMETER:  4.5271E-01  9.3050E-01  2.2845E-02  1.4510E+00  1.8452E-01  1.9292E+00  1.7234E+00  8.1750E-02  1.6226E+00  3.7894E-01
             2.3242E+00
 GRADIENT:   2.6757E+01  3.9385E+01 -4.1004E+01  1.2579E+02 -3.9848E+01  3.4138E+02  5.7924E+01  2.9389E+00  3.5216E+01  1.6540E+01
             6.8009E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -672.421145471404        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.6527E+00  3.2019E+00  4.1713E+00  6.0137E+00  3.4575E+00  4.6857E+00  3.8415E+00  8.8503E-01  1.2942E+01  1.5100E+00
             8.7991E+00
 PARAMETER:  6.0241E-01  1.2637E+00  1.5282E+00  1.8940E+00  1.3405E+00  1.6445E+00  1.4459E+00 -2.2138E-02  2.6604E+00  5.1213E-01
             2.2746E+00
 GRADIENT:   6.5865E+01  3.1131E+01 -1.1497E+01  5.6882E+01  1.4892E-01  1.7240E+02  4.2455E+01 -3.6463E+00  1.0106E+02  7.3028E+00
             8.5202E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -712.924205643926        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      339
 NPARAMETR:  1.3618E+00  5.4770E+00  7.5008E+00  4.7716E+00  3.0062E+00  3.6082E+00  4.3840E+00  4.5650E-01  1.4629E+01  7.4721E-01
             8.9009E+00
 PARAMETER:  4.0879E-01  1.8005E+00  2.1150E+00  1.6627E+00  1.2007E+00  1.3832E+00  1.5780E+00 -6.8417E-01  2.7830E+00 -1.9140E-01
             2.2862E+00
 GRADIENT:   1.9799E+01  2.7952E+01 -5.5136E+00  3.0607E+01 -2.6270E+00  6.7021E+01  1.0232E+01  8.5932E-01 -6.3113E+00  2.2916E+00
             6.0080E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -766.391365284484        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  8.1340E-01  1.3728E+00  1.8794E+01  9.1148E-01  3.4452E+00  2.4725E+00  6.2740E+00  3.3535E+00  2.4279E+00  9.6270E-01
             9.1472E+00
 PARAMETER: -1.0654E-01  4.1682E-01  3.0335E+00  7.3186E-03  1.3370E+00  1.0052E+00  1.9364E+00  1.3100E+00  9.8702E-01  6.1986E-02
             2.3135E+00
 GRADIENT:  -1.0663E+02 -1.8352E+01 -1.5795E-01 -4.8843E+01  7.8273E+00 -6.0255E+01  3.2324E+01  1.0675E-02  2.0746E+01  3.1261E+00
             1.0111E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -808.935369854894        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      693
 NPARAMETR:  1.2126E+00  9.0934E-01  1.9824E+01  1.5179E+00  2.4842E+00  2.9221E+00  7.4595E+00  3.5467E+00  1.7702E+00  4.5595E-01
             8.6288E+00
 PARAMETER:  2.9273E-01  4.9638E-03  3.0869E+00  5.1730E-01  1.0100E+00  1.1723E+00  2.1095E+00  1.3660E+00  6.7110E-01 -6.8537E-01
             2.2551E+00
 GRADIENT:   4.8788E+00  6.5468E-01  1.0007E-01  3.4577E-01 -1.2765E+01  1.3962E+01  8.0041E+00  7.5802E-02  8.3999E-01  1.1355E+00
             8.0176E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -810.362252688322        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      869
 NPARAMETR:  1.1708E+00  9.4223E-01  2.3055E+01  1.4975E+00  2.8054E+00  2.7868E+00  7.1318E+00  2.9014E+00  1.7185E+00  1.6829E-01
             8.6898E+00
 PARAMETER:  2.5772E-01  4.0491E-02  3.2379E+00  5.0377E-01  1.1315E+00  1.1249E+00  2.0646E+00  1.1652E+00  6.4145E-01 -1.6820E+00
             2.2621E+00
 GRADIENT:  -4.5449E+00 -1.2544E+00 -3.2813E-01  6.6161E-01  2.2717E-01 -1.5724E+00  1.3741E+00  3.6264E-02 -1.2306E-01  1.2873E-01
             9.1012E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -810.633983219108        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1049
 NPARAMETR:  1.1870E+00  1.1464E+00  8.2308E+01  1.3917E+00  2.8407E+00  2.8193E+00  6.8054E+00  8.1286E-01  1.6379E+00  1.0000E-02
             8.6224E+00
 PARAMETER:  2.7146E-01  2.3665E-01  4.5105E+00  4.3054E-01  1.1441E+00  1.1365E+00  2.0177E+00 -1.0720E-01  5.9344E-01 -4.5507E+00
             2.2544E+00
 GRADIENT:  -5.3687E-01  7.8257E-01 -1.4045E-03 -2.1825E-01 -1.1421E+00  2.0527E+00  3.3615E+00  1.5863E-04 -3.1245E-01  0.0000E+00
            -1.6510E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -810.660285116172        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1224
 NPARAMETR:  1.1899E+00  1.1068E+00  1.6144E+02  1.4194E+00  2.8897E+00  2.8050E+00  6.7920E+00  6.9927E-01  1.6778E+00  1.0000E-02
             8.6274E+00
 PARAMETER:  2.7384E-01  2.0148E-01  5.1841E+00  4.5024E-01  1.1611E+00  1.1314E+00  2.0157E+00 -2.5772E-01  6.1749E-01 -4.5669E+00
             2.2549E+00
 GRADIENT:  -3.6203E-02  1.3120E-01 -7.5145E-03  4.6419E-01 -6.8602E-02  3.4039E-01  6.8775E-01  3.0658E-05 -9.5598E-02  0.0000E+00
             2.5663E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -810.685327631952        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1404
 NPARAMETR:  1.1907E+00  1.1026E+00  7.7281E+02  1.4187E+00  2.8955E+00  2.8183E+00  6.9096E+00  6.7024E-01  1.6830E+00  1.0000E-02
             8.6231E+00
 PARAMETER:  2.7451E-01  1.9766E-01  6.7500E+00  4.4973E-01  1.1632E+00  1.1361E+00  2.0329E+00 -3.0012E-01  6.2060E-01 -4.5805E+00
             2.2544E+00
 GRADIENT:   1.7734E-01  5.4302E-01 -9.4063E-04 -1.8663E+00 -9.4947E-02  1.8754E+00  4.4511E+00  1.0667E-06  6.0750E-01  0.0000E+00
             5.2575E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -810.703134048131        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1579
 NPARAMETR:  1.1923E+00  1.0720E+00  1.7046E+03  1.4392E+00  2.9065E+00  2.8097E+00  6.8945E+00  6.7024E-01  1.6889E+00  1.0000E-02
             8.6311E+00
 PARAMETER:  2.7588E-01  1.6950E-01  7.5411E+00  4.6411E-01  1.1669E+00  1.1331E+00  2.0307E+00 -3.0012E-01  6.2405E-01 -4.5805E+00
             2.2554E+00
 GRADIENT:   4.4886E-01  8.9246E-02 -7.0490E-04  1.9983E-02  1.3689E-02  9.4833E-01  1.7775E+00  6.0280E-07 -1.7016E-01  0.0000E+00
             6.8040E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -810.718772104374        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1758
 NPARAMETR:  1.1924E+00  1.0533E+00  5.6022E+03  1.4496E+00  2.9065E+00  2.8120E+00  6.9489E+00  6.7225E-01  1.6958E+00  1.0000E-02
             8.6304E+00
 PARAMETER:  2.7600E-01  1.5194E-01  8.7309E+00  4.7126E-01  1.1670E+00  1.1339E+00  2.0386E+00 -2.9712E-01  6.2814E-01 -4.5805E+00
             2.2553E+00
 GRADIENT:   5.1163E-01  7.2227E-02 -2.3598E-04 -1.2816E-01 -1.8555E-02  1.2421E+00  2.3107E+00 -3.5052E-07 -1.8015E-01  0.0000E+00
             6.1535E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -810.728760179906        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1937
 NPARAMETR:  1.1923E+00  1.0404E+00  8.3075E+03  1.4575E+00  2.9072E+00  2.8139E+00  6.9920E+00  6.7286E-01  1.7034E+00  1.0000E-02
             8.6291E+00
 PARAMETER:  2.7590E-01  1.3960E-01  9.1249E+00  4.7673E-01  1.1672E+00  1.1346E+00  2.0448E+00 -2.9622E-01  6.3264E-01 -4.5805E+00
             2.2551E+00
 GRADIENT:   4.3461E-01  1.0404E-01 -1.7826E-04 -3.0182E-01 -1.5329E-03  1.4917E+00  2.8379E+00 -4.5411E-07 -9.7274E-02  0.0000E+00
             5.2899E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -810.735576423485        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2115
 NPARAMETR:  1.1923E+00  1.0286E+00  9.7540E+03  1.4647E+00  2.9071E+00  2.8155E+00  7.0258E+00  6.7297E-01  1.7095E+00  1.0000E-02
             8.6272E+00
 PARAMETER:  2.7584E-01  1.2815E-01  9.2854E+00  4.8165E-01  1.1671E+00  1.1351E+00  2.0496E+00 -2.9605E-01  6.3619E-01 -4.5805E+00
             2.2549E+00
 GRADIENT:   4.7038E-01  9.6691E-02 -1.6546E-04 -3.4285E-01 -1.0983E-02  1.7024E+00  3.1144E+00 -2.0510E-07 -8.2405E-02  0.0000E+00
             2.9979E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -810.739999307268        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2293
 NPARAMETR:  1.1918E+00  1.0207E+00  1.3858E+04  1.4708E+00  2.9077E+00  2.8160E+00  7.0526E+00  6.7306E-01  1.7154E+00  1.0000E-02
             8.6271E+00
 PARAMETER:  2.7547E-01  1.2046E-01  9.6366E+00  4.8581E-01  1.1674E+00  1.1353E+00  2.0534E+00 -2.9592E-01  6.3964E-01 -4.5805E+00
             2.2549E+00
 GRADIENT:   2.8899E-01  1.4341E-01 -1.2613E-04 -2.9856E-01 -8.8523E-03  1.8727E+00  3.3651E+00  1.2125E-07 -5.4126E-02  0.0000E+00
             3.1202E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -810.742912579831        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2471
 NPARAMETR:  1.1920E+00  1.0138E+00  1.8032E+04  1.4758E+00  2.9083E+00  2.8162E+00  7.0719E+00  6.7179E-01  1.7200E+00  1.0000E-02
             8.6266E+00
 PARAMETER:  2.7565E-01  1.1373E-01  9.8999E+00  4.8919E-01  1.1676E+00  1.1354E+00  2.0561E+00 -2.9781E-01  6.4233E-01 -4.5805E+00
             2.2548E+00
 GRADIENT:   8.9624E-02  1.8369E-01 -1.0317E-04 -1.9709E-01 -9.0963E-03  1.9288E+00  3.4893E+00 -1.9465E-07 -3.1619E-02  0.0000E+00
             2.6276E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -810.745783258305        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2659
 NPARAMETR:  1.1908E+00  1.0065E+00  2.4574E+04  1.4793E+00  2.9074E+00  2.8181E+00  7.0948E+00  6.7707E-01  1.7223E+00  1.0000E-02
             8.6231E+00
 PARAMETER:  2.7460E-01  1.0646E-01  1.0209E+01  4.9154E-01  1.1673E+00  1.1361E+00  2.0594E+00 -2.8998E-01  6.4368E-01 -4.5645E+00
             2.2544E+00
 GRADIENT:  -5.0464E-01  1.9279E-01 -7.7053E-05 -2.0511E-01 -3.0089E-02  2.3339E+00  3.7085E+00 -2.7859E-07 -6.4791E-02  0.0000E+00
            -3.4412E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -810.746601870971        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     2799
 NPARAMETR:  1.1908E+00  9.9648E-01  2.4289E+04  1.4825E+00  2.9077E+00  2.8161E+00  7.0949E+00  6.7637E-01  1.7235E+00  1.0000E-02
             8.6265E+00
 PARAMETER:  2.7466E-01  9.6477E-02  1.0198E+01  4.9370E-01  1.1674E+00  1.1353E+00  2.0594E+00 -2.9102E-01  6.4436E-01 -4.5645E+00
             2.2548E+00
 GRADIENT:  -4.3204E-01 -5.6544E-02 -8.1557E-05 -2.1930E-01 -8.3352E-03  2.1185E+00  3.1968E+00  1.6102E-07 -1.3213E-01  0.0000E+00
             4.4506E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -810.749411957163        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2985
 NPARAMETR:  1.1898E+00  9.9135E-01  3.7588E+04  1.4876E+00  2.9084E+00  2.8193E+00  7.1207E+00  6.7045E-01  1.7280E+00  1.0000E-02
             8.6225E+00
 PARAMETER:  2.7378E-01  9.1317E-02  1.0634E+01  4.9720E-01  1.1676E+00  1.1365E+00  2.0630E+00 -2.9981E-01  6.4699E-01 -4.5645E+00
             2.2544E+00
 GRADIENT:  -2.7798E-01  6.3122E-02 -5.6200E-05 -6.2256E-02 -6.6239E-03  2.2428E+00  3.5204E+00  7.7527E-07 -1.4217E-01  0.0000E+00
            -1.8926E-01

0ITERATION NO.:   93    OBJECTIVE VALUE:  -810.750027644256        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:     3088
 NPARAMETR:  1.1907E+00  9.8973E-01  4.2012E+04  1.4936E+00  2.9071E+00  2.8184E+00  7.1340E+00  6.6933E-01  1.7299E+00  1.0000E-02
             8.6222E+00
 PARAMETER:  2.7449E-01  9.0682E-02  1.0759E+01  4.9766E-01  1.1676E+00  1.1361E+00  2.0646E+00 -3.0281E-01  6.4856E-01 -4.5645E+00
             2.2544E+00
 GRADIENT:  -7.5755E-03  8.3217E-02  6.0846E-05 -3.1721E-01  1.3980E-02 -1.1221E-02 -2.1586E-02 -3.5582E-03  1.2162E-02  0.0000E+00
             1.6410E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3088
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4109E-03  4.2639E-02  8.1576E-08 -8.4218E-02 -2.0462E-05
 SE:             2.8897E-02  2.3266E-02  4.7726E-08  1.4970E-02  4.4298E-05
 N:                     100         100         100         100         100

 P VAL.:         7.7100E-01  6.6856E-02  8.7406E-02  1.8508E-08  6.4414E-01

 ETASHRINKSD(%)  3.1920E+00  2.2055E+01  1.0000E+02  4.9849E+01  9.9852E+01
 ETASHRINKVR(%)  6.2821E+00  3.9246E+01  1.0000E+02  7.4849E+01  1.0000E+02
 EBVSHRINKSD(%)  3.5082E+00  1.7577E+01  1.0000E+02  4.9566E+01  9.9790E+01
 EBVSHRINKVR(%)  6.8933E+00  3.2064E+01  1.0000E+02  7.4564E+01  1.0000E+02
 RELATIVEINF(%)  9.1884E+01  3.4412E+01  1.1176E-10  1.2217E+01  7.4092E-05
 EPSSHRINKSD(%)  1.0924E+01
 EPSSHRINKVR(%)  2.0655E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -810.75002764425608     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       291.97621220135102     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    75.41
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -810.750       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.19E+00  9.91E-01  4.26E+04  1.49E+00  2.91E+00  2.82E+00  7.13E+00  6.68E-01  1.73E+00  1.00E-02  8.62E+00
 


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
+        9.01E+01
 
 TH 2
+        8.33E-01  2.74E+01
 
 TH 3
+       -4.78E-07 -9.34E-07  9.06E-13
 
 TH 4
+       -1.29E+00  1.89E+01 -3.26E-08  5.79E+01
 
 TH 5
+       -2.50E-01 -1.41E+00  2.98E-08 -3.36E+00  5.43E+00
 
 TH 6
+       -3.40E-01 -6.84E-02 -1.01E-08  5.98E-01 -1.00E-01  2.13E+01
 
 TH 7
+        1.43E-01  2.31E+00 -2.33E-08 -4.65E+00  1.20E-01 -9.48E-02  2.05E+00
 
 TH 8
+        3.65E-01  5.16E-01  2.17E-06  4.66E-01  9.58E-02  4.21E-02 -3.57E-02  5.05E+00
 
 TH 9
+        2.99E-01 -1.12E+00 -7.47E-09 -1.63E+01  8.89E-01 -4.05E-01  1.16E+00 -2.46E-01  1.25E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  4.86E+01
 
 TH11
+       -3.44E+00 -1.40E+00  1.65E-09 -4.51E+00 -7.33E-02  7.65E-01  1.39E-01  7.73E-03  1.98E+00  0.00E+00  8.25E+00
 
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
 #CPUT: Total CPU Time in Seconds,       87.157
Stop Time:
Thu Sep 30 09:37:19 CDT 2021
