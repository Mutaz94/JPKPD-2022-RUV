Sat Sep 25 02:00:31 CDT 2021
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
$DATA ../../../../data/int/SL3/dat19.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      979
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      879
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -244.571512321992        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.3205E+01  3.6087E+01  6.7408E+01  3.4950E+01  6.5089E+01  1.5373E+01 -1.1488E+02 -1.2814E+02 -2.4366E+01 -5.8038E+01
            -6.8319E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2698.25285231370        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0560E+00  1.1782E+00  1.0344E+00  9.4388E-01  1.1061E+00  8.8995E-01  1.1160E+00  9.6704E-01  6.7545E-01  1.2175E+00
             2.9229E+00
 PARAMETER:  1.5444E-01  2.6395E-01  1.3378E-01  4.2247E-02  2.0084E-01 -1.6590E-02  2.0971E-01  6.6483E-02 -2.9238E-01  2.9684E-01
             1.1726E+00
 GRADIENT:   2.1924E+01  3.4035E+00 -1.0684E+01 -3.1960E+00 -1.0899E+01 -2.1640E+01  1.4779E+01  3.4305E+00  5.4277E+00 -5.1816E+00
             1.1309E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2701.79408897291        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0638E+00  1.3925E+00  1.3321E+00  8.6152E-01  1.3696E+00  9.5547E-01  9.5880E-01  9.6371E-01  6.8466E-01  1.4924E+00
             2.9217E+00
 PARAMETER:  1.6184E-01  4.3111E-01  3.8679E-01 -4.9054E-02  4.1450E-01  5.4449E-02  5.7931E-02  6.3030E-02 -2.7883E-01  5.0038E-01
             1.1722E+00
 GRADIENT:   3.6458E+01  5.7117E+01 -1.9095E+00  6.6010E+01  5.3875E+00  4.5386E+00  9.4456E+00 -4.7646E-01  4.0815E+00  7.8111E+00
             1.0978E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2708.67366113540        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0459E+00  1.5538E+00  1.4316E+00  7.2646E-01  1.5285E+00  9.4792E-01  8.5710E-01  1.9631E+00  6.3166E-01  1.5186E+00
             2.7357E+00
 PARAMETER:  1.4492E-01  5.4070E-01  4.5878E-01 -2.1957E-01  5.2426E-01  4.6515E-02 -5.4206E-02  7.7455E-01 -3.5941E-01  5.1782E-01
             1.1064E+00
 GRADIENT:   4.7308E+00  2.1141E+01 -1.2685E+00  1.4503E+01  8.0875E-01  9.0292E-01  2.3028E+00  2.6038E-01  3.6843E-01 -2.0358E+00
            -2.0172E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2709.64621095228        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  1.0394E+00  1.6109E+00  2.0021E+00  6.7469E-01  1.7343E+00  9.3853E-01  8.2375E-01  2.9911E+00  5.5786E-01  1.6386E+00
             2.7938E+00
 PARAMETER:  1.3867E-01  5.7678E-01  7.9422E-01 -2.9350E-01  6.5059E-01  3.6563E-02 -9.3883E-02  1.1957E+00 -4.8364E-01  5.9382E-01
             1.1274E+00
 GRADIENT:  -1.2179E+01 -2.7898E+01 -5.4868E+00 -1.1539E+01  1.4746E+01 -2.6193E+00 -2.0635E+00  5.8326E-01  7.7819E-01  4.9915E+00
             4.1154E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2709.66108943191        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  1.0400E+00  1.6045E+00  1.9273E+00  6.7893E-01  1.7059E+00  9.3949E-01  8.2617E-01  2.8323E+00  5.6262E-01  1.6259E+00
             2.7889E+00
 PARAMETER:  1.3927E-01  5.7279E-01  7.5611E-01 -2.8724E-01  6.3412E-01  3.7583E-02 -9.0956E-02  1.1411E+00 -4.7516E-01  5.8605E-01
             1.1257E+00
 GRADIENT:  -1.0356E+01 -2.4885E+01 -4.7727E+00 -1.1617E+01  1.0585E+01 -2.1722E+00 -1.9316E+00 -2.9991E-01  5.9553E-01  4.4159E+00
             3.5435E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2709.66386456503        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      449
 NPARAMETR:  1.0404E+00  1.6018E+00  1.8923E+00  6.8080E-01  1.6920E+00  9.4006E-01  8.2791E-01  2.7604E+00  5.6297E-01  1.6193E+00
             2.7856E+00
 PARAMETER:  1.3963E-01  5.7115E-01  7.3780E-01 -2.8448E-01  6.2592E-01  3.8192E-02 -8.8848E-02  1.1154E+00 -4.7453E-01  5.8199E-01
             1.1245E+00
 GRADIENT:  -9.2541E+00 -2.2861E+01 -4.3116E+00 -1.1487E+01  8.3567E+00 -1.9200E+00 -1.8046E+00 -6.4593E-01  4.6756E-01  4.0378E+00
             3.1808E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2709.66408047546        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  1.0406E+00  1.5988E+00  1.8854E+00  6.8286E-01  1.6866E+00  9.4026E-01  8.2957E-01  2.7367E+00  5.6211E-01  1.6161E+00
             2.7844E+00
 PARAMETER:  1.3976E-01  5.6928E-01  7.3416E-01 -2.8146E-01  6.2270E-01  3.8398E-02 -8.6851E-02  1.1068E+00 -4.7606E-01  5.8001E-01
             1.1240E+00
 GRADIENT:  -8.8754E+00 -2.2156E+01 -4.1641E+00 -1.1440E+01  7.6231E+00 -1.8357E+00 -1.7511E+00 -7.6050E-01  4.2061E-01  3.9037E+00
             3.0524E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2709.66433217549        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      598
 NPARAMETR:  1.0407E+00  1.5927E+00  1.8854E+00  6.8701E-01  1.6798E+00  9.4045E-01  8.3268E-01  2.7130E+00  5.6003E-01  1.6115E+00
             2.7832E+00
 PARAMETER:  1.3988E-01  5.6542E-01  7.3414E-01 -2.7540E-01  6.1868E-01  3.8600E-02 -8.3100E-02  1.0980E+00 -4.7976E-01  5.7714E-01
             1.1236E+00
 GRADIENT:  -8.5030E+00 -2.1479E+01 -4.0245E+00 -1.1438E+01  6.8523E+00 -1.7526E+00 -1.7010E+00 -8.9210E-01  3.6910E-01  3.7760E+00
             2.9262E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2709.66488328910        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  1.0409E+00  1.5766E+00  1.8992E+00  6.9784E-01  1.6663E+00  9.4076E-01  8.4069E-01  2.6756E+00  5.5426E-01  1.6012E+00
             2.7811E+00
 PARAMETER:  1.4006E-01  5.5527E-01  7.4143E-01 -2.5977E-01  6.1062E-01  3.8929E-02 -7.3537E-02  1.0842E+00 -4.9012E-01  5.7074E-01
             1.1229E+00
 GRADIENT:  -7.9015E+00 -2.0342E+01 -3.8063E+00 -1.1476E+01  5.4846E+00 -1.6167E+00 -1.6306E+00 -1.1501E+00  2.7335E-01  3.5797E+00
             2.7221E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2709.66626054691        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      745
 NPARAMETR:  1.0415E+00  1.5088E+00  1.9991E+00  7.4440E-01  1.6195E+00  9.4182E-01  8.7642E-01  2.5917E+00  5.2772E-01  1.5608E+00
             2.7741E+00
 PARAMETER:  1.4065E-01  5.1132E-01  7.9270E-01 -1.9518E-01  5.8213E-01  4.0063E-02 -3.1912E-02  1.0523E+00 -5.3920E-01  5.4521E-01
             1.1203E+00
 GRADIENT:  -6.0006E+00 -1.4801E+01 -3.0818E+00 -1.0329E+01  5.2281E-01 -1.1829E+00 -1.5382E+00 -2.2048E+00 -1.3122E-01  3.0201E+00
             2.0606E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2709.69598629842        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      820
 NPARAMETR:  1.0428E+00  1.4144E+00  2.2171E+00  8.1130E-01  1.5627E+00  9.4418E-01  9.3673E-01  2.5936E+00  4.8098E-01  1.5015E+00
             2.7590E+00
 PARAMETER:  1.4189E-01  4.4669E-01  8.9621E-01 -1.0912E-01  5.4642E-01  4.2561E-02  3.4637E-02  1.0530E+00 -6.3192E-01  5.0645E-01
             1.1149E+00
 GRADIENT:  -2.1498E+00 -3.2865E+00 -1.3093E+00 -6.2847E+00 -8.9765E+00 -3.3062E-01 -1.4912E+00 -4.0972E+00 -9.7740E-01  1.9205E+00
             8.2932E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2710.29581510664        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      894
 NPARAMETR:  1.0465E+00  1.3701E+00  2.8390E+00  8.4819E-01  1.5842E+00  9.5120E-01  1.0104E+00  3.2392E+00  3.9690E-01  1.4358E+00
             2.6989E+00
 PARAMETER:  1.4549E-01  4.1490E-01  1.1435E+00 -6.4650E-02  5.6011E-01  4.9970E-02  1.1032E-01  1.2753E+00 -8.2408E-01  4.6172E-01
             1.0928E+00
 GRADIENT:   8.9001E+00  1.1505E+01  2.7888E+00 -3.2630E+00 -1.3841E+01  1.7379E+00  2.1720E+00 -2.2841E+00 -1.5004E+00 -3.4415E+00
            -3.2659E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2710.96692583224        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  1.0439E+00  1.3686E+00  2.9035E+00  8.5336E-01  1.6141E+00  9.4689E-01  9.8048E-01  3.3101E+00  4.5865E-01  1.4647E+00
             2.7352E+00
 PARAMETER:  1.4293E-01  4.1376E-01  1.1659E+00 -5.8572E-02  5.7877E-01  4.5425E-02  8.0290E-02  1.2970E+00 -6.7946E-01  4.8166E-01
             1.1062E+00
 GRADIENT:  -8.9650E+00  3.7752E+00  1.2188E-01  3.5491E+00 -4.4268E+00 -5.9541E-01  3.9258E-02 -1.0862E-01 -6.0246E-01 -1.6471E-01
            -2.3171E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2711.08989809891        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1165
 NPARAMETR:  1.0477E+00  1.4726E+00  2.7448E+00  7.8329E-01  1.6773E+00  9.4782E-01  9.1815E-01  3.3685E+00  4.8015E-01  1.5288E+00
             2.7404E+00
 PARAMETER:  1.4661E-01  4.8702E-01  1.1097E+00 -1.4425E-01  6.1717E-01  4.6412E-02  1.4608E-02  1.3145E+00 -6.3366E-01  5.2446E-01
             1.1081E+00
 GRADIENT:  -4.2349E-01 -3.5897E-01 -7.8223E-02 -1.3583E-01  3.0180E-01 -1.0044E-01  2.0836E-03  5.9013E-02 -8.0861E-02  6.5644E-02
             3.7201E-01

0ITERATION NO.:   74    OBJECTIVE VALUE:  -2711.09061598544        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1292
 NPARAMETR:  1.0479E+00  1.4701E+00  2.7492E+00  7.8516E-01  1.6750E+00  9.4807E-01  9.1813E-01  3.3637E+00  4.8491E-01  1.5265E+00
             2.7398E+00
 PARAMETER:  1.4677E-01  4.8535E-01  1.1113E+00 -1.4186E-01  6.1583E-01  4.6677E-02  1.4582E-02  1.3131E+00 -6.2379E-01  5.2301E-01
             1.1079E+00
 GRADIENT:  -2.0426E-02 -4.2364E-02 -2.4839E-03 -3.9786E-02 -5.7929E-03 -6.0692E-04  6.4156E-03  2.0932E-03  3.8460E-04 -1.2706E-03
            -1.0401E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1292
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8392E-03 -6.0155E-03 -3.4660E-02 -4.0367E-03 -2.6778E-02
 SE:             2.9351E-02  2.5449E-02  1.6517E-02  1.2001E-02  2.3976E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5004E-01  8.1314E-01  3.5868E-02  7.3660E-01  2.6406E-01

 ETASHRINKSD(%)  1.6704E+00  1.4744E+01  4.4665E+01  5.9795E+01  1.9676E+01
 ETASHRINKVR(%)  3.3130E+00  2.7314E+01  6.9380E+01  8.3835E+01  3.5480E+01
 EBVSHRINKSD(%)  1.8015E+00  1.4918E+01  5.0474E+01  6.1671E+01  1.6557E+01
 EBVSHRINKVR(%)  3.5705E+00  2.7610E+01  7.5472E+01  8.5309E+01  3.0373E+01
 RELATIVEINF(%)  9.6317E+01  5.6282E+00  9.5325E+00  1.0617E+00  3.5361E+01
 EPSSHRINKSD(%)  1.6689E+01
 EPSSHRINKVR(%)  3.0592E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          879
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1615.4939413738146     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2711.0906159854362     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1095.5966746116217     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.99
 Elapsed covariance  time in seconds:    13.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2711.091       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.47E+00  2.75E+00  7.85E-01  1.68E+00  9.48E-01  9.18E-01  3.36E+00  4.85E-01  1.53E+00  2.74E+00
 


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
 
         3.19E-02  4.68E-01  8.81E-01  3.06E-01  2.49E-01  7.65E-02  2.27E-01  7.42E-01  1.61E-01  3.29E-01  2.93E-01
 


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
+        1.01E-03
 
 TH 2
+        2.44E-03  2.19E-01
 
 TH 3
+       -7.38E-03 -3.37E-01  7.77E-01
 
 TH 4
+       -1.67E-03 -1.43E-01  2.23E-01  9.39E-02
 
 TH 5
+        1.05E-03  1.10E-01 -1.71E-01 -7.21E-02  6.20E-02
 
 TH 6
+       -2.12E-04 -4.80E-03  1.26E-02  3.24E-03 -2.81E-03  5.85E-03
 
 TH 7
+       -1.56E-03 -9.69E-02  1.45E-01  6.33E-02 -4.99E-02  2.86E-03  5.13E-02
 
 TH 8
+       -5.66E-04  9.04E-02  8.67E-02 -5.75E-02  4.07E-02  1.48E-02 -3.94E-02  5.51E-01
 
 TH 9
+        9.10E-04  7.03E-03 -2.64E-02 -4.85E-03  3.71E-03 -6.97E-04 -7.13E-03 -2.70E-02  2.60E-02
 
 TH10
+        1.29E-03  1.41E-01 -2.35E-01 -9.29E-02  7.11E-02 -5.55E-03 -6.36E-02  2.15E-02  7.28E-03  1.08E-01
 
 TH11
+        3.52E-03  6.54E-02 -1.03E-01 -4.19E-02  2.55E-02 -2.57E-03 -3.14E-02 -1.69E-02 -6.45E-04  3.97E-02  8.61E-02
 
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
+        3.19E-02
 
 TH 2
+        1.64E-01  4.68E-01
 
 TH 3
+       -2.63E-01 -8.16E-01  8.81E-01
 
 TH 4
+       -1.71E-01 -9.95E-01  8.27E-01  3.06E-01
 
 TH 5
+        1.33E-01  9.45E-01 -7.80E-01 -9.45E-01  2.49E-01
 
 TH 6
+       -8.71E-02 -1.34E-01  1.87E-01  1.38E-01 -1.48E-01  7.65E-02
 
 TH 7
+       -2.16E-01 -9.14E-01  7.28E-01  9.11E-01 -8.85E-01  1.65E-01  2.27E-01
 
 TH 8
+       -2.40E-02  2.60E-01  1.32E-01 -2.53E-01  2.20E-01  2.61E-01 -2.34E-01  7.42E-01
 
 TH 9
+        1.77E-01  9.31E-02 -1.85E-01 -9.80E-02  9.23E-02 -5.64E-02 -1.95E-01 -2.26E-01  1.61E-01
 
 TH10
+        1.23E-01  9.15E-01 -8.10E-01 -9.21E-01  8.67E-01 -2.20E-01 -8.52E-01  8.80E-02  1.37E-01  3.29E-01
 
 TH11
+        3.77E-01  4.76E-01 -3.97E-01 -4.66E-01  3.49E-01 -1.14E-01 -4.72E-01 -7.76E-02 -1.36E-02  4.11E-01  2.93E-01
 
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
+        1.39E+03
 
 TH 2
+        1.15E+02  4.86E+02
 
 TH 3
+        3.39E+01  2.58E+00  7.94E+00
 
 TH 4
+        3.76E+01  6.38E+02 -2.46E+01  1.18E+03
 
 TH 5
+       -3.73E+01 -6.54E+01 -4.28E+00  7.39E+01  2.00E+02
 
 TH 6
+        4.94E+01 -5.45E+00  3.80E+00 -2.46E+00  1.57E+00  2.02E+02
 
 TH 7
+        4.27E+01  3.06E+01  3.64E+00 -5.15E+00  4.21E+01 -1.51E+01  1.47E+02
 
 TH 8
+       -2.01E+01 -6.88E+00 -3.87E+00  1.45E+01  8.27E+00 -7.73E+00  2.84E+00  4.99E+00
 
 TH 9
+       -4.99E+01 -1.62E+00 -7.56E-01  4.27E+00  1.55E+01 -8.57E+00  2.61E+01  4.40E+00  4.97E+01
 
 TH10
+        2.61E+01 -9.01E+00 -2.90E-01  7.53E+01  1.90E+01  1.65E+01  1.47E+01  5.91E+00  1.69E+00  7.73E+01
 
 TH11
+       -7.37E+01 -2.73E+01 -3.79E+00  5.34E+00  3.12E+01 -3.72E+00  1.17E+01  5.14E+00  9.57E+00  7.75E+00  2.58E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       37.533
Stop Time:
Sat Sep 25 02:01:10 CDT 2021
