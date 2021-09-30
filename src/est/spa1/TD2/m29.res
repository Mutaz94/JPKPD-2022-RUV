Thu Sep 30 02:00:44 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat29.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2044.10511103661        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7195E+02 -2.5026E+01 -4.7164E+01  3.4420E+01  6.7927E+01  2.2980E+01 -1.8549E+01  7.6808E+00 -1.4524E+01  1.3717E+01
            -8.9936E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2056.65872041235        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9105E-01  1.0533E+00  1.0825E+00  1.0262E+00  1.0166E+00  1.0758E+00  1.0809E+00  9.7217E-01  1.0784E+00  9.3122E-01
             1.1040E+00
 PARAMETER:  9.1005E-02  1.5193E-01  1.7924E-01  1.2585E-01  1.1642E-01  1.7310E-01  1.7777E-01  7.1772E-02  1.7548E-01  2.8738E-02
             1.9892E-01
 GRADIENT:   6.6932E+00  3.1542E-01 -1.8460E+01  2.4268E+01  2.1343E+01  2.1926E+00 -1.0035E+01  2.7678E+00 -1.7464E+00  3.2732E+00
            -4.9437E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2057.76099327415        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.9153E-01  9.9031E-01  1.1242E+00  1.0560E+00  9.9472E-01  1.0435E+00  1.3579E+00  9.1766E-01  1.0074E+00  8.9063E-01
             1.0956E+00
 PARAMETER:  9.1498E-02  9.0258E-02  2.1705E-01  1.5451E-01  9.4706E-02  1.4262E-01  4.0596E-01  1.4072E-02  1.0737E-01 -1.5828E-02
             1.9134E-01
 GRADIENT:   9.6551E+00  6.0436E+00 -4.1428E-01  3.9276E+00  2.6718E+00 -1.0071E+01  4.2105E+00  2.4808E-01 -3.0225E-01 -4.5820E-01
            -1.2250E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2058.31683900903        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.8456E-01  8.3536E-01  1.2472E+00  1.1595E+00  9.7466E-01  1.0723E+00  1.3841E+00  9.6468E-01  9.9026E-01  8.9879E-01
             1.1161E+00
 PARAMETER:  8.4442E-02 -7.9896E-02  3.2092E-01  2.4800E-01  7.4335E-02  1.6980E-01  4.2503E-01  6.4041E-02  9.0210E-02 -6.7078E-03
             2.0988E-01
 GRADIENT:  -1.5132E+00  5.9470E+00  2.4579E+00  6.1435E+00 -5.9166E+00  1.6999E+00 -1.9093E-01 -2.3507E-01  9.2221E-01 -4.2012E-01
             2.5620E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2058.76125492858        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.8115E-01  5.4279E-01  1.5430E+00  1.3599E+00  9.7707E-01  1.0645E+00  1.6070E+00  1.1451E+00  9.1600E-01  9.5376E-01
             1.1126E+00
 PARAMETER:  8.0966E-02 -5.1102E-01  5.3375E-01  4.0739E-01  7.6808E-02  1.6253E-01  5.7438E-01  2.3548E-01  1.2256E-02  5.2660E-02
             2.0668E-01
 GRADIENT:  -1.1429E+00  7.1493E+00  3.1055E+00  1.4615E+01 -6.7332E+00  1.1228E-01 -6.2382E-01 -3.3640E-01 -1.7525E+00  6.7618E-01
            -1.3573E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2059.23638261122        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.7883E-01  2.7940E-01  1.7557E+00  1.5328E+00  9.6278E-01  1.0576E+00  2.0706E+00  1.2940E+00  8.6161E-01  9.5770E-01
             1.1096E+00
 PARAMETER:  7.8607E-02 -1.1751E+00  6.6289E-01  5.2708E-01  6.2070E-02  1.5603E-01  8.2784E-01  3.5771E-01 -4.8958E-02  5.6775E-02
             2.0400E-01
 GRADIENT:   9.4525E-01  3.9661E+00  2.5697E+00  1.5280E+01 -4.8580E+00 -1.4880E+00 -5.0297E-01 -7.3855E-01 -1.6323E-02 -1.9015E-01
            -2.4368E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2059.38869864266        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.7727E-01  1.7755E-01  1.8424E+00  1.5998E+00  9.5811E-01  1.0584E+00  2.5085E+00  1.3585E+00  8.3672E-01  9.5658E-01
             1.1109E+00
 PARAMETER:  7.7011E-02 -1.6285E+00  7.1106E-01  5.6986E-01  5.7207E-02  1.5679E-01  1.0197E+00  4.0638E-01 -7.8261E-02  5.5611E-02
             2.0515E-01
 GRADIENT:   4.0812E-01  2.7063E+00  2.1035E+00  1.4371E+01 -4.3413E+00 -8.5887E-01 -3.9926E-01 -6.9230E-01 -8.4159E-01 -3.0432E-01
            -1.5349E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2059.42384892159        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  9.7624E-01  1.2263E-01  1.9023E+00  1.6380E+00  9.5809E-01  1.0598E+00  2.9528E+00  1.4108E+00  8.2287E-01  9.5922E-01
             1.1122E+00
 PARAMETER:  7.5953E-02 -1.9986E+00  7.4307E-01  5.9347E-01  5.7183E-02  1.5805E-01  1.1827E+00  4.4415E-01 -9.4961E-02  5.8366E-02
             2.0631E-01
 GRADIENT:  -3.5616E-01  2.1401E+00  1.7430E+00  1.6310E+01 -5.1964E+00 -2.1974E-01 -2.8970E-01 -9.4851E-02 -1.5601E+00  2.2540E-01
            -4.0097E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2059.50134050530        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1409
 NPARAMETR:  9.7533E-01  7.0066E-02  1.9568E+00  1.6723E+00  9.5833E-01  1.0608E+00  3.8229E+00  1.4528E+00  8.1124E-01  9.6189E-01
             1.1133E+00
 PARAMETER:  7.5023E-02 -2.5583E+00  7.7133E-01  6.1417E-01  5.7433E-02  1.5906E-01  1.4410E+00  4.7347E-01 -1.0919E-01  6.1140E-02
             2.0729E-01
 GRADIENT:  -8.7831E-01  1.2440E+00  1.1920E+00  1.3712E+01 -4.5752E+00  3.2309E-01 -1.5457E-01  2.5005E-01 -1.5136E+00  5.4683E-01
             5.2656E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2059.62765062670        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1586
 NPARAMETR:  9.7491E-01  2.5835E-02  1.9959E+00  1.6978E+00  9.5815E-01  1.0606E+00  6.1727E+00  1.4842E+00  8.0187E-01  9.6085E-01
             1.1133E+00
 PARAMETER:  7.4587E-02 -3.5560E+00  7.9112E-01  6.2936E-01  5.7251E-02  1.5879E-01  1.9201E+00  4.9488E-01 -1.2081E-01  6.0059E-02
             2.0730E-01
 GRADIENT:  -5.3462E-01  4.1739E-01  1.2846E-01  5.4346E+00 -1.9586E+00  3.5436E-01 -4.7129E-02  5.0918E-01 -1.2384E+00  3.8964E-01
             6.4684E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2059.68896729264        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1768             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7476E-01  1.0206E-02  2.0112E+00  1.7058E+00  9.5868E-01  1.0600E+00  1.0161E+01  1.4882E+00  8.0035E-01  9.6179E-01
             1.1130E+00
 PARAMETER:  7.4440E-02 -4.4848E+00  7.9871E-01  6.3402E-01  5.7804E-02  1.5824E-01  2.4186E+00  4.9759E-01 -1.2271E-01  6.1039E-02
             2.0708E-01
 GRADIENT:   3.5490E+02  1.6867E+00  9.1781E+00  1.1945E+03  4.5523E+00  6.7463E+01  5.8898E-01  1.4998E+00  1.9759E+01  7.9367E-01
             2.4858E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2059.76427682262        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1930
 NPARAMETR:  9.7862E-01  1.0000E-02  2.0060E+00  1.6937E+00  9.5838E-01  1.0614E+00  8.7138E+00  1.4819E+00  8.0173E-01  9.6001E-01
             1.1122E+00
 PARAMETER:  7.8386E-02 -4.5694E+00  7.9616E-01  6.2694E-01  5.7492E-02  1.5955E-01  2.2649E+00  4.9330E-01 -1.2098E-01  5.9188E-02
             2.0631E-01
 GRADIENT:   7.5214E+00  0.0000E+00  4.8672E-01 -2.3069E+01  1.6646E+00  7.3506E-01 -1.4165E-02  6.8667E-02  4.1824E-01 -1.7448E-01
            -1.2786E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2059.78219967631        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:     2052
 NPARAMETR:  9.7621E-01  1.0000E-02  1.9838E+00  1.6950E+00  9.5300E-01  1.0614E+00  8.4320E+00  1.4700E+00  8.0064E-01  9.5623E-01
             1.1124E+00
 PARAMETER:  7.5925E-02 -4.5689E+00  7.8501E-01  6.2769E-01  5.1864E-02  1.5959E-01  2.2320E+00  4.8527E-01 -1.2234E-01  5.5240E-02
             2.0656E-01
 GRADIENT:   3.5831E+02  0.0000E+00  9.1179E+00  1.1643E+03  6.2106E+00  6.8295E+01  3.8584E-01  1.4157E+00  1.9359E+01  4.5612E-01
             2.1862E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2059.79160358837        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     2215
 NPARAMETR:  9.7503E-01  1.0000E-02  1.9585E+00  1.6925E+00  9.4774E-01  1.0591E+00  7.8317E+00  1.4495E+00  8.0263E-01  9.5067E-01
             1.1122E+00
 PARAMETER:  7.4711E-02 -4.5689E+00  7.7220E-01  6.2620E-01  4.6321E-02  1.5747E-01  2.1582E+00  4.7119E-01 -1.1986E-01  4.9414E-02
             2.0634E-01
 GRADIENT:   5.0014E-01  0.0000E+00 -3.2459E-01 -2.0855E+01  1.9107E+00 -1.2180E-02 -1.6482E-02 -1.1601E-01  3.9237E-01 -5.4505E-01
            -1.4051E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2059.79711169897        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2411             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7540E-01  1.0000E-02  1.9569E+00  1.6912E+00  9.4556E-01  1.0600E+00  8.7392E+00  1.4487E+00  8.0184E-01  9.5574E-01
             1.1123E+00
 PARAMETER:  7.5096E-02 -4.5689E+00  7.7135E-01  6.2545E-01  4.4021E-02  1.5825E-01  2.2678E+00  4.7070E-01 -1.2085E-01  5.4726E-02
             2.0640E-01
 GRADIENT:   3.5678E+02  0.0000E+00  9.8061E+00  1.1549E+03  4.1656E+00  6.7346E+01  4.1721E-01  1.3429E+00  1.9444E+01  9.6371E-01
             2.1858E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2059.80376824914        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     2583
 NPARAMETR:  9.7537E-01  1.0000E-02  1.9549E+00  1.6911E+00  9.4618E-01  1.0599E+00  1.1316E+01  1.4473E+00  8.0178E-01  9.5281E-01
             1.1122E+00
 PARAMETER:  7.5065E-02 -4.5689E+00  7.7034E-01  6.2536E-01  4.4673E-02  1.5821E-01  2.5262E+00  4.6973E-01 -1.2092E-01  5.1656E-02
             2.0633E-01
 GRADIENT:   3.5675E+02  0.0000E+00  9.2109E+00  1.1549E+03  5.9604E+00  6.7305E+01  7.5352E-01  1.2219E+00  1.9585E+01  5.0563E-01
             2.0369E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2059.80455430183        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2766             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7556E-01  1.0000E-02  1.9518E+00  1.6909E+00  9.4553E-01  1.0601E+00  1.1850E+01  1.4452E+00  8.0139E-01  9.5258E-01
             1.1122E+00
 PARAMETER:  7.5254E-02 -4.5689E+00  7.6875E-01  6.2525E-01  4.3990E-02  1.5840E-01  2.5723E+00  4.6826E-01 -1.2141E-01  5.1418E-02
             2.0634E-01
 GRADIENT:   3.5714E+02  0.0000E+00  9.1202E+00  1.1546E+03  6.0125E+00  6.7429E+01  8.5873E-01  1.2231E+00  1.9487E+01  5.2673E-01
             2.0482E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2059.80488855252        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2955
 NPARAMETR:  9.7539E-01  1.0000E-02  1.9485E+00  1.6909E+00  9.4501E-01  1.0599E+00  1.1712E+01  1.4437E+00  8.0138E-01  9.5143E-01
             1.1122E+00
 PARAMETER:  7.5078E-02 -4.5689E+00  7.6706E-01  6.2526E-01  4.3443E-02  1.5821E-01  2.5606E+00  4.6722E-01 -1.2142E-01  5.0213E-02
             2.0634E-01
 GRADIENT:   1.4063E+00  0.0000E+00 -2.5902E-01 -2.3051E+01  1.0824E+00  3.3857E-01  1.8406E-02  4.9516E-02  3.5237E-02 -1.4263E-01
            -2.2676E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2059.80560328485        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3147             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7537E-01  1.0000E-02  1.9468E+00  1.6907E+00  9.4370E-01  1.0599E+00  1.1717E+01  1.4411E+00  8.0146E-01  9.5238E-01
             1.1122E+00
 PARAMETER:  7.5060E-02 -4.5689E+00  7.6617E-01  6.2516E-01  4.2050E-02  1.5820E-01  2.5611E+00  4.6542E-01 -1.2133E-01  5.1213E-02
             2.0636E-01
 GRADIENT:   3.5662E+02  0.0000E+00  9.4531E+00  1.1545E+03  5.0267E+00  6.7286E+01  8.2952E-01  1.2009E+00  1.9455E+01  6.7076E-01
             2.0745E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2059.80577002294        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3335
 NPARAMETR:  9.7536E-01  1.0000E-02  1.9454E+00  1.6906E+00  9.4332E-01  1.0599E+00  1.1725E+01  1.4403E+00  8.0149E-01  9.5219E-01
             1.1122E+00
 PARAMETER:  7.5054E-02 -4.5689E+00  7.6545E-01  6.2511E-01  4.1655E-02  1.5819E-01  2.5617E+00  4.6484E-01 -1.2129E-01  5.1008E-02
             2.0637E-01
 GRADIENT:   1.3652E+00  0.0000E+00  2.9688E-01 -2.3050E+01 -3.7151E-01  3.2309E-01  1.8935E-02  3.8595E-02  6.5109E-02  1.1518E-01
             2.6007E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2059.80589518702        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3527             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7536E-01  1.0000E-02  1.9422E+00  1.6905E+00  9.4333E-01  1.0599E+00  1.1716E+01  1.4387E+00  8.0149E-01  9.5024E-01
             1.1122E+00
 PARAMETER:  7.5052E-02 -4.5689E+00  7.6381E-01  6.2505E-01  4.1664E-02  1.5817E-01  2.5609E+00  4.6374E-01 -1.2129E-01  4.8956E-02
             2.0633E-01
 GRADIENT:   3.5664E+02  0.0000E+00  8.9493E+00  1.1541E+03  6.1520E+00  6.7233E+01  8.2990E-01  1.1572E+00  1.9431E+01  4.0047E-01
             2.0072E+00

0ITERATION NO.:  104    OBJECTIVE VALUE:  -2059.80615424009        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     3667
 NPARAMETR:  9.7535E-01  1.0000E-02  1.9403E+00  1.6904E+00  9.4287E-01  1.0599E+00  1.1723E+01  1.4369E+00  8.0154E-01  9.5043E-01
             1.1122E+00
 PARAMETER:  7.5049E-02 -4.5689E+00  7.6402E-01  6.2504E-01  4.1252E-02  1.5818E-01  2.5616E+00  4.6382E-01 -1.2126E-01  5.0159E-02
             2.0636E-01
 GRADIENT:   4.5245E-03  0.0000E+00  1.3269E-01  6.6841E-02  3.0023E-02  2.4779E-03  2.0268E-06  3.1499E-02 -3.7365E-03  3.5538E-02
             1.0498E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3667
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4435E-04 -2.4779E-04 -2.9314E-02 -5.1696E-03 -3.7801E-02
 SE:             2.9828E-02  1.7545E-03  1.7993E-02  2.9417E-02  2.0306E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9346E-01  8.8769E-01  1.0327E-01  8.6050E-01  6.2669E-02

 ETASHRINKSD(%)  7.1857E-02  9.4122E+01  3.9721E+01  1.4507E+00  3.1971E+01
 ETASHRINKVR(%)  1.4366E-01  9.9655E+01  6.3665E+01  2.8803E+00  5.3720E+01
 EBVSHRINKSD(%)  3.8926E-01  9.4230E+01  4.3073E+01  1.7894E+00  2.9399E+01
 EBVSHRINKVR(%)  7.7700E-01  9.9667E+01  6.7593E+01  3.5467E+00  5.0155E+01
 RELATIVEINF(%)  9.5619E+01  8.4365E-03  6.7866E+00  2.8287E+00  8.5340E+00
 EPSSHRINKSD(%)  3.3343E+01
 EPSSHRINKVR(%)  5.5569E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2059.8061542400906     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1140.8676210354179     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    60.67
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:     7.49
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2059.806       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.00E-02  1.94E+00  1.69E+00  9.43E-01  1.06E+00  1.17E+01  1.44E+00  8.02E-01  9.51E-01  1.11E+00
 


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
 
         3.10E-02  0.00E+00  3.60E-01  4.38E-02  9.04E-02  8.80E-02  1.40E+00  3.49E-01  5.29E-02  1.41E-01  5.56E-02
 


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
+        9.64E-04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.22E-03  0.00E+00  1.30E-01
 
 TH 4
+       -1.09E-04  0.00E+00  5.43E-03  1.92E-03
 
 TH 5
+       -3.32E-04  0.00E+00  2.88E-02  1.34E-03  8.18E-03
 
 TH 6
+       -3.43E-04  0.00E+00 -2.02E-04  2.92E-04 -5.67E-04  7.74E-03
 
 TH 7
+       -8.19E-03  0.00E+00  5.74E-02 -6.52E-03  1.99E-02 -1.18E-02  1.95E+00
 
 TH 8
+       -2.43E-03  0.00E+00  8.36E-02  1.83E-03  1.66E-02  2.42E-03  4.20E-02  1.22E-01
 
 TH 9
+        3.04E-04  0.00E+00 -2.11E-03 -5.21E-04 -3.04E-04 -2.11E-04  1.68E-02 -2.90E-03  2.80E-03
 
 TH10
+       -3.75E-04  0.00E+00  1.30E-02  9.97E-04  5.39E-03  1.30E-03  1.05E-02 -8.41E-03 -4.82E-04  1.98E-02
 
 TH11
+        2.80E-04  0.00E+00  5.14E-04 -1.35E-04  1.40E-04 -3.84E-04  1.03E-02 -1.03E-03  3.58E-04 -5.81E-04  3.09E-03
 
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
+        3.10E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.09E-01  0.00E+00  3.60E-01
 
 TH 4
+       -8.03E-02  0.00E+00  3.44E-01  4.38E-02
 
 TH 5
+       -1.18E-01  0.00E+00  8.85E-01  3.39E-01  9.04E-02
 
 TH 6
+       -1.26E-01  0.00E+00 -6.36E-03  7.57E-02 -7.12E-02  8.80E-02
 
 TH 7
+       -1.89E-01  0.00E+00  1.14E-01 -1.06E-01  1.58E-01 -9.58E-02  1.40E+00
 
 TH 8
+       -2.24E-01  0.00E+00  6.65E-01  1.20E-01  5.26E-01  7.89E-02  8.63E-02  3.49E-01
 
 TH 9
+        1.85E-01  0.00E+00 -1.10E-01 -2.25E-01 -6.35E-02 -4.53E-02  2.27E-01 -1.57E-01  5.29E-02
 
 TH10
+       -8.58E-02  0.00E+00  2.56E-01  1.62E-01  4.24E-01  1.05E-01  5.34E-02 -1.71E-01 -6.48E-02  1.41E-01
 
 TH11
+        1.62E-01  0.00E+00  2.56E-02 -5.53E-02  2.79E-02 -7.85E-02  1.32E-01 -5.33E-02  1.22E-01 -7.43E-02  5.56E-02
 
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
+        1.09E+03
 
 TH 2
+       -1.94E-14 -1.19E-29
 
 TH 3
+       -1.75E+01 -2.04E-14  3.49E+01
 
 TH 4
+       -2.35E+01  2.12E-14 -1.24E+01  6.16E+02
 
 TH 5
+        6.33E+01  9.11E-14 -1.56E+02 -6.36E+01  7.36E+02
 
 TH 6
+        5.05E+01  4.54E-14 -1.59E+01 -2.27E+01  6.14E+01  1.32E+02
 
 TH 7
+       -1.50E-02 -6.30E-18  3.68E-03  1.00E-02 -1.27E-02  3.46E-04  8.17E-06
 
 TH 8
+       -1.07E+00 -1.76E-15 -3.02E-01 -3.03E+00 -2.20E-01  2.67E+00 -4.94E-04  4.66E-01
 
 TH 9
+       -1.10E+02 -3.83E-14  1.82E+01  9.29E+01 -5.76E+01  1.07E+00  5.62E-02 -3.78E+00  3.89E+02
 
 TH10
+       -6.31E+00 -1.09E-14  1.62E+01  1.58E+00 -7.89E+01  8.70E-01  6.52E-04  7.37E-01  6.72E-01  9.64E+00
 
 TH11
+       -6.78E+01 -7.99E-14 -6.42E+00  3.14E+01 -2.49E+01  4.01E+00 -2.71E-03  1.14E+01 -2.38E+01  1.81E+01  3.43E+02
 
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
+        1.03E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        1.81E+01  0.00E+00  4.76E+01
 
 TH 4
+       -3.27E+01  0.00E+00 -7.91E+00  5.84E+02
 
 TH 5
+       -1.31E+01  0.00E+00 -1.49E+02 -5.19E+01  7.14E+02
 
 TH 6
+       -5.55E+01  0.00E+00  5.45E+00  3.10E+01 -8.80E+01  2.32E+02
 
 TH 7
+       -4.04E-03  0.00E+00  4.91E-04 -2.41E-02  1.07E-02 -2.81E-03  1.19E-05
 
 TH 8
+       -3.04E+01  0.00E+00 -1.30E+01 -1.98E+01 -1.81E+00  1.04E+01  4.24E-04  1.83E+01
 
 TH 9
+        8.81E+01  0.00E+00  4.06E+00 -7.66E+01  4.45E+01 -1.07E+01  4.10E-02 -7.66E+00  2.34E+02
 
 TH10
+       -2.70E+01  0.00E+00 -2.32E+00 -7.12E-01 -6.26E+01  2.99E+01  1.57E-03  1.16E+01 -1.19E+01  4.93E+01
 
 TH11
+        6.10E+01  0.00E+00 -5.04E+00 -5.62E+01  1.28E+01 -8.55E+00  1.09E-02  7.95E+00  3.70E+01  3.31E+00  3.25E+02
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       68.252
Stop Time:
Thu Sep 30 02:01:53 CDT 2021
