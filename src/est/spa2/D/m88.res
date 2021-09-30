Thu Sep 30 10:00:17 CDT 2021
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
$DATA ../../../../data/spa2/D/dat88.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   32295.1363241093        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.4001E+02  6.4493E+02  2.2301E+01  6.3248E+02  1.8541E+02 -2.9429E+03 -1.4536E+03 -8.0487E+01 -2.0091E+03 -8.0491E+02
            -6.1419E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -419.508836357044        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.3055E+00  1.2739E+00  9.7826E-01  1.3950E+00  1.0018E+00  1.9639E+00  1.3788E+00  9.8945E-01  1.2100E+00  1.0295E+00
             1.4482E+01
 PARAMETER:  3.6662E-01  3.4208E-01  7.8022E-02  4.3287E-01  1.0181E-01  7.7494E-01  4.2125E-01  8.9396E-02  2.9065E-01  1.2908E-01
             2.7729E+00
 GRADIENT:   1.4857E+01  1.7668E+01  2.3554E+00  4.5196E+01 -7.4156E+00  2.8257E+01 -3.4715E+01  2.2848E+00 -2.1965E+01  8.1594E+00
             1.0406E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -508.265067598275        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.2705E+00  1.2621E+00  2.4455E+00  1.4869E+00  1.3646E+00  2.5503E+00  5.5197E+00  6.6447E-01  1.8147E+00  4.5221E-01
             1.3526E+01
 PARAMETER:  3.3939E-01  3.3275E-01  9.9426E-01  4.9668E-01  4.1088E-01  1.0362E+00  1.8083E+00 -3.0877E-01  6.9593E-01 -6.9360E-01
             2.7046E+00
 GRADIENT:   3.0109E+00  8.0970E+00 -7.1881E+00  8.0023E+00 -1.4157E+01  3.8720E+01  4.7411E+01  5.2298E-01  3.0219E+01  2.0089E+00
             1.3392E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -515.242180558780        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.2455E+00  1.3341E+00  2.7099E+00  1.3462E+00  1.4645E+00  2.5204E+00  4.8767E+00  4.9313E-01  1.3852E+00  7.1319E-01
             1.3327E+01
 PARAMETER:  3.1955E-01  3.8826E-01  1.0969E+00  3.9728E-01  4.8150E-01  1.0244E+00  1.6845E+00 -6.0698E-01  4.2584E-01 -2.3800E-01
             2.6898E+00
 GRADIENT:   1.2655E-01  3.5910E+00 -5.6529E-01  5.1156E+00 -1.6938E+01  3.8663E+01  2.0417E+01  2.1461E-01  1.4487E+01  4.8773E+00
             1.2535E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -519.554083271280        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.2283E+00  1.3898E+00  3.2056E+00  1.2453E+00  1.6499E+00  2.4503E+00  4.5989E+00  6.1065E-01  1.1954E+00  7.3890E-01
             1.3014E+01
 PARAMETER:  3.0560E-01  4.2916E-01  1.2649E+00  3.1941E-01  6.0071E-01  9.9622E-01  1.6258E+00 -3.9324E-01  2.7846E-01 -2.0259E-01
             2.6661E+00
 GRADIENT:  -1.4127E+00 -8.1798E-01 -2.2081E-01  3.4386E-01 -9.5020E+00  2.8731E+01  1.0218E+01  2.0508E-01  9.7163E+00  4.8442E+00
             1.1232E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -532.057520612355        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      484
 NPARAMETR:  1.1770E+00  1.1287E+00  7.6320E+00  1.2575E+00  2.2524E+00  2.1572E+00  5.6121E+00  1.0014E+00  6.1646E-01  4.7573E-01
             1.2279E+01
 PARAMETER:  2.6297E-01  2.2105E-01  2.1324E+00  3.2909E-01  9.1198E-01  8.6882E-01  1.8249E+00  1.0138E-01 -3.8375E-01 -6.4290E-01
             2.6079E+00
 GRADIENT:   7.3206E+00 -5.1298E+00 -4.7194E-01 -2.5560E+00 -2.7784E-02  1.7933E+01 -2.4855E+00  8.0380E-02 -1.5942E+00  1.1695E+00
             4.5705E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -533.609639259287        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  1.1469E+00  1.1714E+00  1.0389E+01  1.3178E+00  2.3106E+00  2.0035E+00  5.6994E+00  6.8116E-01  7.7419E-01  2.8584E-01
             1.2280E+01
 PARAMETER:  2.3702E-01  2.5818E-01  2.4408E+00  3.7594E-01  9.3751E-01  7.9489E-01  1.8404E+00 -2.8396E-01 -1.5594E-01 -1.1523E+00
             2.6080E+00
 GRADIENT:  -3.7158E+00  5.3038E-01 -3.1233E-01  2.4573E+00 -1.9775E-01 -1.6966E+00  1.5064E-01  1.9617E-02 -7.1310E-01  4.1184E-01
            -7.9894E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -533.888518198051        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      836
 NPARAMETR:  1.1541E+00  1.0464E+00  1.4654E+01  1.3923E+00  2.3545E+00  1.9967E+00  5.8762E+00  4.0553E-01  9.1467E-01  1.6229E-01
             1.2313E+01
 PARAMETER:  2.4333E-01  1.4534E-01  2.7847E+00  4.3093E-01  9.5634E-01  7.9149E-01  1.8709E+00 -8.0256E-01  1.0804E-02 -1.7184E+00
             2.6106E+00
 GRADIENT:  -1.5391E+00 -1.0315E+00 -2.2533E-01 -8.0135E-01  5.0688E-01 -2.0850E+00 -6.2987E-01  3.6565E-03  5.4103E-01  1.3229E-01
             1.7442E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -534.053100578027        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1019             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1589E+00  1.0929E+00  3.5730E+01  1.3846E+00  2.4153E+00  2.0200E+00  5.9754E+00  3.3986E-01  8.9393E-01  1.0000E-02
             1.2297E+01
 PARAMETER:  2.4743E-01  1.8882E-01  3.6760E+00  4.2540E-01  9.8182E-01  8.0308E-01  1.8876E+00 -9.7922E-01 -1.2132E-02 -5.0914E+00
             2.6094E+00
 GRADIENT:   1.0146E+01  2.1280E+00  1.2823E-03  5.6177E+00  2.2166E-01  1.2825E+01  5.2689E+01  3.7026E-04  3.9542E-01  0.0000E+00
             3.1803E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -534.069853517987        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  1.1586E+00  1.0722E+00  3.7345E+01  1.3843E+00  2.4268E+00  2.0169E+00  5.9291E+00  3.3252E-01  8.8789E-01  1.0000E-02
             1.2291E+01
 PARAMETER:  2.4719E-01  1.6968E-01  3.7202E+00  4.2521E-01  9.8657E-01  8.0154E-01  1.8799E+00 -1.0010E+00 -1.8911E-02 -5.1013E+00
             2.6089E+00
 GRADIENT:   2.0504E-01 -6.7760E-02 -4.9086E-03 -1.1909E+00 -2.1691E-02  3.8365E-01  2.2004E+00  3.1886E-04  4.3874E-02  0.0000E+00
             3.2166E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -534.075702593582        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1356
 NPARAMETR:  1.1592E+00  1.0676E+00  4.1624E+01  1.3911E+00  2.4330E+00  2.0177E+00  5.9474E+00  2.8561E-01  8.9102E-01  1.0000E-02
             1.2302E+01
 PARAMETER:  2.4777E-01  1.6540E-01  3.8287E+00  4.3012E-01  9.8914E-01  8.0196E-01  1.8830E+00 -1.1531E+00 -1.5391E-02 -5.1013E+00
             2.6097E+00
 GRADIENT:   2.2649E-01  9.1320E-02 -1.5591E-03 -3.8061E-01 -1.1870E-01  6.0176E-01  2.2132E+00  1.8953E-04 -1.4762E-01  0.0000E+00
             5.5055E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -534.076218561842        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1531
 NPARAMETR:  1.1600E+00  1.0642E+00  4.6012E+01  1.3970E+00  2.4393E+00  2.0162E+00  5.9392E+00  2.4341E-01  9.0242E-01  1.0000E-02
             1.2310E+01
 PARAMETER:  2.4840E-01  1.6225E-01  3.9289E+00  4.3431E-01  9.9173E-01  8.0120E-01  1.8816E+00 -1.3130E+00 -2.6700E-03 -5.1013E+00
             2.6104E+00
 GRADIENT:   2.8657E-01  5.8521E-02 -1.4715E-03 -5.3760E-02 -6.8295E-02  4.6487E-01  1.6068E+00  1.1264E-04 -7.8005E-02  0.0000E+00
             8.3926E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -534.085994907685        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1714
 NPARAMETR:  1.1591E+00  1.0527E+00  4.8408E+01  1.4000E+00  2.4439E+00  2.0153E+00  5.9916E+00  2.2636E-01  9.0738E-01  1.0000E-02
             1.2297E+01
 PARAMETER:  2.4766E-01  1.5132E-01  3.9797E+00  4.3648E-01  9.9361E-01  8.0078E-01  1.8904E+00 -1.3856E+00  2.8096E-03 -5.1013E+00
             2.6093E+00
 GRADIENT:   2.7437E-01  4.5582E-02 -3.5317E-03 -1.2537E+00  1.5525E-01  2.7395E-01  2.9157E+00  8.8335E-05  1.2087E-01  0.0000E+00
             3.8148E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -534.093869876023        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1896
 NPARAMETR:  1.1586E+00  1.0381E+00  5.3779E+01  1.4108E+00  2.4453E+00  2.0158E+00  6.0233E+00  2.1166E-01  9.1570E-01  1.0000E-02
             1.2295E+01
 PARAMETER:  2.4725E-01  1.3739E-01  4.0849E+00  4.4417E-01  9.9416E-01  8.0101E-01  1.8956E+00 -1.4528E+00  1.1938E-02 -5.1013E+00
             2.6092E+00
 GRADIENT:   5.4212E-02  9.2248E-02 -8.8726E-04 -2.4576E-01 -5.5383E-02  4.3391E-01  2.7612E+00  6.3383E-05 -1.1275E-01  0.0000E+00
            -3.5104E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -534.096487047844        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2082             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1587E+00  1.0321E+00  5.7792E+01  1.4163E+00  2.4476E+00  2.0143E+00  6.0438E+00  1.9595E-01  9.2079E-01  1.0000E-02
             1.2291E+01
 PARAMETER:  2.4733E-01  1.3162E-01  4.1569E+00  4.4803E-01  9.9511E-01  8.0028E-01  1.8990E+00 -1.5299E+00  1.7477E-02 -5.1013E+00
             2.6089E+00
 GRADIENT:   1.0094E+01  8.8458E-01  5.9116E-04  7.6517E+00  5.7900E-01  1.2331E+01  5.2654E+01  4.8329E-05 -8.4104E-02  0.0000E+00
             3.0202E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -534.097720549655        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2270             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1588E+00  1.0283E+00  6.1759E+01  1.4186E+00  2.4477E+00  2.0149E+00  6.0497E+00  1.9056E-01  9.2367E-01  1.0000E-02
             1.2295E+01
 PARAMETER:  2.4737E-01  1.2792E-01  4.2232E+00  4.4967E-01  9.9516E-01  8.0057E-01  1.9000E+00 -1.5578E+00  2.0595E-02 -5.1013E+00
             2.6092E+00
 GRADIENT:   1.0044E+01  8.1386E-01  2.1941E-03  7.6809E+00  5.1260E-01  1.2452E+01  5.2715E+01  4.0287E-05 -8.4592E-02  0.0000E+00
             3.0535E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -534.098776887084        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2446
 NPARAMETR:  1.1592E+00  1.0249E+00  6.1237E+01  1.4199E+00  2.4507E+00  2.0151E+00  6.0542E+00  1.7233E-01  9.2919E-01  1.0000E-02
             1.2302E+01
 PARAMETER:  2.4772E-01  1.2459E-01  4.2147E+00  4.5056E-01  9.9637E-01  8.0068E-01  1.9008E+00 -1.6583E+00  2.6556E-02 -5.1013E+00
             2.6098E+00
 GRADIENT:   1.4009E-01  1.9124E-02 -1.0413E-03 -5.5586E-01  2.5112E-02  4.1322E-01  2.9793E+00  3.2691E-05  1.3760E-02  0.0000E+00
             1.4418E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -534.100019179573        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2636             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1589E+00  1.0216E+00  6.5754E+01  1.4234E+00  2.4505E+00  2.0146E+00  6.0657E+00  1.6566E-01  9.3047E-01  1.0000E-02
             1.2297E+01
 PARAMETER:  2.4747E-01  1.2141E-01  4.2859E+00  4.5305E-01  9.9628E-01  8.0042E-01  1.9027E+00 -1.6978E+00  2.7932E-02 -5.1013E+00
             2.6093E+00
 GRADIENT:   1.0047E+01  7.5865E-01  1.6371E-03  7.7144E+00  5.5324E-01  1.2438E+01  5.3039E+01  2.7076E-05 -3.7298E-02  0.0000E+00
             3.0614E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -534.100461799649        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     2802
 NPARAMETR:  1.1590E+00  1.0194E+00  6.4806E+01  1.4235E+00  2.4523E+00  2.0146E+00  6.0678E+00  1.6703E-01  9.3383E-01  1.0000E-02
             1.2300E+01
 PARAMETER:  2.4760E-01  1.1924E-01  4.2714E+00  4.5311E-01  9.9704E-01  8.0043E-01  1.9030E+00 -1.6896E+00  3.1534E-02 -5.1013E+00
             2.6096E+00
 GRADIENT:   1.0406E-01  1.5606E-02 -7.9402E-04 -5.1899E-01  2.0775E-02  3.6099E-01  3.0569E+00  2.7759E-05  2.1686E-02  0.0000E+00
            -7.9078E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2802
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7793E-02  3.3711E-02  1.4994E-05 -6.4169E-02  2.5266E-06
 SE:             2.7612E-02  2.3565E-02  6.9192E-06  1.1684E-02  4.0835E-05
 N:                     100         100         100         100         100

 P VAL.:         5.1931E-01  1.5256E-01  3.0233E-02  3.9784E-08  9.5066E-01

 ETASHRINKSD(%)  7.4980E+00  2.1054E+01  9.9977E+01  6.0858E+01  9.9863E+01
 ETASHRINKVR(%)  1.4434E+01  3.7675E+01  1.0000E+02  8.4679E+01  1.0000E+02
 EBVSHRINKSD(%)  8.2465E+00  1.6030E+01  9.9956E+01  6.4418E+01  9.9792E+01
 EBVSHRINKVR(%)  1.5813E+01  2.9491E+01  1.0000E+02  8.7339E+01  1.0000E+02
 RELATIVEINF(%)  8.2079E+01  3.5149E+01  3.0937E-06  5.8122E+00  6.3237E-05
 EPSSHRINKSD(%)  5.6738E+00
 EPSSHRINKVR(%)  1.1026E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -534.10046179964888     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       568.62577804595821     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    64.69
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -534.100       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.16E+00  1.02E+00  6.48E+01  1.42E+00  2.45E+00  2.01E+00  6.07E+00  1.67E-01  9.34E-01  1.00E-02  1.23E+01
 


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
+        1.63E+02
 
 TH 2
+       -2.66E+00  2.12E+01
 
 TH 3
+       -4.29E-04  7.03E-05  3.48E-06
 
 TH 4
+       -8.43E+00  2.31E+01  3.77E-04  1.11E+02
 
 TH 5
+       -3.03E-01 -2.46E+00 -3.52E-03 -8.32E+00  7.78E+00
 
 TH 6
+       -3.52E+00 -8.23E-01  1.20E-04  2.92E+00 -1.70E-01  3.28E+01
 
 TH 7
+        7.97E-01  3.01E+00 -6.57E-05 -7.32E+00  2.39E-01 -2.96E-01  2.84E+00
 
 TH 8
+       -2.95E-02  4.18E-02 -7.44E-05 -2.37E-03 -1.67E-04  5.71E-03 -2.72E-04 -4.12E-02
 
 TH 9
+        6.52E-01 -2.38E+00 -1.04E-03 -3.58E+01  3.02E+00 -1.35E+00  2.07E+00  2.03E-01  2.07E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.11E+00 -2.05E+00  1.26E-04 -8.11E+00  4.28E-01  1.32E+00  3.80E-01  5.84E-04  2.10E+00  0.00E+00  4.09E+00
 
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
 #CPUT: Total CPU Time in Seconds,       76.048
Stop Time:
Thu Sep 30 10:01:45 CDT 2021
