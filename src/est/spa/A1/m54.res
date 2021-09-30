Wed Sep 29 12:12:05 CDT 2021
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
$DATA ../../../../data/spa/A1/dat54.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1246.29270258990        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6390E+02  8.1640E+01  4.4637E+01  9.4787E+01  2.8316E+01  1.5612E+01 -5.7134E+00 -2.7287E+01 -2.6986E+01 -1.3733E+01
            -7.8659E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1483.84142447195        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0681E+00  9.7238E-01  9.4580E-01  9.9878E-01  9.6780E-01  1.0306E+00  9.8437E-01  1.0662E+00  1.0373E+00  9.1931E-01
             2.2749E+00
 PARAMETER:  1.6592E-01  7.1991E-02  4.4276E-02  9.8783E-02  6.7268E-02  1.3013E-01  8.4246E-02  1.6414E-01  1.3664E-01  1.5872E-02
             9.2194E-01
 GRADIENT:   1.5750E+02 -1.1392E+01 -2.9387E+00 -1.0031E+01  3.4602E+00 -3.6529E-01  3.6946E+00 -3.7871E+00  2.8876E+00  1.5508E+01
             4.0562E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1487.71303441488        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0623E+00  1.2174E+00  6.9790E-01  8.6095E-01  9.2733E-01  1.1012E+00  1.0728E+00  1.4552E+00  1.0772E+00  5.8191E-01
             2.2023E+00
 PARAMETER:  1.6046E-01  2.9674E-01 -2.5967E-01 -4.9724E-02  2.4556E-02  1.9641E-01  1.7031E-01  4.7517E-01  1.7435E-01 -4.4143E-01
             8.8951E-01
 GRADIENT:   1.4707E+02  4.5700E+01 -4.8915E+00  2.9448E+01 -2.8214E+00  3.0739E+01  1.4752E+01  8.3932E+00  8.1980E+00  7.1591E+00
             2.3748E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1489.58591097540        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      315
 NPARAMETR:  1.0607E+00  1.1967E+00  7.4486E-01  8.5736E-01  9.3641E-01  1.0717E+00  8.7273E-01  1.3853E+00  1.1411E+00  5.2942E-01
             2.1623E+00
 PARAMETER:  1.5897E-01  2.7958E-01 -1.9456E-01 -5.3892E-02  3.4297E-02  1.6928E-01 -3.6127E-02  4.2588E-01  2.3198E-01 -5.3598E-01
             8.7119E-01
 GRADIENT:  -1.3841E+01 -4.1895E+00  3.6848E+00 -4.1780E+00 -1.2172E+01  3.1483E-02 -2.4039E+00  1.6243E+00  2.5215E-01  2.6207E+00
            -3.8162E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1494.22868035560        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      493
 NPARAMETR:  1.0643E+00  1.3172E+00  3.1342E-01  7.3667E-01  7.0849E-01  1.0649E+00  9.3219E-01  8.4808E-01  9.9798E-01  1.3910E-01
             2.1643E+00
 PARAMETER:  1.6234E-01  3.7554E-01 -1.0602E+00 -2.0561E-01 -2.4463E-01  1.6286E-01  2.9780E-02 -6.4776E-02  9.7978E-02 -1.8726E+00
             8.7208E-01
 GRADIENT:  -6.8794E+00  3.5331E+01  4.3153E+00  2.7119E+01 -1.2362E+01 -9.3975E-01  1.8590E+00 -1.5714E+00  2.4419E-01  6.2746E-01
             7.2354E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1526.92470135261        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  1.0585E+00  1.6004E+00  1.8352E-01  4.9964E-01  7.8246E-01  1.0594E+00  6.3693E-01  2.0434E+00  1.1363E+00  1.0000E-02
             1.9094E+00
 PARAMETER:  1.5681E-01  5.7023E-01 -1.5954E+00 -5.9386E-01 -1.4532E-01  1.5769E-01 -3.5109E-01  8.1464E-01  2.2779E-01 -4.6716E+00
             7.4680E-01
 GRADIENT:   9.4093E+00 -2.2484E+01  1.0708E+01 -1.0918E+01 -2.5203E+01  3.7184E+00 -5.1990E+01 -7.0173E+00 -1.3442E+01  0.0000E+00
             3.0523E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1535.19947348816        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      856             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0629E+00  1.6230E+00  2.1549E-01  5.0419E-01  8.3448E-01  1.0516E+00  8.0293E-01  2.6239E+00  1.1754E+00  1.0000E-02
             1.7836E+00
 PARAMETER:  1.6096E-01  5.8425E-01 -1.4348E+00 -5.8480E-01 -8.0948E-02  1.5032E-01 -1.1949E-01  1.0647E+00  2.6160E-01 -4.5616E+00
             6.7866E-01
 GRADIENT:   2.5806E+02  1.7957E+02  1.0931E+01  5.1258E+01  3.7195E+00  2.1230E+01  2.9484E+00  1.9044E+01 -1.2507E+00  0.0000E+00
             2.4389E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1536.40765741922        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      994
 NPARAMETR:  1.0511E+00  1.6460E+00  2.0863E-01  4.9736E-01  8.3421E-01  1.0512E+00  8.0063E-01  2.4193E+00  1.2711E+00  1.0000E-02
             1.6681E+00
 PARAMETER:  1.4984E-01  5.9837E-01 -1.4672E+00 -5.9843E-01 -8.1265E-02  1.4993E-01 -1.2236E-01  9.8350E-01  3.3985E-01 -4.5616E+00
             6.1166E-01
 GRADIENT:   2.5317E+02  2.4763E+02  1.5233E+01  6.2442E+01  7.6468E+00  2.3638E+01  5.3631E+00  9.6548E+00  1.1253E+00  0.0000E+00
            -3.2823E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1536.53954046509        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1066
 NPARAMETR:  1.0451E+00  1.6615E+00  1.8893E-01  4.7838E-01  8.2841E-01  1.0499E+00  7.9105E-01  2.2995E+00  1.3505E+00  1.0000E-02
             1.6105E+00
 PARAMETER:  1.4410E-01  6.0770E-01 -1.5664E+00 -6.3735E-01 -8.8249E-02  1.4867E-01 -1.3439E-01  9.3271E-01  4.0049E-01 -4.5616E+00
             5.7652E-01
 GRADIENT:   2.5301E+02  2.7704E+02  1.5898E+01  6.2544E+01  7.1072E+00  2.4584E+01  5.8185E+00  6.4920E+00  4.0172E+00  0.0000E+00
            -1.7793E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1536.74982388499        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1209
 NPARAMETR:  1.0537E+00  1.6591E+00  1.8877E-01  4.7850E-01  8.2873E-01  1.0574E+00  7.7276E-01  2.2944E+00  1.3580E+00  1.0469E-02
             1.6252E+00
 PARAMETER:  1.5235E-01  6.0630E-01 -1.5672E+00 -6.3711E-01 -8.7865E-02  1.5578E-01 -1.5779E-01  9.3048E-01  4.0603E-01 -4.4593E+00
             5.8563E-01
 GRADIENT:   2.7934E+02  2.6519E+02  1.5413E+01  6.1412E+01  8.2683E+00  2.7694E+01  7.0209E-01  6.2021E+00  4.5364E+00  1.7497E-04
            -1.3937E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1536.98958997922        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1388
 NPARAMETR:  1.0482E+00  1.6701E+00  1.8761E-01  4.7763E-01  8.2872E-01  1.0543E+00  7.8953E-01  2.3009E+00  1.3572E+00  1.1467E-01
             1.6366E+00
 PARAMETER:  1.4712E-01  6.1287E-01 -1.5734E+00 -6.3891E-01 -8.7872E-02  1.5288E-01 -1.3632E-01  9.3328E-01  4.0540E-01 -2.0657E+00
             5.9263E-01
 GRADIENT:  -5.5142E+00  7.2317E+00  4.4820E+00  8.8839E+00 -7.0440E+00 -3.9694E-01  3.1854E+00 -4.4931E+00  4.3725E-01 -1.4920E-01
            -1.3357E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1537.05920919328        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1561
 NPARAMETR:  1.0517E+00  1.6662E+00  1.8824E-01  4.7833E-01  8.2896E-01  1.0551E+00  7.7795E-01  2.2956E+00  1.3556E+00  1.3573E-01
             1.6403E+00
 PARAMETER:  1.5042E-01  6.1053E-01 -1.5700E+00 -6.3746E-01 -8.7586E-02  1.5365E-01 -1.5109E-01  9.3101E-01  4.0425E-01 -1.8971E+00
             5.9486E-01
 GRADIENT:   7.9654E-01  9.3094E-01  4.3531E+00  8.3692E+00 -4.4683E+00  5.2169E-03 -3.7053E-01 -4.7393E+00  1.3360E-01 -1.9576E-01
            -1.2046E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1537.24480315124        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1745
 NPARAMETR:  1.0477E+00  1.6668E+00  1.8617E-01  4.7657E-01  8.2907E-01  1.0524E+00  7.7844E-01  2.3052E+00  1.3522E+00  1.5522E-01
             1.6547E+00
 PARAMETER:  1.4660E-01  6.1092E-01 -1.5811E+00 -6.4114E-01 -8.7455E-02  1.5105E-01 -1.5046E-01  9.3518E-01  4.0173E-01 -1.7629E+00
             6.0361E-01
 GRADIENT:  -6.9073E+00 -2.3802E+00  3.3953E+00  8.1647E+00 -4.6046E+00 -8.9499E-01  2.4740E-01 -3.6267E+00  2.0623E-01 -1.6361E-01
            -7.3891E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1537.34445272861        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1930
 NPARAMETR:  1.0522E+00  1.6667E+00  1.8773E-01  4.7096E-01  8.3142E-01  1.0549E+00  7.7760E-01  2.2919E+00  1.3481E+00  1.8514E-01
             1.6831E+00
 PARAMETER:  1.5093E-01  6.1087E-01 -1.5728E+00 -6.5299E-01 -8.4624E-02  1.5349E-01 -1.5154E-01  9.2938E-01  3.9870E-01 -1.5867E+00
             6.2062E-01
 GRADIENT:   1.0734E+00 -1.1450E+01  6.2812E+00 -6.8351E-01 -1.0031E+01  3.2284E-01  5.5426E-01 -4.2459E+00  1.9827E-01 -4.4525E-02
             6.4473E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1537.68354588183        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2114
 NPARAMETR:  1.0430E+00  1.6946E+00  1.8500E-01  4.5912E-01  8.5255E-01  1.0487E+00  7.6913E-01  2.3024E+00  1.3557E+00  2.7408E-01
             1.6813E+00
 PARAMETER:  1.4213E-01  6.2747E-01 -1.5874E+00 -6.7845E-01 -5.9523E-02  1.4757E-01 -1.6249E-01  9.3397E-01  4.0429E-01 -1.1943E+00
             6.1955E-01
 GRADIENT:  -1.7963E+01 -1.3974E+01  3.9829E+00  4.2098E+00  2.7596E+00 -2.1663E+00 -1.8562E-01 -5.5116E+00 -2.1558E+00  2.9650E-03
             1.2150E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1537.81385470068        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     2287
 NPARAMETR:  1.0516E+00  1.7045E+00  1.8519E-01  4.5966E-01  8.5287E-01  1.0539E+00  7.6492E-01  2.2989E+00  1.3913E+00  2.7288E-01
             1.6845E+00
 PARAMETER:  1.4882E-01  6.3369E-01 -1.5880E+00 -6.7803E-01 -5.9280E-02  1.5195E-01 -1.6872E-01  9.3350E-01  4.3091E-01 -1.1971E+00
             6.2070E-01
 GRADIENT:  -5.5273E+00  4.8346E+00 -2.2596E+01 -5.5760E+01 -2.2054E+02 -3.9902E-01 -3.2751E-01  4.3456E+01  1.4343E-01  3.6378E+01
            -7.1523E+01
 NUMSIGDIG:         1.4         2.6         2.4         2.4         2.3         1.9         1.8         2.4         2.2         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2287
 NO. OF SIG. DIGITS IN FINAL EST.:  1.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.4837E-03 -1.7461E-02 -2.0176E-02  1.2316E-02 -1.7756E-02
 SE:             2.9678E-02  2.6162E-02  1.7306E-02  2.1840E-02  7.4670E-03
 N:                     100         100         100         100         100

 P VAL.:         9.0656E-01  5.0451E-01  2.4369E-01  5.7279E-01  1.7413E-02

 ETASHRINKSD(%)  5.7663E-01  1.2352E+01  4.2022E+01  2.6834E+01  7.4985E+01
 ETASHRINKVR(%)  1.1499E+00  2.3179E+01  6.6385E+01  4.6467E+01  9.3742E+01
 EBVSHRINKSD(%)  1.0363E+00  1.2769E+01  4.4276E+01  2.6302E+01  7.5448E+01
 EBVSHRINKVR(%)  2.0619E+00  2.3908E+01  6.8948E+01  4.5687E+01  9.3972E+01
 RELATIVEINF(%)  9.6874E+01  1.2535E+01  1.0447E+01  9.6749E+00  1.0428E+00
 EPSSHRINKSD(%)  4.1214E+01
 EPSSHRINKVR(%)  6.5442E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1537.8138547006831     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -802.66302813694494     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1537.814       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.71E+00  1.85E-01  4.59E-01  8.53E-01  1.05E+00  7.64E-01  2.30E+00  1.39E+00  2.73E-01  1.68E+00
 


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
+        8.90E+02
 
 TH 2
+       -1.63E+01  5.04E+02
 
 TH 3
+       -2.77E+01  3.27E+02  1.41E+04
 
 TH 4
+       -1.18E+01  3.42E+02 -1.25E+03  1.28E+04
 
 TH 5
+        2.96E+02  1.14E+04 -1.55E+03  4.20E+04  1.52E+05
 
 TH 6
+        5.33E+00 -3.91E+00 -2.88E-01  1.30E+01  2.82E+02  1.73E+02
 
 TH 7
+        1.06E+00  2.83E+00 -2.70E+01  7.49E+00  2.54E+02 -1.10E+00  2.01E+02
 
 TH 8
+        1.82E+00 -4.72E+00 -1.77E+03  3.55E+01  6.38E+01  9.50E-01  1.32E+00  2.55E+02
 
 TH 9
+        1.85E+00 -1.08E+01 -1.29E+01 -8.16E+01 -2.13E+04 -7.57E-01  1.91E+01  1.88E+00  2.87E+01
 
 TH10
+       -4.39E+01 -3.13E+03 -2.13E+01 -1.09E+04 -3.97E+04 -3.62E+01 -2.59E+01  1.65E+00  5.61E+03  1.04E+04
 
 TH11
+       -1.03E+01 -5.34E+00  2.26E+02 -5.06E+01 -1.78E+02  2.39E+00  9.32E+00 -1.29E+01  7.95E+00  6.00E+00  1.13E+03
 
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
 #CPUT: Total CPU Time in Seconds,       37.313
Stop Time:
Wed Sep 29 12:12:43 CDT 2021
