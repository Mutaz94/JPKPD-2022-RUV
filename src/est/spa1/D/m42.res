Thu Sep 30 03:04:09 CDT 2021
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
$DATA ../../../../data/spa1/D/dat42.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16767.6539001001        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.1791E+02  3.1660E+02 -2.9362E+01  4.3380E+02  1.4558E+02 -1.2980E+03 -8.2511E+02 -2.6853E+01 -1.0056E+03 -2.9760E+02
            -3.3875E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -682.091703411482        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.4146E+00  1.1358E+00  1.0077E+00  1.4307E+00  1.0713E+00  1.6582E+00  1.2388E+00  9.8495E-01  1.2680E+00  9.9850E-01
             1.4839E+01
 PARAMETER:  4.4683E-01  2.2736E-01  1.0765E-01  4.5819E-01  1.6886E-01  6.0574E-01  3.1416E-01  8.4836E-02  3.3745E-01  9.8501E-02
             2.7973E+00
 GRADIENT:   1.4387E+01  1.8859E+01 -3.9912E+00  2.7904E+01 -7.1775E+00  3.6464E+01 -6.2522E-01  4.0109E+00  4.2996E+00  5.4659E+00
             2.3420E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -706.363346117377        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.3200E+00  6.3578E-01  1.5863E+00  1.8056E+00  1.5126E+00  1.6314E+00  3.5227E+00  5.1263E-01  1.2839E+00  1.8107E+00
             1.3242E+01
 PARAMETER:  3.7762E-01 -3.5290E-01  5.6137E-01  6.9090E-01  5.1381E-01  5.8942E-01  1.3592E+00 -5.6819E-01  3.4991E-01  6.9370E-01
             2.6834E+00
 GRADIENT:  -5.5877E+00  2.0832E+01  2.5947E+00  4.8009E+01 -1.0990E+01  1.2902E+01  1.1655E+01  5.5120E-01  1.0186E+01  3.7457E+00
             2.0019E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -748.023315299999        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.0494E+00  4.4465E-01  8.1357E-01  1.2923E+00  5.7332E+00  1.4848E+00  2.0994E+00  5.8129E-02  7.4807E-01  5.3664E+00
             1.0057E+01
 PARAMETER:  1.4818E-01 -7.1048E-01 -1.0633E-01  3.5645E-01  1.8463E+00  4.9527E-01  8.4164E-01 -2.7451E+00 -1.9026E-01  1.7802E+00
             2.4082E+00
 GRADIENT:  -1.8350E+01  1.2983E+01  1.5623E+01 -2.8031E+01 -8.8955E+00  2.2222E+01  3.3896E-01 -7.5001E-03 -1.8716E+01  1.8085E+00
            -1.7947E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -789.660845011942        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.3218E-01  2.0704E-01  2.9181E-01  1.1680E+00  2.0051E+01  1.6141E+00  1.5380E+00  1.0000E-02  1.1208E+00  4.4337E+00
             8.6048E+00
 PARAMETER:  2.9772E-02 -1.4748E+00 -1.1317E+00  2.5533E-01  3.0983E+00  5.7876E-01  5.3048E-01 -6.3602E+00  2.1400E-01  1.5892E+00
             2.2523E+00
 GRADIENT:  -3.7872E+00  4.1449E+01 -1.2632E+01  7.0436E+01 -1.9155E+00  6.0052E+00  4.1839E+00  0.0000E+00 -1.9862E+01  1.3369E-01
            -7.7056E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -820.595374999985        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      400
 NPARAMETR:  6.1334E-01  4.4175E-02  5.7856E-02  5.2849E-01  1.0020E+02  1.6807E+00  2.6218E-01  1.0000E-02  7.2359E-01  5.4592E+00
             9.6698E+00
 PARAMETER: -3.8884E-01 -3.0196E+00 -2.7498E+00 -5.3774E-01  4.7072E+00  6.1920E-01 -1.2387E+00 -1.3628E+01 -2.2353E-01  1.7973E+00
             2.3690E+00
 GRADIENT:  -1.1922E+01  1.9696E+01 -1.1830E+02  1.6167E+02 -1.1458E-01  5.4563E+01  3.1392E-01  0.0000E+00 -3.3534E+01  1.2029E-02
             1.8244E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -842.832435599632        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      580
 NPARAMETR:  5.3854E-01  1.5049E-02  3.9931E-02  3.7795E-01  1.5002E+02  1.3254E+00  1.3507E-01  1.0000E-02  8.6920E-01  4.6681E+00
             9.0549E+00
 PARAMETER: -5.1889E-01 -4.0964E+00 -3.1206E+00 -8.7300E-01  5.1107E+00  3.8175E-01 -1.9020E+00 -1.5542E+01 -4.0188E-02  1.6408E+00
             2.3033E+00
 GRADIENT:   1.5411E+01 -1.4579E-01 -1.8056E+01  2.1321E+01  1.4949E-02 -9.9572E+00  3.1082E-05  0.0000E+00 -3.1694E+00 -2.3890E-05
            -5.3392E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -843.661120954100        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      756
 NPARAMETR:  5.0441E-01  1.1751E-02  3.4898E-02  3.3885E-01  1.7831E+02  1.3530E+00  1.0904E-01  1.0000E-02  8.8200E-01  4.7222E+00
             9.0191E+00
 PARAMETER: -5.8436E-01 -4.3438E+00 -3.2553E+00 -9.8221E-01  5.2835E+00  4.0230E-01 -2.1161E+00 -1.6317E+01 -2.5560E-02  1.6523E+00
             2.2993E+00
 GRADIENT:  -5.3054E-01 -3.8106E-02  2.9395E-01 -8.5095E-02  4.2143E-03  5.0582E-01  1.4974E-05  0.0000E+00 -1.9230E-01  7.2197E-06
            -1.0027E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -843.672530749873        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      933            RESET HESSIAN, TYPE II
 NPARAMETR:  5.0551E-01  1.3772E-02  3.4798E-02  3.3851E-01  1.7620E+02  1.3504E+00  4.0708E-02  1.0000E-02  8.8309E-01  4.5861E+00
             9.0182E+00
 PARAMETER: -5.8218E-01 -4.1851E+00 -3.2582E+00 -9.8321E-01  5.2716E+00  4.0041E-01 -3.1013E+00 -1.6317E+01 -2.4325E-02  1.6230E+00
             2.2992E+00
 GRADIENT:   5.8770E+01  4.1806E-02  7.0418E+01  3.6274E+01  4.9145E-03  6.1337E+00  7.8067E-06  0.0000E+00  1.7172E-01  3.7525E-05
             2.2217E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -843.737640359340        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  4.9979E-01  1.4240E-02  3.3835E-02  3.3255E-01  2.9688E+01  1.3460E+00  2.3291E-02  1.0000E-02  8.8421E-01  3.5062E+00
             8.9979E+00
 PARAMETER: -5.9358E-01 -4.1517E+00 -3.2863E+00 -1.0010E+00  3.4907E+00  3.9713E-01 -3.6597E+00 -1.6317E+01 -2.3059E-02  1.3545E+00
             2.2970E+00
 GRADIENT:   5.7052E+01  1.8546E-01  6.9580E+01  4.1606E+01  2.2717E-02  5.4539E+00  3.5529E-06  0.0000E+00  1.7453E-01  7.9659E-04
             2.0766E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -843.738054286148        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  4.9838E-01  1.4051E-02  3.3696E-02  3.3173E-01  2.8790E+01  1.3453E+00  2.3280E-02  1.0000E-02  8.8448E-01  2.0265E+00
             8.9902E+00
 PARAMETER: -5.9639E-01 -4.1651E+00 -3.2904E+00 -1.0034E+00  3.4600E+00  3.9658E-01 -3.6602E+00 -1.6317E+01 -2.2761E-02  8.0630E-01
             2.2961E+00
 GRADIENT:  -2.7604E+00  1.1358E-01 -5.6092E+00  8.1970E+00  1.9878E-02 -7.5143E-01  1.5379E-06  0.0000E+00 -1.4346E-01  8.7535E-05
            -2.8206E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -843.754602483480        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1379
 NPARAMETR:  5.0310E-01  1.2760E-02  3.4267E-02  3.3467E-01  1.3485E+01  1.3492E+00  2.2988E-02  1.0000E-02  8.8253E-01  3.0366E+00
             9.0175E+00
 PARAMETER: -5.8697E-01 -4.2614E+00 -3.2736E+00 -9.9462E-01  2.7016E+00  3.9952E-01 -3.6728E+00 -1.6317E+01 -2.4964E-02  1.2108E+00
             2.2992E+00
 GRADIENT:  -1.0103E+00  4.9852E-02  5.1228E-01 -7.2382E-01  7.0614E-03 -1.8404E-01  1.0563E-06  0.0000E+00  1.1949E-01  1.7412E-03
             1.6483E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -843.759457832246        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1522
 NPARAMETR:  5.0088E-01  1.1521E-02  3.3951E-02  3.3370E-01  1.2957E+01  1.3478E+00  2.2928E-02  1.0000E-02  8.8234E-01  6.0156E-01
             8.9940E+00
 PARAMETER: -5.9138E-01 -4.3636E+00 -3.2828E+00 -9.9752E-01  2.6617E+00  3.9846E-01 -3.6754E+00 -1.6317E+01 -2.5181E-02 -4.0823E-01
             2.2966E+00
 GRADIENT:  -1.0766E+00  7.3326E-03 -4.9359E+00  6.6002E+00  2.8529E-02 -6.3664E-01  7.2273E-07  0.0000E+00 -2.8712E-01  4.3705E-05
            -2.7610E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -843.777624569529        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:     1672
 NPARAMETR:  4.9605E-01  1.0574E-02  3.3487E-02  3.2893E-01  1.1213E+01  1.3443E+00  1.0000E-02  1.0000E-02  8.8539E-01  5.8057E-01
             8.9915E+00
 PARAMETER: -6.0108E-01 -4.4494E+00 -3.2966E+00 -1.0119E+00  2.5171E+00  3.9588E-01 -8.9257E+00 -1.6317E+01 -2.1727E-02 -4.4375E-01
             2.2963E+00
 GRADIENT:  -6.1394E+00  3.5699E-02  2.9014E+00 -1.6779E+00 -1.1486E-01 -1.2196E+00  0.0000E+00  0.0000E+00  6.7460E-01  1.1612E-04
            -1.0337E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -843.803032904836        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1854
 NPARAMETR:  4.9885E-01  1.1124E-02  3.3537E-02  3.2960E-01  1.2008E+01  1.3480E+00  1.0000E-02  1.0000E-02  8.8110E-01  5.7585E-01
             9.0151E+00
 PARAMETER: -5.9544E-01 -4.3987E+00 -3.2951E+00 -1.0099E+00  2.5855E+00  3.9866E-01 -8.3315E+00 -1.6317E+01 -2.6580E-02 -4.5191E-01
             2.2989E+00
 GRADIENT:  -1.0717E+00  6.0774E-03 -2.9397E-01  2.5267E-02 -1.8417E-02 -1.0536E-01  0.0000E+00  0.0000E+00 -1.4533E-01  6.9177E-05
             1.1189E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -843.809747792352        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2039
 NPARAMETR:  4.9945E-01  1.1677E-02  3.3425E-02  3.2907E-01  1.2209E+01  1.3485E+00  1.0000E-02  1.0000E-02  8.8188E-01  5.5523E-01
             9.0158E+00
 PARAMETER: -5.9424E-01 -4.3501E+00 -3.2985E+00 -1.0115E+00  2.6022E+00  3.9896E-01 -8.3315E+00 -1.6317E+01 -2.5705E-02 -4.8836E-01
             2.2990E+00
 GRADIENT:   1.2665E+00  1.0108E-02 -2.3983E+00  1.5941E+00 -4.3306E-03  9.9279E-02  0.0000E+00  0.0000E+00 -8.4870E-02  7.0348E-05
            -4.4244E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -843.813363777279        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2229
 NPARAMETR:  4.9902E-01  1.1502E-02  3.3384E-02  3.2862E-01  1.2344E+01  1.3480E+00  1.0000E-02  1.0000E-02  8.8187E-01  5.1370E-01
             9.0161E+00
 PARAMETER: -5.9511E-01 -4.3652E+00 -3.2997E+00 -1.0128E+00  2.6132E+00  3.9866E-01 -8.3315E+00 -1.6317E+01 -2.5714E-02 -5.6612E-01
             2.2990E+00
 GRADIENT:   9.9538E-01  1.6646E-03 -1.7076E+00  7.9850E-01  2.6619E-03  3.5793E-02  0.0000E+00  0.0000E+00 -5.2395E-02  5.2186E-05
            -2.9162E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -843.823583047272        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     2398
 NPARAMETR:  4.9778E-01  1.0857E-02  3.3211E-02  3.2698E-01  1.2418E+01  1.3479E+00  1.0000E-02  1.0000E-02  8.8286E-01  2.4018E-01
             9.0129E+00
 PARAMETER: -5.9760E-01 -4.4229E+00 -3.3049E+00 -1.0179E+00  2.6191E+00  3.9856E-01 -8.3315E+00 -1.6317E+01 -2.4590E-02 -1.3264E+00
             2.2987E+00
 GRADIENT:   9.7878E-01 -1.9912E-02 -7.6482E-02 -1.4495E+00  9.9731E-03  9.9677E-02  0.0000E+00  0.0000E+00  2.0385E-01  9.4086E-06
            -3.2213E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -843.826069312758        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2590
 NPARAMETR:  4.9742E-01  1.1547E-02  3.3163E-02  3.2654E-01  1.2044E+01  1.3475E+00  1.0000E-02  1.0000E-02  8.8249E-01  2.1982E-01
             9.0188E+00
 PARAMETER: -5.9833E-01 -4.3613E+00 -3.3063E+00 -1.0192E+00  2.5886E+00  3.9826E-01 -8.3315E+00 -1.6317E+01 -2.5011E-02 -1.4150E+00
             2.2993E+00
 GRADIENT:  -4.1592E-02  1.7129E-02  6.0869E-01 -1.9511E+00 -2.7440E-02  1.1478E-01  0.0000E+00  0.0000E+00  1.7974E-01  1.2569E-05
             5.2287E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -843.832710538061        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2775             RESET HESSIAN, TYPE I
 NPARAMETR:  4.9731E-01  1.1193E-02  3.3031E-02  3.2612E-01  1.2522E+01  1.3472E+00  1.0000E-02  1.0000E-02  8.8189E-01  1.9022E-01
             9.0141E+00
 PARAMETER: -5.9854E-01 -4.3924E+00 -3.3103E+00 -1.0205E+00  2.6275E+00  3.9799E-01 -8.3315E+00 -1.6317E+01 -2.5688E-02 -1.5596E+00
             2.2988E+00
 GRADIENT:   6.0892E+01  2.1649E-02  7.3950E+01  3.6566E+01  2.8121E-02  6.0923E+00  0.0000E+00  0.0000E+00  2.0532E-01  1.6409E-05
             2.2383E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -843.834198095606        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     2938
 NPARAMETR:  4.9677E-01  1.1302E-02  3.3039E-02  3.2585E-01  1.2406E+01  1.3469E+00  1.0000E-02  1.0000E-02  8.8203E-01  1.8780E-01
             9.0162E+00
 PARAMETER: -5.9962E-01 -4.3828E+00 -3.3101E+00 -1.0213E+00  2.6182E+00  3.9781E-01 -8.3315E+00 -1.6317E+01 -2.5526E-02 -1.5724E+00
             2.2990E+00
 GRADIENT:   4.0454E-01  4.9968E-04 -6.2367E-01 -5.4944E-01  5.2992E-04  2.3141E-02  0.0000E+00  0.0000E+00  2.6179E-02  6.7927E-06
             5.1018E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2938
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8171E-03  8.5323E-07  7.4711E-05 -1.9134E-02 -6.8884E-06
 SE:             2.8789E-02  1.4244E-06  2.3782E-04  2.4802E-02  2.1323E-05
 N:                     100         100         100         100         100

 P VAL.:         9.4967E-01  5.4916E-01  7.5341E-01  4.4041E-01  7.4666E-01

 ETASHRINKSD(%)  3.5536E+00  9.9995E+01  9.9203E+01  1.6911E+01  9.9929E+01
 ETASHRINKVR(%)  6.9810E+00  1.0000E+02  9.9994E+01  3.0962E+01  1.0000E+02
 EBVSHRINKSD(%)  3.3964E+00  9.9994E+01  9.9260E+01  1.7445E+01  9.9928E+01
 EBVSHRINKVR(%)  6.6774E+00  1.0000E+02  9.9995E+01  3.1847E+01  1.0000E+02
 RELATIVEINF(%)  2.7582E+00  2.1319E-08  4.2862E-05  5.0194E-01  2.8251E-06
 EPSSHRINKSD(%)  1.1591E+01
 EPSSHRINKVR(%)  2.1838E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -843.83419809560633     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       75.104335109066369     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    58.19
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     9.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -843.834       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.97E-01  1.13E-02  3.30E-02  3.26E-01  1.24E+01  1.35E+00  1.00E-02  1.00E-02  8.82E-01  1.88E-01  9.02E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.90E+02
 
 TH 2
+       -1.42E+01  1.06E+00
 
 TH 3
+       -9.01E+03  6.74E+02  4.27E+05
 
 TH 4
+        1.07E+03 -8.02E+01 -5.09E+04  6.05E+03
 
 TH 5
+        8.23E-02 -6.15E-03 -3.90E+00  4.64E-01  3.56E-05
 
 TH 6
+        8.00E-01 -5.98E-02 -3.79E+01  4.51E+00  3.46E-04  3.36E-03
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.41E+01  1.06E+00  6.70E+02 -7.97E+01 -6.11E-03 -5.94E-02  0.00E+00  0.00E+00  1.05E+00
 
 TH10
+       -9.60E-05  7.17E-06  4.55E-03 -5.42E-04 -4.15E-08 -4.03E-07  0.00E+00  0.00E+00  7.13E-06  4.84E-11
 
 TH11
+       -4.87E+00  3.64E-01  2.31E+02 -2.75E+01 -2.11E-03 -2.05E-02  0.00E+00  0.00E+00  3.62E-01  2.46E-06  1.25E-01
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.29E+03
 
 TH 2
+       -7.42E+02  1.91E+03
 
 TH 3
+       -1.12E+04  8.55E+02  5.25E+05
 
 TH 4
+       -1.86E+02  3.46E+02 -6.24E+04  8.64E+03
 
 TH 5
+        4.07E-01 -8.00E-01 -4.80E+00  3.95E-01  1.87E-03
 
 TH 6
+        1.06E+00  2.71E+01 -5.04E+01 -3.26E+01  9.46E-03  9.16E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.05E+00 -1.74E+01  8.22E+02 -1.05E+02 -1.67E-02 -3.55E+00  0.00E+00  0.00E+00  1.11E+02
 
 TH10
+       -7.54E-04 -3.02E-03  6.02E-03  3.42E-03 -2.79E-05 -2.14E-03  0.00E+00  0.00E+00  1.39E-03 -2.64E-03
 
 TH11
+       -2.67E+01  8.07E+00  2.84E+02 -2.14E+01 -6.23E-03  1.69E+00  0.00E+00  0.00E+00  3.96E+00  1.36E-05  6.10E+00
 
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
+        2.32E+03
 
 TH 2
+       -6.71E+02  2.27E+02
 
 TH 3
+       -1.59E+04  1.28E+03  6.45E+05
 
 TH 4
+        2.30E+02  3.40E+02 -7.14E+04  9.23E+03
 
 TH 5
+        4.01E-01 -1.01E-01 -5.67E+00  4.51E-01  9.93E-05
 
 TH 6
+        1.89E+02 -1.08E+01 -3.14E+03  1.98E+02  4.18E-02  1.24E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.17E+02 -3.67E+01 -1.85E+02 -9.93E+01 -6.53E-04  2.55E+01  0.00E+00  0.00E+00  9.06E+01
 
 TH10
+        6.26E-04 -7.39E-05 -7.00E-03  1.25E-04  6.24E-08  3.22E-04  0.00E+00  0.00E+00  1.67E-04  1.64E-09
 
 TH11
+       -7.10E+01  2.75E+01  5.84E+02 -4.24E+01 -1.56E-02  7.27E+00  0.00E+00  0.00E+00 -1.45E+01  1.94E-05  1.04E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       67.948
Stop Time:
Thu Sep 30 03:05:18 CDT 2021
