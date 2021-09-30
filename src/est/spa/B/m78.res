Wed Sep 29 11:31:57 CDT 2021
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
$DATA ../../../../data/spa/B/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1672.11352374190        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5790E+02 -7.7684E+01 -6.2143E+01 -2.4751E+01  4.1245E+01  3.2090E+01 -1.5347E+01  1.8214E+01  1.1754E+01  6.5587E+00
             2.0666E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1690.56767072704        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.9233E-01  1.3475E+00  1.7538E+00  8.7275E-01  1.3947E+00  9.9685E-01  1.2299E+00  7.5991E-01  8.9145E-01  1.0714E+00
             9.4429E-01
 PARAMETER:  9.2303E-02  3.9826E-01  6.6179E-01 -3.6110E-02  4.3265E-01  9.6848E-02  3.0697E-01 -1.7456E-01 -1.4906E-02  1.6894E-01
             4.2676E-02
 GRADIENT:   4.5179E+01  1.1329E+01  7.3731E+00 -7.8876E+00  3.4367E+01 -1.3273E+00  1.9274E+01 -1.8115E+00 -5.1755E+00 -4.8991E+01
            -2.3280E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1696.29701723054        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.6952E-01  1.3641E+00  2.4238E+00  8.3986E-01  1.4811E+00  1.0015E+00  9.7615E-01  5.9869E-01  1.0589E+00  1.3011E+00
             1.0116E+00
 PARAMETER:  6.9043E-02  4.1049E-01  9.8533E-01 -7.4517E-02  4.9280E-01  1.0148E-01  7.5857E-02 -4.1301E-01  1.5719E-01  3.6319E-01
             1.1158E-01
 GRADIENT:  -8.7038E+00 -3.2328E+01  4.8206E+00 -3.6927E+01  4.0642E+00  1.6440E+00  2.3785E+00 -2.6932E-01 -5.7350E+00 -9.9607E+00
             1.1251E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1698.35808425129        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.7384E-01  1.3803E+00  1.8566E+00  8.5341E-01  1.4117E+00  1.0078E+00  9.5782E-01  3.2109E-01  1.1117E+00  1.3231E+00
             9.7612E-01
 PARAMETER:  7.3495E-02  4.2229E-01  7.1875E-01 -5.8514E-02  4.4483E-01  1.0782E-01  5.6906E-02 -1.0360E+00  2.0588E-01  3.8001E-01
             7.5832E-02
 GRADIENT:   2.0843E+00  2.3127E+00  2.9338E+00 -6.5154E-01 -2.9962E+00  4.0423E+00  2.2684E+00 -8.5695E-02 -8.5249E-01 -9.6352E-02
             6.8332E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1698.52404858335        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.7335E-01  1.4216E+00  1.6841E+00  8.2540E-01  1.3954E+00  9.9498E-01  9.0746E-01  2.7055E-01  1.1561E+00  1.3081E+00
             9.7099E-01
 PARAMETER:  7.2985E-02  4.5181E-01  6.2123E-01 -9.1881E-02  4.3321E-01  9.4964E-02  2.8948E-03 -1.2073E+00  2.4505E-01  3.6856E-01
             7.0562E-02
 GRADIENT:   5.9790E-01  2.4746E+00  5.4107E-01  2.9687E+00 -2.2178E+00 -1.0310E+00 -7.4293E-01 -9.6739E-03 -7.8328E-01  1.0564E-01
            -5.0176E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1698.52429761407        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.7349E-01  1.4511E+00  1.6700E+00  8.0642E-01  1.4038E+00  9.9488E-01  8.8988E-01  2.7129E-01  1.1771E+00  1.3124E+00
             9.7111E-01
 PARAMETER:  7.3133E-02  4.7232E-01  6.1281E-01 -1.1515E-01  4.3915E-01  9.4864E-02 -1.6666E-02 -1.2046E+00  2.6305E-01  3.7185E-01
             7.0684E-02
 GRADIENT:   6.2467E-01  3.2528E+00  5.7058E-01  3.6419E+00 -2.3868E+00 -1.1241E+00 -1.2466E+00 -9.0339E-03 -1.1230E+00  1.8830E-01
            -6.3323E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1698.54834049689        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  9.7369E-01  1.4550E+00  1.6582E+00  7.9904E-01  1.4076E+00  9.9773E-01  8.9911E-01  2.9351E-01  1.1897E+00  1.3126E+00
             9.7184E-01
 PARAMETER:  7.3341E-02  4.7500E-01  6.0573E-01 -1.2435E-01  4.4186E-01  9.7727E-02 -6.3465E-03 -1.1258E+00  2.7367E-01  3.7200E-01
             7.1438E-02
 GRADIENT:   1.1282E+00 -3.7589E+00  1.2646E-01 -8.6827E-01 -2.4178E-01  3.5561E-04  3.2102E-01 -1.7559E-03  7.8386E-01  1.2516E-01
             2.6256E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1698.55295650126        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1246
 NPARAMETR:  9.7293E-01  1.4711E+00  1.6369E+00  7.9120E-01  1.4105E+00  9.9783E-01  9.0394E-01  2.9533E-01  1.1802E+00  1.3126E+00
             9.7084E-01
 PARAMETER:  7.2554E-02  4.8602E-01  5.9278E-01 -1.3420E-01  4.4393E-01  9.7830E-02 -9.9137E-04 -1.1197E+00  2.6568E-01  3.7203E-01
             7.0408E-02
 GRADIENT:  -7.9823E-01  5.6747E-02  1.8342E-01  1.6839E+00  4.7444E-02  2.9305E-04 -2.9524E-03 -9.4691E-04 -4.4034E-01 -1.0291E-01
            -1.6085E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1698.55363803888        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1425
 NPARAMETR:  9.7256E-01  1.4984E+00  1.5909E+00  7.7373E-01  1.4138E+00  9.9796E-01  9.0565E-01  2.1685E-01  1.1828E+00  1.3130E+00
             9.7019E-01
 PARAMETER:  7.2180E-02  5.0439E-01  5.6431E-01 -1.5654E-01  4.4632E-01  9.7962E-02  9.0237E-04 -1.4285E+00  2.6786E-01  3.7231E-01
             6.9733E-02
 GRADIENT:  -1.8922E+00  1.1725E+00  1.9235E-01  2.3674E+00  4.6433E-01  2.2313E-03  7.8389E-02 -1.8222E-04 -1.0666E+00 -1.4222E-01
            -2.9412E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1698.56894996426        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1616             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7392E-01  1.5042E+00  1.5759E+00  7.6624E-01  1.4143E+00  9.9811E-01  9.0475E-01  1.9226E-01  1.1958E+00  1.3136E+00
             9.7042E-01
 PARAMETER:  7.3574E-02  5.0827E-01  5.5482E-01 -1.6626E-01  4.4664E-01  9.8112E-02 -1.0091E-04 -1.5489E+00  2.7879E-01  3.7281E-01
             6.9972E-02
 GRADIENT:   4.3532E+02  3.6920E+02  1.1325E+00  6.2570E+01  2.1189E+01  3.8234E+01  6.0139E+00  1.1037E-02  1.2537E+01  4.7206E+00
             8.6809E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1698.57379755435        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1759
 NPARAMETR:  9.7361E-01  1.5044E+00  1.5744E+00  7.6685E-01  1.4126E+00  9.9783E-01  8.9336E-01  1.3909E-01  1.1994E+00  1.3128E+00
             9.7037E-01
 PARAMETER:  7.3258E-02  5.0837E-01  5.5388E-01 -1.6547E-01  4.4543E-01  9.7833E-02 -1.2771E-02 -1.8726E+00  2.8181E-01  3.7214E-01
             6.9919E-02
 GRADIENT:   4.5617E-01 -2.5983E+00 -1.4682E-02  1.1692E-01  1.6657E-01 -3.8220E-02 -3.3929E-01  1.3103E-03 -6.4943E-01  1.3880E-01
            -1.8600E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1698.57686576610        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1920
 NPARAMETR:  9.7404E-01  1.5057E+00  1.5724E+00  7.6582E-01  1.4126E+00  9.9801E-01  8.9524E-01  1.0000E-02  1.2046E+00  1.3122E+00
             9.7068E-01
 PARAMETER:  7.3695E-02  5.0924E-01  5.5261E-01 -1.6681E-01  4.4544E-01  9.8007E-02 -1.0663E-02 -1.4567E+01  2.8617E-01  3.7170E-01
             7.0239E-02
 GRADIENT:   1.4224E+00 -2.8152E+00  1.2409E-01  5.8578E-02  9.4874E-02  2.6463E-02  2.1524E-01  0.0000E+00 -2.9061E-02  5.9430E-02
            -8.4402E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1698.57809538365        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2104             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7466E-01  1.5093E+00  1.5669E+00  7.6182E-01  1.4129E+00  9.9821E-01  8.9192E-01  2.6588E-02  1.2101E+00  1.3120E+00
             9.7093E-01
 PARAMETER:  7.4334E-02  5.1167E-01  5.4909E-01 -1.7205E-01  4.4564E-01  9.8210E-02 -1.4381E-02 -3.5273E+00  2.9072E-01  3.7152E-01
             7.0499E-02
 GRADIENT:   4.3658E+02  3.7426E+02  1.1658E+00  6.3190E+01  2.0479E+01  3.8146E+01  5.3663E+00  4.1696E-04  1.3344E+01  4.6499E+00
             8.7058E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1698.58078206567        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2280             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7430E-01  1.5120E+00  1.5637E+00  7.6174E-01  1.4130E+00  9.9819E-01  8.9081E-01  2.1177E-02  1.2111E+00  1.3117E+00
             9.7084E-01
 PARAMETER:  7.3962E-02  5.1342E-01  5.4708E-01 -1.7215E-01  4.4573E-01  9.8192E-02 -1.5622E-02 -3.7549E+00  2.9149E-01  3.7136E-01
             7.0410E-02
 GRADIENT:   4.3576E+02  3.7844E+02  1.0732E+00  6.4778E+01  2.0519E+01  3.8184E+01  5.2614E+00  2.8845E-04  1.3359E+01  4.6098E+00
             8.1204E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1698.58137275842        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2459
 NPARAMETR:  9.7417E-01  1.5133E+00  1.5612E+00  7.6107E-01  1.4129E+00  9.9819E-01  8.9031E-01  1.3968E-02  1.2121E+00  1.3114E+00
             9.7083E-01
 PARAMETER:  7.3833E-02  5.1426E-01  5.4547E-01 -1.7303E-01  4.4568E-01  9.8186E-02 -1.6184E-02 -4.1710E+00  2.9239E-01  3.7111E-01
             7.0392E-02
 GRADIENT:   1.6468E+00 -2.4866E+00 -2.3061E-02  5.8188E-01 -2.9402E-02  8.3387E-02  2.7519E-02  1.9815E-05  6.8148E-02  2.0430E-02
             1.6334E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1698.58224620827        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2635
 NPARAMETR:  9.7436E-01  1.5145E+00  1.5596E+00  7.5947E-01  1.4132E+00  9.9822E-01  8.9018E-01  1.0756E-02  1.2127E+00  1.3116E+00
             9.7081E-01
 PARAMETER:  7.4024E-02  5.1508E-01  5.4443E-01 -1.7513E-01  4.4587E-01  9.8214E-02 -1.6328E-02 -4.4323E+00  2.9288E-01  3.7126E-01
             7.0372E-02
 GRADIENT:   2.0691E+00 -3.4623E+00  1.9597E-02 -2.0731E-01 -9.9393E-02  9.6897E-02  4.8325E-02  1.2545E-05  2.6205E-02  4.5029E-02
             2.6066E-02

0ITERATION NO.:   79    OBJECTIVE VALUE:  -1698.58273700904        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     2769
 NPARAMETR:  9.7418E-01  1.5185E+00  1.5576E+00  7.5894E-01  1.4130E+00  9.9818E-01  8.8901E-01  1.0000E-02  1.2143E+00  1.3114E+00
             9.7060E-01
 PARAMETER:  7.3849E-02  5.1564E-01  5.4379E-01 -1.7569E-01  4.4597E-01  9.8211E-02 -1.6643E-02 -4.6700E+00  2.9322E-01  3.7125E-01
             7.0357E-02
 GRADIENT:   5.0313E-03 -7.0738E-01  1.2185E-02  2.5685E-02  4.9413E-02  3.1617E-03  2.4708E-02  0.0000E+00 -2.9515E-02  5.9553E-03
             2.3309E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2769
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.6190E-04 -1.7509E-02 -9.8295E-05  6.0281E-03 -3.7153E-02
 SE:             2.9766E-02  2.0214E-02  5.5465E-05  2.1876E-02  2.3907E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8494E-01  3.8639E-01  7.6363E-02  7.8288E-01  1.2017E-01

 ETASHRINKSD(%)  2.7869E-01  3.2282E+01  9.9814E+01  2.6714E+01  1.9909E+01
 ETASHRINKVR(%)  5.5660E-01  5.4142E+01  1.0000E+02  4.6291E+01  3.5854E+01
 EBVSHRINKSD(%)  4.2307E-01  3.0574E+01  9.9821E+01  2.8341E+01  1.6678E+01
 EBVSHRINKVR(%)  8.4434E-01  5.1801E+01  1.0000E+02  4.8650E+01  3.0575E+01
 RELATIVEINF(%)  9.8798E+01  9.1137E-01  6.3243E-05  9.9108E-01  1.6157E+01
 EPSSHRINKSD(%)  4.1073E+01
 EPSSHRINKVR(%)  6.5276E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1698.5827370090408     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -963.43191044530261     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1698.583       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.74E-01  1.52E+00  1.56E+00  7.59E-01  1.41E+00  9.98E-01  8.90E-01  1.00E-02  1.21E+00  1.31E+00  9.71E-01
 


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
+       -7.56E+00  3.03E+02
 
 TH 3
+        9.18E-01  1.36E+01  1.65E+01
 
 TH 4
+       -6.08E+00  4.03E+02 -2.44E+01  6.76E+02
 
 TH 5
+       -2.45E+00 -6.28E+01 -4.21E+01  4.63E+01  2.26E+02
 
 TH 6
+        9.54E-01 -2.03E+00  2.49E-01 -3.07E+00 -1.82E+00  1.96E+02
 
 TH 7
+        1.58E+00  2.17E-01  7.73E+00 -9.27E+00 -9.92E+00 -3.98E-02  6.42E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.06E+00 -8.72E+00 -4.27E+00  2.31E+01  3.83E+00 -3.39E-01  3.30E+01  0.00E+00  4.35E+01
 
 TH10
+       -1.87E-01 -3.96E+00 -3.81E+00 -2.86E+00 -3.41E+01  5.80E-01 -5.11E-01  0.00E+00  1.73E+00  6.12E+01
 
 TH11
+       -1.02E+01 -1.53E+01 -1.14E+01 -3.65E+00  6.98E+00  1.49E+00  6.59E+00  0.00E+00  3.73E+00  1.71E+01  2.59E+02
 
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
 #CPUT: Total CPU Time in Seconds,       44.208
Stop Time:
Wed Sep 29 11:32:43 CDT 2021
