Wed Sep 29 09:58:27 CDT 2021
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
$DATA ../../../../data/int/D/dat89.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   50575.5177837134        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.6753E+02  7.9085E+02  7.0288E+01  7.9803E+02  2.0807E+02 -3.6398E+03 -1.8877E+03 -1.3339E+02 -2.4797E+03 -1.0011E+03
            -9.8458E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -553.824060414552        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0230E+00  2.2535E+00  8.0341E-01  2.0897E+00  9.4077E-01  4.7669E+00  5.0120E+00  9.9309E-01  1.9851E+00  1.4372E+00
             1.2741E+01
 PARAMETER:  1.2274E-01  9.1246E-01 -1.1889E-01  8.3701E-01  3.8946E-02  1.6617E+00  1.7118E+00  9.3066E-02  7.8568E-01  4.6271E-01
             2.6448E+00
 GRADIENT:  -2.4674E+01  4.9884E+01 -4.8142E+01  1.9209E+02  4.5092E+00  1.9470E+02  5.5153E+01  4.2264E+00 -1.9222E+01  2.2837E+01
             5.0373E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -636.153998983154        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  8.8127E-01  2.1504E+00  3.2470E+01  2.6754E+00  2.7895E+00  4.2730E+00  7.8807E+00  7.2763E-01  2.2902E+00  1.3710E+00
             1.2699E+01
 PARAMETER: -2.6392E-02  8.6568E-01  3.5803E+00  1.0841E+00  1.1259E+00  1.5523E+00  2.1644E+00 -2.1796E-01  9.2864E-01  4.1557E-01
             2.6415E+00
 GRADIENT:  -3.8361E+01  4.3823E+01 -4.9721E+00  1.3228E+02  2.1429E+01  1.5185E+02  5.7810E+01  2.7298E-02 -8.5448E+00  2.2012E+01
             8.7280E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -736.225349582245        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.2566E+00  1.4759E+00  1.5304E+01  1.1741E+00  2.4440E+00  1.9442E+00  4.8151E+00  2.2829E+00  1.2544E+00  7.8392E-01
             1.3261E+01
 PARAMETER:  3.2842E-01  4.8930E-01  2.8281E+00  2.6047E-01  9.9364E-01  7.6487E-01  1.6718E+00  9.2544E-01  3.2662E-01 -1.4345E-01
             2.6848E+00
 GRADIENT:   4.1063E+01 -9.7831E+00 -1.9726E+00 -2.4143E+01  1.1253E+01 -5.5584E+00  4.4313E+01  3.3982E-01  4.0002E+00  7.3881E+00
             1.7138E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -738.628812343847        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.1194E+00  1.6352E+00  1.3288E+01  1.1517E+00  2.3775E+00  2.0006E+00  4.2506E+00  1.5192E+00  1.4736E+00  6.2617E-01
             1.2685E+01
 PARAMETER:  2.1275E-01  5.9175E-01  2.6869E+00  2.4121E-01  9.6607E-01  7.9344E-01  1.5471E+00  5.1821E-01  4.8771E-01 -3.6813E-01
             2.6404E+00
 GRADIENT:  -8.0166E+00 -4.6447E+00 -1.7560E+00 -1.1665E+00  4.4534E+00  3.3384E+00  2.4151E+00  1.7903E-01  4.8664E+00  4.7706E+00
             1.2594E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -739.390118438632        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0886E+00  1.7134E+00  1.3072E+01  1.0748E+00  2.3704E+00  2.0063E+00  4.1030E+00  9.0305E-01  1.4921E+00  5.0654E-01
             1.2144E+01
 PARAMETER:  1.8487E-01  6.3849E-01  2.6705E+00  1.7217E-01  9.6305E-01  7.9629E-01  1.5117E+00 -1.9720E-03  5.0019E-01 -5.8015E-01
             2.5969E+00
 GRADIENT:  -1.1781E+01 -2.5044E+00 -1.1839E+00 -1.6053E-01  1.9364E+00  2.0376E+00 -2.7629E+00  5.1074E-02  4.0922E+00  3.1025E+00
             6.5582E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -739.681037214396        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  1.0860E+00  1.7461E+00  1.3879E+01  1.0261E+00  2.3866E+00  2.0022E+00  4.0701E+00  5.9146E-01  1.4259E+00  3.7378E-01
             1.1796E+01
 PARAMETER:  1.8254E-01  6.5736E-01  2.7304E+00  1.2580E-01  9.6988E-01  7.9423E-01  1.5037E+00 -4.2516E-01  4.5477E-01 -8.8408E-01
             2.5678E+00
 GRADIENT:  -6.4416E+00 -1.0657E+00 -7.5813E-01 -7.4392E-01  1.0958E+00  3.2311E-01 -1.4655E+00  1.6300E-02  2.1850E+00  1.6443E+00
             1.7909E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -740.234421270164        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      514
 NPARAMETR:  1.0985E+00  1.7397E+00  1.7800E+01  1.0099E+00  2.4299E+00  1.9970E+00  4.1059E+00  3.5486E-01  1.2627E+00  1.6070E-01
             1.1604E+01
 PARAMETER:  1.9390E-01  6.5369E-01  2.9792E+00  1.0981E-01  9.8784E-01  7.9163E-01  1.5124E+00 -9.3602E-01  3.3327E-01 -1.7282E+00
             2.5513E+00
 GRADIENT:   2.8368E+00  8.1968E-02 -2.8419E-01  5.6937E-01  4.5899E-01 -7.9621E-01  9.4378E-01  3.2305E-03 -1.7704E+00  2.8825E-01
            -1.6132E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -754.583649871154        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  1.1311E+00  1.1972E+00  7.0278E+01  1.5342E+00  2.5725E+00  2.0182E+00  5.9273E+00  1.8671E-01  1.4722E+00  2.1216E-02
             1.2144E+01
 PARAMETER:  2.2318E-01  2.7998E-01  4.3525E+00  5.2804E-01  1.0449E+00  8.0220E-01  1.8796E+00 -1.5782E+00  4.8674E-01 -3.7530E+00
             2.5969E+00
 GRADIENT:  -1.5293E-01  8.0927E+00 -4.1713E-01  8.8705E+00  6.5287E-01 -3.4477E+00 -4.9645E+00  1.9872E-04 -2.6231E+00  4.3385E-03
            -5.6624E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -758.646421741176        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      848
 NPARAMETR:  1.1397E+00  6.2070E-01  1.5632E+02  1.8205E+00  2.5873E+00  2.0421E+00  7.1811E+00  1.8234E-01  1.6591E+00  3.1934E-02
             1.2161E+01
 PARAMETER:  2.3073E-01 -3.7691E-01  5.1519E+00  6.9910E-01  1.0506E+00  8.1397E-01  2.0715E+00 -1.6019E+00  6.0626E-01 -3.3441E+00
             2.5983E+00
 GRADIENT:   3.3967E+00  1.1678E+00 -3.0027E-01  1.9176E+00  1.4731E+00  5.3600E-01  3.2578E+00  7.6837E-05  7.6562E-01  9.8909E-03
             3.6072E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -759.258191098746        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1007
 NPARAMETR:  1.1259E+00  5.7509E-01  7.0201E+02  1.8158E+00  2.5784E+00  2.0381E+00  7.5581E+00  1.8138E-01  1.6405E+00  1.0469E-02
             1.2116E+01
 PARAMETER:  2.1862E-01 -4.5323E-01  6.6540E+00  6.9655E-01  1.0472E+00  8.1201E-01  2.1226E+00 -1.6072E+00  5.9502E-01 -4.4594E+00
             2.5945E+00
 GRADIENT:  -5.9140E-01  1.4563E+00 -4.7873E-02 -5.2034E+00 -2.0941E+00  1.8934E-01  1.1873E+01  3.8529E-06  4.4004E-01  1.0525E-03
             3.2569E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -759.754271417540        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1171
 NPARAMETR:  1.1272E+00  4.9108E-01  6.0431E+04  1.8763E+00  2.6078E+00  2.0360E+00  7.9487E+00  1.8255E-01  1.6738E+00  1.0000E-02
             1.2095E+01
 PARAMETER:  2.1974E-01 -6.1115E-01  1.1109E+01  7.2932E-01  1.0585E+00  8.1096E-01  2.1730E+00 -1.6007E+00  6.1510E-01 -3.1123E+01
             2.5928E+00
 GRADIENT:   5.2870E+00  4.2531E+00 -6.7211E-04  1.5931E+01  2.8512E+00  1.0683E+01  1.6627E+02  2.5473E-06  3.9119E+00  0.0000E+00
             4.8556E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -759.907105792418        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1356
 NPARAMETR:  1.1320E+00  4.3129E-01  7.9335E+04  1.9198E+00  2.6055E+00  2.0293E+00  7.9569E+00  3.8656E-01  1.6939E+00  1.0000E-02
             1.2093E+01
 PARAMETER:  2.2402E-01 -7.4097E-01  1.1381E+01  7.5222E-01  1.0576E+00  8.0770E-01  2.1740E+00 -8.5047E-01  6.2705E-01 -3.1123E+01
             2.5926E+00
 GRADIENT:  -1.8596E+00 -1.5476E-01 -5.5345E-04 -5.9274E-01 -3.5109E-01  1.3073E+00  9.9111E+00  3.6133E-06 -1.3000E-01  0.0000E+00
            -1.8808E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -760.122169935337        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     1528
 NPARAMETR:  1.1270E+00  3.7134E-01  9.4394E+04  1.9487E+00  2.6073E+00  2.0345E+00  8.2332E+00  3.6554E-01  1.6918E+00  1.0000E-02
             1.2088E+01
 PARAMETER:  2.1956E-01 -8.9062E-01  1.1555E+01  7.6716E-01  1.0583E+00  8.1026E-01  2.2082E+00 -9.0639E-01  6.2578E-01 -3.1123E+01
             2.5922E+00
 GRADIENT:  -3.2981E+01  6.0887E-03 -4.7735E-04  1.2174E+01 -1.5880E-01  1.5618E+01  1.0697E+01 -1.2570E-05 -1.0016E+01  0.0000E+00
             8.3125E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -760.287683250163        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1718
 NPARAMETR:  1.1234E+00  3.3807E-01  1.1387E+06  1.9733E+00  2.6073E+00  2.0342E+00  8.4633E+00  3.6772E-01  1.7074E+00  1.0000E-02
             1.2091E+01
 PARAMETER:  2.1635E-01 -9.8449E-01  1.4045E+01  7.7972E-01  1.0583E+00  8.1011E-01  2.2357E+00 -9.0044E-01  6.3497E-01 -3.1123E+01
             2.5924E+00
 GRADIENT:  -6.7005E+01  8.1927E-01 -4.1246E-05  1.7214E+01 -4.8225E-01  6.7895E+00  1.5941E+01  4.7796E-05 -3.4853E+00  0.0000E+00
             2.4321E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -760.348995807474        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1883
 NPARAMETR:  1.1206E+00  3.5372E-01  1.0688E+06  1.9798E+00  2.6078E+00  2.0316E+00  8.6656E+00  3.6602E-01  1.7198E+00  1.0000E-02
             1.2090E+01
 PARAMETER:  2.1390E-01 -9.3926E-01  1.3982E+01  7.8301E-01  1.0585E+00  8.0883E-01  2.2594E+00 -9.0506E-01  6.4220E-01 -3.1123E+01
             2.5924E+00
 GRADIENT:  -3.8173E+01  9.5518E-01 -4.8186E-05  8.0187E-01  4.8274E-01  1.0685E+01  2.0536E+01  1.3291E-05 -3.7976E+00  0.0000E+00
             1.7646E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -760.506824658007        NO. OF FUNC. EVALS.: 202
 CUMULATIVE NO. OF FUNC. EVALS.:     2085             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1254E+00  2.9584E-01  1.6644E+06  2.0025E+00  2.6080E+00  2.0316E+00  8.8505E+00  3.5212E-01  1.7237E+00  1.0000E-02
             1.2083E+01
 PARAMETER:  2.1810E-01 -1.1179E+00  1.4425E+01  7.9440E-01  1.0586E+00  8.0881E-01  2.2805E+00 -9.4378E-01  6.4447E-01 -3.1123E+01
             2.5918E+00
 GRADIENT:  -4.2697E+01  4.6127E+00 -2.2065E-05  2.1768E+01  4.9263E+00  4.4268E+00  2.1317E+02  2.0219E-05  2.3827E+00  0.0000E+00
             6.4725E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -760.606056683630        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2275             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1207E+00  2.5862E-01  2.6644E+06  2.0226E+00  2.6135E+00  2.0326E+00  9.1110E+00  3.6025E-01  1.7349E+00  1.0000E-02
             1.2101E+01
 PARAMETER:  2.1396E-01 -1.2524E+00  1.4895E+01  8.0436E-01  1.0607E+00  8.0932E-01  2.3095E+00 -9.2095E-01  6.5092E-01 -3.1123E+01
             2.5933E+00
 GRADIENT:  -9.2692E+01  5.8751E+00 -2.7010E-05  3.8600E+01  4.3367E-01  1.8092E+01  2.2406E+02 -2.7770E-05  7.7300E+00  0.0000E+00
             8.1603E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -760.661891349943        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2458
 NPARAMETR:  1.1227E+00  2.6198E-01  1.9804E+06  2.0398E+00  2.6075E+00  2.0319E+00  9.2514E+00  3.5953E-01  1.7286E+00  1.0000E-02
             1.2086E+01
 PARAMETER:  2.1576E-01 -1.2395E+00  1.4599E+01  8.1285E-01  1.0584E+00  8.0896E-01  2.3248E+00 -9.2296E-01  6.4729E-01 -3.1123E+01
             2.5921E+00
 GRADIENT:  -8.4497E+02  1.2159E+02 -4.7570E-04 -4.4839E+02  7.7110E+01 -2.5860E+02  4.2135E+02  8.2665E-04 -1.2938E+02  0.0000E+00
             1.6301E+03

0ITERATION NO.:   95    OBJECTIVE VALUE:  -760.678109920950        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:     2663             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1206E+00  2.5612E-01  2.3621E+06  2.0384E+00  2.6114E+00  2.0312E+00  9.3875E+00  3.5850E-01  1.7294E+00  1.0000E-02
             1.2087E+01
 PARAMETER:  2.1383E-01 -1.2621E+00  1.4775E+01  8.1216E-01  1.0599E+00  8.0862E-01  2.3394E+00 -9.2582E-01  6.4775E-01 -3.1123E+01
             2.5922E+00
 GRADIENT:  -1.6570E+02  4.7311E+00 -3.5544E-05  6.3291E+01  3.7213E+00  7.2800E+01  2.2993E+02  3.7365E-04 -9.2246E+00  0.0000E+00
             9.6463E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -760.724998535658        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     2836
 NPARAMETR:  1.1268E+00  2.4034E-01  2.3101E+06  2.0424E+00  2.6121E+00  2.0268E+00  9.4307E+00  2.0303E-01  1.7318E+00  1.0000E-02
             1.2084E+01
 PARAMETER:  2.1936E-01 -1.3257E+00  1.4753E+01  8.1414E-01  1.0601E+00  8.0644E-01  2.3440E+00 -1.4944E+00  6.4918E-01 -3.1123E+01
             2.5919E+00
 GRADIENT:  -4.8199E+01  8.9861E+00 -4.2392E-06  3.6632E+01  2.3418E+00 -9.1414E+00  2.4171E+02  1.3278E-04  5.8769E+00  0.0000E+00
             5.8475E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -760.749902201729        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     3030             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1249E+00  2.3111E-01  1.5786E+06  2.0454E+00  2.6126E+00  2.0284E+00  9.5679E+00  2.0762E-01  1.7284E+00  1.0000E-02
             1.2090E+01
 PARAMETER:  2.1765E-01 -1.3649E+00  1.4372E+01  8.1560E-01  1.0604E+00  8.0723E-01  2.3584E+00 -1.4720E+00  6.4720E-01 -3.1123E+01
             2.5923E+00
 GRADIENT:  -3.5141E+01  6.6860E+00 -5.0114E-05  1.4769E+01  4.5744E+00  5.2431E-01  2.5355E+02 -2.3088E-04  4.8045E-01  0.0000E+00
             7.2415E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -760.780622269300        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3210
 NPARAMETR:  1.1217E+00  2.2099E-01  1.3678E+06  2.0509E+00  2.6137E+00  2.0286E+00  9.6545E+00  2.1495E-01  1.7315E+00  1.0000E-02
             1.2092E+01
 PARAMETER:  2.1483E-01 -1.4096E+00  1.4229E+01  8.1828E-01  1.0608E+00  8.0733E-01  2.3674E+00 -1.4374E+00  6.4899E-01 -3.1123E+01
             2.5925E+00
 GRADIENT:  -4.5298E+01  1.8325E+00 -4.7134E-05  1.3539E+01 -8.9031E-02 -2.2122E+01  2.7350E+01 -4.5395E-05  3.2185E+00  0.0000E+00
             7.5512E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -760.799395915253        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3392
 NPARAMETR:  1.1249E+00  2.1315E-01  1.1743E+06  2.0537E+00  2.6139E+00  2.0282E+00  9.7571E+00  2.1594E-01  1.7256E+00  1.0000E-02
             1.2092E+01
 PARAMETER:  2.1770E-01 -1.4458E+00  1.4076E+01  8.1965E-01  1.0608E+00  8.0713E-01  2.3780E+00 -1.4328E+00  6.4555E-01 -3.1123E+01
             2.5925E+00
 GRADIENT:  -3.8582E+01  4.1337E+00 -1.0615E-04 -3.4302E+01  3.3400E+00  5.8491E+00  4.0143E+01  2.0065E-05 -1.6269E+01  0.0000E+00
             2.0165E+01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -760.835745957094        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     3583
 NPARAMETR:  1.1273E+00  1.9971E-01  1.2965E+06  2.0633E+00  2.6117E+00  2.0280E+00  9.8761E+00  2.1911E-01  1.7409E+00  1.0000E-02
             1.2092E+01
 PARAMETER:  2.1981E-01 -1.5109E+00  1.4175E+01  8.2433E-01  1.0600E+00  8.0703E-01  2.3901E+00 -1.4182E+00  6.5439E-01 -3.1123E+01
             2.5925E+00
 GRADIENT:  -3.2489E+01  3.0459E+01 -3.5392E-04 -2.2191E+02  1.7063E+01  1.6239E+01  1.1704E+02 -2.6063E-05 -9.3242E+01  0.0000E+00
             7.3992E+01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -760.854112019828        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     3778             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1249E+00  1.9556E-01  1.2484E+06  2.0627E+00  2.6146E+00  2.0271E+00  9.9380E+00  2.2332E-01  1.7274E+00  1.0000E-02
             1.2092E+01
 PARAMETER:  2.1768E-01 -1.5319E+00  1.4137E+01  8.2402E-01  1.0611E+00  8.0661E-01  2.3964E+00 -1.3991E+00  6.4661E-01 -3.1123E+01
             2.5926E+00
 GRADIENT:  -3.5076E+01  6.6244E+00 -2.5629E-05  3.8255E+01  1.7709E+00  1.4757E+00  2.6304E+02  9.4200E-05  9.6065E+00  0.0000E+00
             5.5072E+01

0ITERATION NO.:  129    OBJECTIVE VALUE:  -760.856088167802        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     3923
 NPARAMETR:  1.1249E+00  1.9332E-01  1.0510E+06  2.0615E+00  2.6166E+00  2.0279E+00  1.0023E+01  2.3252E-01  1.7195E+00  1.0000E-02
             1.2097E+01
 PARAMETER:  2.1700E-01 -1.5590E+00  1.3996E+01  8.2535E-01  1.0607E+00  8.0565E-01  2.3976E+00 -1.3711E+00  6.4729E-01 -3.1123E+01
             2.5922E+00
 GRADIENT:  -4.4043E-01 -2.7083E-01  4.0542E-04  2.4339E+00 -4.1362E-01 -4.1689E-01 -3.2969E+00 -7.1101E-03  8.3556E-01  0.0000E+00
            -1.6941E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3923
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4445E-02  7.7536E-02  6.4085E-09 -8.3411E-02 -3.2897E-05
 SE:             2.8177E-02  1.9662E-02  4.1288E-09  1.8661E-02  9.5837E-05
 N:                     100         100         100         100         100

 P VAL.:         6.0820E-01  8.0380E-05  1.2062E-01  7.8338E-06  7.3140E-01

 ETASHRINKSD(%)  5.6035E+00  3.4129E+01  1.0000E+02  3.7484E+01  9.9679E+01
 ETASHRINKVR(%)  1.0893E+01  5.6610E+01  1.0000E+02  6.0918E+01  9.9999E+01
 EBVSHRINKSD(%)  7.5481E+00  4.2546E+01  1.0000E+02  2.3985E+01  9.9608E+01
 EBVSHRINKVR(%)  1.4526E+01  6.6990E+01  1.0000E+02  4.2217E+01  9.9998E+01
 RELATIVEINF(%)  8.5222E+01  2.3116E+01  0.0000E+00  3.9177E+01  3.3477E-04
 EPSSHRINKSD(%)  4.9230E+00
 EPSSHRINKVR(%)  9.6037E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -760.85608816780177     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       893.23327160060899     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   157.27
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    20.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -760.856       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.12E+00  1.90E-01  1.08E+06  2.07E+00  2.61E+00  2.03E+00  9.95E+00  2.30E-01  1.73E+00  1.00E-02  1.21E+01
 


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
+        1.82E+02
 
 TH 2
+        1.99E+01  4.15E+02
 
 TH 3
+       -1.43E-07  9.04E-08  6.74E-15
 
 TH 4
+       -1.32E+01 -1.29E+02 -3.94E-08  1.34E+02
 
 TH 5
+       -3.06E-01  1.00E+01  9.87E-09 -1.33E+01  2.35E+01
 
 TH 6
+       -4.18E+00  1.43E+01  1.60E-08 -7.40E+00  9.01E-01  3.26E+01
 
 TH 7
+        1.72E+00  2.63E+01  6.84E-11 -1.04E+01  7.70E-01  1.03E+00  1.98E+00
 
 TH 8
+        1.09E+00 -5.02E-01 -6.18E-08  1.16E+00  2.02E-01  1.27E+00 -1.45E-01 -5.37E-01
 
 TH 9
+       -9.48E+00 -7.54E+01 -8.54E-09  2.58E+01 -1.14E+00 -6.21E+00 -4.38E+00  1.31E+00  4.63E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.73E+00  1.07E+01 -3.30E-10 -1.12E+01  8.52E-01  2.55E+00  8.21E-01  2.22E-02 -1.06E+00  0.00E+00  6.62E+00
 
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
 #CPUT: Total CPU Time in Seconds,      177.567
Stop Time:
Wed Sep 29 10:01:26 CDT 2021
