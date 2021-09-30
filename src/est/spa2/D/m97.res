Thu Sep 30 10:09:34 CDT 2021
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
$DATA ../../../../data/spa2/D/dat97.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   37242.5442295627        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1334E+03  8.4790E+02  2.3654E+01  8.0274E+02  9.2443E+01 -3.2636E+03 -1.5720E+03 -6.2821E+01 -2.1941E+03 -8.7042E+02
            -7.0445E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -300.507334893055        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0809E+00  1.1142E+00  9.5969E-01  1.2270E+00  1.1822E+00  2.0500E+00  1.3079E+00  9.7542E-01  1.1131E+00  9.4792E-01
             1.4457E+01
 PARAMETER:  1.7776E-01  2.0817E-01  5.8858E-02  3.0455E-01  2.6740E-01  8.1786E-01  3.6840E-01  7.5111E-02  2.0719E-01  4.6519E-02
             2.7712E+00
 GRADIENT:  -1.7533E+01 -1.8052E+01 -1.5323E+01  1.0529E+01  2.0170E+01  4.3339E+01 -1.7382E+01  3.4663E+00 -7.8844E+00  9.8477E+00
            -3.8177E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -425.340075014746        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0410E+00  2.1010E+00  1.2308E+00  7.3487E-01  1.6184E+00  2.3151E+00  5.4081E+00  3.2179E-01  1.2239E+00  2.4109E-02
             1.5481E+01
 PARAMETER:  1.4022E-01  8.4240E-01  3.0769E-01 -2.0806E-01  5.8142E-01  9.3943E-01  1.7879E+00 -1.0339E+00  3.0205E-01 -3.6252E+00
             2.8396E+00
 GRADIENT:  -1.7347E+01 -9.2457E+00 -7.4064E+00 -5.6139E+01 -7.0527E+00  5.7774E+01  6.8429E+01  7.6067E-02  1.0519E+01  5.5833E-03
             2.2124E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -495.676883333907        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0343E+00  9.4834E-01  5.4414E+00  1.4961E+00  2.6410E+00  1.6690E+00  6.0762E+00  2.3396E-01  6.6765E-01  6.9306E-01
             1.4070E+01
 PARAMETER:  1.3368E-01  4.6962E-02  1.7940E+00  5.0287E-01  1.0712E+00  6.1221E-01  1.9044E+00 -1.3526E+00 -3.0399E-01 -2.6663E-01
             2.7440E+00
 GRADIENT:   1.0034E+01 -4.3137E+00 -2.8214E+00  3.0939E+01  3.3067E+00  5.8661E-01  5.5703E+00  1.3290E-02  6.1482E-01  1.5785E+00
             6.7894E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -514.480609826815        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  9.6097E-01  7.3920E-01  9.5310E+00  1.5800E+00  2.9131E+00  1.5920E+00  8.5316E+00  2.8164E-01  6.9989E-01  5.8609E-01
             1.3677E+01
 PARAMETER:  6.0189E-02 -2.0219E-01  2.3546E+00  5.5742E-01  1.1692E+00  5.6498E-01  2.2438E+00 -1.1671E+00 -2.5683E-01 -4.3428E-01
             2.7158E+00
 GRADIENT:  -5.8417E+00  3.8265E+00 -1.1161E+00  1.4821E+00  2.1804E+00  1.2696E+00  3.6351E+00  5.2062E-03  7.8304E-01  9.5920E-01
             1.1595E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -516.175305316682        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  9.5686E-01  5.3028E-01  1.4571E+01  1.6096E+00  2.7912E+00  1.5113E+00  8.7724E+00  2.9581E-01  7.3772E-01  2.9341E-01
             1.3626E+01
 PARAMETER:  5.5902E-02 -5.3435E-01  2.7790E+00  5.7600E-01  1.1265E+00  5.1299E-01  2.2716E+00 -1.1180E+00 -2.0419E-01 -1.1262E+00
             2.7120E+00
 GRADIENT:   6.3257E-01 -1.0427E+00 -5.1698E-01 -4.3956E-01 -1.8249E-01  1.1320E+00 -6.3215E-02  2.5455E-03  6.7890E-01  2.7614E-01
             4.3938E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -516.587179427108        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      763             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5600E-01  6.0233E-01  4.8720E+01  1.6078E+00  2.7824E+00  1.5086E+00  1.0461E+01  3.3684E-01  7.1055E-01  9.1279E-02
             1.3566E+01
 PARAMETER:  5.5004E-02 -4.0694E-01  3.9861E+00  5.7487E-01  1.1233E+00  5.1119E-01  2.4477E+00 -9.8815E-01 -2.4171E-01 -2.2938E+00
             2.7076E+00
 GRADIENT:   8.1274E+00  7.0673E+00 -3.6895E-02 -8.5289E+00 -1.2813E+00  6.5467E+00  2.0672E+02  2.7454E-04  7.8721E-01  2.9754E-02
             3.5094E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -517.996237766818        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      922             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4929E-01  4.3301E-01  1.6811E+02  1.6393E+00  2.8810E+00  1.5024E+00  9.6788E+00  3.3295E-01  7.0302E-01  4.0694E-02
             1.3541E+01
 PARAMETER:  4.7956E-02 -7.3699E-01  5.2246E+00  5.9424E-01  1.1581E+00  5.0707E-01  2.3699E+00 -9.9977E-01 -2.5237E-01 -3.1017E+00
             2.7057E+00
 GRADIENT:   7.5744E-01  2.0845E+00 -1.3114E-02  1.2937E+01 -5.7816E-01  5.1256E+00  1.9462E+02  2.1413E-05 -2.3083E+00  5.7490E-03
             3.0261E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -518.331808677494        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1100             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5580E-01  4.3840E-01  2.1934E+03  1.6552E+00  2.9323E+00  1.5054E+00  1.0342E+01  3.3320E-01  7.6810E-01  1.8353E-02
             1.3575E+01
 PARAMETER:  5.4793E-02 -7.2463E-01  7.7932E+00  6.0395E-01  1.1758E+00  5.0905E-01  2.4363E+00 -9.9902E-01 -1.6383E-01 -3.8980E+00
             2.7082E+00
 GRADIENT:   5.3160E+00  4.1437E+00 -1.3112E-03 -6.4639E-01  1.1044E+00  5.2693E+00  2.1644E+02  1.0229E-06  1.0334E+00  1.1855E-03
             3.6958E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -518.484666174704        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1260
 NPARAMETR:  9.5335E-01  3.9887E-01  3.9333E+03  1.6772E+00  2.9297E+00  1.5039E+00  1.0155E+01  3.3276E-01  7.6694E-01  1.3580E-02
             1.3557E+01
 PARAMETER:  5.2226E-02 -8.1911E-01  8.3772E+00  6.1712E-01  1.1749E+00  5.0807E-01  2.4180E+00 -1.0003E+00 -1.6534E-01 -4.1992E+00
             2.7069E+00
 GRADIENT:  -9.9035E-01  7.1586E-02 -8.5951E-04  8.2551E-01 -2.2097E-01 -1.7427E-02  2.5122E+01  2.7282E-07 -1.0672E+00  6.0911E-04
            -1.4208E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -518.570198356524        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1439
 NPARAMETR:  9.5567E-01  3.8985E-01  5.6297E+03  1.6872E+00  2.9303E+00  1.5065E+00  1.0428E+01  3.3248E-01  7.8621E-01  1.1606E-02
             1.3586E+01
 PARAMETER:  5.4660E-02 -8.4198E-01  8.7358E+00  6.2305E-01  1.1751E+00  5.0981E-01  2.4445E+00 -1.0012E+00 -1.4053E-01 -4.3562E+00
             2.7090E+00
 GRADIENT:  -6.6541E-02  5.6210E-01 -6.3498E-04 -2.1525E+00 -1.7913E-02  4.0033E-01  3.0307E+01  7.4360E-07 -2.5857E-01  4.4638E-04
             1.5797E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -518.610790421371        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1625             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5577E-01  3.6594E-01  2.0896E+04  1.6970E+00  2.9288E+00  1.5051E+00  1.0602E+01  3.4256E-01  8.0298E-01  1.0000E-02
             1.3556E+01
 PARAMETER:  5.4761E-02 -9.0530E-01  1.0047E+01  6.2886E-01  1.1746E+00  5.0885E-01  2.4611E+00 -9.7132E-01 -1.1942E-01 -4.7632E+00
             2.7068E+00
 GRADIENT:   2.5488E+00  3.4750E+00 -1.7325E-04  1.1084E+01  5.6325E-01  4.7434E+00  2.3035E+02  8.8656E-07  2.0044E-01  0.0000E+00
             3.2832E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -518.616838659430        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1794             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5525E-01  3.6124E-01  1.8375E+04  1.7008E+00  2.9271E+00  1.5061E+00  1.0629E+01  3.4879E-01  8.0680E-01  1.0000E-02
             1.3559E+01
 PARAMETER:  5.4220E-02 -9.1821E-01  9.9188E+00  6.3108E-01  1.1740E+00  5.0949E-01  2.4635E+00 -9.5329E-01 -1.1468E-01 -4.7353E+00
             2.7070E+00
 GRADIENT:   2.9898E+00  3.1055E+00 -2.0575E-04  1.1812E+01  5.2583E-01  4.6050E+00  2.3205E+02  4.4094E-06  3.6224E-01  0.0000E+00
             3.3183E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -518.621320263100        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1982
 NPARAMETR:  9.5524E-01  3.5395E-01  1.9046E+04  1.7032E+00  2.9271E+00  1.5065E+00  1.0629E+01  3.5148E-01  8.0817E-01  1.0000E-02
             1.3560E+01
 PARAMETER:  5.4209E-02 -9.3861E-01  9.9546E+00  6.3249E-01  1.1740E+00  5.0979E-01  2.4636E+00 -9.4559E-01 -1.1299E-01 -4.7353E+00
             2.7072E+00
 GRADIENT:  -1.3996E+00 -1.1493E-01 -2.0165E-04 -7.6842E-01  8.4213E-03  6.3686E-01  3.2324E+01 -3.8660E-06 -2.8133E-02  0.0000E+00
             6.5894E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -518.623124320857        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2173             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5555E-01  3.4933E-01  2.7899E+04  1.7067E+00  2.9293E+00  1.5072E+00  1.0664E+01  3.7150E-01  8.1418E-01  1.0000E-02
             1.3565E+01
 PARAMETER:  5.4527E-02 -9.5175E-01  1.0336E+01  6.3455E-01  1.1748E+00  5.1027E-01  2.4668E+00 -8.9020E-01 -1.0557E-01 -4.7353E+00
             2.7075E+00
 GRADIENT:  -1.6965E-01  3.1345E+00 -1.3836E-04  1.2872E+01  5.4834E-01  5.4109E+00  2.3399E+02  4.8478E-06  1.7797E-01  0.0000E+00
             3.4176E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -518.623434720273        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2358
 NPARAMETR:  9.5558E-01  3.4846E-01  7.2389E+04  1.7074E+00  2.9296E+00  1.5071E+00  1.0672E+01  3.7814E-01  8.1520E-01  1.0000E-02
             1.3565E+01
 PARAMETER:  5.4566E-02 -9.5422E-01  1.1290E+01  6.3495E-01  1.1749E+00  5.1021E-01  2.4676E+00 -8.7248E-01 -1.0433E-01 -4.7353E+00
             2.7075E+00
 GRADIENT:  -1.5413E+01 -7.6270E-01 -4.9958E-05  1.1320E+00  8.7474E-02  3.5378E+00  3.2786E+01 -2.0725E-08 -3.3246E-01  0.0000E+00
             9.3766E+00

0ITERATION NO.:   76    OBJECTIVE VALUE:  -518.623434720273        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     2387
 NPARAMETR:  9.5573E-01  3.4852E-01  7.2165E+04  1.7077E+00  2.9287E+00  1.5072E+00  1.0673E+01  3.7815E-01  8.1435E-01  1.0000E-02
             1.3564E+01
 PARAMETER:  5.4566E-02 -9.5422E-01  1.1290E+01  6.3495E-01  1.1749E+00  5.1021E-01  2.4676E+00 -8.7248E-01 -1.0433E-01 -4.7353E+00
             2.7075E+00
 GRADIENT:  -1.6438E-01 -4.0776E-03  1.9040E-04 -2.9808E-01  2.7755E-02 -1.2335E-02 -3.2525E-02 -7.4166E-05  1.4539E-01  0.0000E+00
             1.3401E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2387
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0596E-02  6.4026E-02  5.2979E-08 -7.5262E-02  1.8774E-05
 SE:             2.6271E-02  2.1035E-02  1.4132E-08  1.2523E-02  2.8840E-05
 N:                     100         100         100         100         100

 P VAL.:         2.4418E-01  2.3368E-03  1.7760E-04  1.8644E-09  5.1507E-01

 ETASHRINKSD(%)  1.1987E+01  2.9529E+01  1.0000E+02  5.8046E+01  9.9903E+01
 ETASHRINKVR(%)  2.2537E+01  5.0339E+01  1.0000E+02  8.2398E+01  1.0000E+02
 EBVSHRINKSD(%)  1.6117E+01  2.4609E+01  1.0000E+02  5.7178E+01  9.9837E+01
 EBVSHRINKVR(%)  2.9637E+01  4.3161E+01  1.0000E+02  8.1662E+01  1.0000E+02
 RELATIVEINF(%)  6.6959E+01  3.8571E+01  0.0000E+00  8.8828E+00  4.7460E-05
 EPSSHRINKSD(%)  3.6131E+00
 EPSSHRINKVR(%)  7.0957E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -518.62343472027339     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       584.10280512533370     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    60.83
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -518.623       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.56E-01  3.48E-01  7.24E+04  1.71E+00  2.93E+00  1.51E+00  1.07E+01  3.78E-01  8.15E-01  1.00E-02  1.36E+01
 


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
+        3.72E+02
 
 TH 2
+       -4.06E+00  5.72E+01
 
 TH 3
+       -8.12E-07 -4.83E-07  6.33E-13
 
 TH 4
+       -3.72E+01  1.97E+01 -4.95E-07  1.67E+02
 
 TH 5
+        1.37E+00 -1.18E+00 -3.79E-08 -7.30E+00  3.52E+00
 
 TH 6
+        8.53E+00 -3.90E+00  3.16E-07 -7.57E+00 -6.47E-01  6.14E+01
 
 TH 7
+        9.68E-01  3.57E+00  1.54E-08 -3.57E+00  9.74E-02  1.26E-01  8.38E-01
 
 TH 8
+        3.92E+00  1.91E+00 -7.99E-07  1.03E+00 -8.91E-02  8.37E-01 -3.81E-02  2.13E+00
 
 TH 9
+       -1.69E+01  8.21E+00  5.45E-07 -4.78E+01  3.01E+00 -4.00E+00  5.57E-01  6.06E-01  2.34E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -9.52E+00 -2.12E+00  9.11E-09 -9.84E+00  3.78E-01  1.24E+00  1.01E-01  8.99E-03  3.69E+00  0.00E+00  3.52E+00
 
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
 #CPUT: Total CPU Time in Seconds,       73.425
Stop Time:
Thu Sep 30 10:10:49 CDT 2021
