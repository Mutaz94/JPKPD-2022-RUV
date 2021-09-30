Thu Sep 30 09:31:29 CDT 2021
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
$DATA ../../../../data/spa2/D/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   18938.1202131279        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.1369E+02  3.9812E+02  3.2280E+01  1.6496E+02  4.1766E+02 -2.7166E+03 -1.0800E+03 -7.2497E+01 -1.8838E+03 -1.0576E+03
            -3.5708E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -621.874945869523        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.3835E+00  1.2268E+00  8.0457E-01  2.3824E+00  1.0595E+00  3.0722E+00  2.1422E+00  9.5735E-01  2.6182E+00  1.3873E+00
             1.1876E+01
 PARAMETER:  4.2459E-01  3.0440E-01 -1.1745E-01  9.6810E-01  1.5783E-01  1.2224E+00  8.6182E-01  5.6417E-02  1.0625E+00  4.2735E-01
             2.5745E+00
 GRADIENT:  -5.4407E+00  1.2165E+01 -5.4728E+01  1.2292E+02  3.8263E+01  6.9804E+01 -2.4537E+01  5.0326E+00 -4.2318E+01  1.3661E+01
             1.2396E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -751.231945304803        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.4940E+00  5.7419E+00  2.4725E+00  2.0157E+00  6.2034E+00  2.9600E+00  3.5096E+00  7.4712E-01  1.6723E+01  1.7610E+00
             1.0176E+01
 PARAMETER:  5.0143E-01  1.8478E+00  1.0052E+00  8.0098E-01  1.9251E+00  1.1852E+00  1.3555E+00 -1.9153E-01  2.9168E+00  6.6590E-01
             2.4200E+00
 GRADIENT:   6.5400E+01  1.2021E+02  4.0440E+01  2.5880E+01  7.8818E+00  7.4345E+01  1.2657E+01 -3.7588E+00  3.1745E+01  5.2432E-01
             1.3503E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -762.585824256047        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.4189E+00  4.6488E+00  2.0691E+00  1.8504E+00  4.6031E+00  2.7582E+00  3.4193E+00  8.1921E-01  1.2827E+01  2.7414E+00
             9.8825E+00
 PARAMETER:  4.4985E-01  1.6366E+00  8.2711E-01  7.1542E-01  1.6267E+00  1.1146E+00  1.3294E+00 -9.9418E-02  2.6515E+00  1.1085E+00
             2.3908E+00
 GRADIENT:   5.8393E+01  8.6832E+01  1.6973E+01  3.0374E+01 -1.3519E+01  5.7685E+01  4.1331E+00 -5.9720E-01  5.1883E+00  2.3680E+01
             1.1715E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -769.585441114025        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      330
 NPARAMETR:  1.1935E+00  3.2737E+00  1.7592E+00  1.5485E+00  3.0317E+00  2.5009E+00  3.4778E+00  4.9825E-01  1.0620E+01  1.8045E+00
             9.6479E+00
 PARAMETER:  2.7688E-01  1.2859E+00  6.6485E-01  5.3730E-01  1.2091E+00  1.0166E+00  1.3464E+00 -5.9666E-01  2.4628E+00  6.9031E-01
             2.3667E+00
 GRADIENT:   8.0913E+00 -3.1242E+01 -4.3279E-01  2.6372E+01  7.3837E-01  3.0252E+01 -3.0845E+01 -1.2572E+00 -2.7520E+01  1.7981E+01
             8.6895E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -786.556699155615        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  1.1938E+00  3.2546E+00  1.7593E+00  1.5439E+00  3.0339E+00  2.4906E+00  3.2819E+00  1.6352E+00  1.0709E+01  2.3539E-01
             8.3510E+00
 PARAMETER:  2.7713E-01  1.2801E+00  6.6493E-01  5.3428E-01  1.2098E+00  1.0125E+00  1.2884E+00  5.9175E-01  2.4711E+00 -1.3465E+00
             2.2224E+00
 GRADIENT:   3.1942E+01  3.6418E+01  3.4379E+00  3.3317E+01  2.1658E+01  3.8250E+01 -7.9452E+00  3.3908E+00 -8.9909E+00  2.4137E-01
             5.5349E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -787.618250498576        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      500
 NPARAMETR:  1.1937E+00  3.2442E+00  1.7592E+00  1.5420E+00  3.0307E+00  2.4861E+00  3.3552E+00  1.3804E+00  1.0737E+01  1.8970E-01
             8.2252E+00
 PARAMETER:  2.7710E-01  1.2769E+00  6.6484E-01  5.3306E-01  1.2088E+00  1.0107E+00  1.3105E+00  4.2241E-01  2.4737E+00 -1.5623E+00
             2.2072E+00
 GRADIENT:   3.4503E+01  3.9595E+01  5.9862E+00  3.3384E+01  2.1971E+01  3.7180E+01  3.2833E+00 -8.0627E-01 -8.1030E+00  1.5938E-01
            -1.3875E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -787.631886149577        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      570
 NPARAMETR:  1.1935E+00  3.2281E+00  1.7586E+00  1.5394E+00  3.0229E+00  2.4789E+00  3.3285E+00  1.4348E+00  1.0774E+01  1.9029E-01
             8.2976E+00
 PARAMETER:  2.7692E-01  1.2719E+00  6.6452E-01  5.3140E-01  1.2062E+00  1.0078E+00  1.3025E+00  4.6100E-01  2.4771E+00 -1.5592E+00
             2.2160E+00
 GRADIENT:   3.3363E+01  3.6318E+01  5.4337E+00  3.3203E+01  2.1652E+01  3.6380E+01 -1.0327E+00  2.5042E-01 -7.3427E+00  1.6016E-01
            -2.3498E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -791.480012444001        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  1.1933E+00  3.2255E+00  1.7580E+00  1.5387E+00  3.0170E+00  2.4762E+00  3.6323E+00  1.5731E+00  1.0817E+01  2.2543E-01
             8.5120E+00
 PARAMETER:  2.7673E-01  1.2711E+00  6.6418E-01  5.3092E-01  1.2043E+00  1.0067E+00  1.3899E+00  5.5303E-01  2.4811E+00 -1.3897E+00
             2.2415E+00
 GRADIENT:   1.1559E+01 -2.8371E+01 -2.3856E+00  2.8140E+01  1.7989E+01  4.0490E+00 -2.3695E+01  2.0547E+00 -3.6887E+01  2.1530E-01
             1.5115E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -801.864892539271        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      871
 NPARAMETR:  1.1923E+00  3.4182E+00  1.7718E+00  1.5008E+00  2.9132E+00  2.4659E+00  3.7777E+00  1.5683E+00  1.5225E+01  1.0000E-02
             8.9855E+00
 PARAMETER:  2.7592E-01  1.3291E+00  6.7201E-01  5.0598E-01  1.1692E+00  1.0025E+00  1.4291E+00  5.5000E-01  2.8229E+00 -6.5408E+00
             2.2956E+00
 GRADIENT:   2.8907E+01 -4.1519E+01  2.7433E+01  1.9908E+01  1.0742E+01  3.1262E+01 -2.6358E+00 -1.1441E+01 -1.6659E+01  0.0000E+00
             6.2820E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -847.841158259070        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1053             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0956E+00  4.5092E+00  1.3234E+00  1.0023E-01  2.6376E+00  2.2310E+00  3.8200E+00  1.9378E+00  3.1365E+01  1.0000E-02
             8.3673E+00
 PARAMETER:  1.9134E-01  1.6061E+00  3.8021E-01 -2.2003E+00  1.0699E+00  9.0246E-01  1.4403E+00  7.6155E-01  3.5457E+00 -1.8229E+01
             2.2243E+00
 GRADIENT:   3.0621E+01  1.0718E+02  2.8376E+01  5.6540E+00  5.6806E+00  3.3528E+01  1.1197E+02  5.3446E+00  5.6607E+01  0.0000E+00
             1.4543E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -854.981837559552        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     1148
 NPARAMETR:  1.0371E+00  4.3363E+00  9.2618E-01  5.2128E-02  2.7093E+00  2.1159E+00  3.5199E+00  2.2812E+00  2.8612E+01  1.0000E-02
             8.3415E+00
 PARAMETER:  1.3641E-01  1.5670E+00  2.3313E-02 -2.8541E+00  1.0967E+00  8.4947E-01  1.3584E+00  9.2469E-01  3.4538E+00 -1.8229E+01
             2.2212E+00
 GRADIENT:  -5.6408E+00  1.2004E+01  1.3123E+00  1.0855E+00  8.5095E+00  1.9553E+00 -1.2566E+01  9.3089E+00 -2.1460E+01  0.0000E+00
            -1.0555E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -857.836382151959        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1333
 NPARAMETR:  1.0480E+00  4.2284E+00  9.1964E-01  3.3992E-02  2.7463E+00  2.0533E+00  3.5954E+00  1.7753E+00  3.5520E+01  1.0000E-02
             8.4150E+00
 PARAMETER:  1.4687E-01  1.5418E+00  1.6228E-02 -3.2816E+00  1.1102E+00  8.1946E-01  1.3797E+00  6.7394E-01  3.6701E+00 -1.8229E+01
             2.2300E+00
 GRADIENT:  -5.6790E-01  7.6398E+00 -4.1121E-01 -1.1487E+00  7.0016E+00 -4.1832E+00  1.2979E+00 -1.3254E-01 -2.8723E+01  0.0000E+00
            -1.8830E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -857.899382844535        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1511
 NPARAMETR:  1.0525E+00  4.2113E+00  9.1867E-01  3.3902E-02  2.7450E+00  2.0554E+00  3.5966E+00  1.8882E+00  3.6224E+01  1.0000E-02
             8.3998E+00
 PARAMETER:  1.5117E-01  1.5378E+00  1.5172E-02 -3.2843E+00  1.1098E+00  8.2048E-01  1.3800E+00  7.3560E-01  3.6897E+00 -1.8229E+01
             2.2282E+00
 GRADIENT:   1.3313E+00  6.4281E+00  8.5986E-01 -6.1457E-01  7.0860E+00 -4.2315E+00 -8.4714E-01  2.4124E+00 -2.7243E+01  0.0000E+00
            -3.4768E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -857.914107239663        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1688
 NPARAMETR:  1.0548E+00  4.1908E+00  8.8267E-01  3.3680E-02  2.7432E+00  2.0582E+00  3.6062E+00  1.9596E+00  3.7116E+01  1.0000E-02
             8.3808E+00
 PARAMETER:  1.5337E-01  1.5329E+00 -2.4805E-02 -3.2909E+00  1.1091E+00  8.2185E-01  1.3827E+00  7.7272E-01  3.7140E+00 -1.8229E+01
             2.2259E+00
 GRADIENT:   2.2807E+00  5.6063E+00 -8.9650E-01 -5.7498E-01  6.8620E+00 -4.0191E+00 -1.2817E-01  3.3207E+00 -2.8073E+01  0.0000E+00
            -5.2720E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -857.916024709307        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1868
 NPARAMETR:  1.0552E+00  4.1802E+00  8.6879E-01  3.3524E-02  2.7422E+00  2.0598E+00  3.6135E+00  1.9808E+00  3.7614E+01  1.0000E-02
             8.3717E+00
 PARAMETER:  1.5374E-01  1.5304E+00 -4.0652E-02 -3.2955E+00  1.1088E+00  8.2262E-01  1.3847E+00  7.8348E-01  3.7274E+00 -1.8229E+01
             2.2249E+00
 GRADIENT:   2.4885E+00  5.1221E+00 -1.5560E+00 -5.1968E-01  6.8019E+00 -3.9108E+00  4.1989E-01  3.4875E+00 -2.8371E+01  0.0000E+00
            -6.0716E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -857.916457177628        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2052
 NPARAMETR:  1.0553E+00  4.1723E+00  8.6059E-01  3.3393E-02  2.7414E+00  2.0611E+00  3.6197E+00  1.9930E+00  3.7997E+01  1.0000E-02
             8.3656E+00
 PARAMETER:  1.5387E-01  1.5285E+00 -5.0133E-02 -3.2994E+00  1.1085E+00  8.2322E-01  1.3864E+00  7.8966E-01  3.7375E+00 -1.8229E+01
             2.2241E+00
 GRADIENT:   2.5913E+00  4.7190E+00 -1.9151E+00 -4.6109E-01  6.7758E+00 -3.8279E+00  8.3800E-01  3.5668E+00 -2.8519E+01  0.0000E+00
            -6.5780E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -858.787002499212        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     2199
 NPARAMETR:  1.0477E+00  3.8959E+00  8.9545E-01  3.6366E-02  2.4949E+00  2.0799E+00  3.5510E+00  1.8229E+00  3.6245E+01  1.0000E-02
             8.3748E+00
 PARAMETER:  1.4660E-01  1.4599E+00 -1.0434E-02 -3.2141E+00  1.0142E+00  8.3230E-01  1.3672E+00  7.0041E-01  3.6903E+00 -1.8229E+01
             2.2252E+00
 GRADIENT:   9.9674E+00  7.9546E+01  9.8830E-01  2.9164E+00  2.3738E+00  1.3830E+01  8.9164E+01  4.4169E+00  5.4362E+01  0.0000E+00
             2.8891E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -858.835943243146        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2385
 NPARAMETR:  1.0472E+00  3.9634E+00  9.0541E-01  3.3626E-02  2.4758E+00  2.1101E+00  3.5843E+00  1.7552E+00  3.6768E+01  1.0000E-02
             8.3552E+00
 PARAMETER:  1.4612E-01  1.4771E+00  6.2901E-04 -3.2925E+00  1.0066E+00  8.4672E-01  1.3766E+00  6.6259E-01  3.7046E+00 -1.8229E+01
             2.2229E+00
 GRADIENT:   5.7851E-01 -1.3113E+00 -8.9625E-01  2.0484E-02  1.0336E-01  1.6528E+00 -5.2772E-01 -5.1182E-01 -2.3988E+01  0.0000E+00
            -9.8435E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -858.850983799818        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     2543
 NPARAMETR:  1.0452E+00  3.9417E+00  9.1060E-01  3.2229E-02  2.4736E+00  2.0955E+00  3.5616E+00  1.7758E+00  3.6969E+01  1.0000E-02
             8.3941E+00
 PARAMETER:  1.4392E-01  1.4681E+00  6.1441E-03 -3.3429E+00  1.0072E+00  8.4182E-01  1.3735E+00  6.8108E-01  3.7013E+00 -1.8229E+01
             2.2228E+00
 GRADIENT:  -2.7620E-01 -2.0090E+02 -3.3603E-02 -8.7328E+01  4.0134E-01  3.4169E+02  4.0369E+02  3.0896E-01 -1.7935E+02  0.0000E+00
            -1.1708E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2543
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.4861E-03 -2.9071E-02  1.3162E-02  2.9686E-02 -1.5279E-04
 SE:             2.7390E-02  2.6727E-02  1.0160E-02  1.0279E-02  4.6574E-05
 N:                     100         100         100         100         100

 P VAL.:         8.6990E-01  2.7673E-01  1.9517E-01  3.8786E-03  1.0358E-03

 ETASHRINKSD(%)  8.2389E+00  1.0460E+01  6.5962E+01  6.5563E+01  9.9844E+01
 ETASHRINKVR(%)  1.5799E+01  1.9826E+01  8.8414E+01  8.8141E+01  1.0000E+02
 EBVSHRINKSD(%)  9.4990E+00  2.4319E+00  6.7689E+01  6.0831E+01  9.9746E+01
 EBVSHRINKVR(%)  1.8096E+01  4.8046E+00  8.9560E+01  8.4658E+01  9.9999E+01
 RELATIVEINF(%)  8.1417E+01  8.7918E+01  1.0330E+01  1.4349E+01  6.3055E-04
 EPSSHRINKSD(%)  1.2596E+01
 EPSSHRINKVR(%)  2.3606E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -858.85098379981844     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       243.87525604578866     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    69.08
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -858.851       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  3.93E+00  9.10E-01  3.20E-02  2.48E+00  2.10E+00  3.57E+00  1.79E+00  3.66E+01  1.00E-02  8.35E+00
 


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
+        3.23E+05
 
 TH 2
+       -1.60E+04  5.10E+02
 
 TH 3
+        3.87E+00  8.00E+00  4.81E+01
 
 TH 4
+        2.16E+03  2.24E+04  9.18E+02  1.46E+06
 
 TH 5
+        1.95E+04 -3.37E+00 -1.06E+00 -1.42E+02  1.06E+03
 
 TH 6
+       -2.42E+04 -5.97E+02 -4.25E+01  1.70E+03  1.46E+03  4.28E+03
 
 TH 7
+       -8.93E+03 -4.87E+02 -4.37E+01 -1.40E+04  1.14E+03  1.56E+03  3.22E+02
 
 TH 8
+       -5.41E-01  3.27E-01  4.98E+00  2.37E+01 -4.09E-01 -1.48E+00 -1.15E+00  6.77E+00
 
 TH 9
+        3.23E+00  1.99E+01  1.48E+00  6.40E+02 -4.47E-02  2.41E+00 -1.42E+01  8.81E-02  9.25E-01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.60E+03  5.85E+01  4.99E+00  6.94E+03 -1.11E+00 -2.07E+02 -1.52E+02  1.49E-01  6.21E+00  0.00E+00  5.43E+01
 
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
 #CPUT: Total CPU Time in Seconds,       83.712
Stop Time:
Thu Sep 30 09:32:54 CDT 2021
