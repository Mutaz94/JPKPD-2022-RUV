Sat Sep 25 02:05:48 CDT 2021
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
$DATA ../../../../data/int/SL3/dat26.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      969
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

 TOT. NO. OF OBS RECS:      869
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   225.735943311988        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1125E+01  9.4121E+00  1.8742E+02  7.6550E+01  1.3466E+02  2.6126E+01 -1.6593E+02 -4.2318E+02 -1.4233E+02 -6.8746E+01
            -7.1172E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2623.90452311059        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0255E+00  1.2494E+00  9.6966E-01  8.8449E-01  1.0989E+00  8.6716E-01  1.1081E+00  9.8983E-01  8.9586E-01  1.0561E+00
             2.8494E+00
 PARAMETER:  1.2521E-01  3.2265E-01  6.9187E-02 -2.2747E-02  1.9430E-01 -4.2529E-02  2.0262E-01  8.9782E-02 -9.9719E-03  1.5460E-01
             1.1471E+00
 GRADIENT:   2.9410E+01  1.6914E+00 -6.1451E+00  6.7491E+00 -1.3689E+01 -2.1901E+01  1.0997E+01  5.7283E+00 -2.9883E+00 -1.0195E+01
             7.9390E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2630.37621831311        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0272E+00  1.5432E+00  1.3722E+00  7.6249E-01  1.3852E+00  8.7727E-01  9.5142E-01  7.3233E-01  9.5138E-01  1.4611E+00
             2.8557E+00
 PARAMETER:  1.2683E-01  5.3388E-01  4.1639E-01 -1.7117E-01  4.2586E-01 -3.0940E-02  5.0203E-02 -2.1152E-01  5.0163E-02  4.7920E-01
             1.1493E+00
 GRADIENT:   3.1527E+01  7.4674E+01  4.9685E+00  5.9355E+01 -2.7257E+01 -1.7268E+01  1.0250E+01 -2.8597E-01 -7.3416E+00  1.1520E+01
             9.2033E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2643.14163984610        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0116E+00  1.7514E+00  1.2360E+00  5.6514E-01  1.5971E+00  9.0772E-01  7.4046E-01  4.2793E-01  1.3406E+00  1.5359E+00
             2.6848E+00
 PARAMETER:  1.1150E-01  6.6039E-01  3.1185E-01 -4.7068E-01  5.6822E-01  3.1765E-03 -2.0048E-01 -7.4881E-01  3.9309E-01  5.2911E-01
             1.0876E+00
 GRADIENT:  -5.7976E-01 -1.3526E+01  2.5844E-01  8.9598E+00  5.2762E+00 -4.8476E+00 -1.8260E+00 -3.4030E-01 -4.9192E-01  5.1317E-01
            -1.5348E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2650.91285989946        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  1.0124E+00  2.1187E+00  8.5782E-01  3.2024E-01  1.7765E+00  9.2906E-01  6.7677E-01  2.4321E-01  1.9433E+00  1.6601E+00
             2.6772E+00
 PARAMETER:  1.1228E-01  8.5081E-01 -5.3355E-02 -1.0387E+00  6.7462E-01  2.6418E-02 -2.9043E-01 -1.3138E+00  7.6441E-01  6.0686E-01
             1.0848E+00
 GRADIENT:   1.5778E+00  1.0914E+01  1.6421E+00  4.3038E+00 -6.3654E+00  3.4315E+00 -4.4654E-01 -4.1868E-02 -1.7044E-01 -5.2761E-01
            -2.2345E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2654.97763185988        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      464
 NPARAMETR:  1.0157E+00  2.4022E+00  4.2693E-01  1.3763E-01  1.9766E+00  9.3657E-01  6.5850E-01  7.7215E-02  2.8502E+00  1.7961E+00
             2.6668E+00
 PARAMETER:  1.1561E-01  9.7639E-01 -7.5115E-01 -1.8832E+00  7.8137E-01  3.4467E-02 -3.1779E-01 -2.4612E+00  1.1474E+00  6.8560E-01
             1.0809E+00
 GRADIENT:   1.9247E+00 -4.0714E+00 -1.9995E+00  2.2454E+00  4.5527E+00  5.2933E+00 -8.2236E-01  5.0244E-03 -4.1010E-02  1.0423E+00
            -4.7619E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2656.68119637435        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      641
 NPARAMETR:  1.0149E+00  2.4963E+00  3.5209E-01  7.1089E-02  2.0458E+00  9.4481E-01  6.5805E-01  4.7854E-02  3.5969E+00  1.8690E+00
             2.6523E+00
 PARAMETER:  1.1481E-01  1.0148E+00 -9.4386E-01 -2.5438E+00  8.1581E-01  4.3224E-02 -3.1848E-01 -2.9396E+00  1.3801E+00  7.2542E-01
             1.0754E+00
 GRADIENT:   6.9680E-01 -7.4956E+00 -1.4491E+00  7.3604E-01  3.9054E+00  8.4090E+00 -1.0419E+00  1.2954E-03 -2.7830E+00  2.2396E+00
            -1.1312E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2659.57922075880        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      819
 NPARAMETR:  1.0158E+00  2.5813E+00  3.2142E-01  1.5747E-02  2.0783E+00  9.0952E-01  6.5123E-01  1.3358E-02  8.0072E+00  1.8641E+00
             2.6771E+00
 PARAMETER:  1.1563E-01  1.0483E+00 -1.0350E+00 -4.0511E+00  8.3155E-01  5.1662E-03 -3.2889E-01 -4.2157E+00  2.1803E+00  7.2278E-01
             1.0847E+00
 GRADIENT:   1.1705E+00  7.8686E+00 -1.2678E+00 -1.2022E+00  5.5118E+00 -5.3845E+00  3.8736E+00  8.0410E-05 -4.0448E+00 -8.3226E-01
             1.0141E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2659.71996796641        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:      987
 NPARAMETR:  1.0148E+00  2.5809E+00  3.2267E-01  1.5515E-02  2.0568E+00  9.2385E-01  6.4348E-01  1.3079E-02  8.0842E+00  1.8681E+00
             2.6660E+00
 PARAMETER:  1.1466E-01  1.0481E+00 -1.0311E+00 -4.0659E+00  8.2116E-01  2.0797E-02 -3.4086E-01 -4.2367E+00  2.1899E+00  7.2491E-01
             1.0806E+00
 GRADIENT:  -6.6590E-01  8.4160E+00 -1.2521E+00 -1.1164E+00 -1.3043E-01  3.6504E-01 -6.3910E-02  7.8524E-05 -3.8770E+00 -1.2892E-02
             7.7640E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2659.73332168135        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1170
 NPARAMETR:  1.0150E+00  2.5777E+00  3.2283E-01  1.5541E-02  2.0600E+00  9.1913E-01  6.4361E-01  1.0000E-02  8.1036E+00  1.8628E+00
             2.6672E+00
 PARAMETER:  1.1486E-01  1.0469E+00 -1.0306E+00 -4.0643E+00  8.2269E-01  1.5669E-02 -3.4066E-01 -2.3365E+01  2.1923E+00  7.2209E-01
             1.0810E+00
 GRADIENT:  -1.7003E-01  2.3967E+00 -1.2750E+00 -1.0956E+00  9.6358E-01 -1.4819E+00 -3.4378E-02  0.0000E+00 -3.7347E+00 -7.8903E-01
             1.7852E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2659.75258862090        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1345
 NPARAMETR:  1.0149E+00  2.5724E+00  3.2310E-01  1.5586E-02  2.0562E+00  9.2283E-01  6.4462E-01  1.0000E-02  8.1368E+00  1.8682E+00
             2.6644E+00
 PARAMETER:  1.1478E-01  1.0448E+00 -1.0298E+00 -4.0614E+00  8.2085E-01  1.9688E-02 -3.3910E-01 -5.5737E+01  2.1964E+00  7.2496E-01
             1.0800E+00
 GRADIENT:  -1.8493E-03 -7.6091E+00 -1.3415E+00 -8.4627E-01  3.9008E-02  1.9522E-02  2.6067E-02  0.0000E+00 -3.1251E+00 -4.0737E-02
            -5.5616E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2659.76882242326        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1526
 NPARAMETR:  1.0150E+00  2.5753E+00  3.2324E-01  1.5641E-02  2.0550E+00  9.2202E-01  6.4728E-01  1.0000E-02  8.1707E+00  1.8753E+00
             2.6658E+00
 PARAMETER:  1.1488E-01  1.0460E+00 -1.0294E+00 -4.0579E+00  8.2028E-01  1.8814E-02 -3.3498E-01 -6.5251E+01  2.2006E+00  7.2878E-01
             1.0805E+00
 GRADIENT:   1.2786E-01 -2.8280E+00 -1.3730E+00 -5.0471E-01 -5.8721E-01 -3.1947E-01  9.7358E-01  0.0000E+00 -2.5679E+00  9.6837E-01
             1.2312E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2659.77575795755        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1702
 NPARAMETR:  1.0150E+00  2.5783E+00  3.2330E-01  1.5676E-02  2.0580E+00  9.2294E-01  6.4877E-01  1.0000E-02  8.1919E+00  1.8701E+00
             2.6650E+00
 PARAMETER:  1.1485E-01  1.0471E+00 -1.0292E+00 -4.0556E+00  8.2172E-01  1.9809E-02 -3.3268E-01 -6.8011E+01  2.2031E+00  7.2599E-01
             1.0802E+00
 GRADIENT:   2.6858E-03 -4.1660E-01 -1.5716E+00  8.0784E-01 -5.1070E-03  4.4864E-04 -3.9910E-02  0.0000E+00 -3.3948E-01 -2.2658E-03
            -3.3119E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2659.79953626769        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1794
 NPARAMETR:  1.0124E+00  2.5806E+00  3.3895E-01  1.5613E-02  2.0478E+00  9.2137E-01  6.4898E-01  1.0000E-02  8.1742E+00  1.8623E+00
             2.6628E+00
 PARAMETER:  1.1229E-01  1.0480E+00 -9.8191E-01 -4.0596E+00  8.1675E-01  1.8105E-02 -3.3235E-01 -6.8011E+01  2.2010E+00  7.2182E-01
             1.0794E+00
 GRADIENT:   7.9594E-01  5.8325E+01 -9.9195E-01 -7.8872E-01  3.7186E-01  7.7188E-02  3.0261E+00  0.0000E+00 -2.2916E+00  2.0771E-01
             2.7170E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2659.80888318788        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1892
 NPARAMETR:  1.0114E+00  2.5801E+00  3.4758E-01  1.5615E-02  2.0435E+00  9.2080E-01  6.4898E-01  1.0000E-02  8.1746E+00  1.8585E+00
             2.6618E+00
 PARAMETER:  1.1134E-01  1.0478E+00 -9.5675E-01 -4.0595E+00  8.1467E-01  1.7491E-02 -3.3235E-01 -6.8011E+01  2.2010E+00  7.1975E-01
             1.0790E+00
 GRADIENT:  -8.8404E+00  9.0231E+00 -8.4022E-01 -1.0446E+00 -3.2564E+00 -9.3882E-01  2.3079E+00  0.0000E+00 -3.6707E+00 -9.6135E-01
            -2.3456E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2659.83201457400        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     2033
 NPARAMETR:  1.0149E+00  2.5787E+00  3.4742E-01  1.5647E-02  2.0427E+00  9.2290E-01  6.4887E-01  1.0000E-02  8.1837E+00  1.8578E+00
             2.6604E+00
 PARAMETER:  1.1477E-01  1.0473E+00 -9.5723E-01 -4.0575E+00  8.1426E-01  1.9770E-02 -3.3252E-01 -6.8011E+01  2.2021E+00  7.1939E-01
             1.0785E+00
 GRADIENT:  -3.2989E-02  6.2950E+00 -8.6130E-01 -1.0048E+00 -3.4587E+00 -3.4702E-02  2.1937E+00  0.0000E+00 -3.5512E+00 -1.0704E+00
            -3.8171E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2659.91060241193        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:     2149
 NPARAMETR:  1.0149E+00  2.5776E+00  3.6334E-01  1.5637E-02  2.0467E+00  9.2314E-01  6.4842E-01  1.0000E-02  8.2724E+00  1.8603E+00
             2.6617E+00
 PARAMETER:  1.1474E-01  1.0469E+00 -9.1242E-01 -4.0581E+00  8.1623E-01  2.0029E-02 -3.3322E-01 -6.8011E+01  2.2129E+00  7.2074E-01
             1.0790E+00
 GRADIENT:   7.3248E+00  5.1509E+01 -6.7654E-01  1.3331E-01  2.2753E-01  8.3347E-01  1.3194E+00  0.0000E+00 -4.5787E-01  1.2141E-02
            -6.1297E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2659.92113080118        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:     2299
 NPARAMETR:  1.0149E+00  2.5773E+00  3.6566E-01  1.5603E-02  2.0493E+00  9.2308E-01  6.4841E-01  1.0000E-02  8.2895E+00  1.8605E+00
             2.6617E+00
 PARAMETER:  1.1474E-01  1.0467E+00 -9.0604E-01 -4.0603E+00  8.1748E-01  1.9958E-02 -3.3324E-01 -6.8011E+01  2.2150E+00  7.2084E-01
             1.0790E+00
 GRADIENT:  -9.9541E-02 -1.4281E+00 -8.3533E-01  1.2047E+00 -1.8760E+00  4.0283E-02 -1.1914E+00  0.0000E+00  4.2233E-01 -9.1551E-01
            -2.6963E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2659.92445835730        NO. OF FUNC. EVALS.: 206
 CUMULATIVE NO. OF FUNC. EVALS.:     2505             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0149E+00  2.5772E+00  3.6564E-01  1.5606E-02  2.0561E+00  9.2306E-01  6.4841E-01  1.0000E-02  8.2905E+00  1.8619E+00
             2.6615E+00
 PARAMETER:  1.1481E-01  1.0467E+00 -9.0610E-01 -4.0601E+00  8.2081E-01  1.9939E-02 -3.3323E-01 -6.8011E+01  2.2151E+00  7.2161E-01
             1.0789E+00
 GRADIENT:   7.5296E+00  4.2697E+01 -1.1088E+00  3.3436E+00  2.3430E+00  7.8657E-01 -3.1646E+00  0.0000E+00  5.2117E+00 -3.2788E-01
            -1.2166E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2659.92580531441        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     2661
 NPARAMETR:  1.0149E+00  2.5772E+00  3.6564E-01  1.5606E-02  2.0569E+00  9.2305E-01  6.4835E-01  1.0000E-02  8.2905E+00  1.8674E+00
             2.6615E+00
 PARAMETER:  1.1483E-01  1.0467E+00 -9.0610E-01 -4.0601E+00  8.2119E-01  1.9932E-02 -3.3333E-01 -6.8011E+01  2.2151E+00  7.2454E-01
             1.0789E+00
 GRADIENT:   6.4285E-02 -4.8575E-01 -7.1802E-01  5.2032E-01  1.7935E-02  2.3175E-02 -2.7899E-01  0.0000E+00 -7.6602E-01 -2.5214E-02
            -2.4931E+00

0ITERATION NO.:   96    OBJECTIVE VALUE:  -2659.92580531441        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     2683
 NPARAMETR:  1.0149E+00  2.5772E+00  3.6564E-01  1.5606E-02  2.0569E+00  9.2305E-01  6.4835E-01  1.0000E-02  8.2905E+00  1.8674E+00
             2.6615E+00
 PARAMETER:  1.1483E-01  1.0467E+00 -9.0610E-01 -4.0601E+00  8.2119E-01  1.9932E-02 -3.3333E-01 -6.8011E+01  2.2151E+00  7.2454E-01
             1.0789E+00
 GRADIENT:   6.4285E-02 -4.8575E-01 -7.1802E-01  5.2032E-01  1.7935E-02  2.3175E-02 -2.7899E-01  0.0000E+00 -7.6602E-01 -2.5214E-02
            -2.4931E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2683
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0820E-04 -1.0797E-02  4.3496E-06  7.6787E-03 -1.7537E-02
 SE:             2.9317E-02  2.7838E-02  4.4427E-06  6.2233E-03  2.6626E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8345E-01  6.9811E-01  3.2756E-01  2.1725E-01  5.1012E-01

 ETASHRINKSD(%)  1.7843E+00  6.7404E+00  9.9985E+01  7.9151E+01  1.0799E+01
 ETASHRINKVR(%)  3.5368E+00  1.3026E+01  1.0000E+02  9.5653E+01  2.0431E+01
 EBVSHRINKSD(%)  1.7992E+00  6.2774E+00  9.9949E+01  8.5220E+01  8.8655E+00
 EBVSHRINKVR(%)  3.5660E+00  1.2161E+01  1.0000E+02  9.7816E+01  1.6945E+01
 RELATIVEINF(%)  9.6352E+01  4.0464E+01  2.6092E-05  1.0142E+00  8.2241E+01
 EPSSHRINKSD(%)  1.5860E+01
 EPSSHRINKVR(%)  2.9205E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          869
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1597.1151707097210     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2659.9258053144131     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1062.8106346046920     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    71.83
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2659.926       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  2.58E+00  3.66E-01  1.56E-02  2.06E+00  9.23E-01  6.48E-01  1.00E-02  8.29E+00  1.87E+00  2.66E+00
 


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
+        1.23E+03
 
 TH 2
+        1.33E+02  2.31E+05
 
 TH 3
+        1.15E+03  1.91E+01  1.40E+07
 
 TH 4
+       -6.40E+03  4.54E+04 -1.29E+03  4.19E+08
 
 TH 5
+       -3.77E+00 -2.48E+02 -1.18E+03  1.01E+04  6.62E+01
 
 TH 6
+       -5.55E+00  3.54E+01  3.27E+02 -1.72E+03 -8.72E-01  2.28E+02
 
 TH 7
+        2.00E+03 -7.54E+02  2.14E+07  3.02E+04  4.20E+06  5.19E+02  3.61E+07
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.20E+01  1.56E+02 -3.75E+00  1.29E+05  3.48E+01 -5.85E+00  1.05E+02  0.00E+00  4.99E+03
 
 TH10
+       -1.35E+00 -4.68E+02  3.42E+06  1.99E+04 -6.86E-01  8.11E-01  5.24E+06  0.00E+00  6.86E+01  8.36E+05
 
 TH11
+        1.16E+02 -1.00E+03  1.61E+06  4.19E+04 -1.34E+02  4.18E+01  2.71E+06  0.00E+00  1.45E+02  3.93E+05  2.04E+05
 
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
 #CPUT: Total CPU Time in Seconds,       86.715
Stop Time:
Sat Sep 25 02:07:17 CDT 2021
