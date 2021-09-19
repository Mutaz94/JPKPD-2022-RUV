Sat Sep 18 07:45:20 CDT 2021
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
$DATA ../../../../data/int/D/dat93.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m93.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   54350.8073057763        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.1748E+02  8.4033E+02  7.1080E+01  7.9265E+02  4.8977E+01 -3.8440E+03 -1.9426E+03 -1.0874E+02 -2.7622E+03 -1.0445E+03
            -1.0572E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -395.964467661994        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.8980E-01  2.2422E+00  7.8069E-01  2.1282E+00  1.0916E+00  6.5883E+00  4.6157E+00  9.8048E-01  1.6093E+00  1.2607E+00
             1.2758E+01
 PARAMETER:  8.9750E-02  9.0747E-01 -1.4758E-01  8.5529E-01  1.8765E-01  1.9853E+00  1.6295E+00  8.0285E-02  5.7582E-01  3.3169E-01
             2.6462E+00
 GRADIENT:  -1.2309E+01  4.7250E+01 -6.9168E+01  2.8518E+02  5.3357E+01  1.5159E+02 -1.3309E+01  3.9057E+00 -7.5226E+01  1.7494E+01
            -1.1715E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -489.291600499590        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  6.4250E-01  4.3279E+00  4.8923E+01  3.7014E+00  2.7998E+00  3.6793E+00  1.0765E+01  5.7364E-01  2.3787E+00  1.3440E+00
             1.3047E+01
 PARAMETER: -3.4239E-01  1.5651E+00  3.9903E+00  1.4087E+00  1.1295E+00  1.4027E+00  2.4763E+00 -4.5576E-01  9.6656E-01  3.9566E-01
             2.6686E+00
 GRADIENT:  -7.9934E+01  4.4296E+01 -2.6207E+00  1.8956E+02 -1.8104E+00  7.6097E+01 -4.4792E+01  1.1200E-02 -7.6568E+01  2.2076E+01
            -2.0744E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -659.972395798255        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.1639E+00  1.5726E+00  3.3438E+01  1.4193E+00  3.0944E+00  2.7594E+00  6.9518E+00  1.3701E+00  7.5544E-01  7.2906E-01
             1.3882E+01
 PARAMETER:  2.5177E-01  5.5273E-01  3.6097E+00  4.5015E-01  1.2296E+00  1.1150E+00  2.0390E+00  4.1487E-01 -1.8045E-01 -2.1600E-01
             2.7306E+00
 GRADIENT:   1.8492E+00  1.5381E+01 -2.2311E+00 -1.3071E+01  3.7765E+01  7.7088E+01  3.1981E+01  4.7228E-02  3.0048E+00  5.5094E+00
             1.2977E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -695.902710882571        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0194E+00  8.1490E-01  5.6278E+01  1.5370E+00  2.7523E+00  1.6910E+00  6.8738E+00  1.3132E+00  8.9183E-01  4.7322E-01
             1.3084E+01
 PARAMETER:  1.1926E-01 -1.0468E-01  4.1303E+00  5.2981E-01  1.1124E+00  6.2531E-01  2.0277E+00  3.7250E-01 -1.4484E-02 -6.4820E-01
             2.6714E+00
 GRADIENT:  -2.6210E+01  1.9790E+00 -9.0827E-01  4.6288E+00  1.5398E+01 -1.4474E+01  1.8297E+00  1.9901E-02 -2.9477E+00  2.4596E+00
             4.4243E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -699.071031392148        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0459E+00  6.2212E-01  5.9815E+01  1.6034E+00  2.5457E+00  1.7372E+00  7.2328E+00  1.0761E+00  1.0611E+00  3.4532E-01
             1.2579E+01
 PARAMETER:  1.4492E-01 -3.7462E-01  4.1913E+00  5.7214E-01  1.0344E+00  6.5225E-01  2.0786E+00  1.7337E-01  1.5929E-01 -9.6329E-01
             2.6320E+00
 GRADIENT:   2.7008E+00  4.3244E-01 -4.0766E-01  6.3850E+00 -4.7149E+00 -3.1935E+00  4.8587E-01  1.1357E-02 -1.9260E+00  1.2375E+00
            -2.4492E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -700.611841701666        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  1.0513E+00  5.1837E-01  1.0776E+02  1.6659E+00  2.6289E+00  1.7541E+00  7.5265E+00  1.0651E+00  1.1186E+00  7.5887E-02
             1.2796E+01
 PARAMETER:  1.5004E-01 -5.5706E-01  4.7799E+00  6.1037E-01  1.0666E+00  6.6197E-01  2.1184E+00  1.6307E-01  2.1211E-01 -2.4785E+00
             2.6491E+00
 GRADIENT:   5.5941E-01 -1.6258E+00 -3.4003E-01 -1.2604E+00  1.0269E+00 -4.5197E-01  8.1764E-01  3.6659E-03 -4.6674E-01  5.8892E-02
             7.4100E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -702.363833384236        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      624
 NPARAMETR:  1.0545E+00  2.7353E-01  5.6219E+02  1.8581E+00  2.6688E+00  1.7641E+00  9.5471E+00  9.4495E-01  1.2638E+00  1.0000E-02
             1.2837E+01
 PARAMETER:  1.5310E-01 -1.1963E+00  6.4318E+00  7.1955E-01  1.0816E+00  6.6764E-01  2.3562E+00  4.3375E-02  3.3415E-01 -6.0284E+00
             2.6524E+00
 GRADIENT:  -2.1139E-01 -2.4858E-01 -9.8504E-02 -2.5209E+00  4.2301E-02 -2.8375E-01  4.7772E-01  1.3128E-04 -1.9957E-01  0.0000E+00
             3.7001E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -702.407954466991        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      804
 NPARAMETR:  1.0554E+00  2.7765E-01  1.0061E+03  1.8561E+00  2.6738E+00  1.7681E+00  9.5198E+00  9.3157E-01  1.2621E+00  1.0000E-02
             1.2828E+01
 PARAMETER:  1.5393E-01 -1.1814E+00  7.0138E+00  7.1846E-01  1.0835E+00  6.6990E-01  2.3534E+00  2.9121E-02  3.3280E-01 -6.8337E+00
             2.6516E+00
 GRADIENT:   1.2069E+00 -1.7077E-01 -5.3682E-02 -2.2967E+00 -1.9414E-01  1.4763E+00  4.6214E-01  4.3967E-05 -3.4115E-01  0.0000E+00
             1.7355E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -702.468124985582        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      980
 NPARAMETR:  1.0549E+00  2.8138E-01  4.4051E+04  1.8594E+00  2.6822E+00  1.7647E+00  9.4938E+00  8.3661E-01  1.2692E+00  1.0000E-02
             1.2816E+01
 PARAMETER:  1.5341E-01 -1.1680E+00  1.0793E+01  7.2024E-01  1.0866E+00  6.6797E-01  2.3506E+00 -7.8399E-02  3.3836E-01 -1.2287E+01
             2.6507E+00
 GRADIENT:  -1.2523E+01 -1.0608E+00 -1.2737E-03  1.1638E+00  5.3363E-01  7.6918E+00  3.5692E-01 -9.9374E-06 -1.0719E-01  0.0000E+00
             5.9557E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -702.468507041848        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1177
 NPARAMETR:  1.0548E+00  2.8217E-01  4.4694E+04  1.8608E+00  2.6832E+00  1.7649E+00  9.4986E+00  8.3650E-01  1.2697E+00  1.0000E-02
             1.2818E+01
 PARAMETER:  1.5332E-01 -1.1653E+00  1.0808E+01  7.2103E-01  1.0870E+00  6.6809E-01  2.3511E+00 -7.8527E-02  3.3875E-01 -1.2292E+01
             2.6508E+00
 GRADIENT:  -7.0956E+00 -6.7754E-02 -1.2052E-03  2.2322E+00  2.2893E-02  3.0655E+00  1.6015E-01  8.9991E-05 -1.3972E+00  0.0000E+00
             2.3983E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -702.468821171778        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1368
 NPARAMETR:  1.0546E+00  2.8149E-01  4.4997E+04  1.8615E+00  2.6830E+00  1.7647E+00  9.5030E+00  8.3655E-01  1.2707E+00  1.0000E-02
             1.2815E+01
 PARAMETER:  1.5315E-01 -1.1677E+00  1.0814E+01  7.2136E-01  1.0869E+00  6.6800E-01  2.3516E+00 -7.8465E-02  3.3958E-01 -1.2292E+01
             2.6506E+00
 GRADIENT:  -4.2308E+00 -3.6287E-02 -1.1926E-03  6.3379E-01  8.8615E-03  2.3152E+00  1.6263E-01 -8.7173E-05 -3.9450E-01  0.0000E+00
             2.5404E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -702.469284562428        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1559
 NPARAMETR:  1.0539E+00  2.8003E-01  5.0342E+04  1.8629E+00  2.6823E+00  1.7657E+00  9.5129E+00  8.3417E-01  1.2718E+00  1.0000E-02
             1.2814E+01
 PARAMETER:  1.5246E-01 -1.1729E+00  1.0927E+01  7.2212E-01  1.0867E+00  6.6857E-01  2.3527E+00 -8.1317E-02  3.4040E-01 -1.2292E+01
             2.6506E+00
 GRADIENT:  -5.3389E+00 -1.7341E-01 -1.0770E-03  8.6395E-01 -9.4694E-02  2.1757E+00  8.3463E-02 -1.6934E-04 -3.1668E-01  0.0000E+00
             2.0317E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -702.469432474944        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1730
 NPARAMETR:  1.0537E+00  2.7975E-01  5.1482E+04  1.8623E+00  2.6829E+00  1.7653E+00  9.5200E+00  8.3413E-01  1.2727E+00  1.0000E-02
             1.2815E+01
 PARAMETER:  1.5234E-01 -1.1738E+00  1.0949E+01  7.2180E-01  1.0869E+00  6.6833E-01  2.3534E+00 -8.1371E-02  3.4115E-01 -1.2292E+01
             2.6506E+00
 GRADIENT:  -5.6265E+00  5.4799E-01 -1.0295E-03  4.3486E+00  2.2512E-01  2.3118E+00  2.1935E+01  3.9575E-04 -7.6577E-01  0.0000E+00
             6.6783E+00

0ITERATION NO.:   69    OBJECTIVE VALUE:  -702.469488433706        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     1875
 NPARAMETR:  1.0539E+00  2.7944E-01  5.1636E+04  1.8626E+00  2.6837E+00  1.7653E+00  9.5247E+00  8.3436E-01  1.2722E+00  1.0000E-02
             1.2817E+01
 PARAMETER:  1.5238E-01 -1.1741E+00  1.0947E+01  7.2214E-01  1.0868E+00  6.6848E-01  2.3533E+00 -8.1033E-02  3.4088E-01 -1.2292E+01
             2.6505E+00
 GRADIENT:  -4.1749E-01  3.5148E-02 -3.1151E-03  7.9068E-02 -6.6200E-02  8.3453E-02 -5.6580E-02  4.8521E-01  8.0672E-02  0.0000E+00
            -1.9731E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1875
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1395E-02  7.0841E-02  2.9235E-07 -9.2528E-02 -7.2202E-07
 SE:             2.7761E-02  2.0496E-02  2.9552E-07  1.6340E-02  8.9604E-05
 N:                     100         100         100         100         100

 P VAL.:         4.4089E-01  5.4772E-04  3.2254E-01  1.4950E-08  9.9357E-01

 ETASHRINKSD(%)  6.9964E+00  3.1335E+01  9.9999E+01  4.5258E+01  9.9700E+01
 ETASHRINKVR(%)  1.3503E+01  5.2851E+01  1.0000E+02  7.0034E+01  9.9999E+01
 EBVSHRINKSD(%)  1.0568E+01  3.1261E+01  9.9999E+01  3.6832E+01  9.9624E+01
 EBVSHRINKVR(%)  2.0019E+01  5.2750E+01  1.0000E+02  6.0099E+01  9.9999E+01
 RELATIVEINF(%)  7.9642E+01  2.9912E+01  4.3497E-09  2.4258E+01  2.7752E-04
 EPSSHRINKSD(%)  3.8640E+00
 EPSSHRINKVR(%)  7.5786E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -702.46948843370649     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       951.61987133470427     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    62.89
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.47
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -702.469       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  2.80E-01  5.14E+04  1.86E+00  2.68E+00  1.77E+00  9.52E+00  8.34E-01  1.27E+00  1.00E-02  1.28E+01
 


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
+        1.14E+03
 
 TH 2
+       -2.72E+02  2.14E+02
 
 TH 3
+        2.48E-05  1.05E-04  9.50E-11
 
 TH 4
+       -1.25E+02  7.19E+01 -5.32E-06  1.22E+02
 
 TH 5
+       -2.36E+01  1.19E+00  2.46E-06 -9.42E+00  2.14E+01
 
 TH 6
+        6.33E+01 -1.90E+01 -2.51E-05 -1.14E+01 -2.53E+00  1.03E+02
 
 TH 7
+       -4.41E-01  5.82E+00 -9.95E-07 -4.37E+00  3.52E-01 -6.53E-01  1.15E+00
 
 TH 8
+       -2.39E+02 -4.55E+02 -2.53E-04  4.18E+01 -4.65E+01  1.41E+02 -4.41E+00  2.63E+03
 
 TH 9
+        1.42E+02 -5.66E+01  7.07E-05 -2.11E+01 -1.12E-01  4.76E+01 -4.15E-01  3.98E+02  2.18E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.12E+01 -1.95E+00  6.43E-08 -8.62E+00  1.30E+00 -3.34E+00  1.92E-01 -2.73E-01  1.83E+00  0.00E+00  5.80E+00
 
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
 #CPUT: Total CPU Time in Seconds,       80.504
Stop Time:
Sat Sep 18 07:46:42 CDT 2021
