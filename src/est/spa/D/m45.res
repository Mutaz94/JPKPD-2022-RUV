Wed Sep 29 20:03:04 CDT 2021
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
$DATA ../../../../data/spa/D/dat45.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12601.7802022585        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.9728E+02  2.8103E+02 -6.0024E+01  2.7452E+02  1.7319E+02 -1.5815E+03 -7.8067E+02 -4.3926E+01 -1.2037E+03 -3.1159E+02
            -2.4132E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -554.356238855480        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0898E+00  1.1024E+00  1.0076E+00  1.4323E+00  1.2478E+00  1.7069E+00  1.3448E+00  9.7460E-01  1.4184E+00  1.0589E+00
             1.4550E+01
 PARAMETER:  1.8600E-01  1.9746E-01  1.0756E-01  4.5926E-01  3.2140E-01  6.3470E-01  3.9622E-01  7.4270E-02  4.4955E-01  1.5721E-01
             2.7776E+00
 GRADIENT:  -7.9321E+01  1.7235E+01 -5.5775E+00  3.4380E+01 -9.8572E+00  1.6655E+01  1.9018E+00  3.0047E+00  1.5600E+01  3.0328E+00
             1.8872E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -569.219214101315        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.1595E+00  1.0124E+00  1.4555E+00  1.5981E+00  2.6778E+00  1.5128E+00  2.8290E+00  6.3531E-01  1.0925E+00  6.5889E+00
             1.3553E+01
 PARAMETER:  2.4796E-01  1.1236E-01  4.7536E-01  5.6880E-01  1.0850E+00  5.1394E-01  1.1399E+00 -3.5364E-01  1.8846E-01  1.9854E+00
             2.7066E+00
 GRADIENT:  -1.9354E+01  3.1497E+01  2.1388E+00  5.7991E+01 -1.4153E+01 -1.9096E+01  8.5088E+00  2.4891E-01  5.9539E+00  7.5328E+00
             1.1991E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -597.803771659089        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.8919E-01  4.4791E-01  7.8084E-01  1.4419E+00  4.0271E+00  1.6806E+00  2.8182E+00  7.8817E-02  7.4450E-01  7.5345E+00
             1.0045E+01
 PARAMETER:  8.9128E-02 -7.0316E-01 -1.4739E-01  4.6598E-01  1.4930E+00  6.1915E-01  1.1361E+00 -2.4406E+00 -1.9504E-01  2.1195E+00
             2.4071E+00
 GRADIENT:  -6.8604E+00  1.7521E+01 -9.0208E+00  3.5913E+01 -6.5219E+00  2.7372E+01  1.3721E+00 -2.4853E-04 -1.0967E+01  2.2599E+00
            -4.2832E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -610.573879605241        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  8.8554E-01  1.8937E-01  2.7298E-01  1.2289E+00  1.3805E+01  1.2070E+00  4.9847E-01  1.0000E-02  6.3370E-01  6.6874E+00
             1.0738E+01
 PARAMETER: -2.1561E-02 -1.5640E+00 -1.1984E+00  3.0614E-01  2.7251E+00  2.8814E-01 -5.9621E-01 -7.2126E+00 -3.5619E-01  2.0002E+00
             2.4738E+00
 GRADIENT:  -1.9031E+01  5.9456E+00 -1.9470E+01  1.1827E+02 -3.2048E+00 -1.1674E+02  5.3388E-01  0.0000E+00  2.3457E+00  1.6103E+01
            -2.2680E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -649.594888963750        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  7.7148E-01  1.2191E-01  1.7218E-01  1.0121E+00  2.7263E+01  1.7160E+00  1.0649E-01  1.0000E-02  6.5833E-01  6.8453E+00
             9.5560E+00
 PARAMETER: -1.5944E-01 -2.0045E+00 -1.6592E+00  1.1198E-01  3.4055E+00  6.4001E-01 -2.1397E+00 -1.1709E+01 -3.1805E-01  2.0236E+00
             2.3572E+00
 GRADIENT:  -1.3758E+01  1.4303E+01  6.8067E+00  3.2640E+01  3.1708E+00  1.0262E+01  5.0420E-02  0.0000E+00 -5.2913E+00 -7.7416E+00
            -2.3768E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -654.041728059238        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      497             RESET HESSIAN, TYPE I
 NPARAMETR:  7.9472E-01  9.2155E-02  1.3515E-01  9.0209E-01  3.5587E+01  1.7054E+00  1.0000E-02  1.0000E-02  5.9703E-01  7.8187E+00
             1.0037E+01
 PARAMETER: -1.2976E-01 -2.2843E+00 -1.9014E+00 -3.0452E-03  3.6720E+00  6.3378E-01 -4.9053E+00 -1.3455E+01 -4.1578E-01  2.1565E+00
             2.4063E+00
 GRADIENT:   4.1032E+01  6.6888E+00 -4.1138E+00  2.3484E+01  2.8787E+00  1.3557E+01  0.0000E+00  0.0000E+00 -4.9823E+00 -5.2425E+00
             6.1217E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -655.647552097032        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      660             RESET HESSIAN, TYPE I
 NPARAMETR:  7.6186E-01  9.1392E-02  1.3600E-01  8.9336E-01  3.2608E+01  1.6808E+00  1.0000E-02  1.0000E-02  6.9086E-01  8.1484E+00
             1.0274E+01
 PARAMETER: -1.7200E-01 -2.2926E+00 -1.8951E+00 -1.2761E-02  3.5846E+00  6.1924E-01 -5.1922E+00 -1.3455E+01 -2.6982E-01  2.1978E+00
             2.4296E+00
 GRADIENT:   7.9537E+00  9.2198E+00  1.1232E+01  7.2254E+00  1.4492E+00  1.0786E+01  0.0000E+00  0.0000E+00 -4.1614E+00 -3.2747E+00
             4.2748E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -656.640339299234        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  7.5996E-01  9.0658E-02  1.3621E-01  8.8963E-01  3.2815E+01  1.6593E+00  1.0000E-02  1.0000E-02  7.8373E-01  8.1167E+00
             1.0247E+01
 PARAMETER: -1.7449E-01 -2.3007E+00 -1.8936E+00 -1.6952E-02  3.5909E+00  6.0642E-01 -5.1922E+00 -1.3455E+01 -1.4369E-01  2.1939E+00
             2.4270E+00
 GRADIENT:   5.1402E+00  8.8076E+00  1.6990E+01  3.6154E+00  6.5161E-01  4.5475E+00  0.0000E+00  0.0000E+00 -4.3679E+00 -2.3891E+00
             5.1014E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -656.969465713970        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      861
 NPARAMETR:  7.6085E-01  9.1469E-02  1.3355E-01  8.8198E-01  3.1892E+01  1.6727E+00  1.0000E-02  1.0000E-02  7.8295E-01  8.2676E+00
             9.7652E+00
 PARAMETER: -1.7332E-01 -2.2918E+00 -1.9133E+00 -2.5588E-02  3.5624E+00  6.1444E-01 -5.1922E+00 -1.3455E+01 -1.4469E-01  2.2123E+00
             2.3788E+00
 GRADIENT:   1.6474E+01  1.2895E+01  1.0168E+01  3.7411E+00  1.1347E+00  5.7399E+00  0.0000E+00  0.0000E+00 -5.0656E+00 -2.7852E+00
             1.6825E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -656.974728436246        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      932
 NPARAMETR:  7.6074E-01  9.0839E-02  1.3210E-01  8.7532E-01  3.1637E+01  1.6618E+00  1.0000E-02  1.0000E-02  7.8282E-01  8.3184E+00
             9.6667E+00
 PARAMETER: -1.7347E-01 -2.2987E+00 -1.9242E+00 -3.3167E-02  3.5543E+00  6.0790E-01 -5.1922E+00 -1.3455E+01 -1.4485E-01  2.2185E+00
             2.3687E+00
 GRADIENT:   2.0375E+01  1.3650E+01  9.6238E+00  2.5509E+00  1.0447E+00  3.3099E+00  0.0000E+00  0.0000E+00 -5.8675E+00 -2.6055E+00
             8.9250E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -657.122853560892        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:     1049
 NPARAMETR:  7.6030E-01  8.9087E-02  1.2979E-01  8.6478E-01  3.1201E+01  1.6737E+00  1.0000E-02  1.0000E-02  7.8270E-01  8.4065E+00
             9.5666E+00
 PARAMETER: -1.7405E-01 -2.3181E+00 -1.9419E+00 -4.5277E-02  3.5405E+00  6.1501E-01 -5.1922E+00 -1.3455E+01 -1.4500E-01  2.2290E+00
             2.3583E+00
 GRADIENT:   1.4769E+01  1.2236E+01 -5.7737E-01 -3.1424E+00  1.0225E+00 -5.0141E+00  0.0000E+00  0.0000E+00 -7.1632E+00 -1.5227E+00
            -1.7425E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -657.774871131497        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1190
 NPARAMETR:  7.4230E-01  8.5340E-02  1.3015E-01  8.6427E-01  3.1746E+01  1.6770E+00  1.0000E-02  1.0000E-02  7.8344E-01  8.3166E+00
             9.6580E+00
 PARAMETER: -1.9800E-01 -2.3611E+00 -1.9390E+00 -4.5871E-02  3.5578E+00  6.1701E-01 -5.1922E+00 -1.3455E+01 -1.4406E-01  2.2183E+00
             2.3678E+00
 GRADIENT:  -2.9256E+00  8.6242E+00  6.4305E+00 -1.0708E+00  1.3575E+00 -3.7156E+00  0.0000E+00  0.0000E+00 -5.8859E+00 -1.8903E+00
            -5.2815E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -661.740125305558        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1333
 NPARAMETR:  7.2738E-01  8.1294E-02  1.2730E-01  8.5262E-01  3.1257E+01  1.6747E+00  1.0000E-02  1.0000E-02  1.0525E+00  8.3980E+00
             9.3794E+00
 PARAMETER: -2.1831E-01 -2.4097E+00 -1.9612E+00 -5.9445E-02  3.5422E+00  6.1565E-01 -5.1922E+00 -1.3455E+01  1.5113E-01  2.2280E+00
             2.3385E+00
 GRADIENT:   4.5110E+00  8.5580E+00  2.3478E+01  3.4327E+00 -3.7857E-01  1.5319E+00  0.0000E+00  0.0000E+00  2.1420E+00  1.8475E-01
             2.4076E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -662.059529562937        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1468
 NPARAMETR:  7.3315E-01  8.0904E-02  1.2653E-01  8.5065E-01  3.1270E+01  1.7316E+00  1.0000E-02  1.0000E-02  1.0446E+00  8.3979E+00
             9.2660E+00
 PARAMETER: -2.1040E-01 -2.4145E+00 -1.9673E+00 -6.1749E-02  3.5427E+00  6.4903E-01 -5.1922E+00 -1.3455E+01  1.4365E-01  2.2280E+00
             2.3264E+00
 GRADIENT:  -4.5798E-01  7.2173E+00  8.4981E+00  1.6541E+00 -3.9638E-01  5.6993E-01  0.0000E+00  0.0000E+00 -9.1595E-01  1.7650E-01
            -3.8650E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -669.415985433034        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1647
 NPARAMETR:  5.5649E-01  3.2359E-02  5.2494E-02  4.9318E-01  3.4647E+01  1.6122E+00  1.0000E-02  1.0000E-02  1.0014E+00  8.2916E+00
             9.0810E+00
 PARAMETER: -4.8610E-01 -3.3309E+00 -2.8471E+00 -6.0688E-01  3.6452E+00  5.7763E-01 -5.1922E+00 -1.3455E+01  1.0136E-01  2.2152E+00
             2.3062E+00
 GRADIENT:  -6.2879E-01  5.0388E+00  1.0961E+00 -7.5590E+00 -2.1129E-01  6.6250E-01  0.0000E+00  0.0000E+00 -2.1574E+00  5.8376E-02
             1.1703E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -669.912888740824        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:     1763
 NPARAMETR:  5.3601E-01  2.3338E-02  4.8659E-02  4.6768E-01  9.5166E+01  1.5998E+00  1.0000E-02  1.0000E-02  1.0194E+00  2.9504E+00
             8.9919E+00
 PARAMETER: -5.2361E-01 -3.6577E+00 -2.9229E+00 -6.5998E-01  4.6556E+00  5.6990E-01 -5.1922E+00 -1.3455E+01  1.1922E-01  1.1819E+00
             2.2963E+00
 GRADIENT:   3.6173E+01  4.4409E-01  4.2779E+01  1.3992E+01 -1.3511E-03  8.9089E+00  0.0000E+00  0.0000E+00  1.3463E+00 -4.2314E-05
             1.5403E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -669.914446355859        NO. OF FUNC. EVALS.:  82
 CUMULATIVE NO. OF FUNC. EVALS.:     1845
 NPARAMETR:  5.3357E-01  2.2279E-02  4.8261E-02  4.6514E-01  1.1980E+02  1.5991E+00  1.0000E-02  1.0000E-02  1.0195E+00  2.2992E+00
             8.9830E+00
 PARAMETER: -5.2816E-01 -3.7041E+00 -2.9311E+00 -6.6541E-01  4.8858E+00  5.6943E-01 -5.1922E+00 -1.3455E+01  1.1935E-01  9.3256E-01
             2.2953E+00
 GRADIENT:   3.5811E+01  2.9496E-01  4.3139E+01  1.4695E+01 -4.7259E-04  8.8579E+00  0.0000E+00  0.0000E+00  1.3266E+00 -2.3534E-05
             1.4934E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -669.956506247960        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2027
 NPARAMETR:  5.3849E-01  2.1877E-02  4.8305E-02  4.6669E-01  1.0080E+02  1.6032E+00  1.0000E-02  1.0000E-02  1.0131E+00  2.7506E+00
             9.0476E+00
 PARAMETER: -5.1899E-01 -3.7223E+00 -2.9302E+00 -6.6209E-01  4.7131E+00  5.7199E-01 -5.1922E+00 -1.3455E+01  1.1300E-01  1.1118E+00
             2.3025E+00
 GRADIENT:   1.2681E+00  9.1418E-02 -3.6735E-01 -1.0370E+00  3.6805E-03 -1.6227E-01  0.0000E+00  0.0000E+00 -4.6322E-02 -6.3889E-05
            -1.4224E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -669.966287511059        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2216
 NPARAMETR:  5.3644E-01  1.9823E-02  4.8528E-02  4.6820E-01  2.9072E+01  1.6048E+00  1.0000E-02  1.0000E-02  1.0126E+00  1.0112E+01
             9.0427E+00
 PARAMETER: -5.2280E-01 -3.8209E+00 -2.9256E+00 -6.5886E-01  3.4698E+00  5.7302E-01 -5.1922E+00 -1.3455E+01  1.1255E-01  2.4138E+00
             2.3020E+00
 GRADIENT:  -3.0474E+00  2.9607E-02  1.5116E+00 -1.5019E+00  1.6512E-02 -2.0097E-01  0.0000E+00  0.0000E+00  1.0584E-01 -1.0350E-02
             4.6148E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -669.983182065863        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2397
 NPARAMETR:  5.3782E-01  1.9614E-02  4.8500E-02  4.6872E-01  1.7486E+01  1.6059E+00  1.0000E-02  1.0000E-02  1.0118E+00  8.8940E+00
             9.0459E+00
 PARAMETER: -5.2024E-01 -3.8315E+00 -2.9262E+00 -6.5775E-01  2.9614E+00  5.7369E-01 -5.1922E+00 -1.3455E+01  1.1178E-01  2.2854E+00
             2.3023E+00
 GRADIENT:  -1.8471E+00  8.5117E-02  2.6510E-01 -4.5335E-01 -4.7693E-03 -7.9026E-02  0.0000E+00  0.0000E+00  2.3235E-02  6.8705E-03
             4.1314E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -669.989224064319        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2559
 NPARAMETR:  5.3906E-01  1.7829E-02  4.8568E-02  4.6956E-01  1.4360E+01  1.6072E+00  1.0000E-02  1.0000E-02  1.0108E+00  7.5257E+00
             9.0458E+00
 PARAMETER: -5.1794E-01 -3.9269E+00 -2.9248E+00 -6.5597E-01  2.7644E+00  5.7449E-01 -5.1922E+00 -1.3455E+01  1.1078E-01  2.1183E+00
             2.3023E+00
 GRADIENT:  -2.1039E-01  2.4843E-03 -5.1987E-01 -1.5673E-02 -1.4273E-02  9.2141E-04  0.0000E+00  0.0000E+00 -7.2975E-02  1.4761E-02
            -1.0119E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2559
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3958E-03 -3.7415E-07  6.9697E-05 -2.3213E-02 -8.7312E-04
 SE:             2.8700E-02  1.5022E-06  1.7965E-04  2.4005E-02  8.3632E-04
 N:                     100         100         100         100         100

 P VAL.:         9.3347E-01  8.0331E-01  6.9805E-01  3.3352E-01  2.9648E-01

 ETASHRINKSD(%)  3.8526E+00  9.9995E+01  9.9398E+01  1.9582E+01  9.7198E+01
 ETASHRINKVR(%)  7.5567E+00  1.0000E+02  9.9996E+01  3.5329E+01  9.9922E+01
 EBVSHRINKSD(%)  3.3703E+00  9.9994E+01  9.9423E+01  1.9718E+01  9.7735E+01
 EBVSHRINKVR(%)  6.6270E+00  1.0000E+02  9.9997E+01  3.5548E+01  9.9949E+01
 RELATIVEINF(%)  5.9452E+00  5.8924E-08  5.5988E-05  1.0069E+00  6.2813E-03
 EPSSHRINKSD(%)  1.4068E+01
 EPSSHRINKVR(%)  2.6157E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -669.98922406431905     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       65.161602499419132     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.55
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -669.989       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.39E-01  1.78E-02  4.86E-02  4.70E-01  1.44E+01  1.61E+00  1.00E-02  1.00E-02  1.01E+00  7.53E+00  9.05E+00
 


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
+        1.37E+03
 
 TH 2
+       -3.66E+02  1.61E+03
 
 TH 3
+       -4.38E+03 -2.06E+02  1.48E+05
 
 TH 4
+       -1.01E+02  1.44E+02 -1.94E+04  3.00E+03
 
 TH 5
+        2.91E-01  1.02E+00  4.70E+02  9.52E-02  1.67E+00
 
 TH 6
+        2.64E+00  1.24E+01 -4.15E+01 -2.15E+01  4.52E-02  6.41E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -9.31E-01 -2.43E+01  5.84E+02 -8.20E+01  1.32E-01 -8.45E-01  0.00E+00  0.00E+00  8.48E+01
 
 TH10
+        2.24E-01  2.41E+03 -1.20E+03  3.02E-01 -4.24E+00  3.59E-03  0.00E+00  0.00E+00 -2.30E-01  1.07E+01
 
 TH11
+       -1.76E+01  4.20E+00  1.23E+02 -1.15E+01  3.20E+00  1.53E+00  0.00E+00  0.00E+00  2.96E+00 -8.13E+00  4.58E+00
 
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
 #CPUT: Total CPU Time in Seconds,       45.540
Stop Time:
Wed Sep 29 20:03:51 CDT 2021
