Wed Sep 29 12:38:14 CDT 2021
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
$DATA ../../../../data/spa/A2/dat19.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -935.276966426347        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6549E+02  3.7294E+01  8.9280E+01 -2.4796E+00  6.9247E+01  5.8456E+01 -8.4027E+00 -3.2982E+01 -3.0566E+00 -8.8269E+01
            -1.3603E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1406.18451293571        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0419E+00  9.9369E-01  8.9269E-01  1.0624E+00  9.1849E-01  8.9120E-01  8.9903E-01  9.5186E-01  8.8909E-01  1.1947E+00
             2.0387E+00
 PARAMETER:  1.4104E-01  9.3666E-02 -1.3519E-02  1.6050E-01  1.4979E-02 -1.5191E-02 -6.4343E-03  5.0660E-02 -1.7552E-02  2.7792E-01
             8.1232E-01
 GRADIENT:   1.4828E+02  2.7097E+01  1.1307E+01  3.3157E+01 -5.9047E+00 -5.8591E+00 -1.2170E-01  3.1688E+00 -9.5957E+00  5.2897E-01
            -1.6549E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1416.76715035196        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0362E+00  9.4114E-01  4.3896E-01  1.0379E+00  6.0188E-01  9.0519E-01  4.0206E-01  1.5392E-01  9.9355E-01  7.6466E-01
             2.1665E+00
 PARAMETER:  1.3553E-01  3.9337E-02 -7.2334E-01  1.3720E-01 -4.0769E-01  3.8792E-04 -8.1115E-01 -1.7713E+00  9.3527E-02 -1.6832E-01
             8.7312E-01
 GRADIENT:   8.8791E+01  3.0779E+01  2.3969E+00  6.1467E+01  3.7162E+01 -5.8806E+00 -4.3294E+00 -1.9081E-02 -1.7929E+00 -2.4212E+01
            -1.2380E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1430.71269483942        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      320
 NPARAMETR:  1.0612E+00  9.0502E-01  3.5546E-01  1.0368E+00  5.2598E-01  9.4545E-01  6.0367E-01  2.7867E-02  8.9230E-01  7.4974E-01
             2.5604E+00
 PARAMETER:  1.5938E-01  2.0300E-04 -9.3434E-01  1.3610E-01 -5.4249E-01  4.3911E-02 -4.0474E-01 -3.4803E+00 -1.3957E-02 -1.8803E-01
             1.0402E+00
 GRADIENT:   2.0038E+01  1.5009E+01 -1.9189E+01  5.3491E+01  1.9368E+01  3.9535E+00 -2.3405E+00  5.6196E-03 -4.2039E+00  6.9101E+00
            -9.2513E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1434.84451232466        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      498
 NPARAMETR:  1.0422E+00  6.8472E-01  3.1714E-01  1.0899E+00  4.1801E-01  9.2868E-01  1.0525E+00  5.1123E-02  8.2478E-01  5.4000E-01
             2.5635E+00
 PARAMETER:  1.4129E-01 -2.7874E-01 -1.0484E+00  1.8611E-01 -7.7224E-01  2.6004E-02  1.5112E-01 -2.8735E+00 -9.2639E-02 -5.1618E-01
             1.0414E+00
 GRADIENT:  -2.1505E+01  5.5429E+00  9.4041E-01  5.5478E+00  3.3689E+00 -2.1257E+00  1.8221E+00 -3.3132E-03  1.2605E+00 -6.5065E-01
             6.6527E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1436.17749208312        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  1.0456E+00  6.3265E-01  2.3250E-01  1.0354E+00  3.4630E-01  9.4445E-01  1.0066E+00  5.3632E-02  8.7997E-01  4.7349E-01
             2.4466E+00
 PARAMETER:  1.4456E-01 -3.5784E-01 -1.3589E+00  1.3476E-01 -9.6046E-01  4.2843E-02  1.0661E-01 -2.8256E+00 -2.7871E-02 -6.4762E-01
             9.9469E-01
 GRADIENT:   4.0631E-01 -1.1833E-01  7.7367E-02 -1.5320E+00  1.2799E+00  1.8218E-01  2.1321E-02 -3.6229E-02  5.1273E-01 -1.2497E+00
            -3.4230E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1436.22487098107        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      848
 NPARAMETR:  1.0467E+00  6.1157E-01  2.3710E-01  1.0490E+00  3.4275E-01  9.4407E-01  1.0263E+00  5.5915E-02  8.7313E-01  4.9421E-01
             2.4435E+00
 PARAMETER:  1.4560E-01 -3.9173E-01 -1.3393E+00  1.4784E-01 -9.7075E-01  4.2449E-02  1.2600E-01 -2.7839E+00 -3.5665E-02 -6.0480E-01
             9.9342E-01
 GRADIENT:   1.0612E-01 -1.4932E-03  1.0493E-01  9.5222E-03 -7.4552E-02  3.9837E-02  2.2831E-02 -3.4969E-02  5.7579E-02 -5.9238E-02
            -2.1827E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1437.64369303476        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      946
 NPARAMETR:  1.0519E+00  6.1496E-01  2.3714E-01  1.0489E+00  3.4199E-01  9.4710E-01  9.8495E-01  6.6232E-01  8.6309E-01  4.2650E-01
             2.3888E+00
 PARAMETER:  1.5059E-01 -3.8620E-01 -1.3391E+00  1.4775E-01 -9.7297E-01  4.5652E-02  8.4840E-02 -3.1201E-01 -4.7231E-02 -7.5214E-01
             9.7078E-01
 GRADIENT:   1.0506E+02  1.0995E+01  2.7411E+01  2.5927E+01  6.2531E+01  6.3407E+00 -3.5653E+00 -1.4973E+00  2.6842E+00  8.8886E-01
             3.4238E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1438.86876557082        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     1022
 NPARAMETR:  1.0452E+00  6.7890E-01  2.2059E-01  1.0075E+00  3.5022E-01  9.3917E-01  1.0271E+00  1.0464E+00  8.2340E-01  2.3417E-01
             2.3700E+00
 PARAMETER:  1.4424E-01 -2.8728E-01 -1.4114E+00  1.0747E-01 -9.4919E-01  3.7240E-02  1.2677E-01  1.4535E-01 -9.4313E-02 -1.3517E+00
             9.6288E-01
 GRADIENT:   9.8215E+01  1.7301E+01  2.6112E+01  2.9798E+01  6.0234E+01  5.2170E+00  8.6355E-01 -4.9925E-01 -3.0088E+00 -1.8153E-01
             9.7175E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1439.17825988766        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:     1136
 NPARAMETR:  1.0410E+00  7.2687E-01  2.0593E-01  9.6950E-01  3.5700E-01  9.3675E-01  9.5883E-01  1.1318E+00  8.6159E-01  1.6918E-01
             2.3609E+00
 PARAMETER:  1.4019E-01 -2.1901E-01 -1.4802E+00  6.9030E-02 -9.3002E-01  3.4659E-02  5.7960E-02  2.2384E-01 -4.8973E-02 -1.6768E+00
             9.5904E-01
 GRADIENT:   1.0006E+01  3.5455E+00  1.4611E+00  7.4744E+00  4.3863E+00  5.7852E-01 -2.6344E-01 -1.3952E-02 -3.1569E-01 -5.2296E-01
             2.0990E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1439.48899344025        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1313
 NPARAMETR:  1.0394E+00  7.2658E-01  2.0402E-01  9.6588E-01  3.5570E-01  9.3709E-01  9.2320E-01  1.0892E+00  8.7635E-01  2.9054E-01
             2.3327E+00
 PARAMETER:  1.3861E-01 -2.1941E-01 -1.4895E+00  6.5284E-02 -9.3368E-01  3.5025E-02  2.0086E-02  1.8548E-01 -3.1987E-02 -1.1360E+00
             9.4703E-01
 GRADIENT:   7.8028E+00  1.0828E+00  9.8145E-01  1.8611E+00  6.2867E+00  3.5598E-01 -3.4869E-01 -8.7466E-02  7.6608E-02  8.4697E-02
             9.6459E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1439.53541928248        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1448
 NPARAMETR:  1.0296E+00  7.2197E-01  2.0147E-01  9.6406E-01  3.5235E-01  9.3354E-01  9.3086E-01  1.0945E+00  8.7464E-01  2.7781E-01
             2.3168E+00
 PARAMETER:  1.2912E-01 -2.2577E-01 -1.5021E+00  6.3403E-02 -9.4314E-01  3.1225E-02  2.8357E-02  1.9034E-01 -3.3945E-02 -1.1808E+00
             9.4020E-01
 GRADIENT:   6.9569E+01  1.3343E+01  2.2808E+01  1.4166E+01  6.5320E+01  4.0028E+00 -5.7826E-02  1.8868E-02 -8.3274E-02  4.7380E-01
             6.7414E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1439.53900516717        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:     1526
 NPARAMETR:  1.0258E+00  7.1668E-01  1.9755E-01  9.6112E-01  3.4806E-01  9.3256E-01  9.3196E-01  1.0954E+00  8.7972E-01  2.7232E-01
             2.3061E+00
 PARAMETER:  1.2552E-01 -2.3313E-01 -1.5217E+00  6.0348E-02 -9.5537E-01  3.0181E-02  2.9539E-02  1.9113E-01 -2.8150E-02 -1.2008E+00
             9.3557E-01
 GRADIENT:   6.0225E+01  1.4983E+01  2.2418E+01  1.3164E+01  6.5003E+01  3.3840E+00 -9.0663E-02  3.4241E-02 -4.0623E-01  4.2650E-01
             6.2105E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1439.68847207964        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1705
 NPARAMETR:  1.0342E+00  7.0680E-01  1.9338E-01  9.6004E-01  3.4226E-01  9.3620E-01  9.3364E-01  1.0964E+00  8.9325E-01  2.7615E-01
             2.3042E+00
 PARAMETER:  1.3366E-01 -2.4701E-01 -1.5431E+00  5.9215E-02 -9.7219E-01  3.4079E-02  3.1335E-02  1.9200E-01 -1.2888E-02 -1.1868E+00
             9.3472E-01
 GRADIENT:  -5.6292E-01  4.5637E+00 -2.4478E+00 -1.1610E+00 -8.1733E-01 -5.3325E-01  8.8739E-02  3.7055E-01 -6.1866E-03  4.3380E-02
             4.1644E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1439.99621806517        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1881
 NPARAMETR:  1.0374E+00  6.3431E-01  1.9566E-01  9.8866E-01  3.2159E-01  9.4137E-01  1.0039E+00  1.0539E+00  8.9216E-01  2.8173E-01
             2.2781E+00
 PARAMETER:  1.3669E-01 -3.5521E-01 -1.5314E+00  8.8592E-02 -1.0345E+00  3.9585E-02  1.0387E-01  1.5252E-01 -1.4112E-02 -1.1668E+00
             9.2332E-01
 GRADIENT:  -2.9600E-01  8.8325E-02  1.2631E-01 -1.5909E-01 -3.1595E-01  4.7131E-02 -6.2069E-02  3.5039E-02  1.0200E-01  5.0060E-02
            -1.9658E-02

0ITERATION NO.:   72    OBJECTIVE VALUE:  -1439.99625258056        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1938
 NPARAMETR:  1.0375E+00  6.3475E-01  1.9573E-01  9.8860E-01  3.2180E-01  9.4118E-01  1.0054E+00  1.0551E+00  8.9138E-01  2.7914E-01
             2.2786E+00
 PARAMETER:  1.3679E-01 -3.5453E-01 -1.5310E+00  8.8531E-02 -1.0338E+00  3.9383E-02  1.0536E-01  1.5361E-01 -1.4986E-02 -1.1760E+00
             9.2356E-01
 GRADIENT:  -5.1351E-02  3.3068E-02  2.7480E-02 -9.9656E-02 -8.3437E-02 -9.2753E-03 -1.4378E-03  3.4891E-02  4.5250E-02  1.2789E-02
            -1.2184E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1938
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8760E-04  3.8445E-03 -8.0467E-03 -9.7922E-03  3.5048E-03
 SE:             2.9261E-02  2.0772E-02  1.6330E-02  2.4886E-02  9.3973E-03
 N:                     100         100         100         100         100

 P VAL.:         9.9488E-01  8.5316E-01  6.2220E-01  6.9397E-01  7.0918E-01

 ETASHRINKSD(%)  1.9734E+00  3.0412E+01  4.5291E+01  1.6628E+01  6.8518E+01
 ETASHRINKVR(%)  3.9078E+00  5.1576E+01  7.0069E+01  3.0490E+01  9.0089E+01
 EBVSHRINKSD(%)  2.0971E+00  3.0231E+01  4.5438E+01  1.6197E+01  6.8977E+01
 EBVSHRINKVR(%)  4.1503E+00  5.1322E+01  7.0229E+01  2.9770E+01  9.0375E+01
 RELATIVEINF(%)  9.1821E+01  2.5022E+00  5.5654E+00  2.3293E+01  3.8665E-01
 EPSSHRINKSD(%)  3.7718E+01
 EPSSHRINKVR(%)  6.1210E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1439.9962525805593     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -704.84542601682108     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.78
 Elapsed covariance  time in seconds:     6.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1439.996       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  6.35E-01  1.96E-01  9.89E-01  3.22E-01  9.41E-01  1.01E+00  1.06E+00  8.91E-01  2.79E-01  2.28E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.06E-02  1.95E-01  2.86E-02  6.16E-02  6.91E-02  6.70E-02  2.01E-01  2.74E-01  1.41E-01  1.60E-01  1.88E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.38E-04
 
 TH 2
+       -1.86E-03  3.80E-02
 
 TH 3
+       -1.72E-04  1.54E-03  8.16E-04
 
 TH 4
+        5.32E-05 -9.19E-03  4.00E-04  3.80E-03
 
 TH 5
+       -6.33E-04  1.30E-02  9.94E-04 -2.75E-03  4.78E-03
 
 TH 6
+       -2.83E-04 -1.18E-04 -1.91E-04 -7.70E-05 -2.55E-04  4.49E-03
 
 TH 7
+        8.13E-04 -2.08E-02  2.01E-03  8.62E-03 -5.36E-03 -2.48E-03  4.03E-02
 
 TH 8
+       -3.42E-04  1.42E-02 -5.44E-03 -9.27E-03  7.43E-04  2.70E-03 -2.94E-02  7.53E-02
 
 TH 9
+        8.48E-04 -1.52E-02 -2.90E-03  1.74E-03 -6.88E-03  2.09E-03 -1.53E-03  1.82E-02  1.99E-02
 
 TH10
+        6.05E-04 -8.10E-03  3.92E-04  1.71E-03 -2.46E-03  3.64E-04 -1.69E-03 -6.04E-03  4.84E-03  2.55E-02
 
 TH11
+        2.54E-03  5.65E-03  1.47E-03 -2.17E-03  3.03E-03  3.00E-03 -1.75E-04 -1.01E-02 -3.76E-03  5.22E-03  3.52E-02
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.06E-02
 
 TH 2
+       -3.11E-01  1.95E-01
 
 TH 3
+       -1.97E-01  2.76E-01  2.86E-02
 
 TH 4
+        2.82E-02 -7.65E-01  2.27E-01  6.16E-02
 
 TH 5
+       -2.99E-01  9.62E-01  5.04E-01 -6.45E-01  6.91E-02
 
 TH 6
+       -1.38E-01 -9.06E-03 -9.99E-02 -1.86E-02 -5.50E-02  6.70E-02
 
 TH 7
+        1.32E-01 -5.31E-01  3.51E-01  6.97E-01 -3.86E-01 -1.84E-01  2.01E-01
 
 TH 8
+       -4.07E-02  2.65E-01 -6.94E-01 -5.48E-01  3.92E-02  1.47E-01 -5.34E-01  2.74E-01
 
 TH 9
+        1.96E-01 -5.54E-01 -7.19E-01  2.00E-01 -7.06E-01  2.21E-01 -5.39E-02  4.71E-01  1.41E-01
 
 TH10
+        1.24E-01 -2.60E-01  8.59E-02  1.74E-01 -2.23E-01  3.40E-02 -5.28E-02 -1.38E-01  2.15E-01  1.60E-01
 
 TH11
+        4.42E-01  1.54E-01  2.74E-01 -1.88E-01  2.34E-01  2.38E-01 -4.65E-03 -1.97E-01 -1.42E-01  1.74E-01  1.88E-01
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.21E+03
 
 TH 2
+        1.78E+02  2.19E+03
 
 TH 3
+        5.42E+02  3.23E+03  1.86E+04
 
 TH 4
+        1.94E+02 -1.55E+02 -3.46E+03  2.03E+03
 
 TH 5
+        8.37E+01 -6.90E+03 -1.50E+04  2.06E+03  2.47E+04
 
 TH 6
+        3.05E+02 -1.13E+02 -3.49E+02 -4.35E+01  5.01E+02  3.34E+02
 
 TH 7
+       -5.54E+00  5.03E+00 -4.01E+02 -4.19E+01  1.68E+02  3.91E+01  7.56E+01
 
 TH 8
+       -2.23E+01 -8.78E+01  2.35E+02  4.83E+01  1.06E+02 -1.42E+01 -9.81E+00  6.29E+01
 
 TH 9
+        7.35E+01 -1.24E+02  1.51E+02 -1.40E+02  6.30E+02 -1.22E+01  9.27E+00 -5.35E+01  2.54E+02
 
 TH10
+        1.43E+01 -2.10E+01 -4.75E+02  7.82E+01  2.39E+02  2.49E+01  2.93E+01 -2.17E+00 -2.57E+01  6.50E+01
 
 TH11
+       -2.32E+02  6.05E+01 -6.89E+01  9.39E+01 -2.49E+02 -7.25E+01 -9.88E+00  1.36E+01 -3.80E+01 -8.65E+00  7.29E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       29.081
Stop Time:
Wed Sep 29 12:38:45 CDT 2021
