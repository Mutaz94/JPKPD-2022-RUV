Wed Sep 29 19:36:44 CDT 2021
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
$DATA ../../../../data/spa/D/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1424.31365708769        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0439E+02  5.5894E+01 -5.1830E+01  1.2134E+02  9.7760E+01 -5.5628E+01 -8.6840E+01 -7.9887E+00 -1.6960E+02 -1.8303E+01
            -5.3047E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1487.58581424266        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8039E-01  8.9110E-01  1.3238E+00  1.1110E+00  9.8922E-01  1.0165E+00  1.4875E+00  1.1270E+00  1.3725E+00  9.8588E-01
             1.1103E+00
 PARAMETER:  8.0193E-02 -1.5304E-02  3.8051E-01  2.0530E-01  8.9162E-02  1.1635E-01  4.9708E-01  2.1956E-01  4.1664E-01  8.5780E-02
             2.0461E-01
 GRADIENT:   3.8101E+02  6.3282E+01  9.4302E+00  2.1295E+02 -2.4769E+01 -4.4770E+01  9.9459E+00 -1.1068E+01  4.6371E+01 -2.8418E+00
             2.9159E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1496.29683442256        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  9.7395E-01  8.6518E-01  1.5014E+00  1.1546E+00  1.1453E+00  1.3794E+00  2.8581E+00  1.4816E+00  1.0291E+00  9.9004E-01
             1.1906E+00
 PARAMETER:  7.3607E-02 -4.4822E-02  5.0638E-01  2.4374E-01  2.3564E-01  4.2162E-01  1.1501E+00  4.9315E-01  1.2872E-01  8.9995E-02
             2.7448E-01
 GRADIENT:  -1.0128E+01  3.0840E+01 -2.9965E+01  6.5463E+01  9.6455E+01  3.0964E+01  2.8764E+01 -2.0511E+00 -1.3622E+01 -1.8095E+01
             2.1152E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1508.94590207951        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      402
 NPARAMETR:  9.7832E-01  6.6253E-01  1.4148E+00  1.1791E+00  9.8305E-01  1.3075E+00  2.6304E+00  1.1553E+00  1.1091E+00  9.4050E-01
             1.1023E+00
 PARAMETER:  7.8077E-02 -3.1169E-01  4.4701E-01  2.6471E-01  8.2902E-02  3.6809E-01  1.0672E+00  2.4439E-01  2.0356E-01  3.8656E-02
             1.9739E-01
 GRADIENT:  -2.4090E+00  2.0415E+00  1.2783E+01 -1.0814E+01 -3.0126E+00  1.1833E+01  7.6887E-01 -7.5807E+00 -9.3406E+00 -4.9190E+00
            -1.9358E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1510.04832938058        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  9.7883E-01  5.9157E-01  1.4978E+00  1.2421E+00  9.8317E-01  1.2689E+00  2.7804E+00  1.3987E+00  1.1638E+00  9.5684E-01
             1.0932E+00
 PARAMETER:  7.8600E-02 -4.2497E-01  5.0403E-01  3.1677E-01  8.3026E-02  3.3816E-01  1.1226E+00  4.3556E-01  2.5170E-01  5.5880E-02
             1.8911E-01
 GRADIENT:  -7.7418E-01  2.4397E+00 -7.6071E-01  3.1695E+00  6.7032E-01  5.4902E-02  1.1352E+00 -1.4219E-01  2.6661E+00 -1.0281E-01
             2.3985E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1510.13918059147        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  9.7499E-01  4.4853E-01  1.6878E+00  1.3482E+00  9.8593E-01  1.2736E+00  3.2623E+00  1.5626E+00  1.1276E+00  9.7666E-01
             1.0938E+00
 PARAMETER:  7.4677E-02 -7.0178E-01  6.2342E-01  3.9877E-01  8.5834E-02  3.4187E-01  1.2824E+00  5.4635E-01  2.2010E-01  7.6386E-02
             1.8970E-01
 GRADIENT:  -3.7724E+00  4.7678E+00  3.0910E-02  9.5419E+00 -3.2511E+00  2.0522E+00  2.3631E+00  6.8577E-01  3.0015E+00 -5.7450E-02
             3.4776E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1510.16971015482        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      929
 NPARAMETR:  9.7404E-01  3.7379E-01  1.8208E+00  1.4023E+00  9.9395E-01  1.2749E+00  3.5699E+00  1.6735E+00  1.0980E+00  9.9504E-01
             1.0943E+00
 PARAMETER:  7.3694E-02 -8.8405E-01  6.9930E-01  4.3813E-01  9.3928E-02  3.4286E-01  1.3725E+00  6.1492E-01  1.9352E-01  9.5028E-02
             1.9013E-01
 GRADIENT:  -3.9953E+00  4.8819E+00  1.1078E+00  1.1989E+01 -6.2137E+00  2.6947E+00  1.8032E+00  1.0243E+00  5.1515E-01  4.1917E-01
             2.2134E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1510.23667835457        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1104
 NPARAMETR:  9.7461E-01  3.0582E-01  1.9590E+00  1.4452E+00  1.0057E+00  1.2726E+00  3.8810E+00  1.7751E+00  1.0718E+00  1.0134E+00
             1.0949E+00
 PARAMETER:  7.4283E-02 -1.0848E+00  7.7241E-01  4.6824E-01  1.0568E-01  3.4106E-01  1.4561E+00  6.7387E-01  1.6932E-01  1.1334E-01
             1.9066E-01
 GRADIENT:  -2.3198E+00  3.1767E+00  2.0380E+00  8.9767E+00 -6.5453E+00  2.1556E+00  3.5676E-01  7.7701E-01 -1.9626E+00  6.6702E-01
            -3.5219E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1510.29754500400        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1281
 NPARAMETR:  9.7596E-01  2.4278E-01  2.0522E+00  1.4797E+00  1.0129E+00  1.2660E+00  4.2519E+00  1.8333E+00  1.0573E+00  1.0234E+00
             1.0952E+00
 PARAMETER:  7.5666E-02 -1.3156E+00  8.1892E-01  4.9183E-01  1.1285E-01  3.3587E-01  1.5474E+00  7.0610E-01  1.5570E-01  1.2309E-01
             1.9093E-01
 GRADIENT:   3.8307E-01  6.7896E-01  1.3725E+00  1.8381E+00 -2.1949E+00  1.7083E-01 -6.0698E-01  9.4775E-02 -1.7999E+00  1.5857E-01
            -2.5578E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1510.43050191403        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1467             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7673E-01  2.4319E-01  2.0355E+00  1.4658E+00  1.0100E+00  1.2718E+00  4.3459E+00  1.8129E+00  1.0615E+00  1.0227E+00
             1.0949E+00
 PARAMETER:  7.6458E-02 -1.3139E+00  8.1074E-01  4.8238E-01  1.0995E-01  3.4047E-01  1.5692E+00  6.9490E-01  1.5969E-01  1.2248E-01
             1.9065E-01
 GRADIENT:   3.9621E+02  6.7407E+01  1.0322E+01  8.0469E+02  3.3789E+00  2.0029E+02  1.8003E+02  2.8456E+00  3.7074E+01  1.0239E+00
             1.4206E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1510.44169897006        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1648             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7669E-01  2.4225E-01  2.0206E+00  1.4653E+00  1.0071E+00  1.2718E+00  4.3405E+00  1.7960E+00  1.0641E+00  1.0208E+00
             1.0947E+00
 PARAMETER:  7.6412E-02 -1.3178E+00  8.0338E-01  4.8207E-01  1.0712E-01  3.4044E-01  1.5680E+00  6.8555E-01  1.6214E-01  1.2058E-01
             1.9045E-01
 GRADIENT:   3.9618E+02  6.7031E+01  1.0402E+01  8.0386E+02  3.3542E+00  2.0021E+02  1.7868E+02  2.5126E+00  3.8495E+01  9.0075E-01
             1.3708E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1510.45186647293        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1835
 NPARAMETR:  9.7668E-01  2.4240E-01  2.0067E+00  1.4648E+00  1.0067E+00  1.2718E+00  4.3393E+00  1.7946E+00  1.0642E+00  1.0198E+00
             1.0948E+00
 PARAMETER:  7.6405E-02 -1.3172E+00  7.9648E-01  4.8175E-01  1.0667E-01  3.4044E-01  1.5677E+00  6.8481E-01  1.6223E-01  1.1963E-01
             1.9059E-01
 GRADIENT:   1.6516E+00  4.1166E-01  7.9004E-01 -1.3637E+01 -9.3739E-02  2.0356E+00  3.2973E+00  1.9668E-01  1.4891E-01  1.8872E-01
             4.0526E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1510.45832123785        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2023             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7667E-01  2.4299E-01  1.9919E+00  1.4638E+00  1.0052E+00  1.2718E+00  4.3331E+00  1.7851E+00  1.0645E+00  1.0163E+00
             1.0946E+00
 PARAMETER:  7.6391E-02 -1.3147E+00  7.8911E-01  4.8107E-01  1.0518E-01  3.4042E-01  1.5663E+00  6.7945E-01  1.6248E-01  1.1614E-01
             1.9043E-01
 GRADIENT:   3.9589E+02  6.7054E+01  8.4443E+00  8.0112E+02  6.6777E+00  1.9996E+02  1.7956E+02  2.9786E+00  3.8650E+01  5.0912E-01
             1.4997E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1510.46107397636        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2210
 NPARAMETR:  9.7666E-01  2.4352E-01  1.9825E+00  1.4631E+00  1.0042E+00  1.2718E+00  4.3272E+00  1.7775E+00  1.0647E+00  1.0139E+00
             1.0945E+00
 PARAMETER:  7.6386E-02 -1.3126E+00  7.8435E-01  4.8055E-01  1.0418E-01  3.4042E-01  1.5649E+00  6.7523E-01  1.6270E-01  1.1384E-01
             1.9031E-01
 GRADIENT:   1.6430E+00  2.1998E-01 -1.8688E-01 -1.3517E+01  2.1433E+00  2.0239E+00  3.2424E+00  1.9323E-01  2.4860E-01 -3.5561E-01
            -9.4605E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1510.46538978616        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2398             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7666E-01  2.4494E-01  1.9766E+00  1.4621E+00  1.0008E+00  1.2718E+00  4.3214E+00  1.7693E+00  1.0647E+00  1.0175E+00
             1.0947E+00
 PARAMETER:  7.6382E-02 -1.3067E+00  7.8138E-01  4.7986E-01  1.0078E-01  3.4042E-01  1.5636E+00  6.7060E-01  1.6270E-01  1.1733E-01
             1.9044E-01
 GRADIENT:   3.9574E+02  6.7265E+01  9.1597E+00  7.9768E+02  3.9716E+00  1.9984E+02  1.7880E+02  2.8853E+00  3.8560E+01  1.0302E+00
             1.6500E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1510.46852606312        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2585
 NPARAMETR:  9.7666E-01  2.4543E-01  1.9709E+00  1.4615E+00  1.0006E+00  1.2718E+00  4.3164E+00  1.7650E+00  1.0649E+00  1.0157E+00
             1.0946E+00
 PARAMETER:  7.6382E-02 -1.3048E+00  7.7850E-01  4.7948E-01  1.0057E-01  3.4042E-01  1.5624E+00  6.6816E-01  1.6287E-01  1.1558E-01
             1.9036E-01
 GRADIENT:   1.6380E+00  3.8525E-01  5.0748E-01 -1.3562E+01 -2.6949E-01  2.0287E+00  3.3366E+00  1.6921E-01  1.4641E-01  1.6659E-01
             6.2213E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1510.47162427413        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2773             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7666E-01  2.4629E-01  1.9625E+00  1.4606E+00  9.9985E-01  1.2718E+00  4.3079E+00  1.7585E+00  1.0652E+00  1.0133E+00
             1.0944E+00
 PARAMETER:  7.6385E-02 -1.3013E+00  7.7421E-01  4.7887E-01  9.9848E-02  3.4042E-01  1.5604E+00  6.6445E-01  1.6314E-01  1.1324E-01
             1.9022E-01
 GRADIENT:   3.9577E+02  6.7327E+01  8.4822E+00  7.9523E+02  5.8031E+00  1.9981E+02  1.7859E+02  2.7707E+00  3.8711E+01  5.4592E-01
             1.4986E+00

0ITERATION NO.:   82    OBJECTIVE VALUE:  -1510.47162427413        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:     2834
 NPARAMETR:  9.7666E-01  2.4629E-01  1.9625E+00  1.4606E+00  9.9985E-01  1.2718E+00  4.3079E+00  1.7585E+00  1.0652E+00  1.0133E+00
             1.0944E+00
 PARAMETER:  7.6385E-02 -1.3013E+00  7.7421E-01  4.7887E-01  9.9848E-02  3.4042E-01  1.5604E+00  6.6445E-01  1.6314E-01  1.1324E-01
             1.9022E-01
 GRADIENT:  -1.6352E-03 -1.0468E-01  2.1843E-01  4.0964E-01  3.7231E-01 -1.0286E-03  1.0669E-01  1.0348E-01 -9.4691E-03 -4.5411E-02
            -2.1659E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2834
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.9587E-04  4.0446E-02 -4.5791E-02 -3.5154E-02 -3.4423E-02
 SE:             2.9902E-02  1.7859E-02  1.7412E-02  2.4488E-02  1.9223E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8143E-01  2.3533E-02  8.5407E-03  1.5112E-01  7.3334E-02

 ETASHRINKSD(%)  1.0000E-10  4.0169E+01  4.1669E+01  1.7964E+01  3.5602E+01
 ETASHRINKVR(%)  1.0000E-10  6.4202E+01  6.5975E+01  3.2700E+01  5.8529E+01
 EBVSHRINKSD(%)  3.2444E-01  4.8237E+01  4.6103E+01  1.3157E+01  3.0188E+01
 EBVSHRINKVR(%)  6.4782E-01  7.3205E+01  7.0951E+01  2.4582E+01  5.1263E+01
 RELATIVEINF(%)  9.9190E+01  8.1073E+00  8.8376E+00  2.4475E+01  1.4793E+01
 EPSSHRINKSD(%)  4.5463E+01
 EPSSHRINKVR(%)  7.0257E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1510.4716242741315     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -775.32079771039332     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    47.36
 Elapsed covariance  time in seconds:     8.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1510.472       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.77E-01  2.46E-01  1.96E+00  1.46E+00  1.00E+00  1.27E+00  4.31E+00  1.76E+00  1.07E+00  1.01E+00  1.09E+00
 


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
 
         3.77E-02  2.27E-01  6.36E-01  1.62E-01  1.15E-01  9.59E-02  1.83E+00  5.99E-01  1.08E-01  1.77E-01  1.00E-01
 


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
+        1.42E-03
 
 TH 2
+        3.52E-04  5.17E-02
 
 TH 3
+        7.15E-03 -8.09E-02  4.05E-01
 
 TH 4
+        8.30E-04 -3.43E-02  6.63E-02  2.62E-02
 
 TH 5
+        1.82E-03 -6.02E-03  6.30E-02  6.56E-03  1.32E-02
 
 TH 6
+       -1.53E-04 -5.63E-03  8.43E-03  4.21E-03  2.78E-04  9.19E-03
 
 TH 7
+        3.28E-03 -4.05E-01  6.72E-01  2.78E-01  5.59E-02  5.99E-02  3.35E+00
 
 TH 8
+        4.34E-03 -6.30E-02  3.34E-01  5.00E-02  5.17E-02  6.95E-03  5.23E-01  3.59E-01
 
 TH 9
+        5.16E-04  1.50E-02 -2.69E-02 -1.03E-02 -2.45E-03 -2.53E-04 -1.00E-01 -2.03E-02  1.16E-02
 
 TH10
+        1.12E-03 -1.23E-02  4.64E-02  1.09E-02  9.82E-03  1.93E-03  1.08E-01  2.42E-02 -4.87E-03  3.14E-02
 
 TH11
+       -3.04E-06  6.29E-04  1.65E-02  1.19E-04  3.48E-03 -2.02E-04 -5.37E-03  9.52E-03  6.81E-05  8.33E-04  1.00E-02
 
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
+        3.77E-02
 
 TH 2
+        4.11E-02  2.27E-01
 
 TH 3
+        2.98E-01 -5.59E-01  6.36E-01
 
 TH 4
+        1.36E-01 -9.33E-01  6.45E-01  1.62E-01
 
 TH 5
+        4.21E-01 -2.31E-01  8.63E-01  3.53E-01  1.15E-01
 
 TH 6
+       -4.23E-02 -2.58E-01  1.38E-01  2.71E-01  2.52E-02  9.59E-02
 
 TH 7
+        4.76E-02 -9.74E-01  5.77E-01  9.38E-01  2.66E-01  3.41E-01  1.83E+00
 
 TH 8
+        1.92E-01 -4.63E-01  8.77E-01  5.16E-01  7.53E-01  1.21E-01  4.77E-01  5.99E-01
 
 TH 9
+        1.27E-01  6.12E-01 -3.92E-01 -5.90E-01 -1.98E-01 -2.44E-02 -5.08E-01 -3.15E-01  1.08E-01
 
 TH10
+        1.69E-01 -3.06E-01  4.12E-01  3.82E-01  4.83E-01  1.14E-01  3.34E-01  2.28E-01 -2.54E-01  1.77E-01
 
 TH11
+       -8.08E-04  2.76E-02  2.59E-01  7.37E-03  3.03E-01 -2.11E-02 -2.93E-02  1.59E-01  6.31E-03  4.70E-02  1.00E-01
 
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
+        1.24E+03
 
 TH 2
+       -3.71E+02  8.58E+02
 
 TH 3
+       -6.31E+00  1.88E+01  3.46E+01
 
 TH 4
+       -2.37E+02  9.43E+01 -2.62E+01  4.86E+02
 
 TH 5
+       -2.32E+02 -8.61E+01 -9.69E+01  6.50E+01  6.31E+02
 
 TH 6
+        7.62E+01 -1.32E+02 -3.72E+00  7.87E-01  2.38E+01  1.50E+02
 
 TH 7
+       -2.91E+01  9.13E+01  1.08E+00 -2.23E+01 -2.19E+00 -1.81E+01  1.28E+01
 
 TH 8
+        2.62E+01 -4.12E+00 -1.26E+01  5.28E+00 -1.19E+01 -3.66E-02 -8.37E-01  1.47E+01
 
 TH 9
+       -3.54E+01 -2.12E+02  9.04E-01  8.08E+01  1.38E+01  1.07E+01 -2.70E+01 -4.32E-01  2.09E+02
 
 TH10
+        4.26E+01 -2.04E+01  2.45E+00 -2.04E+01 -8.36E+01 -2.02E+00 -2.80E+00  9.09E+00  9.38E+00  5.59E+01
 
 TH11
+        7.52E+01 -2.79E+00 -1.20E+01 -7.05E+00 -3.75E+01 -3.82E-01  1.21E+00  9.80E+00 -9.98E+00  1.16E+01  1.23E+02
 
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
 #CPUT: Total CPU Time in Seconds,       55.443
Stop Time:
Wed Sep 29 19:37:42 CDT 2021
