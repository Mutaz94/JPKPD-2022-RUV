Thu Sep 30 03:07:29 CDT 2021
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
$DATA ../../../../data/spa1/D/dat46.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   15158.0690932402        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.3849E+02  3.0199E+02 -1.8022E+01  1.7373E+02  1.2500E+02 -1.6515E+03 -7.9411E+02 -4.7398E+01 -1.2652E+03 -4.1801E+02
            -2.9910E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -611.057829717050        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0557E+00  1.0080E+00  1.0125E+00  1.7870E+00  1.2745E+00  2.3502E+00  1.4530E+00  9.3691E-01  1.8213E+00  1.1328E+00
             1.3510E+01
 PARAMETER:  1.5418E-01  1.0793E-01  1.1241E-01  6.8054E-01  3.4252E-01  9.5450E-01  4.7361E-01  3.4830E-02  6.9955E-01  2.2474E-01
             2.7034E+00
 GRADIENT:  -7.0198E+01  1.8782E+01 -1.5102E+01  3.6972E+01 -4.4646E+00  6.3173E+01  2.0872E+00  5.1789E+00  4.3544E+00  3.1687E+00
             2.6516E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -638.262739893446        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.1109E+00  9.1643E-01  1.7667E+00  2.0820E+00  4.3919E+00  1.9825E+00  5.2037E+00  6.8187E-01  2.0995E+00  5.0803E+00
             1.2208E+01
 PARAMETER:  2.0520E-01  1.2731E-02  6.6912E-01  8.3333E-01  1.5798E+00  7.8437E-01  1.7494E+00 -2.8292E-01  8.4170E-01  1.7254E+00
             2.6021E+00
 GRADIENT:  -4.5056E+01  2.7038E+01  6.5847E+00  5.4020E+01 -4.8583E+00 -4.6844E-01  1.7168E+01 -9.1525E-01  2.6007E+01 -2.8081E-01
             2.2645E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -663.398075669073        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      272             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1107E+00  9.0582E-01  1.8037E+00  2.0767E+00  4.8652E+00  1.9709E+00  4.7961E+00  6.6408E-01  2.0995E+00  8.6504E+00
             1.0102E+01
 PARAMETER:  2.0496E-01  1.0862E-03  6.8984E-01  8.3080E-01  1.6821E+00  7.7847E-01  1.6678E+00 -3.0936E-01  8.4170E-01  2.2576E+00
             2.4128E+00
 GRADIENT:  -7.3811E+00  3.1558E+01  4.3402E+00  8.6476E+01 -1.0282E+01 -1.2404E+01  1.1297E+01 -1.3926E+00  1.1216E+01  2.3308E+01
             8.9889E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -670.411502961588        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      390             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1097E+00  9.0625E-01  1.7956E+00  2.0584E+00  5.0396E+00  2.4028E+00  4.4553E+00  6.7870E-01  2.0872E+00  8.1857E+00
             9.1469E+00
 PARAMETER:  2.0407E-01  1.5592E-03  6.8534E-01  8.2193E-01  1.7173E+00  9.7663E-01  1.5941E+00 -2.8758E-01  8.3584E-01  2.2024E+00
             2.3134E+00
 GRADIENT:   1.6874E+01  3.6272E+01  1.2774E+00  1.0514E+02 -2.5875E+01  5.7659E+01  5.9194E+00 -1.2032E+00 -1.0268E+01  4.5411E+01
             8.1962E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -676.072673080202        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      463
 NPARAMETR:  1.1086E+00  9.0634E-01  1.7863E+00  1.9694E+00  5.3169E+00  1.9543E+00  4.3065E+00  7.1364E-01  2.0857E+00  7.3910E+00
             8.8310E+00
 PARAMETER:  2.0308E-01  1.6614E-03  6.8017E-01  7.7771E-01  1.7709E+00  7.7005E-01  1.5601E+00 -2.3738E-01  8.3511E-01  2.1003E+00
             2.2783E+00
 GRADIENT:   1.2960E+01  3.6220E+01  6.8934E+00  9.5668E+01 -3.3529E+00 -1.8887E+01  8.9251E+00 -1.8860E+00 -3.8638E+00  4.2357E+00
            -7.6778E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -690.216846203728        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      536
 NPARAMETR:  1.0998E+00  9.0585E-01  1.6926E+00  1.2652E+00  8.0105E+00  2.0164E+00  2.9664E+00  1.3563E+00  2.0992E+00  3.5079E+00
             8.2297E+00
 PARAMETER:  1.9514E-01  1.1151E-03  6.2625E-01  3.3524E-01  2.1808E+00  8.0132E-01  1.1873E+00  4.0473E-01  8.4155E-01  1.3550E+00
             2.2077E+00
 GRADIENT:   2.9290E+01  1.0448E+01  1.9950E+01 -4.0552E+01 -2.7511E+00  3.3834E+00  1.9258E+01 -7.6343E+00 -1.7245E+01 -3.1978E-01
            -2.7863E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -701.366483693626        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      609
 NPARAMETR:  1.0966E+00  9.0561E-01  1.6158E+00  1.3423E+00  8.5282E+00  1.8900E+00  2.3301E+00  3.1013E+00  2.1369E+00  3.5901E+00
             8.5813E+00
 PARAMETER:  1.9220E-01  8.4945E-04  5.7985E-01  3.9438E-01  2.2434E+00  7.3660E-01  9.4590E-01  1.2318E+00  8.5937E-01  1.3782E+00
             2.2496E+00
 GRADIENT:   2.5635E+01  2.6186E+01  2.6100E+00 -7.6214E+00 -3.9362E+00 -7.6409E+00  4.5048E+00  4.4883E+00  9.8141E+00 -7.4259E-01
            -8.9873E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -703.483297416989        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      684
 NPARAMETR:  1.0865E+00  9.0403E-01  1.5412E+00  1.3347E+00  1.1437E+01  1.9446E+00  1.7140E+00  2.3042E+00  2.0679E+00  3.4224E+00
             8.5602E+00
 PARAMETER:  1.8296E-01 -8.9358E-04  5.3254E-01  3.8867E-01  2.5368E+00  7.6504E-01  6.3886E-01  9.3473E-01  8.2652E-01  1.3303E+00
             2.2471E+00
 GRADIENT:   1.9849E+01  2.6493E+01  9.3179E+00 -6.4503E+00 -2.2861E+00 -8.3416E-01  3.0190E+00 -2.2324E+00 -1.0897E+01 -3.2719E-01
            -3.8183E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -705.589821662412        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      755
 NPARAMETR:  1.0345E+00  8.9567E-01  1.1252E+00  1.2958E+00  5.3279E+01  1.9276E+00  2.7894E-01  1.8965E+00  1.9169E+00  2.4815E+00
             8.7855E+00
 PARAMETER:  1.3388E-01 -1.0185E-02  2.1800E-01  3.5910E-01  4.0755E+00  7.5626E-01 -1.1768E+00  7.4000E-01  7.5073E-01  1.0089E+00
             2.2731E+00
 GRADIENT:  -6.8747E+00  4.6407E+01  3.4070E+00  7.6327E-01 -5.8293E-01  7.4604E-01 -9.3034E-02  4.3704E+00 -3.8711E+01 -1.0049E-02
             3.7238E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -713.924245109668        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      830
 NPARAMETR:  9.9841E-01  8.8523E-01  8.5608E-01  1.2387E+00  1.7774E+02  1.9666E+00  8.1884E-02  3.4836E-01  2.3197E+00  2.0513E+00
             8.4897E+00
 PARAMETER:  9.8411E-02 -2.1906E-02 -5.5391E-02  3.1405E-01  5.2803E+00  7.7633E-01 -2.4024E+00 -9.5453E-01  9.4146E-01  8.1848E-01
             2.2389E+00
 GRADIENT:  -1.7259E+01  2.1853E+01  1.8968E+00 -1.1243E-01 -1.5898E-01 -8.2478E-01  4.4925E-02  3.8119E-01  1.7930E+00  5.0692E-05
             1.7443E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -713.994150644843        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  9.9521E-01  8.8427E-01  8.3743E-01  1.2358E+00  2.0033E+02  1.9713E+00  7.2848E-02  2.8623E-01  2.3216E+00  2.0234E+00
             8.4747E+00
 PARAMETER:  9.5195E-02 -2.2996E-02 -7.7417E-02  3.1174E-01  5.4000E+00  7.7870E-01 -2.5194E+00 -1.1510E+00  9.4224E-01  8.0477E-01
             2.2371E+00
 GRADIENT:  -1.7987E+01  2.2689E+01  4.3473E-01  1.1158E-01 -1.4749E-01  3.6990E-01  3.6770E-02  3.0204E-01  2.4617E+00  5.9701E-05
            -1.2681E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -715.213594063292        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      980
 NPARAMETR:  9.8918E-01  8.4028E-01  7.4573E-01  1.2322E+00  2.4359E+03  1.9929E+00  4.5804E-02  1.0000E-02  2.2547E+00  4.4569E+00
             8.4765E+00
 PARAMETER:  8.9126E-02 -7.4022E-02 -1.9340E-01  3.0880E-01  7.8981E+00  7.8958E-01 -2.9834E+00 -2.0188E+01  9.1304E-01  1.5944E+00
             2.2373E+00
 GRADIENT:  -1.6024E+01  2.5233E+01 -5.5080E+00 -3.2311E+00 -1.4496E-02  7.7587E+00  1.6077E-02  0.0000E+00  4.4135E+00  5.0998E-06
            -8.6127E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -726.287954686362        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  1.0493E+00  4.2614E-01  3.7712E-01  1.3403E+00  3.2716E+16  2.0086E+00  1.1228E-02  1.0000E-02  1.4716E+00  6.4381E+05
             8.8036E+00
 PARAMETER:  1.4812E-01 -7.5299E-01 -8.7519E-01  3.9286E-01  3.8127E+01  7.9745E-01 -4.3894E+00 -2.8241E+02  4.8637E-01  1.3475E+01
             2.2752E+00
 GRADIENT:   5.7248E+01  6.8625E+01 -5.1874E+01  2.5925E+01  0.0000E+00  2.9431E+01  1.1527E-03  0.0000E+00 -1.2274E+01  0.0000E+00
            -1.2707E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -760.412728993941        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  7.5879E-01  1.1680E-01  1.5736E-01  1.0721E+00  3.7183E+16  1.6186E+00  1.0000E-02  1.0000E-02  1.1036E+00  6.1509E+05
             7.7038E+00
 PARAMETER: -1.7603E-01 -2.0473E+00 -1.7492E+00  1.6964E-01  3.8255E+01  5.8158E-01 -5.1226E+00 -2.8241E+02  1.9857E-01  1.3430E+01
             2.1417E+00
 GRADIENT:   1.4868E+01  3.0302E+01 -4.2305E+01  1.5928E+02  0.0000E+00 -2.6194E+01  0.0000E+00  0.0000E+00 -2.4808E+01  0.0000E+00
            -1.3224E+02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -771.754514874230        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:     1294
 NPARAMETR:  6.7280E-01  6.8232E-02  9.3112E-02  7.9734E-01  3.7409E+16  1.4780E+00  1.0000E-02  1.0000E-02  9.9167E-01  6.1376E+05
             7.4852E+00
 PARAMETER: -2.9631E-01 -2.5848E+00 -2.2740E+00 -1.2648E-01  3.8261E+01  4.9072E-01 -5.1575E+00 -2.8241E+02  9.1634E-02  1.3427E+01
             2.1129E+00
 GRADIENT:   2.7248E+01  2.2823E+01 -8.9183E+01  1.5866E+02  0.0000E+00 -5.7138E+01  0.0000E+00  0.0000E+00 -2.4635E+01  0.0000E+00
            -1.8271E+02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -807.919583791608        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1470
 NPARAMETR:  4.8118E-01  1.8880E-02  3.5932E-02  3.9246E-01  3.7002E+16  1.5213E+00  1.0000E-02  1.0000E-02  9.1073E-01  6.1616E+05
             8.3712E+00
 PARAMETER: -6.3151E-01 -3.8697E+00 -3.2261E+00 -8.3533E-01  3.8250E+01  5.1959E-01 -5.0947E+00 -2.8241E+02  6.4867E-03  1.3431E+01
             2.2248E+00
 GRADIENT:  -1.7118E+00  9.3659E-01 -1.4715E+01  2.1904E+01  0.0000E+00 -9.4578E-01  0.0000E+00  0.0000E+00 -7.0894E-01  0.0000E+00
            -3.4279E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -808.729632512341        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1652
 NPARAMETR:  4.4205E-01  1.1749E-02  2.8700E-02  3.2658E-01  1.7775E+16  1.4981E+00  1.0000E-02  1.0000E-02  9.0639E-01  5.0442E+05
             8.3750E+00
 PARAMETER: -7.1408E-01 -4.3266E+00 -3.4569E+00 -1.0130E+00  3.7483E+01  5.0486E-01 -5.0833E+00 -2.8241E+02  8.8491E-04  1.3365E+01
             2.2229E+00
 GRADIENT:   2.2173E+00  1.7586E-02 -7.1485E+00  5.5726E+00 -2.3302E-05  2.4903E-01  0.0000E+00  0.0000E+00 -1.1211E-01  1.1457E-04
            -2.4652E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1652
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0467E-03  9.7564E-07  1.0993E-04 -2.4027E-02 -1.5236E-13
 SE:             2.8970E-02  1.2684E-06  2.3355E-04  2.4618E-02  1.3556E-13
 N:                     100         100         100         100         100

 P VAL.:         9.7118E-01  4.4179E-01  6.3788E-01  3.2907E-01  2.6104E-01

 ETASHRINKSD(%)  2.9481E+00  9.9996E+01  9.9218E+01  1.7527E+01  1.0000E+02
 ETASHRINKVR(%)  5.8093E+00  1.0000E+02  9.9994E+01  3.1981E+01  1.0000E+02
 EBVSHRINKSD(%)  2.5656E+00  9.9994E+01  9.9278E+01  1.7981E+01  1.0000E+02
 EBVSHRINKVR(%)  5.0654E+00  1.0000E+02  9.9995E+01  3.2729E+01  1.0000E+02
 RELATIVEINF(%)  5.1291E+00  3.8685E-08  5.1095E-05  7.3256E-01  0.0000E+00
 EPSSHRINKSD(%)  1.2083E+01
 EPSSHRINKVR(%)  2.2707E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -808.72963251234069     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       110.20890069233201     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.93
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    18.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -808.730       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.43E-01  1.20E-02  2.85E-02  3.29E-01  1.72E+16  1.50E+00  1.00E-02  1.00E-02  9.06E-01  5.77E+05  8.36E+00
 


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
+        2.05E+02
 
 TH 2
+       -4.90E+00  1.17E-01
 
 TH 3
+       -1.08E+04  2.58E+02  5.66E+05
 
 TH 4
+        1.09E+03 -2.61E+01 -5.73E+04  5.81E+03
 
 TH 5
+       -1.51E-20  3.61E-22  7.94E-19 -8.05E-20  1.11E-42
 
 TH 6
+        2.69E+00 -6.43E-02 -1.41E+02  1.43E+01 -1.98E-22  3.53E-02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.73E+01  4.13E-01  9.07E+02 -9.19E+01  1.27E-21 -2.27E-01  0.00E+00  0.00E+00  1.45E+00
 
 TH10
+        4.60E-10 -1.10E-11 -2.42E-08  2.45E-09 -3.39E-32  6.04E-12  0.00E+00  0.00E+00 -3.88E-11  1.03E-21
 
 TH11
+       -4.94E+00  1.18E-01  2.59E+02 -2.63E+01  3.64E-22 -6.48E-02  0.00E+00  0.00E+00  4.16E-01 -1.11E-11  1.19E-01
 
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
+        2.37E+03
 
 TH 2
+       -6.76E+02  3.38E+03
 
 TH 3
+       -1.32E+04  3.22E+02  6.87E+05
 
 TH 4
+       -4.78E+01  2.56E+02 -6.96E+04  8.07E+03
 
 TH 5
+        1.06E-19  1.22E-18  9.36E-19 -4.12E-19  3.66E-37
 
 TH 6
+        7.31E+00  2.37E+01 -1.75E+02 -2.36E+01  1.89E-19  7.75E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.03E+00 -2.31E+01  1.10E+03 -1.34E+02  1.26E-18 -1.39E+00  0.00E+00  0.00E+00  1.14E+02
 
 TH10
+       -1.03E-08 -2.74E-07 -2.97E-08 -2.10E-09  1.54E-26  5.67E-09  0.00E+00  0.00E+00 -2.95E-08  3.92E-08
 
 TH11
+       -2.50E+01  6.89E+00  3.15E+02 -2.34E+01 -3.54E-21  2.49E+00  0.00E+00  0.00E+00  4.02E+00  5.40E-10  7.13E+00
 
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
+        2.39E+03
 
 TH 2
+       -5.65E+02  1.67E+02
 
 TH 3
+       -1.85E+04  8.90E+02  8.35E+05
 
 TH 4
+        3.37E+02  2.76E+02 -7.87E+04  8.51E+03
 
 TH 5
+       -6.55E-20  1.56E-20  1.94E-19  2.42E-20  3.31E-41
 
 TH 6
+        2.54E+02 -1.83E+01 -3.61E+03  1.62E+02 -9.60E-21  1.24E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.43E+02 -5.04E+01  7.01E+02 -1.94E+02  4.51E-21  1.99E+01  0.00E+00  0.00E+00  1.28E+02
 
 TH10
+       -5.46E-10  1.55E-10 -1.08E-09  4.82E-10  2.54E-31 -6.08E-11  0.00E+00  0.00E+00 -1.07E-11  2.38E-21
 
 TH11
+       -1.82E+02  5.90E+01  1.15E+03 -2.07E+01  2.54E-21  9.62E+00  0.00E+00  0.00E+00 -1.77E+01  1.11E-10  1.46E+02
 
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
 #CPUT: Total CPU Time in Seconds,       48.090
Stop Time:
Thu Sep 30 03:08:18 CDT 2021
