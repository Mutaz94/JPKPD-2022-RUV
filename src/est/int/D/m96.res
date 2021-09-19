Sat Sep 18 07:48:55 CDT 2021
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
$DATA ../../../../data/int/D/dat96.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   54829.5101187313        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.0033E+02  9.0878E+02  4.4408E+01  8.7204E+02 -9.0295E+01 -4.0407E+03 -2.2674E+03 -1.0345E+02 -2.6683E+03 -1.2640E+03
            -1.0615E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -486.306449731701        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1608E+00  1.9168E+00  8.8414E-01  1.7154E+00  9.6740E-01  5.0609E+00  3.6503E+00  9.9236E-01  1.8112E+00  1.2367E+00
             1.2806E+01
 PARAMETER:  2.4915E-01  7.5063E-01 -2.3145E-02  6.3967E-01  6.6858E-02  1.7215E+00  1.3948E+00  9.2329E-02  6.9399E-01  3.1247E-01
             2.6499E+00
 GRADIENT:  -1.4178E+01  3.6959E+01 -2.2145E+01  1.9128E+02 -2.8689E+01  8.2350E+01 -2.0483E+02  3.9830E+00 -4.3145E+01  2.4214E+01
            -4.3586E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -597.180372483444        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.2965E+00  5.6788E+00  1.8533E+01  4.0134E+00  2.9188E+00  2.9609E+00  5.3988E+00  6.9112E-01  7.5217E+00  5.5077E-01
             1.3250E+01
 PARAMETER:  3.5971E-01  1.8367E+00  3.0196E+00  1.4896E+00  1.1712E+00  1.1855E+00  1.7862E+00 -2.6944E-01  2.1178E+00 -4.9644E-01
             2.6840E+00
 GRADIENT:   3.6908E+01  5.6596E+01 -9.0419E+00  6.0140E+01  3.4800E+01  6.4928E+01 -4.7394E+01  1.6768E-01 -2.0731E+01  3.8004E+00
             2.7830E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -607.933570888061        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      239
 NPARAMETR:  1.2674E+00  4.7231E+00  1.8539E+01  4.0762E+00  2.7393E+00  2.9961E+00  5.5505E+00  4.9513E-01  7.4652E+00  5.3793E-01
             1.3194E+01
 PARAMETER:  3.3700E-01  1.6525E+00  3.0199E+00  1.5052E+00  1.1077E+00  1.1973E+00  1.8139E+00 -6.0294E-01  2.1102E+00 -5.2002E-01
             2.6797E+00
 GRADIENT:   3.2086E+01  4.3802E+01 -1.0320E+01  5.9747E+01  2.6104E+01  6.8748E+01 -3.3261E+01  1.0398E-01 -1.0070E+01  3.8539E+00
             3.6960E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -706.747539138189        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      316
 NPARAMETR:  7.8686E-01  4.8561E-01  4.9930E+01  1.8529E+00  2.2835E+00  2.8850E+00  6.4887E+00  5.1692E-01  2.3690E+00  9.3157E-01
             1.2339E+01
 PARAMETER: -1.3971E-01 -6.2235E-01  4.0106E+00  7.1674E-01  9.2570E-01  1.1595E+00  1.9701E+00 -5.5986E-01  9.6247E-01  2.9121E-02
             2.6128E+00
 GRADIENT:  -6.6664E+01 -1.1390E+01 -9.9311E-01  1.5353E+01 -3.5084E+00  5.0260E+01 -7.8974E+01  7.9643E-03  7.6796E+00  1.3911E+01
             4.6559E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -749.386950612801        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      388
 NPARAMETR:  1.0047E+00  2.6315E-01  7.2161E+01  2.0236E+00  2.3380E+00  2.0982E+00  9.1812E+00  1.5808E+00  2.2518E+00  4.0184E-01
             1.2301E+01
 PARAMETER:  1.0466E-01 -1.2350E+00  4.3789E+00  8.0488E-01  9.4928E-01  8.4110E-01  2.3172E+00  5.5793E-01  9.1172E-01 -8.1170E-01
             2.6097E+00
 GRADIENT:  -1.0660E+01 -3.3994E+00 -1.4493E+00  9.9130E+00 -2.6682E+00 -2.0522E+00 -1.5097E+00  2.7210E-02  3.9554E+00  2.2185E+00
            -9.0677E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -756.301806766940        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  1.0444E+00  4.3366E-01  7.0235E+01  1.6757E+00  2.3147E+00  2.1524E+00  9.0028E+00  9.4468E-01  1.7211E+00  1.3129E-01
             1.2877E+01
 PARAMETER:  1.4344E-01 -7.3550E-01  4.3518E+00  6.1624E-01  9.3926E-01  8.6658E-01  2.2975E+00  4.3096E-02  6.4296E-01 -1.9304E+00
             2.6555E+00
 GRADIENT:  -1.9042E+00 -2.8345E+00  2.3530E-02  3.3917E-01 -9.8047E+00  8.6492E+00  4.0225E+00  4.9086E-03  5.3557E+00  2.3438E-01
             4.7649E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -756.510184497086        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  1.0471E+00  5.6758E-01  6.5518E+01  1.6059E+00  2.3587E+00  2.1222E+00  8.5969E+00  8.5495E-01  1.6539E+00  1.5027E-01
             1.2799E+01
 PARAMETER:  1.4598E-01 -4.6637E-01  4.2823E+00  5.7367E-01  9.5810E-01  8.5247E-01  2.2514E+00 -5.6711E-02  6.0314E-01 -1.7953E+00
             2.6494E+00
 GRADIENT:   2.2152E-01 -8.0168E-01 -1.4396E-01 -2.8925E-01 -4.1717E+00  4.2299E+00  4.2847E+00  3.8922E-03  4.4731E+00  3.1105E-01
             3.6761E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -758.357269719818        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      713
 NPARAMETR:  1.0391E+00  3.7866E-01  9.8153E+01  1.7186E+00  2.4248E+00  2.0913E+00  1.0350E+01  2.0532E+00  1.6616E+00  1.1668E-01
             1.2607E+01
 PARAMETER:  1.3837E-01 -8.7111E-01  4.6865E+00  6.4149E-01  9.8574E-01  8.3777E-01  2.4370E+00  8.1938E-01  6.0777E-01 -2.0483E+00
             2.6343E+00
 GRADIENT:  -3.7544E-02  1.9320E-01 -2.7694E-01 -6.4000E-01  1.5369E+00 -1.1028E+00  1.9358E+00  1.1782E-02  4.0381E-01  1.8663E-01
             4.9297E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -758.436848789844        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  1.0347E+00  3.5401E-01  1.2872E+02  1.7133E+00  2.4373E+00  2.0933E+00  1.0289E+01  2.2975E+00  1.6416E+00  7.9901E-02
             1.2566E+01
 PARAMETER:  1.3411E-01 -9.3842E-01  4.9576E+00  6.3842E-01  9.9090E-01  8.3874E-01  2.4311E+00  9.3180E-01  5.9564E-01 -2.4270E+00
             2.6310E+00
 GRADIENT:  -1.1919E+00 -9.2469E-01 -1.8365E-01 -4.0562E-01  1.5316E+00 -9.3099E-01 -1.4037E+00  8.9033E-03 -1.0097E+00  8.7300E-02
             1.0128E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -758.636306008176        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  1.0380E+00  3.8490E-01  1.6919E+03  1.7133E+00  2.4696E+00  2.1021E+00  1.0228E+01  5.4763E+00  1.6602E+00  1.0000E-02
             1.2560E+01
 PARAMETER:  1.3732E-01 -8.5477E-01  7.5336E+00  6.3844E-01  1.0040E+00  8.4294E-01  2.4251E+00  1.8004E+00  6.0691E-01 -6.6548E+00
             2.6305E+00
 GRADIENT:   1.3472E-01  5.1334E-03 -6.8928E-03  5.1402E-02  2.9021E-02 -5.3846E-02  6.3094E-02  3.0392E-04  3.1350E-04  0.0000E+00
            -8.8447E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -758.642371092448        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1241
 NPARAMETR:  1.0379E+00  3.8467E-01  3.2135E+04  1.7130E+00  2.4737E+00  2.1022E+00  1.0226E+01  1.4522E+01  1.6596E+00  1.0000E-02
             1.2558E+01
 PARAMETER:  1.3721E-01 -8.5538E-01  1.0478E+01  6.3823E-01  1.0057E+00  8.4300E-01  2.4250E+00  2.7757E+00  6.0655E-01 -1.1476E+01
             2.6303E+00
 GRADIENT:  -9.0441E-01  1.5789E-02 -3.8916E-04  3.8933E-01  1.4518E-01 -1.0096E-01 -2.2092E-03  5.2767E-05 -8.1915E-02  0.0000E+00
             2.6216E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -758.642848255266        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1428
 NPARAMETR:  1.0381E+00  3.8448E-01  2.8410E+04  1.7142E+00  2.4727E+00  2.1018E+00  1.0248E+01  1.4514E+01  1.6603E+00  1.0000E-02
             1.2560E+01
 PARAMETER:  1.3740E-01 -8.5587E-01  1.0354E+01  6.3895E-01  1.0053E+00  8.4281E-01  2.4271E+00  2.7751E+00  6.0701E-01 -1.1476E+01
             2.6305E+00
 GRADIENT:  -3.2913E-01  6.0543E-02 -3.3473E-04  1.9766E-01  5.2121E-03  1.1988E-01  3.3901E-01 -3.5309E-05 -1.0740E-01  0.0000E+00
             5.8916E-02

0ITERATION NO.:   61    OBJECTIVE VALUE:  -758.642848255266        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1457
 NPARAMETR:  1.0380E+00  3.8415E-01  2.8441E+04  1.7141E+00  2.4728E+00  2.1022E+00  1.0250E+01  1.4531E+01  1.6603E+00  1.0000E-02
             1.2561E+01
 PARAMETER:  1.3740E-01 -8.5587E-01  1.0354E+01  6.3895E-01  1.0053E+00  8.4281E-01  2.4271E+00  2.7751E+00  6.0701E-01 -1.1476E+01
             2.6305E+00
 GRADIENT:   4.4791E-01  4.1492E-02 -1.0505E-03  4.2624E-02 -1.1640E-02 -4.5507E-02 -4.0988E-02 -9.6063E-03  7.0118E-04  0.0000E+00
            -8.1811E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1457
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5295E-02  6.4316E-02  3.0494E-06 -9.8409E-02 -1.1145E-05
 SE:             2.7836E-02  2.1189E-02  7.9157E-06  1.6459E-02  8.8803E-05
 N:                     100         100         100         100         100

 P VAL.:         5.8269E-01  2.4034E-03  7.0006E-01  2.2513E-09  9.0013E-01

 ETASHRINKSD(%)  6.7453E+00  2.9013E+01  9.9973E+01  4.4861E+01  9.9702E+01
 ETASHRINKVR(%)  1.3036E+01  4.9608E+01  1.0000E+02  6.9597E+01  9.9999E+01
 EBVSHRINKSD(%)  8.5736E+00  2.6893E+01  9.9960E+01  3.7609E+01  9.9603E+01
 EBVSHRINKVR(%)  1.6412E+01  4.6553E+01  1.0000E+02  6.1073E+01  9.9998E+01
 RELATIVEINF(%)  8.3339E+01  3.3965E+01  3.6223E-06  2.4191E+01  3.5144E-04
 EPSSHRINKSD(%)  4.4833E+00
 EPSSHRINKVR(%)  8.7655E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -758.64284825526556     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       895.44651151314520     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    43.94
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -758.643       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  3.84E-01  2.84E+04  1.71E+00  2.47E+00  2.10E+00  1.02E+01  1.45E+01  1.66E+00  1.00E-02  1.26E+01
 


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
+        1.06E+03
 
 TH 2
+       -1.90E+02  1.64E+02
 
 TH 3
+        1.19E-04 -2.71E-05  2.18E-10
 
 TH 4
+       -1.12E+01  2.72E+01  2.75E-05  9.16E+01
 
 TH 5
+        9.58E+00 -8.43E-01  4.06E-06 -6.16E+00  3.03E+01
 
 TH 6
+       -1.32E+01  1.11E+01  5.37E-06 -1.32E+00  3.92E+00  3.55E+01
 
 TH 7
+        1.52E-02  3.97E+00  7.70E-07 -2.36E+00  1.44E-01  3.85E-01  9.60E-01
 
 TH 8
+       -9.84E-01  4.60E-01 -4.67E-07  3.61E-02 -4.18E-02 -1.38E-02  6.80E-03  1.27E-02
 
 TH 9
+        1.72E+01 -4.24E+00  7.88E-06 -2.63E+01  4.36E-01 -1.18E+00  7.38E-01  6.26E-01  3.42E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.68E+00 -1.53E+00  1.12E-06 -5.72E+00  2.98E-01  1.59E+00  6.46E-02  5.13E-03  2.05E+00  0.00E+00  5.11E+00
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       61.460
Stop Time:
Sat Sep 18 07:49:58 CDT 2021
