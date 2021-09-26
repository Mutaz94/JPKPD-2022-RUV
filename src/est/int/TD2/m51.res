Sat Sep 25 04:54:09 CDT 2021
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
$DATA ../../../../data/int/TD2/dat51.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m51.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3414.48823822085        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4616E+02  9.1293E+01  5.5204E+01 -1.1356E+01  6.6112E+01  2.3710E+01 -3.1345E+01 -5.1508E+02 -1.4317E+02 -5.2689E+01
            -7.6129E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3675.97705768280        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.2044E-01  9.0834E-01  1.1578E+00  1.0772E+00  1.0108E+00  9.2692E-01  1.1527E+00  2.6427E+00  9.0583E-01  1.0886E+00
             1.0080E+00
 PARAMETER:  1.7091E-02  3.8620E-03  2.4650E-01  1.7440E-01  1.1072E-01  2.4108E-02  2.4214E-01  1.0718E+00  1.0940E-03  1.8491E-01
             1.0795E-01
 GRADIENT:  -5.1318E+01  4.0961E+01 -2.0898E+01  7.5062E+01  1.9433E+01 -3.0284E+00 -4.9360E+00 -5.4826E+01 -8.0293E+00 -9.7615E+00
             3.7855E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3678.26473880586        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.2411E-01  9.0362E-01  1.2623E+00  1.0830E+00  9.9294E-01  9.1738E-01  1.1887E+00  2.6712E+00  9.1253E-01  1.1143E+00
             1.0062E+00
 PARAMETER:  2.1072E-02 -1.3484E-03  3.3294E-01  1.7975E-01  9.2917E-02  1.3766E-02  2.7289E-01  1.0825E+00  8.4675E-03  2.0819E-01
             1.0615E-01
 GRADIENT:  -4.2516E+01  4.6912E+01  9.5893E-01  7.3135E+01 -1.4755E+01 -6.6085E+00  1.8589E-01 -5.2797E+01 -3.1753E+00 -7.0510E-01
             3.6286E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3687.16106231365        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.3791E-01  8.6657E-01  1.3118E+00  1.0687E+00  1.0067E+00  9.2636E-01  1.1507E+00  3.1220E+00  8.8797E-01  1.1187E+00
             9.8903E-01
 PARAMETER:  3.5893E-02 -4.3218E-02  3.7136E-01  1.6645E-01  1.0663E-01  2.3508E-02  2.4039E-01  1.2385E+00 -1.8823E-02  2.1219E-01
             8.8964E-02
 GRADIENT:  -8.4510E-01  7.0089E+00  2.6205E-02  9.1894E+00  1.4413E+00 -7.4550E-01 -7.4617E+00 -7.8042E+00 -3.8385E+00  6.1317E-01
             1.0350E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3687.23645323982        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      399
 NPARAMETR:  9.3793E-01  8.6648E-01  1.3116E+00  1.0688E+00  1.0067E+00  9.4080E-01  1.1509E+00  3.1234E+00  8.8803E-01  1.1186E+00
             9.8426E-01
 PARAMETER:  3.5918E-02 -4.3316E-02  3.7128E-01  1.6649E-01  1.0663E-01  3.8980E-02  2.4057E-01  1.2389E+00 -1.8753E-02  2.1208E-01
             8.4131E-02
 GRADIENT:  -4.2545E+01  1.3096E+00 -2.8964E+00 -9.6475E+00 -5.5047E+00  4.8506E-01 -9.2325E+00 -1.4451E+01 -4.6725E+00 -8.5919E+00
             7.7870E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3687.35685527187        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      540
 NPARAMETR:  9.3840E-01  8.6596E-01  1.3123E+00  1.0686E+00  1.0067E+00  9.3880E-01  1.1559E+00  3.1330E+00  8.8979E-01  1.1197E+00
             9.8389E-01
 PARAMETER:  3.6420E-02 -4.3916E-02  3.7177E-01  1.6639E-01  1.0667E-01  3.6849E-02  2.4486E-01  1.2420E+00 -1.6773E-02  2.1304E-01
             8.3755E-02
 GRADIENT:   2.1850E+00  6.8609E+00  1.9016E-01  8.8999E+00  2.0375E+00  4.6700E+00 -6.8157E+00 -6.9440E+00 -3.1664E+00  9.2120E-01
             5.5268E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3687.42890119662        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      667             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3840E-01  8.6492E-01  1.3125E+00  1.0685E+00  1.0064E+00  9.3923E-01  1.1613E+00  3.1426E+00  8.9058E-01  1.1186E+00
             9.8355E-01
 PARAMETER:  3.6424E-02 -4.5118E-02  3.7196E-01  1.6630E-01  1.0640E-01  3.7305E-02  2.4951E-01  1.2450E+00 -1.5880E-02  2.1204E-01
             8.3417E-02
 GRADIENT:   2.2761E+00  6.1478E+00  1.4434E-01  7.8520E+00  2.0296E+00  4.8520E+00 -6.0562E+00 -6.0645E+00 -2.7692E+00  8.5006E-01
             2.5121E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3687.93762723647        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      848
 NPARAMETR:  9.5302E-01  8.6496E-01  1.3128E+00  1.0685E+00  1.0064E+00  9.3955E-01  1.1614E+00  3.1406E+00  8.9903E-01  1.1584E+00
             9.8352E-01
 PARAMETER:  5.1886E-02 -4.5068E-02  3.7215E-01  1.6621E-01  1.0635E-01  3.7651E-02  2.4963E-01  1.2444E+00 -6.4344E-03  2.4700E-01
             8.3380E-02
 GRADIENT:  -2.6349E+00 -3.1945E+00 -2.5309E+00 -9.2349E+00 -7.3947E+00  7.0605E-01 -6.9056E+00 -1.3150E+01 -5.2692E-01 -6.9296E-01
             9.1183E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3687.96558247861        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  9.5384E-01  8.6513E-01  1.3138E+00  1.0684E+00  1.0068E+00  9.3773E-01  1.1648E+00  3.1405E+00  9.0060E-01  1.1614E+00
             9.8302E-01
 PARAMETER:  5.2739E-02 -4.4876E-02  3.7296E-01  1.6621E-01  1.0676E-01  3.5703E-02  2.5252E-01  1.2444E+00 -4.6943E-03  2.4959E-01
             8.2875E-02
 GRADIENT:  -4.9314E-01 -3.3342E+00 -2.4046E+00 -8.7947E+00 -7.2627E+00 -5.9737E-02 -6.3823E+00 -1.3228E+01  3.5977E-02 -1.1170E-01
            -1.8167E-05

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3688.04630197099        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     1160
 NPARAMETR:  9.5387E-01  8.6525E-01  1.3146E+00  1.0689E+00  1.0073E+00  9.3784E-01  1.1680E+00  3.1542E+00  9.0066E-01  1.1609E+00
             9.8301E-01
 PARAMETER:  5.2776E-02 -4.4741E-02  3.7354E-01  1.6664E-01  1.0730E-01  3.5819E-02  2.5525E-01  1.2487E+00 -4.6325E-03  2.4920E-01
             8.2867E-02
 GRADIENT:   4.3340E+01  3.2223E+00  5.8608E-01  1.1244E+01  6.7987E-01  4.7536E+00 -4.2027E+00 -5.3063E+00  1.0447E+00  1.0636E+01
             3.3984E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3688.06510887365        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     1308
 NPARAMETR:  9.5393E-01  8.6515E-01  1.3143E+00  1.0691E+00  1.0075E+00  9.3785E-01  1.1680E+00  3.1589E+00  9.0070E-01  1.1605E+00
             9.8280E-01
 PARAMETER:  5.2834E-02 -4.4854E-02  3.7328E-01  1.6681E-01  1.0745E-01  3.5831E-02  2.5525E-01  1.2502E+00 -4.5796E-03  2.4885E-01
             8.2651E-02
 GRADIENT:  -2.5668E-01 -2.8250E+00 -2.6034E+00 -7.6175E+00 -6.6703E+00 -1.4706E-02 -6.0130E+00 -1.1543E+01  2.7656E-01 -3.0668E-01
            -2.1679E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3688.51447410627        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1468
 NPARAMETR:  9.5283E-01  8.6755E-01  1.3346E+00  1.0725E+00  1.0146E+00  9.3700E-01  1.2179E+00  3.2853E+00  8.9798E-01  1.1580E+00
             9.8284E-01
 PARAMETER:  5.1685E-02 -4.2077E-02  3.8866E-01  1.6997E-01  1.1445E-01  3.4925E-02  2.9710E-01  1.2895E+00 -7.6072E-03  2.4670E-01
             8.2686E-02
 GRADIENT:  -3.1902E+00  1.9784E+00 -2.0746E+00  4.3543E-02 -3.8937E+00 -4.1650E-01  1.8302E-01 -9.0644E-01  1.5566E+00 -1.0461E+00
             1.2204E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3688.51614070554        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1655
 NPARAMETR:  9.5286E-01  8.6753E-01  1.3348E+00  1.0725E+00  1.0146E+00  9.3791E-01  1.2179E+00  3.2860E+00  8.9790E-01  1.1580E+00
             9.8224E-01
 PARAMETER:  5.1712E-02 -4.2102E-02  3.8880E-01  1.6998E-01  1.1452E-01  3.5901E-02  2.9709E-01  1.2897E+00 -7.6947E-03  2.4674E-01
             8.2085E-02
 GRADIENT:  -3.1082E+00  1.9873E+00 -2.0470E+00  4.1901E-02 -3.7739E+00 -3.3611E-02  1.5467E-01 -8.7248E-01  1.5363E+00 -1.0618E+00
             3.7758E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3688.51902675212        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1829
 NPARAMETR:  9.5304E-01  8.6725E-01  1.3348E+00  1.0725E+00  1.0146E+00  9.3801E-01  1.2179E+00  3.2860E+00  8.9726E-01  1.1589E+00
             9.8223E-01
 PARAMETER:  5.1897E-02 -4.2431E-02  3.8880E-01  1.6998E-01  1.1451E-01  3.6003E-02  2.9709E-01  1.2897E+00 -8.4074E-03  2.4748E-01
             8.2073E-02
 GRADIENT:  -2.2507E+00  2.0550E+01  1.0527E+05 -2.4079E+05  2.4217E+01  7.7280E-03  6.8875E+04 -1.5887E+04 -9.7046E-01 -4.1364E+01
             5.1172E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1829
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4359E-03 -3.2831E-02 -2.8727E-02  1.5143E-02 -5.3235E-02
 SE:             2.9903E-02  2.2140E-02  2.6714E-02  2.6560E-02  2.3665E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6170E-01  1.3811E-01  2.8223E-01  5.6857E-01  2.4482E-02

 ETASHRINKSD(%)  1.0000E-10  2.5829E+01  1.0503E+01  1.1022E+01  2.0718E+01
 ETASHRINKVR(%)  1.0000E-10  4.4986E+01  1.9903E+01  2.0829E+01  3.7143E+01
 EBVSHRINKSD(%)  2.7662E-01  2.4482E+01  1.0841E+01  1.2001E+01  2.1082E+01
 EBVSHRINKVR(%)  5.5247E-01  4.2970E+01  2.0506E+01  2.2561E+01  3.7719E+01
 RELATIVEINF(%)  9.9446E+01  3.7043E+01  7.3600E+01  5.7169E+01  4.5192E+01
 EPSSHRINKSD(%)  2.3786E+01
 EPSSHRINKVR(%)  4.1914E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3688.5190267521193     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2034.4296669837086     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    53.32
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3688.519       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.53E-01  8.67E-01  1.33E+00  1.07E+00  1.01E+00  9.38E-01  1.22E+00  3.29E+00  8.97E-01  1.16E+00  9.82E-01
 


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
+        2.25E+09
 
 TH 2
+       -1.20E+01  1.36E+09
 
 TH 3
+        2.07E+08 -2.27E+08  3.80E+07
 
 TH 4
+       -5.94E+02 -2.95E+04 -1.66E+01  6.16E+08
 
 TH 5
+       -9.24E+08  1.02E+09  1.70E+08 -3.77E+04  7.58E+08
 
 TH 6
+       -9.61E+00 -5.22E+00  3.67E+02 -1.04E+03  1.52E+01  2.28E+02
 
 TH 7
+        2.94E+02 -3.26E+08 -9.46E+03  6.61E+00  2.43E+08  5.24E+02  7.81E+07
 
 TH 8
+       -2.46E+01 -1.28E+03 -4.88E+00 -4.42E+00 -2.08E+07 -4.51E+01 -1.37E+00  1.14E+06
 
 TH 9
+        1.36E+00  1.31E+09  2.20E+08  3.64E+03 -9.82E+08  3.75E+00  3.15E+08  1.58E+02  1.27E+09
 
 TH10
+       -3.74E+08 -4.11E+08 -6.87E+07  4.79E+04  3.07E+08 -7.83E+00 -9.86E+07  2.06E+03  1.19E+01  2.49E+08
 
 TH11
+        1.16E+00 -4.10E+01 -1.04E+04  2.96E+04 -4.40E+01 -3.50E+00 -1.49E+04  1.28E+03 -1.10E+01  1.06E+01  1.06E+09
 
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
 #CPUT: Total CPU Time in Seconds,       67.638
Stop Time:
Sat Sep 25 04:55:18 CDT 2021
