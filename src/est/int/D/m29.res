Sat Sep 18 06:47:17 CDT 2021
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
$DATA ../../../../data/int/D/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   20744.4686485815        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.8346E+01  1.5615E+02 -8.0126E+01  7.7429E+00  1.9880E+02 -1.5364E+03 -7.2355E+02 -5.0892E+01 -1.3666E+03 -4.4659E+02
            -4.5115E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1036.27927107082        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.6845E+00  1.7686E+00  1.0375E+00  2.2619E+00  9.4114E-01  5.3765E+00  4.3120E+00  9.9025E-01  3.5873E+00  1.7359E+00
             1.2206E+01
 PARAMETER:  6.2148E-01  6.7019E-01  1.3680E-01  9.1620E-01  3.9341E-02  1.7820E+00  1.5614E+00  9.0199E-02  1.3774E+00  6.5151E-01
             2.6019E+00
 GRADIENT:   1.5154E+01  1.0016E+01 -4.0391E+01  7.1887E+01 -1.9375E+01  1.6148E+02  6.9373E+01  4.0722E+00  7.6918E+01  4.2527E+01
             5.8835E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1128.50045467463        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.2520E+00  1.2483E+00  3.8221E+01  4.4995E+00  1.5603E+00  3.2629E+00  1.1306E+01  8.5699E-01  3.0327E+00  2.2450E+00
             1.1743E+01
 PARAMETER:  3.2475E-01  3.2175E-01  3.7434E+00  1.6040E+00  5.4487E-01  1.2826E+00  2.5254E+00 -5.4332E-02  1.2095E+00  9.0873E-01
             2.5633E+00
 GRADIENT:  -2.4312E+01  2.6307E+01  5.4221E+00  1.2644E+02 -8.0412E+01  9.5503E+01  1.1684E+01  9.4297E-02 -1.1729E+01  5.7554E+01
             5.6663E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1449.89395481417        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0418E+00  1.7851E+00  1.2157E+01  8.3531E-01  2.5911E+00  2.4545E+00  3.6413E+00  4.4123E-01  2.7415E+00  5.7501E-01
             7.6695E+00
 PARAMETER:  1.4093E-01  6.7949E-01  2.5979E+00 -7.9956E-02  1.0521E+00  9.9793E-01  1.3923E+00 -7.1820E-01  1.1085E+00 -4.5337E-01
             2.1373E+00
 GRADIENT:  -5.6354E+01 -1.7569E+01 -4.8467E+00 -1.9603E+01  1.0387E+02  1.5696E+01  4.5693E+01 -2.5455E-03  2.1532E+01  4.4107E+00
             1.9633E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1486.53881880548        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      295
 NPARAMETR:  1.1992E+00  2.0452E+00  3.7315E+00  7.0213E-01  1.9041E+00  2.2789E+00  2.8662E+00  2.1814E-01  2.4467E+00  4.0158E-01
             6.8956E+00
 PARAMETER:  2.8165E-01  8.1551E-01  1.4168E+00 -2.5364E-01  7.4403E-01  9.2371E-01  1.1530E+00 -1.4226E+00  9.9475E-01 -8.1234E-01
             2.0309E+00
 GRADIENT:   4.3569E+00  9.5274E-01  9.7643E-01 -5.0711E+00 -2.4994E+00 -3.0883E+00 -2.3490E+00  6.0259E-03  8.2717E+00  1.9858E+00
            -3.3460E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1495.48160464732        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      367
 NPARAMETR:  1.1630E+00  1.4687E+00  3.4975E+00  9.4717E-01  1.7200E+00  2.3068E+00  3.6622E+00  5.2058E-01  8.0378E-01  2.2201E-01
             6.9241E+00
 PARAMETER:  2.5100E-01  4.8437E-01  1.3521E+00  4.5719E-02  6.4231E-01  9.3585E-01  1.3981E+00 -5.5281E-01 -1.1843E-01 -1.4050E+00
             2.0350E+00
 GRADIENT:  -6.9268E+00 -1.7963E+00 -9.7126E-01 -8.5789E-01 -5.7929E+00  4.3137E+00 -1.1665E+00  1.1132E-01  4.9119E-01  1.8160E-01
            -1.0244E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1496.04438154397        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      437
 NPARAMETR:  1.1873E+00  1.4722E+00  4.3137E+00  9.5777E-01  1.8157E+00  2.2758E+00  3.7221E+00  3.9899E-01  7.7245E-01  2.2127E-01
             6.9742E+00
 PARAMETER:  2.7170E-01  4.8678E-01  1.5618E+00  5.6855E-02  6.9645E-01  9.2235E-01  1.4143E+00 -8.1883E-01 -1.5819E-01 -1.4084E+00
             2.0422E+00
 GRADIENT:   5.5943E-01 -8.2183E-01 -1.1634E-01  1.7890E-01  1.0105E+00 -6.4903E-01  1.1862E-01  4.4044E-02  3.5266E-01  1.4339E-01
             5.9328E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1496.05329048288        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      507
 NPARAMETR:  1.1872E+00  1.5137E+00  4.2390E+00  9.3382E-01  1.8247E+00  2.2764E+00  3.6719E+00  3.0046E-01  7.1850E-01  2.1011E-01
             6.9767E+00
 PARAMETER:  2.7161E-01  5.1456E-01  1.5443E+00  3.1524E-02  7.0143E-01  9.2258E-01  1.4007E+00 -1.1024E+00 -2.3059E-01 -1.4601E+00
             2.0426E+00
 GRADIENT:   5.2785E-01 -2.9711E-01  1.4419E-01 -1.0266E-02  2.6146E-01 -5.7934E-01 -1.9316E-01  2.8252E-02  7.4337E-02  1.0987E-01
             5.4550E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1496.05333742679        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  1.1869E+00  1.5190E+00  4.2080E+00  9.3068E-01  1.8240E+00  2.2770E+00  3.6656E+00  2.8423E-01  7.1171E-01  2.0817E-01
             6.9764E+00
 PARAMETER:  2.7136E-01  5.1805E-01  1.5370E+00  2.8158E-02  7.0106E-01  9.2287E-01  1.3990E+00 -1.1580E+00 -2.4009E-01 -1.4694E+00
             2.0425E+00
 GRADIENT:   4.3337E-01 -2.2569E-01  1.2305E-01 -1.5537E-02  1.8764E-01 -4.7276E-01 -1.7353E-01  2.5957E-02  4.1664E-02  1.0620E-01
             4.1573E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1496.05337461083        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      650
 NPARAMETR:  1.1866E+00  1.5241E+00  4.1763E+00  9.2765E-01  1.8232E+00  2.2777E+00  3.6597E+00  2.6322E-01  7.0511E-01  2.0582E-01
             6.9761E+00
 PARAMETER:  2.7109E-01  5.2139E-01  1.5294E+00  2.4901E-02  7.0062E-01  9.2316E-01  1.3974E+00 -1.2348E+00 -2.4940E-01 -1.4808E+00
             2.0425E+00
 GRADIENT:   3.3465E-01 -1.5708E-01  9.7883E-02 -1.8421E-02  1.2173E-01 -3.6468E-01 -1.4877E-01  2.2833E-02  1.0452E-02  1.0198E-01
             2.8666E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1496.05342991739        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  1.1863E+00  1.5282E+00  4.1497E+00  9.2521E-01  1.8225E+00  2.2782E+00  3.6549E+00  2.4170E-01  6.9978E-01  2.0340E-01
             6.9759E+00
 PARAMETER:  2.7086E-01  5.2407E-01  1.5230E+00  2.2268E-02  7.0021E-01  9.2340E-01  1.3961E+00 -1.3200E+00 -2.5700E-01 -1.4926E+00
             2.0425E+00
 GRADIENT:   2.4931E-01 -1.0089E-01  7.5077E-02 -1.8987E-02  6.9871E-02 -2.7282E-01 -1.2537E-01  1.9641E-02 -1.4514E-02  9.7907E-02
             1.7785E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1496.05347164097        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      792
 NPARAMETR:  1.1860E+00  1.5322E+00  4.1229E+00  9.2284E-01  1.8217E+00  2.2788E+00  3.6503E+00  2.1656E-01  6.9459E-01  2.0043E-01
             6.9758E+00
 PARAMETER:  2.7062E-01  5.2668E-01  1.5166E+00  1.9706E-02  6.9978E-01  9.2365E-01  1.3948E+00 -1.4299E+00 -2.6443E-01 -1.5073E+00
             2.0424E+00
 GRADIENT:   1.6174E-01 -4.5510E-02  5.0851E-02 -1.8249E-02  2.0545E-02 -1.7953E-01 -9.9879E-02  1.6062E-02 -3.8594E-02  9.3217E-02
             6.8754E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1496.05351387417        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      862
 NPARAMETR:  1.1857E+00  1.5368E+00  4.0911E+00  9.2012E-01  1.8207E+00  2.2795E+00  3.6449E+00  1.8136E-01  6.8864E-01  1.9585E-01
             6.9756E+00
 PARAMETER:  2.7033E-01  5.2968E-01  1.5088E+00  1.6744E-02  6.9925E-01  9.2395E-01  1.3933E+00 -1.6073E+00 -2.7304E-01 -1.5304E+00
             2.0424E+00
 GRADIENT:   5.5089E-02  1.9640E-02  2.0822E-02 -1.6006E-02 -3.5734E-02 -6.7581E-02 -6.7739E-02  1.1494E-02 -6.6223E-02  8.6564E-02
            -6.0198E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1496.05368047224        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      932
 NPARAMETR:  1.1853E+00  1.5416E+00  4.0563E+00  9.1725E-01  1.8196E+00  2.2802E+00  3.6392E+00  1.3509E-01  6.8245E-01  1.8866E-01
             6.9756E+00
 PARAMETER:  2.7002E-01  5.3283E-01  1.5003E+00  1.3620E-02  6.9864E-01  9.2427E-01  1.3918E+00 -1.9018E+00 -2.8207E-01 -1.5678E+00
             2.0424E+00
 GRADIENT:  -6.1678E-02  8.9494E-02 -1.3083E-02 -1.3076E-02 -9.4480E-02  5.6321E-02 -3.1132E-02  6.4939E-03 -9.5137E-02  7.7262E-02
            -1.9757E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1496.05458904531        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1003
 NPARAMETR:  1.1849E+00  1.5472E+00  4.0145E+00  9.1397E-01  1.8183E+00  2.2810E+00  3.6325E+00  6.7251E-02  6.7569E-01  1.7314E-01
             6.9760E+00
 PARAMETER:  2.6964E-01  5.3647E-01  1.4899E+00  1.0039E-02  6.9792E-01  9.2463E-01  1.3899E+00 -2.5993E+00 -2.9201E-01 -1.6536E+00
             2.0425E+00
 GRADIENT:  -1.9896E-01  1.7146E-01 -5.3635E-02 -1.1072E-02 -1.6159E-01  2.0503E-01  1.0917E-02  1.6296E-03 -1.2795E-01  6.0615E-02
            -3.4478E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1496.07002005942        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  1.1845E+00  1.5506E+00  3.9565E+00  9.1206E-01  1.8163E+00  2.2812E+00  3.6266E+00  1.0000E-02  6.7574E-01  5.9117E-02
             6.9809E+00
 PARAMETER:  2.6932E-01  5.3864E-01  1.4753E+00  7.9461E-03  6.9680E-01  9.2469E-01  1.3883E+00 -1.1778E+01 -2.9194E-01 -2.7282E+00
             2.0432E+00
 GRADIENT:  -3.2052E-01  1.8976E-01 -1.3625E-01 -3.7704E-02 -4.5152E-02  3.0621E-01  6.3083E-02  0.0000E+00 -1.2614E-01  5.2112E-03
            -1.2922E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1496.24456432007        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1239
 NPARAMETR:  1.1944E+00  1.4668E+00  4.4152E+00  9.6819E-01  1.8224E+00  2.3106E+00  3.8346E+00  1.0000E-02  7.4612E-01  1.0000E-02
             6.9946E+00
 PARAMETER:  2.7765E-01  4.8307E-01  1.5851E+00  6.7677E-02  7.0017E-01  9.3751E-01  1.4441E+00 -4.8055E+01 -1.9286E-01 -6.9164E+00
             2.0451E+00
 GRADIENT:   6.4672E-02  3.2110E-02  9.9430E-03 -3.5681E-02 -1.1355E-02  1.2904E-01  3.6023E-02  0.0000E+00  3.2935E-02  0.0000E+00
             2.3050E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1496.24465588861        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1415
 NPARAMETR:  1.1942E+00  1.4674E+00  4.4088E+00  9.6741E-01  1.8224E+00  2.3099E+00  3.8333E+00  1.0000E-02  7.4312E-01  1.0000E-02
             6.9939E+00
 PARAMETER:  2.7747E-01  4.8347E-01  1.5836E+00  6.6863E-02  7.0015E-01  9.3719E-01  1.4437E+00 -4.7266E+01 -1.9689E-01 -6.8264E+00
             2.0450E+00
 GRADIENT:   3.9331E-03  2.3707E-03  1.9654E-04 -2.6763E-03  2.0314E-03  1.2853E-02  1.3943E-03  0.0000E+00  1.2058E-03  0.0000E+00
            -5.4398E-03

0ITERATION NO.:   86    OBJECTIVE VALUE:  -1496.24465588861        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1437
 NPARAMETR:  1.1942E+00  1.4674E+00  4.4088E+00  9.6741E-01  1.8224E+00  2.3099E+00  3.8333E+00  1.0000E-02  7.4312E-01  1.0000E-02
             6.9939E+00
 PARAMETER:  2.7747E-01  4.8347E-01  1.5836E+00  6.6863E-02  7.0015E-01  9.3719E-01  1.4437E+00 -4.7266E+01 -1.9689E-01 -6.8264E+00
             2.0450E+00
 GRADIENT:   3.9331E-03  2.3707E-03  1.9654E-04 -2.6763E-03  2.0314E-03  1.2853E-02  1.3943E-03  0.0000E+00  1.2058E-03  0.0000E+00
            -5.4398E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1437
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.4788E-03  1.3656E-02 -3.6394E-05 -4.3397E-02  3.1819E-05
 SE:             2.9230E-02  2.6877E-02  3.6873E-05  1.0497E-02  1.9497E-04
 N:                     100         100         100         100         100

 P VAL.:         7.9806E-01  6.1140E-01  3.2363E-01  3.5661E-05  8.7036E-01

 ETASHRINKSD(%)  2.0759E+00  9.9575E+00  9.9876E+01  6.4832E+01  9.9347E+01
 ETASHRINKVR(%)  4.1087E+00  1.8923E+01  1.0000E+02  8.7632E+01  9.9996E+01
 EBVSHRINKSD(%)  1.8698E+00  6.2191E+00  9.9870E+01  7.0661E+01  9.9361E+01
 EBVSHRINKVR(%)  3.7046E+00  1.2051E+01  1.0000E+02  9.1392E+01  9.9996E+01
 RELATIVEINF(%)  9.6170E+01  3.6972E+01  2.9928E-05  3.1699E+00  7.9935E-04
 EPSSHRINKSD(%)  7.8817E+00
 EPSSHRINKVR(%)  1.5142E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1496.2446558886115     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       157.84470387979923     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.03
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    15.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1496.245       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.19E+00  1.47E+00  4.41E+00  9.67E-01  1.82E+00  2.31E+00  3.83E+00  1.00E-02  7.43E-01  1.00E-02  6.99E+00
 


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
+        5.97E+01
 
 TH 2
+        1.16E+01  5.84E+00
 
 TH 3
+        3.83E+00 -3.95E-01  6.09E-01
 
 TH 4
+        1.07E+00  2.70E+01 -8.48E+00  2.01E+02
 
 TH 5
+       -4.47E+01  2.62E+00 -6.48E+00  8.41E+01  6.93E+01
 
 TH 6
+        4.89E+00  7.34E-01  3.83E-01 -1.55E+00 -4.36E+00  4.20E-01
 
 TH 7
+       -1.28E+00 -2.27E+00  5.64E-01 -1.52E+01 -5.45E+00  1.84E-02  1.17E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.35E+00 -4.78E+00  9.76E-01 -2.96E+01 -9.21E+00 -1.17E-01  2.32E+00  0.00E+00  4.65E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.61E+00 -1.74E+00  2.31E-01 -9.30E+00 -2.04E+00 -8.15E-02  7.52E-01  0.00E+00  1.53E+00  0.00E+00  1.04E+00
 
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
+        1.39E+02
 
 TH 2
+       -6.30E-01  3.10E+01
 
 TH 3
+        3.34E-01  8.95E-01  1.21E+00
 
 TH 4
+       -3.46E+00  4.14E+01 -5.27E+00  2.38E+02
 
 TH 5
+       -3.40E+00 -1.28E+01 -1.08E+01  4.34E+01  1.23E+02
 
 TH 6
+        1.21E+00 -6.04E-02  5.48E-02  8.24E-01 -1.50E+00  3.32E+01
 
 TH 7
+        2.77E-01  2.96E+00 -4.52E-01 -2.17E+01  3.31E+00 -3.05E-01  9.33E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.84E-01 -3.63E+00 -3.14E-01 -4.05E+01  5.16E+00 -1.09E-01  4.43E+00  0.00E+00  1.87E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.57E+00 -3.42E+00 -9.04E-02 -1.18E+01  2.75E-01  1.31E+00  1.22E+00  0.00E+00  2.21E+00  0.00E+00  2.29E+01
 
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
+        1.41E+02
 
 TH 2
+        5.30E+01  3.10E+01
 
 TH 3
+        8.68E-01  5.40E-01  9.19E-01
 
 TH 4
+        5.93E+01  4.76E+01 -4.55E+00  2.40E+02
 
 TH 5
+       -1.23E+01 -9.17E+00 -1.02E+01  4.91E+01  1.35E+02
 
 TH 6
+        2.97E+01  1.25E+01 -1.73E-01 -4.92E+00 -3.80E-01  3.34E+01
 
 TH 7
+        8.71E+00  2.21E+00 -3.87E-01 -2.36E+01  3.05E+00  1.14E+01  9.18E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.17E+00 -3.05E+00 -1.16E-01 -3.39E+01  3.75E+00 -5.08E-01  4.25E+00  0.00E+00  1.52E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.59E+01 -3.09E+01 -7.68E-01 -1.64E+01  1.66E+01  3.13E+00 -3.69E+00  0.00E+00 -9.77E+00  0.00E+00  9.74E+02
 
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
 #CPUT: Total CPU Time in Seconds,       44.337
Stop Time:
Sat Sep 18 06:48:03 CDT 2021
