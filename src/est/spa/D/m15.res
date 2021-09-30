Wed Sep 29 19:45:56 CDT 2021
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
$DATA ../../../../data/spa/D/dat15.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -906.476238893792        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1713E+02 -2.2044E+02 -3.1016E+01 -2.2422E+02  2.0989E+02 -3.6093E+02 -4.2794E+02 -5.9091E+01 -4.4933E+02 -1.2201E+02
            -1.0511E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1388.37440898194        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0805E+00  1.7061E+00  5.7216E-01  6.1111E-01  1.3205E+00  1.4343E+00  1.8437E+00  1.9666E+00  2.0003E+00  1.2170E+00
             1.2880E+00
 PARAMETER:  1.7741E-01  6.3424E-01 -4.5833E-01 -3.9248E-01  3.7801E-01  4.6069E-01  7.1178E-01  7.7630E-01  7.9330E-01  2.9637E-01
             3.5306E-01
 GRADIENT:   3.9666E+02  3.0376E+02 -4.3229E+01  1.0172E+02  9.1067E+01  2.7640E+02  1.5625E+02  8.6909E+00  4.7642E+01  9.1298E+00
             5.7213E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1415.07934889390        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      268
 NPARAMETR:  1.0708E+00  1.9630E+00  6.6366E-01  5.6238E-01  1.2948E+00  1.4141E+00  2.0525E+00  1.2564E+00  1.9176E+00  1.1313E+00
             1.0762E+00
 PARAMETER:  1.6839E-01  7.7450E-01 -3.0998E-01 -4.7558E-01  3.5837E-01  4.4652E-01  8.1906E-01  3.2826E-01  7.5110E-01  2.2335E-01
             1.7345E-01
 GRADIENT:  -7.8745E+01 -2.3278E+01  1.6181E-01  1.3561E+01  4.7376E+01 -1.2760E+02 -4.9589E+01 -5.8433E+00 -4.7204E-01 -5.6768E+00
            -1.1130E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1435.47347916582        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.1452E+00  2.3136E+00  4.1905E-01  3.9672E-01  1.2625E+00  1.7148E+00  2.1150E+00  7.2643E-01  2.1353E+00  1.1093E+00
             1.1105E+00
 PARAMETER:  2.3560E-01  9.3879E-01 -7.6976E-01 -8.2453E-01  3.3310E-01  6.3929E-01  8.4905E-01 -2.1961E-01  8.5861E-01  2.0375E-01
             2.0483E-01
 GRADIENT:   1.2733E+00 -2.2057E+00  1.8811E+00 -6.2445E+00  6.4130E-01 -7.1817E+00 -5.3427E+00 -5.9095E-02 -6.3813E+00 -3.2711E-01
             2.6681E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1435.99122962634        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      625             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1621E+00  2.1532E+00  4.0225E-01  3.9413E-01  1.2637E+00  1.7372E+00  2.1380E+00  3.4828E-01  2.2380E+00  1.1121E+00
             1.1037E+00
 PARAMETER:  2.5021E-01  8.6696E-01 -8.1068E-01 -8.3108E-01  3.3403E-01  6.5226E-01  8.5989E-01 -9.5474E-01  9.0559E-01  2.0628E-01
             1.9870E-01
 GRADIENT:   9.1214E+02  9.2507E+02 -8.9223E-01  9.3403E+01  2.6386E+01  7.0377E+02  4.0206E+02  2.0462E-01  3.9103E+01  2.1650E+00
             2.8827E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1436.67214368758        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      787
 NPARAMETR:  1.1563E+00  2.2058E+00  4.1652E-01  3.9843E-01  1.2494E+00  1.9079E+00  2.1505E+00  1.4841E-01  2.2262E+00  1.1064E+00
             1.1011E+00
 PARAMETER:  2.4523E-01  8.9108E-01 -7.7583E-01 -8.2022E-01  3.2270E-01  7.4603E-01  8.6568E-01 -1.8078E+00  9.0031E-01  2.0110E-01
             1.9631E-01
 GRADIENT:   8.9466E+02  9.6218E+02  6.7290E+00  9.2882E+01  1.1356E+01  8.5855E+02  4.0156E+02  1.7328E-02  3.8950E+01  2.4234E+00
             1.5390E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1437.21115736526        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      947             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1540E+00  2.2116E+00  4.1458E-01  4.0947E-01  1.2438E+00  1.8304E+00  2.1846E+00  5.4775E-02  2.1591E+00  1.0937E+00
             1.1018E+00
 PARAMETER:  2.4322E-01  8.9372E-01 -7.8049E-01 -7.9289E-01  3.1820E-01  7.0455E-01  8.8142E-01 -2.8045E+00  8.6970E-01  1.8960E-01
             1.9694E-01
 GRADIENT:   8.8352E+02  9.5329E+02  4.8961E+00  9.4527E+01  1.7438E+01  7.9308E+02  4.0925E+02  4.6118E-03  3.4623E+01  1.2252E+00
             1.5606E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1437.32100221797        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1128
 NPARAMETR:  1.1526E+00  2.1968E+00  4.1475E-01  4.2291E-01  1.2354E+00  1.8132E+00  2.1841E+00  2.2346E-02  2.1616E+00  1.0921E+00
             1.1018E+00
 PARAMETER:  2.4202E-01  8.8699E-01 -7.8009E-01 -7.6059E-01  3.1141E-01  6.9508E-01  8.8122E-01 -3.7011E+00  8.7087E-01  1.8809E-01
             1.9693E-01
 GRADIENT:   6.9263E+00 -1.4824E+01 -2.0956E+00 -2.8723E-01  3.7182E+00  2.0376E+01 -2.6600E+00  4.1374E-04  3.1507E+00  8.9438E-01
             1.2722E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1437.43694171315        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1311
 NPARAMETR:  1.1544E+00  2.1768E+00  4.1656E-01  4.3323E-01  1.2287E+00  1.8273E+00  2.1979E+00  1.0622E-02  2.1252E+00  1.0806E+00
             1.1024E+00
 PARAMETER:  2.4357E-01  8.7784E-01 -7.7573E-01 -7.3648E-01  3.0599E-01  7.0283E-01  8.8751E-01 -4.4448E+00  8.5386E-01  1.7752E-01
             1.9750E-01
 GRADIENT:   8.0242E+00 -1.5148E+01 -3.5238E+00  1.2252E+00  7.2651E+00  2.3873E+01 -2.0719E+00  1.1075E-04  3.1937E+00  2.1162E-01
             2.6038E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1437.52160968353        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1495
 NPARAMETR:  1.1549E+00  2.1491E+00  4.2067E-01  4.4367E-01  1.2215E+00  1.8385E+00  2.2094E+00  1.0000E-02  2.0934E+00  1.0703E+00
             1.1023E+00
 PARAMETER:  2.4401E-01  8.6503E-01 -7.6591E-01 -7.1267E-01  3.0006E-01  7.0893E-01  8.9272E-01 -7.4046E+00  8.3880E-01  1.6794E-01
             1.9741E-01
 GRADIENT:   8.3247E+00 -1.6606E+01 -4.2634E+00  2.0369E+00  9.6292E+00  2.6663E+01 -2.4200E+00  0.0000E+00  3.3824E+00 -2.7477E-01
             1.3380E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1437.62647063921        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1683
 NPARAMETR:  1.1544E+00  2.1513E+00  4.2568E-01  4.4759E-01  1.2142E+00  1.8289E+00  2.2168E+00  1.0000E-02  2.0673E+00  1.0668E+00
             1.1023E+00
 PARAMETER:  2.4359E-01  8.6608E-01 -7.5407E-01 -7.0387E-01  2.9408E-01  7.0372E-01  8.9604E-01 -7.4046E+00  8.2622E-01  1.6470E-01
             1.9740E-01
 GRADIENT:   8.0278E+00 -1.3692E+01 -1.9487E+00 -7.8169E-02  5.2111E+00  2.4227E+01 -2.5458E+00  0.0000E+00  2.0251E+00 -8.0106E-03
            -7.7610E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1437.68090526035        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1877
 NPARAMETR:  1.1546E+00  2.1417E+00  4.2839E-01  4.5278E-01  1.2090E+00  1.8302E+00  2.2235E+00  1.0000E-02  2.0484E+00  1.0624E+00
             1.1024E+00
 PARAMETER:  2.4372E-01  8.6161E-01 -7.4773E-01 -6.9234E-01  2.8977E-01  7.0441E-01  8.9907E-01 -7.4046E+00  8.1704E-01  1.6048E-01
             1.9745E-01
 GRADIENT:   8.1133E+00 -1.3348E+01 -1.6082E+00 -3.3623E-01  4.5640E+00  2.4518E+01 -2.6678E+00  0.0000E+00  1.7857E+00 -1.8014E-02
            -4.6142E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1437.73402855280        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2065
 NPARAMETR:  1.1539E+00  2.1293E+00  4.2995E-01  4.6010E-01  1.2029E+00  1.8273E+00  2.2337E+00  1.0000E-02  2.0241E+00  1.0564E+00
             1.1026E+00
 PARAMETER:  2.4314E-01  8.5578E-01 -7.4409E-01 -6.7631E-01  2.8473E-01  7.0282E-01  9.0367E-01 -7.4046E+00  8.0513E-01  1.5488E-01
             1.9765E-01
 GRADIENT:   7.6590E+00 -1.3125E+01 -2.2545E+00  3.9227E-01  5.5161E+00  2.3778E+01 -2.2968E+00  0.0000E+00  1.7653E+00 -1.0790E-01
             7.2554E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1437.78387449721        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2258             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1546E+00  2.1147E+00  4.3821E-01  4.6410E-01  1.1904E+00  1.8306E+00  2.2430E+00  1.0000E-02  1.9864E+00  1.0530E+00
             1.1027E+00
 PARAMETER:  2.4379E-01  8.4890E-01 -7.2506E-01 -6.6766E-01  2.7431E-01  7.0462E-01  9.0781E-01 -7.4046E+00  7.8631E-01  1.5169E-01
             1.9781E-01
 GRADIENT:   8.8093E+02  8.3682E+02  6.5208E+00  9.8119E+01  1.0043E+01  7.8924E+02  4.1285E+02  0.0000E+00  3.1538E+01  1.6141E+00
             1.5269E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1437.83398462453        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2446             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1543E+00  2.1083E+00  4.3798E-01  4.6967E-01  1.1886E+00  1.8289E+00  2.2486E+00  1.0000E-02  1.9798E+00  1.0486E+00
             1.1028E+00
 PARAMETER:  2.4348E-01  8.4588E-01 -7.2558E-01 -6.5573E-01  2.7281E-01  7.0372E-01  9.1033E-01 -7.4046E+00  7.8298E-01  1.4748E-01
             1.9789E-01
 GRADIENT:   8.7893E+02  8.2910E+02  4.9661E+00  9.9802E+01  1.3196E+01  7.8766E+02  4.1372E+02  0.0000E+00  3.2108E+01  1.2640E+00
             1.5405E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1437.85612447806        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2635
 NPARAMETR:  1.1546E+00  2.0986E+00  4.4165E-01  4.7421E-01  1.1829E+00  1.8306E+00  2.2555E+00  1.0000E-02  1.9615E+00  1.0457E+00
             1.1030E+00
 PARAMETER:  2.4379E-01  8.4128E-01 -7.1725E-01 -6.4611E-01  2.6793E-01  7.0463E-01  9.1338E-01 -7.4046E+00  7.7371E-01  1.4470E-01
             1.9805E-01
 GRADIENT:   8.1598E+00 -1.1615E+01  1.7605E+00 -4.1723E+00 -2.9378E+00  2.4489E+01 -3.2594E+00  0.0000E+00 -1.4903E-01  5.8686E-01
            -5.7381E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1437.89820598227        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2828             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1546E+00  2.0888E+00  4.4048E-01  4.8225E-01  1.1814E+00  1.8306E+00  2.2637E+00  1.0000E-02  1.9529E+00  1.0389E+00
             1.1032E+00
 PARAMETER:  2.4377E-01  8.3661E-01 -7.1990E-01 -6.2930E-01  2.6670E-01  7.0463E-01  9.1699E-01 -7.4046E+00  7.6932E-01  1.3815E-01
             1.9819E-01
 GRADIENT:   8.7844E+02  8.0662E+02  3.3021E+00  1.0189E+02  1.6787E+01  7.8776E+02  4.1553E+02  0.0000E+00  3.2119E+01  7.5265E-01
             1.5572E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1437.91186926275        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3016
 NPARAMETR:  1.1546E+00  2.0809E+00  4.4111E-01  4.8715E-01  1.1785E+00  1.8306E+00  2.2698E+00  1.0000E-02  1.9408E+00  1.0349E+00
             1.1033E+00
 PARAMETER:  2.4377E-01  8.3282E-01 -7.1845E-01 -6.1917E-01  2.6423E-01  7.0463E-01  9.1970E-01 -7.4046E+00  7.6312E-01  1.3428E-01
             1.9830E-01
 GRADIENT:   8.0674E+00 -1.1859E+01 -1.9792E+00  9.1018E-01  4.8136E+00  2.4444E+01 -2.3682E+00  0.0000E+00  1.3182E+00 -2.6695E-01
            -2.2058E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1437.93988984259        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3209             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1546E+00  2.0738E+00  4.4514E-01  4.8884E-01  1.1720E+00  1.8306E+00  2.2748E+00  1.0000E-02  1.9250E+00  1.0340E+00
             1.1034E+00
 PARAMETER:  2.4379E-01  8.2939E-01 -7.0937E-01 -6.1572E-01  2.5869E-01  7.0463E-01  9.2188E-01 -7.4046E+00  7.5490E-01  1.3343E-01
             1.9842E-01
 GRADIENT:   8.7794E+02  7.8956E+02  4.9814E+00  1.0066E+02  1.2708E+01  7.8700E+02  4.1644E+02  0.0000E+00  3.0838E+01  1.0646E+00
             1.5621E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1437.95191250819        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     3403
 NPARAMETR:  1.1546E+00  2.0678E+00  4.4696E-01  4.9205E-01  1.1686E+00  1.8306E+00  2.2795E+00  1.0000E-02  1.9151E+00  1.0320E+00
             1.1035E+00
 PARAMETER:  2.4379E-01  8.2649E-01 -7.0529E-01 -6.0917E-01  2.5584E-01  7.0463E-01  9.2395E-01 -7.4046E+00  7.4975E-01  1.3152E-01
             1.9852E-01
 GRADIENT:   8.1075E+00 -1.1138E+01  8.4020E-01 -2.5099E+00 -1.4686E+00  2.4408E+01 -2.9362E+00  0.0000E+00  2.3784E-01  3.6415E-01
            -2.7171E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1437.96914285191        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3596             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1546E+00  2.0604E+00  4.4676E-01  4.9714E-01  1.1669E+00  1.8306E+00  2.2856E+00  1.0000E-02  1.9071E+00  1.0270E+00
             1.1036E+00
 PARAMETER:  2.4378E-01  8.2290E-01 -7.0574E-01 -5.9889E-01  2.5436E-01  7.0463E-01  9.2663E-01 -7.4046E+00  7.4559E-01  1.2669E-01
             1.9861E-01
 GRADIENT:   8.7686E+02  7.7455E+02  3.9392E+00  1.0184E+02  1.4991E+01  7.8622E+02  4.1781E+02  0.0000E+00  3.0740E+01  6.8759E-01
             1.5566E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1437.97763523731        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3784
 NPARAMETR:  1.1546E+00  2.0558E+00  4.4770E-01  4.9958E-01  1.1646E+00  1.8306E+00  2.2891E+00  1.0000E-02  1.8998E+00  1.0252E+00
             1.1037E+00
 PARAMETER:  2.4378E-01  8.2067E-01 -7.0363E-01 -5.9398E-01  2.5238E-01  7.0463E-01  9.2818E-01 -7.4046E+00  7.4173E-01  1.2491E-01
             1.9869E-01
 GRADIENT:   8.0671E+00 -1.1222E+01 -6.0290E-01 -5.2844E-01  1.5927E+00  2.4375E+01 -2.5362E+00  0.0000E+00  7.0254E-01 -4.4533E-02
            -2.7872E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1437.98636115390        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     3977             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1547E+00  2.0484E+00  4.5036E-01  5.0213E-01  1.1595E+00  1.8306E+00  2.2945E+00  1.0000E-02  1.8854E+00  1.0241E+00
             1.1039E+00
 PARAMETER:  2.4380E-01  8.1706E-01 -6.9771E-01 -5.8889E-01  2.4799E-01  7.0464E-01  9.3052E-01 -7.4046E+00  7.3414E-01  1.2377E-01
             1.9885E-01
 GRADIENT:   8.7622E+02  7.6123E+02  5.3775E+00  1.0061E+02  1.1460E+01  7.8556E+02  4.1857E+02  0.0000E+00  2.9660E+01  1.0197E+00
             1.5895E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1437.99155283228        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     4168
 NPARAMETR:  1.1547E+00  2.0448E+00  4.5113E-01  5.0424E-01  1.1576E+00  1.8306E+00  2.2974E+00  1.0000E-02  1.8798E+00  1.0226E+00
             1.1040E+00
 PARAMETER:  2.4380E-01  8.1530E-01 -6.9600E-01 -5.8471E-01  2.4632E-01  7.0464E-01  9.3176E-01 -7.4046E+00  7.3115E-01  1.2230E-01
             1.9892E-01
 GRADIENT:   8.0855E+00 -1.0862E+01  7.9391E-01 -2.2901E+00 -1.6897E+00  2.4343E+01 -2.7694E+00  0.0000E+00  7.7653E-02  3.1180E-01
             2.8672E-04

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1438.00010364523        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     4361             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1546E+00  2.0405E+00  4.5087E-01  5.0818E-01  1.1571E+00  1.8306E+00  2.3012E+00  1.0000E-02  1.8756E+00  1.0188E+00
             1.1040E+00
 PARAMETER:  2.4378E-01  8.1318E-01 -6.9657E-01 -5.7693E-01  2.4589E-01  7.0463E-01  9.3343E-01 -7.4046E+00  7.2895E-01  1.1866E-01
             1.9897E-01
 GRADIENT:   8.7542E+02  7.5218E+02  4.1178E+00  1.0195E+02  1.4421E+01  7.8501E+02  4.1923E+02  0.0000E+00  2.9819E+01  5.9168E-01
             1.5633E+00

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1438.00358910754        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     4549
 NPARAMETR:  1.1546E+00  2.0377E+00  4.5155E-01  5.0950E-01  1.1555E+00  1.8306E+00  2.3034E+00  1.0000E-02  1.8713E+00  1.0180E+00
             1.1041E+00
 PARAMETER:  2.4379E-01  8.1184E-01 -6.9506E-01 -5.7433E-01  2.4456E-01  7.0463E-01  9.3438E-01 -7.4046E+00  7.2664E-01  1.1785E-01
             1.9903E-01
 GRADIENT:   8.0499E+00 -1.0824E+01 -3.9472E-01 -4.8620E-01  1.0426E+00  2.4319E+01 -2.5177E+00  0.0000E+00  4.9899E-01 -4.2719E-02
            -1.8486E-02

0ITERATION NO.:  126    OBJECTIVE VALUE:  -1438.00358910754        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     4577
 NPARAMETR:  1.1547E+00  2.0344E+00  4.5346E-01  5.0966E-01  1.1527E+00  1.8306E+00  2.3055E+00  1.0000E-02  1.8633E+00  1.0189E+00
             1.1042E+00
 PARAMETER:  2.4379E-01  8.1184E-01 -6.9506E-01 -5.7433E-01  2.4456E-01  7.0463E-01  9.3438E-01 -7.4046E+00  7.2664E-01  1.1785E-01
             1.9903E-01
 GRADIENT:  -4.6694E-03  2.7423E-01 -3.7698E-01 -4.2137E-02  1.0860E+00 -7.7454E-04 -1.3135E-01  0.0000E+00  2.8621E-01 -4.3757E-02
            -1.6868E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4577
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.6287E-05 -1.0814E-02 -3.4292E-04  2.9087E-02 -4.0127E-02
 SE:             2.9927E-02  2.8675E-02  1.1384E-04  2.0760E-02  1.9311E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9877E-01  7.0607E-01  2.5927E-03  1.6118E-01  3.7711E-02

 ETASHRINKSD(%)  1.0000E-10  3.9341E+00  9.9619E+01  3.0451E+01  3.5307E+01
 ETASHRINKVR(%)  1.0000E-10  7.7135E+00  9.9999E+01  5.1629E+01  5.8148E+01
 EBVSHRINKSD(%)  1.6427E-01  3.8785E+00  9.9753E+01  3.3087E+01  3.2700E+01
 EBVSHRINKVR(%)  3.2828E-01  7.6066E+00  9.9999E+01  5.5227E+01  5.4707E+01
 RELATIVEINF(%)  9.9662E+01  4.2419E+01  1.8830E-04  9.9315E+00  1.6332E+01
 EPSSHRINKSD(%)  4.5511E+01
 EPSSHRINKVR(%)  7.0310E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1438.0035891075411     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -702.85276254380290     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    79.14
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1438.004       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.15E+00  2.04E+00  4.52E-01  5.09E-01  1.16E+00  1.83E+00  2.30E+00  1.00E-02  1.87E+00  1.02E+00  1.10E+00
 


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
+        2.48E+02
 
 TH 2
+       -2.70E-01  4.74E+01
 
 TH 3
+        2.54E+00  3.27E+01  5.18E+02
 
 TH 4
+       -1.95E+00  3.27E+01 -3.44E+02  6.07E+02
 
 TH 5
+       -1.89E+00 -2.46E+01 -2.16E+02  2.72E+02  3.92E+02
 
 TH 6
+        3.89E-03 -3.80E-02 -8.60E-01 -3.60E-01 -5.88E-01  5.93E+01
 
 TH 7
+       -2.02E-02  2.52E+00 -3.50E+01 -1.15E+01  1.15E+01 -1.82E-01  3.16E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.99E-01 -3.00E+00 -2.59E+01  5.31E+01 -1.30E+01 -9.79E-02  9.65E-01  0.00E+00  2.24E+01
 
 TH10
+       -6.95E-02 -2.62E+00 -1.21E+01 -6.15E+00 -6.44E+01 -5.34E-02  1.23E+00  0.00E+00  5.25E+00  5.57E+01
 
 TH11
+       -2.18E+00 -3.00E+00 -2.38E+01  3.26E-01 -1.02E+01  7.64E-01  3.52E+00  0.00E+00  2.46E+00  1.41E+01  1.68E+02
 
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
 #CPUT: Total CPU Time in Seconds,       86.481
Stop Time:
Wed Sep 29 19:47:25 CDT 2021
