Thu Sep 30 03:12:58 CDT 2021
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
$DATA ../../../../data/spa1/D/dat53.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   19445.0020286165        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.2377E+02  3.9340E+02 -3.0506E+01  2.2451E+02  5.7052E+01 -1.7717E+03 -7.9178E+02 -4.9969E+01 -1.3006E+03 -5.1929E+02
            -3.8214E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -611.121360729995        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.2527E+00  1.0850E+00  9.9587E-01  1.5546E+00  1.1738E+00  1.9591E+00  1.1973E+00  9.7938E-01  1.3272E+00  1.0615E+00
             1.4385E+01
 PARAMETER:  3.2527E-01  1.8160E-01  9.5862E-02  5.4120E-01  2.6027E-01  7.7249E-01  2.8009E-01  7.9168E-02  3.8307E-01  1.5969E-01
             2.7662E+00
 GRADIENT:  -2.1006E+01  1.4669E+01 -4.0412E+00  1.0272E+01 -8.1198E+00  4.2156E+01 -1.4526E+00  4.9689E+00 -2.0371E+00  4.1093E+00
             1.9635E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -644.051023937947        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.2530E+00  9.1929E-01  2.1792E+00  2.0452E+00  1.4887E+01  1.8010E+00  5.6274E+00  5.1301E-01  1.7894E+00  2.7761E+00
             1.2527E+01
 PARAMETER:  3.2557E-01  1.5851E-02  8.7898E-01  8.1550E-01  2.8005E+00  6.8835E-01  1.8276E+00 -5.6746E-01  6.8185E-01  1.1211E+00
             2.6278E+00
 GRADIENT:   1.7927E+01  2.5741E+01  8.2345E+00  5.0285E+01 -4.4355E-02 -1.7945E+01  2.5787E+01 -1.2082E-01  2.0272E+01 -1.8004E-02
             1.3724E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -665.876570427903        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  1.0375E+00  9.2016E-01  2.2257E+00  1.9985E+00  1.2698E+03  1.8341E+00  5.8644E+00  5.5624E-01  1.7583E+00  3.0973E+00
             1.0917E+01
 PARAMETER:  1.3681E-01  1.6790E-02  9.0008E-01  7.9237E-01  7.2466E+00  7.0658E-01  1.8689E+00 -4.8655E-01  6.6436E-01  1.2305E+00
             2.4903E+00
 GRADIENT:  -3.8117E+01  2.3250E+01  4.4198E+00  8.0097E+01  7.9528E-04  1.1053E+01  1.8087E+01 -5.7453E-02  1.3600E+01 -1.1566E-06
             6.7733E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -668.905834792448        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      482
 NPARAMETR:  1.1190E+00  9.2056E-01  2.2340E+00  1.9791E+00  2.3445E+03  1.8402E+00  5.9062E+00  5.9862E-01  1.7518E+00  3.1606E+00
             1.0479E+01
 PARAMETER:  2.1246E-01  1.7223E-02  9.0377E-01  7.8263E-01  7.8598E+00  7.0989E-01  1.8760E+00 -4.1312E-01  6.6063E-01  1.2508E+00
             2.4494E+00
 GRADIENT:   5.0641E+00  2.3860E+01  4.2490E+00  6.1843E+01  2.2948E-03  5.6246E-01 -1.1653E+01 -9.8948E-02  6.0851E+00 -2.9979E-07
            -6.5110E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -672.412741102390        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  1.0855E+00  9.2268E-01  2.2747E+00  1.7160E+00  9.2470E+02  1.8736E+00  6.4467E+00  3.4655E+00  1.7051E+00  3.5619E+00
             1.0152E+01
 PARAMETER:  1.8200E-01  1.9522E-02  9.2183E-01  6.3998E-01  6.9295E+00  7.2784E-01  1.9636E+00  1.3429E+00  6.3361E-01  1.3703E+00
             2.4176E+00
 GRADIENT:  -8.7578E-01  1.9010E+01  3.1656E+00  3.6514E+00 -2.9484E-03  5.9728E+00  1.5945E+01 -1.9727E-02  2.1546E+01 -5.5806E-06
            -6.8905E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -673.176356885252        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      838
 NPARAMETR:  1.1095E+00  9.2165E-01  2.2523E+00  1.7757E+00  3.0529E+04  1.8518E+00  6.0468E+00  4.6577E+00  1.6984E+00  3.3977E+00
             1.0401E+01
 PARAMETER:  2.0391E-01  1.8416E-02  9.1196E-01  6.7421E-01  1.0426E+01  7.1615E-01  1.8995E+00  1.6385E+00  6.2967E-01  1.3231E+00
             2.4419E+00
 GRADIENT:  -2.4587E+00  2.5951E+01  3.8478E+00  1.7183E+01 -1.8931E-04 -6.2488E+00  1.2499E+01 -2.5113E+00  1.7282E+01 -8.4550E-09
             1.3325E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -677.390884499484        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      987
 NPARAMETR:  1.0789E+00  9.1024E-01  2.2639E+00  1.6466E+00  2.0272E+08  1.8515E+00  5.1145E+00  4.5822E+00  1.5919E+00  3.2093E+00
             9.8945E+00
 PARAMETER:  1.7596E-01  5.9512E-03  9.1707E-01  5.9868E-01  1.9227E+01  7.1602E-01  1.7321E+00  1.6222E+00  5.6490E-01  1.2661E+00
             2.3920E+00
 GRADIENT:  -1.7092E+01  2.6375E+01  6.9890E+00  3.2857E+01 -4.3568E-08 -1.8415E+01  2.4535E+01  5.1656E+00  1.5779E+01  0.0000E+00
             1.4537E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -682.979657507082        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:     1101
 NPARAMETR:  1.0887E+00  9.0215E-01  2.2325E+00  1.5590E+00  2.0347E+08  1.8356E+00  4.2417E+00  4.5374E+00  1.5428E+00  3.2093E+00
             9.9182E+00
 PARAMETER:  1.8495E-01 -2.9691E-03  9.0314E-01  5.4405E-01  1.9231E+01  7.0738E-01  1.5450E+00  1.6124E+00  5.3361E-01  1.2661E+00
             2.3944E+00
 GRADIENT:  -2.2244E+01  3.5116E+01  1.2891E+01  6.8648E+00 -7.6348E-08 -2.5193E+01 -3.8934E+00 -1.8883E+01  7.4167E+00  0.0000E+00
            -4.9770E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -684.374198611624        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:     1221
 NPARAMETR:  1.1474E+00  8.9800E-01  2.2129E+00  1.5718E+00  1.7854E+08  1.8434E+00  4.2856E+00  4.6883E+00  1.5483E+00  3.0291E+00
             1.0153E+01
 PARAMETER:  2.3749E-01 -7.5903E-03  8.9429E-01  5.5220E-01  1.9100E+01  7.1159E-01  1.5553E+00  1.6451E+00  5.3717E-01  1.2083E+00
             2.4178E+00
 GRADIENT:   1.9262E+01  3.7379E+01  1.4653E+01  1.1485E+01 -7.4467E-08 -3.1303E+00  1.2708E+01  4.5814E+00  7.6586E+00  0.0000E+00
             3.1615E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -687.473132486947        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1410             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1364E+00  8.9783E-01  2.2087E+00  1.5594E+00  1.4181E+08  1.8810E+00  4.3258E+00  4.7486E+00  1.5336E+00  2.9907E+00
             9.9819E+00
 PARAMETER:  2.2787E-01 -7.7778E-03  8.9238E-01  5.4429E-01  1.8870E+01  7.3182E-01  1.5646E+00  1.6579E+00  5.2762E-01  1.1955E+00
             2.4008E+00
 GRADIENT:   2.2568E+01  2.8201E+01  1.0074E+01  7.6668E+00 -8.8794E-08  1.0561E+01  1.1839E+01  4.1823E+00  4.4729E+00  0.0000E+00
             1.4821E+01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -687.473132486947        NO. OF FUNC. EVALS.:  69
 CUMULATIVE NO. OF FUNC. EVALS.:     1479
 NPARAMETR:  1.1370E+00  8.9738E-01  2.2133E+00  1.5572E+00  1.7127E+08  1.8850E+00  4.3167E+00  4.7322E+00  1.5315E+00  3.0086E+00
             1.0015E+01
 PARAMETER:  2.2787E-01 -7.7778E-03  8.9238E-01  5.4429E-01  1.8870E+01  7.3182E-01  1.5646E+00  1.6579E+00  5.2762E-01  1.1955E+00
             2.4008E+00
 GRADIENT:  -9.6900E+01  3.0107E+01 -1.1517E+02  1.0641E+02 -5.8478E-07 -8.8622E+01  1.7135E+02  5.3440E+01  1.1245E+02 -2.5575E-05
            -1.2332E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1479
 NO. OF SIG. DIGITS UNREPORTABLE

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.5060E-03 -1.2816E-02 -1.4991E-02 -6.4487E-02  2.9495E-11
 SE:             2.9225E-02  2.1196E-02  1.0823E-02  1.5202E-02  5.3847E-11
 N:                     100         100         100         100         100

 P VAL.:         8.5056E-01  5.4540E-01  1.6603E-01  2.2177E-05  5.8386E-01

 ETASHRINKSD(%)  2.0919E+00  2.8991E+01  6.3740E+01  4.9071E+01  1.0000E+02
 ETASHRINKVR(%)  4.1401E+00  4.9577E+01  8.6852E+01  7.4062E+01  1.0000E+02
 EBVSHRINKSD(%)  6.1501E+00  3.7001E+01  7.5273E+01  3.8760E+01  1.0000E+02
 EBVSHRINKVR(%)  1.1922E+01  6.0311E+01  9.3886E+01  6.2497E+01  1.0000E+02
 RELATIVEINF(%)  8.0370E+01  1.2169E+01  3.9041E+00  9.4099E+00  0.0000E+00
 EPSSHRINKSD(%)  1.0795E+01
 EPSSHRINKVR(%)  2.0425E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -687.47313248694684     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       231.46540071772586     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.67
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -687.473       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.14E+00  8.98E-01  2.21E+00  1.56E+00  1.42E+08  1.88E+00  4.33E+00  4.75E+00  1.53E+00  2.99E+00  9.98E+00
 


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
+        9.18E+03
 
 TH 2
+        2.58E+04  7.35E+04
 
 TH 3
+        5.91E+01  1.21E+04  7.00E+02
 
 TH 4
+       -1.21E+02  3.29E+02 -1.27E+03  3.93E+03
 
 TH 5
+        1.33E-12  6.83E-13  3.39E-14  7.44E-13  6.00E-23
 
 TH 6
+        6.78E+01 -1.84E+02  4.36E+03 -1.02E+04  1.16E-13  6.62E+03
 
 TH 7
+        1.47E+00 -5.99E+01 -9.89E+00  3.01E+00  7.22E-16 -1.59E+01  2.69E+02
 
 TH 8
+       -1.52E+01  2.58E+01 -1.69E+02  1.75E+03  6.68E-14 -1.08E+03  9.41E+00  1.96E+02
 
 TH 9
+       -1.32E+02  2.86E+02 -1.34E+03  3.15E+03 -5.81E-13 -1.07E+04  3.77E+00  1.84E+03  4.28E+03
 
 TH10
+        2.60E-04  7.15E-04 -4.91E-05 -3.33E-05 -3.84E-14  9.75E-05  4.76E-06  9.65E-06 -3.48E-04 -6.62E-05
 
 TH11
+       -7.67E+00 -1.01E+03 -4.40E+01  1.02E+02  1.55E-14 -5.98E+01  1.48E-01  9.83E+00  1.15E+02  1.01E-05  2.74E+01
 
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
 #CPUT: Total CPU Time in Seconds,       44.190
Stop Time:
Thu Sep 30 03:13:44 CDT 2021
