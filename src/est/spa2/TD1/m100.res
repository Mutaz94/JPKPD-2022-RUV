Thu Sep 30 07:48:49 CDT 2021
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
$DATA ../../../../data/spa2/TD1/dat100.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2155.87300185779        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0893E+02 -5.4097E+00  8.2427E+01  1.1873E+01  4.6908E+01  3.1115E+01  6.2482E+00 -4.2797E+02 -7.1350E+01 -4.1856E+00
            -7.2965E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2334.86674110195        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      179
 NPARAMETR:  1.3523E+00  1.0626E+00  9.2676E-01  1.0471E+00  9.7174E-01  1.1693E+00  9.8957E-01  1.7183E+00  1.0550E+00  9.9999E-01
             1.0468E+00
 PARAMETER:  4.0183E-01  1.6069E-01  2.3939E-02  1.4604E-01  7.1336E-02  2.5642E-01  8.9518E-02  6.4134E-01  1.5358E-01  9.9987E-02
             1.4578E-01
 GRADIENT:   5.7656E+02  4.3268E+01 -4.1540E+01  1.0332E+02 -2.0144E+01 -1.4441E+02  9.3730E+00  3.2968E+01  2.9500E+01  2.4331E+00
             3.2612E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2368.48296952392        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      358
 NPARAMETR:  1.3540E+00  1.0150E+00  9.2631E-01  1.0478E+00  9.7221E-01  1.5491E+00  9.8909E-01  7.9123E-01  9.1319E-01  1.0438E+00
             1.0476E+00
 PARAMETER:  4.0308E-01  1.1490E-01  2.3457E-02  1.4673E-01  7.1820E-02  5.3765E-01  8.9030E-02 -1.3417E-01  9.1839E-03  1.4283E-01
             1.4649E-01
 GRADIENT:   3.4064E+02  8.0028E+00 -4.2894E+01  9.9494E+01  3.4146E+01  1.9926E+00  9.8649E-01 -3.1888E-01  7.0923E-01  1.3348E+00
             1.4013E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2393.75267605604        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      540
 NPARAMETR:  1.1648E+00  1.0206E+00  9.2757E-01  1.0419E+00  9.7104E-01  1.7205E+00  9.8919E-01  1.1895E+00  8.3451E-01  9.6612E-01
             1.0465E+00
 PARAMETER:  2.5258E-01  1.2036E-01  2.4812E-02  1.4104E-01  7.0615E-02  6.4262E-01  8.9127E-02  2.7354E-01 -8.0911E-02  6.5529E-02
             1.4549E-01
 GRADIENT:   1.5351E+02  1.8176E+01 -5.5556E+01  9.5692E+01  2.4857E+01  1.0152E+02 -3.7040E+00  1.3439E+01 -1.7942E+01 -6.0056E+00
             1.7660E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2448.91732892605        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  9.5531E-01  1.0155E+00  9.2901E-01  1.0344E+00  9.6973E-01  1.0769E+00  9.8907E-01  8.1023E-01  9.0911E-01  1.0441E+00
             1.0455E+00
 PARAMETER:  5.4277E-02  1.1541E-01  2.6365E-02  1.3386E-01  6.9264E-02  1.7407E-01  8.9005E-02 -1.1044E-01  4.7145E-03  1.4312E-01
             1.4453E-01
 GRADIENT:  -1.0043E+01  2.8804E+00 -3.9603E+01  7.5971E+01  2.4411E+01  1.5895E+01  7.2874E-01  2.9560E-01  3.3771E-01  2.0101E+00
             1.5314E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2449.26921816199        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      892
 NPARAMETR:  9.5864E-01  1.0124E+00  9.2898E-01  1.0345E+00  9.6976E-01  1.0325E+00  9.8905E-01  8.1365E-01  9.0661E-01  1.0328E+00
             1.0456E+00
 PARAMETER:  5.7758E-02  1.1230E-01  2.6333E-02  1.3395E-01  6.9296E-02  1.3198E-01  8.8990E-02 -1.0622E-01  1.9535E-03  1.3231E-01
             1.4456E-01
 GRADIENT:  -3.8942E+00  2.5970E-02 -4.0439E+01  7.4079E+01  2.7672E+01 -1.9888E-01  1.0610E-01  2.4447E-01 -9.0015E-02 -6.3013E-02
             1.5198E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2451.37778488091        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  9.5976E-01  1.0399E+00  9.4025E-01  9.8887E-01  9.6141E-01  1.0322E+00  9.8999E-01  8.2810E-01  9.1077E-01  1.0345E+00
             1.0350E+00
 PARAMETER:  5.8928E-02  1.3913E-01  3.8389E-02  8.8805E-02  6.0643E-02  1.3170E-01  8.9937E-02 -8.8621E-02  6.5351E-03  1.3391E-01
             1.3436E-01
 GRADIENT:  -7.9032E-01  5.3496E-01 -1.3518E+01  9.2953E-01 -2.6395E+01 -1.8532E-01  1.4873E+00  1.7270E-01  9.0668E-03  1.1461E-02
             4.5930E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2453.39978745251        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1213
 NPARAMETR:  9.4738E-01  1.0465E+00  1.0958E+00  9.8647E-01  1.0498E+00  1.0253E+00  9.6111E-01  9.6188E-01  9.1005E-01  1.0566E+00
             1.0232E+00
 PARAMETER:  4.5950E-02  1.4548E-01  1.9152E-01  8.6373E-02  1.4863E-01  1.2497E-01  6.0336E-02  6.1136E-02  5.7494E-03  1.5509E-01
             1.2292E-01
 GRADIENT:   3.8152E+02  5.4016E+01  7.3020E+00  7.3851E+01  3.4884E+01  6.6217E+01  8.5548E-01 -2.6231E+00  7.8014E+00 -1.9848E+00
            -6.1644E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2454.17143455467        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1393
 NPARAMETR:  9.4849E-01  1.0879E+00  1.1192E+00  9.6985E-01  1.0760E+00  1.0190E+00  9.9947E-01  1.0850E+00  8.8133E-01  1.1012E+00
             1.0336E+00
 PARAMETER:  4.7120E-02  1.8421E-01  2.1260E-01  6.9391E-02  1.7329E-01  1.1887E-01  9.9466E-02  1.8156E-01 -2.6323E-02  1.9642E-01
             1.3308E-01
 GRADIENT:  -2.5066E+01 -4.1437E+00 -2.7111E+00  3.5292E+00 -7.2185E+00 -5.5840E+00  1.9375E+00  1.3961E+00 -2.2222E+00  2.1520E+00
             4.1969E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2454.51966415129        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1557
 NPARAMETR:  9.6515E-01  1.0998E+00  1.1199E+00  9.6237E-01  1.0868E+00  1.0406E+00  9.5883E-01  1.0846E+00  9.0572E-01  1.0879E+00
             1.0328E+00
 PARAMETER:  6.4531E-02  1.9516E-01  2.1326E-01  6.1643E-02  1.8328E-01  1.3983E-01  5.7957E-02  1.8124E-01  9.7755E-04  1.8421E-01
             1.3225E-01
 GRADIENT:   1.1494E+01 -5.4811E+00 -2.1879E+00  3.3725E+00 -1.2032E+00  3.0645E+00 -7.6305E-01  6.7188E-01  1.4602E-01 -2.5014E+00
             3.2466E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2454.58426023483        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1733            RESET HESSIAN, TYPE II
 NPARAMETR:  9.6023E-01  1.1051E+00  1.1199E+00  9.6237E-01  1.0869E+00  1.0333E+00  9.6517E-01  1.0846E+00  9.0322E-01  1.1007E+00
             1.0328E+00
 PARAMETER:  5.9413E-02  1.9998E-01  2.1325E-01  6.1647E-02  1.8329E-01  1.3273E-01  6.4553E-02  1.8125E-01 -1.7858E-03  1.9596E-01
             1.3225E-01
 GRADIENT:   4.0286E+02  1.0876E+02  5.6757E+00  7.8559E+01  3.2411E+01  7.0039E+01  4.4638E+00  1.5595E+00  5.4185E+00  4.4148E+00
             4.7619E+00

0ITERATION NO.:   53    OBJECTIVE VALUE:  -2454.58442306020        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     1828
 NPARAMETR:  9.6023E-01  1.1052E+00  1.1199E+00  9.6237E-01  1.0869E+00  1.0326E+00  9.6366E-01  1.0846E+00  9.0348E-01  1.1016E+00
             1.0328E+00
 PARAMETER:  5.9414E-02  2.0006E-01  2.1325E-01  6.1647E-02  1.8329E-01  1.3209E-01  6.2979E-02  1.8125E-01 -1.4958E-03  1.9675E-01
             1.3225E-01
 GRADIENT:   1.0053E+00 -7.1968E-03 -1.7515E+00  8.0835E+00 -4.9444E+00  3.7004E-02  4.0451E-02  8.7396E-01 -8.4753E-02 -3.9272E-02
             3.1927E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1828
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.3182E-04 -1.8835E-02 -2.6310E-02  7.3113E-03 -2.2514E-02
 SE:             2.9895E-02  2.0592E-02  1.4330E-02  2.5116E-02  2.4297E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8314E-01  3.6037E-01  6.6366E-02  7.7097E-01  3.5413E-01

 ETASHRINKSD(%)  1.0000E-10  3.1014E+01  5.1991E+01  1.5859E+01  1.8600E+01
 ETASHRINKVR(%)  1.0000E-10  5.2409E+01  7.6951E+01  2.9202E+01  3.3741E+01
 EBVSHRINKSD(%)  3.2876E-01  3.0865E+01  5.5346E+01  1.7452E+01  1.6617E+01
 EBVSHRINKVR(%)  6.5643E-01  5.2203E+01  8.0060E+01  3.1859E+01  3.0473E+01
 RELATIVEINF(%)  9.9325E+01  1.1659E+01  9.9078E+00  1.9187E+01  2.2598E+01
 EPSSHRINKSD(%)  3.0087E+01
 EPSSHRINKVR(%)  5.1121E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2454.5844230601992     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1351.8581832145921     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.49
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2454.584       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.60E-01  1.11E+00  1.12E+00  9.62E-01  1.09E+00  1.03E+00  9.64E-01  1.08E+00  9.03E-01  1.10E+00  1.03E+00
 


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
+        7.45E+07
 
 TH 2
+        3.88E+03  4.76E+02
 
 TH 3
+       -2.17E+04 -6.51E+06  1.42E+02
 
 TH 4
+       -3.53E+02 -1.61E+07 -2.17E+04  3.71E+07
 
 TH 5
+        1.32E+02  7.80E+06 -7.22E+06  1.62E+02  8.65E+06
 
 TH 6
+       -3.18E+02 -1.57E+00 -1.27E+02 -1.47E+00 -1.54E+02  1.85E+02
 
 TH 7
+       -7.91E+02  1.11E+01 -3.25E+02 -3.70E+07 -1.79E+07 -3.18E+02  3.70E+07
 
 TH 8
+        1.82E+07  7.92E+06 -2.23E+01  1.82E+07  8.79E+06 -7.60E-01  3.86E+00  8.94E+06
 
 TH 9
+        1.90E+04 -1.02E+01  7.63E+03  3.56E+01  9.16E+03 -3.38E+02  1.89E+04  1.24E+02  4.21E+07
 
 TH10
+        1.65E+07 -7.16E+06 -6.64E+06 -1.65E+07  2.60E+01 -1.41E+02 -3.39E+02  1.61E+07  8.41E+03  7.31E+06
 
 TH11
+       -2.62E+07 -1.88E+01 -1.20E+01 -2.61E+07 -2.52E+07  1.47E+00  4.40E-01 -1.28E+07  2.31E+01 -1.16E+07  1.84E+07
 
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
 #CPUT: Total CPU Time in Seconds,       44.681
Stop Time:
Thu Sep 30 07:49:35 CDT 2021
