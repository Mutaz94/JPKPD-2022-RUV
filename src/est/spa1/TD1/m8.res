Thu Sep 30 01:08:30 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat8.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m8.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2114.83156230815        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7147E+02  2.2310E+01  1.3212E+01  6.4676E+01  1.2595E+01  3.4039E+01 -3.9614E-01 -6.4423E+00  1.5994E+01  1.8885E+00
             3.2550E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2116.90670021829        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.8029E-01  1.0804E+00  8.7343E-01  9.6520E-01  9.5483E-01  1.0683E+00  1.0669E+00  1.0927E+00  9.1439E-01  9.6729E-01
             9.6554E-01
 PARAMETER:  8.0089E-02  1.7731E-01 -3.5331E-02  6.4581E-02  5.3775E-02  1.6608E-01  1.6473E-01  1.8863E-01  1.0498E-02  6.6747E-02
             6.4934E-02
 GRADIENT:  -7.4367E+00  1.0113E+01 -6.4717E-01  1.2870E+01 -1.7942E+01  1.5771E+01 -2.2719E+00  3.0767E+00 -6.7049E+00  1.1800E+01
             6.3195E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2117.41810225233        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.8066E-01  1.1417E+00  8.1131E-01  9.2903E-01  9.4360E-01  1.0335E+00  1.1007E+00  1.1493E+00  9.1979E-01  8.5496E-01
             9.6148E-01
 PARAMETER:  8.0467E-02  2.3256E-01 -1.0911E-01  2.6388E-02  4.1951E-02  1.3297E-01  1.9595E-01  2.3917E-01  1.6391E-02 -5.6704E-02
             6.0716E-02
 GRADIENT:  -9.0809E+00  2.1882E+01 -1.6729E+00  1.7403E+01 -1.6982E+01  2.6124E+00  4.0092E+00  4.6393E+00 -6.5880E+00  3.4974E+00
             9.9176E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2118.64679391749        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.8607E-01  1.2311E+00  7.2961E-01  8.5762E-01  9.6399E-01  1.0252E+00  9.8033E-01  8.7217E-01  1.0385E+00  8.8698E-01
             9.5919E-01
 PARAMETER:  8.5974E-02  3.0789E-01 -2.1524E-01 -5.3591E-02  6.3325E-02  1.2487E-01  8.0135E-02 -3.6773E-02  1.3782E-01 -1.9932E-02
             5.8332E-02
 GRADIENT:   8.9785E-01  3.3255E+00  7.9956E-01  6.2236E+00 -1.6001E+00 -8.8792E-01 -5.8133E-01  2.5842E-01  3.1323E-02  1.8122E-01
            -7.5262E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2119.14729177039        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  9.8665E-01  1.4865E+00  5.4958E-01  6.8595E-01  1.0172E+00  1.0321E+00  8.6811E-01  4.8693E-01  1.2200E+00  9.3112E-01
             9.5773E-01
 PARAMETER:  8.6558E-02  4.9643E-01 -4.9859E-01 -2.7695E-01  1.1704E-01  1.3160E-01 -4.1436E-02 -6.1963E-01  2.9881E-01  2.8628E-02
             5.6806E-02
 GRADIENT:  -1.0173E+00 -5.1716E-02 -1.9683E+00  1.4656E+00 -2.8294E-01  9.1397E-01  1.7922E-01  6.4468E-01 -3.1865E-01  6.3029E-01
             6.8327E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2119.15227528574        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  9.8701E-01  1.5336E+00  5.2739E-01  6.5565E-01  1.0354E+00  1.0314E+00  8.5126E-01  4.1872E-01  1.2624E+00  9.4442E-01
             9.5723E-01
 PARAMETER:  8.6929E-02  5.2765E-01 -5.3981E-01 -3.2213E-01  1.3479E-01  1.3090E-01 -6.1035E-02 -7.7056E-01  3.3299E-01  4.2816E-02
             5.6289E-02
 GRADIENT:  -5.3164E-01  4.1640E-01 -1.5231E+00  1.1610E+00 -6.9277E-02  5.5184E-01  1.2300E-01  5.2829E-01  1.5660E-02  4.0384E-01
             4.2045E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2119.15267687129        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  9.8714E-01  1.5521E+00  5.1859E-01  6.4335E-01  1.0430E+00  1.0310E+00  8.4521E-01  3.8340E-01  1.2797E+00  9.5024E-01
             9.5703E-01
 PARAMETER:  8.7056E-02  5.3958E-01 -5.5665E-01 -3.4107E-01  1.4206E-01  1.3056E-01 -6.8168E-02 -8.5868E-01  3.4664E-01  4.8960E-02
             5.6077E-02
 GRADIENT:  -3.6368E-01 -2.8216E-02 -1.2639E+00  5.5439E-01  5.0720E-02  3.8812E-01  1.2346E-01  4.5893E-01  7.2419E-02  3.3833E-01
             4.5345E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2119.15294316274        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1230
 NPARAMETR:  9.8727E-01  1.5717E+00  5.0906E-01  6.3015E-01  1.0512E+00  1.0307E+00  8.3899E-01  3.3708E-01  1.2989E+00  9.5656E-01
             9.5680E-01
 PARAMETER:  8.7184E-02  5.5218E-01 -5.7519E-01 -3.6180E-01  1.4989E-01  1.3020E-01 -7.5562E-02 -9.8744E-01  3.6149E-01  5.5590E-02
             5.5842E-02
 GRADIENT:  -1.9197E-01 -5.2583E-01 -9.4900E-01 -9.3684E-02  1.9958E-01  2.1675E-01  1.1805E-01  3.6775E-01  1.1582E-01  2.5574E-01
             4.6036E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2119.15407628086        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1405
 NPARAMETR:  9.8736E-01  1.5863E+00  5.0185E-01  6.2034E-01  1.0572E+00  1.0304E+00  8.3455E-01  2.9408E-01  1.3135E+00  9.6126E-01
             9.5663E-01
 PARAMETER:  8.7278E-02  5.6140E-01 -5.8945E-01 -3.7749E-01  1.5565E-01  1.2992E-01 -8.0863E-02 -1.1239E+00  3.7269E-01  6.0490E-02
             5.5662E-02
 GRADIENT:  -6.4836E-02 -9.2739E-01 -7.0570E-01 -5.7934E-01  2.9728E-01  8.7223E-02  1.1725E-01  2.8817E-01  1.4360E-01  1.9519E-01
             4.7299E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2119.15736650359        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1580
 NPARAMETR:  9.8746E-01  1.6039E+00  4.9328E-01  6.0844E-01  1.0649E+00  1.0300E+00  8.2932E-01  2.2557E-01  1.3318E+00  9.6713E-01
             9.5640E-01
 PARAMETER:  8.7383E-02  5.7243E-01 -6.0668E-01 -3.9686E-01  1.6287E-01  1.2958E-01 -8.7147E-02 -1.3891E+00  3.8650E-01  6.6574E-02
             5.5423E-02
 GRADIENT:   8.1841E-02 -1.4633E+00 -3.1237E-01 -1.2051E+00  5.2964E-01 -6.4391E-02  8.8404E-02  1.7420E-01  1.5018E-01  7.8512E-02
             3.8304E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2119.18052969911        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1758
 NPARAMETR:  9.8755E-01  1.6179E+00  4.8502E-01  5.9884E-01  1.0699E+00  1.0297E+00  8.2573E-01  9.9471E-02  1.3469E+00  9.7117E-01
             9.5622E-01
 PARAMETER:  8.7474E-02  5.8111E-01 -6.2357E-01 -4.1276E-01  1.6755E-01  1.2928E-01 -9.1489E-02 -2.2079E+00  3.9779E-01  7.0744E-02
             5.5237E-02
 GRADIENT:   1.3116E-01 -1.8448E+00 -7.1205E-02 -1.6483E+00  4.3156E-01 -2.1762E-01  1.4536E-01  3.8804E-02  1.9334E-01  4.8834E-02
             5.0615E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2119.34424295567        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1940
 NPARAMETR:  9.8622E-01  1.5499E+00  5.0878E-01  6.4522E-01  1.0341E+00  1.0286E+00  8.4730E-01  1.0000E-02  1.2891E+00  9.4325E-01
             9.5732E-01
 PARAMETER:  8.6129E-02  5.3822E-01 -5.7574E-01 -3.3816E-01  1.3351E-01  1.2821E-01 -6.5705E-02 -6.4645E+00  3.5391E-01  4.1574E-02
             5.6377E-02
 GRADIENT:  -2.3913E+00  3.6547E+00  4.3481E-01  2.7835E+00 -9.8001E-01 -5.6955E-01 -3.1968E-01  0.0000E+00  8.6521E-01 -2.7245E-01
             7.0014E-02

0ITERATION NO.:   59    OBJECTIVE VALUE:  -2119.36170805032        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     2081
 NPARAMETR:  9.8847E-01  1.5413E+00  5.0850E-01  6.4256E-01  1.0345E+00  1.0306E+00  8.4912E-01  1.0000E-02  1.2842E+00  9.4541E-01
             9.5731E-01
 PARAMETER:  8.7972E-02  5.3703E-01 -5.7645E-01 -3.3991E-01  1.3397E-01  1.3015E-01 -6.2552E-02 -6.4645E+00  3.4826E-01  4.3102E-02
             5.6196E-02
 GRADIENT:  -3.2704E-01  3.3368E+00 -2.3646E-02  7.6797E-01  1.2739E-02  1.1577E-03  6.3651E-02  0.0000E+00 -1.2361E-01 -3.7697E-02
            -5.1656E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2081
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.8548E-04 -2.6035E-02 -2.9760E-04  2.0839E-02 -3.0062E-02
 SE:             2.9900E-02  2.3816E-02  1.3577E-04  2.4226E-02  2.2308E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8705E-01  2.7431E-01  2.8384E-02  3.8970E-01  1.7779E-01

 ETASHRINKSD(%)  1.0000E-10  2.0215E+01  9.9545E+01  1.8838E+01  2.5265E+01
 ETASHRINKVR(%)  1.0000E-10  3.6343E+01  9.9998E+01  3.4128E+01  4.4147E+01
 EBVSHRINKSD(%)  2.9146E-01  2.0005E+01  9.9604E+01  1.9760E+01  2.4005E+01
 EBVSHRINKVR(%)  5.8208E-01  3.6008E+01  9.9998E+01  3.5616E+01  4.2247E+01
 RELATIVEINF(%)  9.9361E+01  5.9884E+00  3.2466E-04  6.4136E+00  1.1276E+01
 EPSSHRINKSD(%)  3.4201E+01
 EPSSHRINKVR(%)  5.6705E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2119.3617080503213     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1200.4231748456486     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.75
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2119.362       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  1.55E+00  5.08E-01  6.44E-01  1.03E+00  1.03E+00  8.50E-01  1.00E-02  1.28E+00  9.45E-01  9.57E-01
 


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
+        1.07E+03
 
 TH 2
+       -3.97E+00  4.30E+02
 
 TH 3
+        7.52E+00  2.51E+02  8.03E+02
 
 TH 4
+       -4.42E+00  3.11E+02 -4.26E+02  1.06E+03
 
 TH 5
+        1.54E+00 -2.43E+02 -4.71E+02  3.40E+02  6.75E+02
 
 TH 6
+        7.71E-01 -7.13E-01  3.62E+00 -2.00E+00  3.23E-01  1.86E+02
 
 TH 7
+        7.58E-01  9.16E+00 -4.84E+01  1.60E+00 -4.71E+00  1.41E-01  1.21E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.00E-01 -2.03E+01 -2.94E+01  5.34E+01 -1.43E+00 -3.50E-01  1.70E+01  0.00E+00  5.59E+01
 
 TH10
+        1.27E+00 -1.66E+01 -4.89E+01 -7.30E+00 -6.49E+01  4.59E-01  2.14E+01  0.00E+00  9.66E+00  7.62E+01
 
 TH11
+       -7.25E+00 -1.43E+01 -1.16E+01 -9.39E+00  2.68E-01  1.22E+00  5.39E+00  0.00E+00  8.05E+00  1.70E+01  4.37E+02
 
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
 #CPUT: Total CPU Time in Seconds,       39.120
Stop Time:
Thu Sep 30 01:09:10 CDT 2021
