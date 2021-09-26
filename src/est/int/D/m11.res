Sat Sep 25 05:22:21 CDT 2021
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
$DATA ../../../../data/int/D/dat11.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m11.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3439.94359367313        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.7935E+01 -1.6357E+02 -7.9512E+01 -1.5355E+02  1.2821E+02 -2.9804E+02 -2.5308E+02 -1.8587E+01 -3.9862E+02 -5.7512E+01
             3.9630E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3638.85310437027        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.7864E-01  1.1968E+00  1.3006E+00  6.5795E-01  1.6096E+00  1.2188E+00  2.6670E+00  1.2391E+00  1.1676E+00  1.0661E+00
             1.0926E+00
 PARAMETER:  7.8404E-02  2.7969E-01  3.6286E-01 -3.1863E-01  5.7598E-01  2.9786E-01  1.0809E+00  3.1441E-01  2.5499E-01  1.6404E-01
             1.8852E-01
 GRADIENT:  -7.4437E+00 -6.8218E+01 -3.3199E-01 -7.2680E+01  1.9330E+02 -1.1823E+02  1.0553E+02 -2.1410E+01 -1.6158E+01 -3.1480E+01
             1.0377E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3651.93164369224        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      167
 NPARAMETR:  1.0328E+00  1.2890E+00  1.3946E+00  6.0575E-01  1.6797E+00  1.3949E+00  2.6020E+00  2.0705E+00  1.0801E+00  1.2303E+00
             1.1043E+00
 PARAMETER:  1.3232E-01  3.5387E-01  4.3263E-01 -4.0129E-01  6.1864E-01  4.3281E-01  1.0563E+00  8.2780E-01  1.7709E-01  3.0728E-01
             1.9919E-01
 GRADIENT:   7.3734E+01 -5.4358E+01  2.0980E+00 -1.0028E+02  1.2575E+02 -2.6422E+01  9.7791E+01  1.4733E+00 -1.1825E+01 -1.5622E+00
             1.2125E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3683.54538332886        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      241
 NPARAMETR:  9.9371E-01  1.2990E+00  1.3345E+00  7.4048E-01  1.4776E+00  1.4561E+00  1.9610E+00  1.3570E+00  1.1893E+00  1.1503E+00
             1.0436E+00
 PARAMETER:  9.3686E-02  3.6159E-01  3.8855E-01 -2.0045E-01  4.9044E-01  4.7575E-01  7.7345E-01  4.0526E-01  2.7339E-01  2.4005E-01
             1.4267E-01
 GRADIENT:   2.4346E+01 -6.2768E+01 -5.7880E+00 -5.1900E+01  1.1132E+02  6.8171E+00 -2.4412E+01 -1.3527E+01 -2.6965E+01  7.0744E+00
             4.9230E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3687.04630860526        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:      346
 NPARAMETR:  9.8946E-01  1.2901E+00  1.3428E+00  7.6184E-01  1.4415E+00  1.4513E+00  1.9529E+00  1.3651E+00  1.2033E+00  1.1179E+00
             1.0381E+00
 PARAMETER:  8.9405E-02  3.5475E-01  3.9474E-01 -1.7202E-01  4.6566E-01  4.7249E-01  7.6934E-01  4.1121E-01  2.8503E-01  2.1143E-01
             1.3739E-01
 GRADIENT:  -1.7082E+01 -9.8343E+01 -8.8902E+00 -5.4966E+01  8.1576E+01 -4.4365E+01 -7.3695E+01 -1.1855E+01 -2.8658E+01  4.0135E+00
             4.2750E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3687.69679168027        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  9.8946E-01  1.2901E+00  1.3428E+00  7.6184E-01  1.4414E+00  1.4637E+00  1.9529E+00  1.3651E+00  1.2164E+00  1.1179E+00
             1.0381E+00
 PARAMETER:  8.9407E-02  3.5475E-01  3.9474E-01 -1.7202E-01  4.6565E-01  4.8100E-01  7.6933E-01  4.1120E-01  2.9591E-01  2.1143E-01
             1.3739E-01
 GRADIENT:  -1.6655E+01 -9.9353E+01 -8.0741E+00 -5.1248E+01  7.9680E+01 -4.0225E+01 -7.4283E+01 -1.0531E+01 -2.3217E+01  2.9215E+00
             4.1078E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3689.50587073947        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      714
 NPARAMETR:  9.8945E-01  1.2902E+00  1.3427E+00  7.6185E-01  1.4415E+00  1.6120E+00  1.9531E+00  1.3651E+00  1.2165E+00  1.1175E+00
             1.0381E+00
 PARAMETER:  8.9395E-02  3.5479E-01  3.9470E-01 -1.7200E-01  4.6570E-01  5.7746E-01  7.6942E-01  4.1125E-01  2.9597E-01  2.1114E-01
             1.3738E-01
 GRADIENT:  -1.3734E+01 -9.9315E+01 -8.1489E+00 -5.1261E+01  7.9803E+01  1.8059E+00 -7.4276E+01 -1.0550E+01 -2.3315E+01  2.9063E+00
             4.1209E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3689.57845565786        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      897            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8952E-01  1.2907E+00  1.3427E+00  7.6187E-01  1.4414E+00  1.6047E+00  1.9535E+00  1.3651E+00  1.2165E+00  1.1111E+00
             1.0381E+00
 PARAMETER:  8.9462E-02  3.5515E-01  3.9470E-01 -1.7197E-01  4.6565E-01  5.7292E-01  7.6961E-01  4.1125E-01  2.9597E-01  2.0531E-01
             1.3735E-01
 GRADIENT:   2.3312E+01 -6.0253E+01 -5.5049E+00 -3.7313E+01  9.4480E+01  6.0513E+01 -2.7616E+01 -8.7788E+00 -1.5430E+01  1.0151E+00
             3.7949E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3692.41836249452        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1027
 NPARAMETR:  9.9536E-01  1.2906E+00  1.3428E+00  7.6596E-01  1.4261E+00  1.5900E+00  2.0085E+00  1.3651E+00  1.2188E+00  1.1104E+00
             1.0325E+00
 PARAMETER:  9.5352E-02  3.5511E-01  3.9475E-01 -1.6662E-01  4.5492E-01  5.6376E-01  7.9740E-01  4.1125E-01  2.9787E-01  2.0468E-01
             1.3195E-01
 GRADIENT:   2.8691E+01 -4.6220E+01 -1.0952E+01 -5.4465E+01  9.8106E+01  5.6246E+01 -7.0734E+00 -1.3686E+01 -3.3324E+01  7.7975E+00
             4.0924E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3692.95018995440        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     1227
 NPARAMETR:  9.9538E-01  1.2906E+00  1.3428E+00  7.6609E-01  1.4259E+00  1.5853E+00  2.0090E+00  1.3651E+00  1.2359E+00  1.1082E+00
             1.0324E+00
 PARAMETER:  9.5373E-02  3.5511E-01  3.9475E-01 -1.6646E-01  4.5479E-01  5.6078E-01  7.9766E-01  4.1125E-01  3.1182E-01  2.0277E-01
             1.3191E-01
 GRADIENT:   2.8658E+01 -4.5877E+01 -1.1560E+01 -5.3692E+01  9.9743E+01  5.4714E+01 -5.9719E+00 -1.3393E+01 -3.1115E+01  7.4770E+00
             4.1925E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3693.02750554018        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     1358
 NPARAMETR:  9.9533E-01  1.2906E+00  1.3425E+00  7.6615E-01  1.4262E+00  1.6070E+00  2.0098E+00  1.3651E+00  1.2373E+00  1.1078E+00
             1.0324E+00
 PARAMETER:  9.5322E-02  3.5511E-01  3.9455E-01 -1.6637E-01  4.5502E-01  5.7437E-01  7.9806E-01  4.1125E-01  3.1292E-01  2.0235E-01
             1.3185E-01
 GRADIENT:   2.8860E+01 -4.5808E+01 -1.1732E+01 -5.3418E+01  1.0028E+02  6.1650E+01 -5.6967E+00 -1.3371E+01 -3.0880E+01  7.3440E+00
             4.1831E+01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -3693.02750554018        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1422
 NPARAMETR:  9.9533E-01  1.2906E+00  1.3425E+00  7.6615E-01  1.4262E+00  1.6070E+00  2.0098E+00  1.3651E+00  1.2373E+00  1.1078E+00
             1.0324E+00
 PARAMETER:  9.5322E-02  3.5511E-01  3.9455E-01 -1.6637E-01  4.5502E-01  5.7437E-01  7.9806E-01  4.1125E-01  3.1292E-01  2.0235E-01
             1.3185E-01
 GRADIENT:  -4.2481E+05  2.3917E+05 -1.0767E+05  5.1059E+05  9.3441E+04 -1.4792E+05  1.0640E+05 -2.6587E+02 -1.3580E+05 -4.1985E+05
            -6.4437E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1422
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.4691E-03  6.2109E-02 -1.3835E-02  5.9435E-02 -7.3082E-02
 SE:             2.9909E-02  3.1043E-02  1.7029E-02  2.6695E-02  2.3129E-02
 N:                     100         100         100         100         100

 P VAL.:         8.2876E-01  4.5422E-02  4.1653E-01  2.5984E-02  1.5790E-03

 ETASHRINKSD(%)  1.0000E-10  1.0000E-10  4.2952E+01  1.0568E+01  2.2516E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E-10  6.7455E+01  2.0020E+01  3.9962E+01
 EBVSHRINKSD(%)  1.0644E-01  7.1322E+00  5.2282E+01  2.4907E+01  1.7445E+01
 EBVSHRINKVR(%)  2.1278E-01  1.3756E+01  7.7230E+01  4.3610E+01  3.1846E+01
 RELATIVEINF(%)  9.9787E+01  6.4992E+01  1.8662E+01  3.6043E+01  4.3982E+01
 EPSSHRINKSD(%)  2.2745E+01
 EPSSHRINKVR(%)  4.0317E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3693.0275055401803     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2038.9381457717695     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.07
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3693.028       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  1.29E+00  1.34E+00  7.66E-01  1.43E+00  1.61E+00  2.01E+00  1.37E+00  1.24E+00  1.11E+00  1.03E+00
 


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
+        2.14E+09
 
 TH 2
+        1.91E+02  1.01E+08
 
 TH 3
+       -1.65E+02  2.26E+04  1.51E+08
 
 TH 4
+        8.56E+04 -1.85E+04  1.60E+04  1.31E+09
 
 TH 5
+        7.70E+03 -1.72E+03  1.40E+03 -5.80E+03  5.04E+07
 
 TH 6
+       -9.93E+02  2.17E+02 -1.88E+02  7.80E+02  1.53E+02  2.49E+07
 
 TH 7
+        3.02E+03 -6.20E+02  5.66E+02 -2.34E+03 -4.68E+02  6.19E+01  8.26E+06
 
 TH 8
+       -3.80E+08  8.26E+07 -7.15E+07  2.97E+08  5.83E+07 -1.64E-01  2.36E+07  1.35E+08
 
 TH 9
+        9.85E+04 -2.14E+04  1.85E+04 -7.69E+04 -1.51E+04 -2.56E+02 -6.10E+03 -9.78E+07  1.42E+08
 
 TH10
+       -1.28E+03  2.79E+02 -2.29E+02  9.91E+02  1.39E+02 -4.43E+02  7.80E+01 -1.75E+00  4.38E+04  4.23E+08
 
 TH11
+        8.31E+04 -1.81E+04  1.56E+04 -6.48E+04 -1.28E+04 -7.30E+02  2.20E+03 -2.78E+08  7.21E+04 -9.07E+02  1.15E+09
 
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
 #CPUT: Total CPU Time in Seconds,       60.136
Stop Time:
Sat Sep 25 05:23:22 CDT 2021
