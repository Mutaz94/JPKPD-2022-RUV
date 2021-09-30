Wed Sep 29 23:15:58 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat33.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m33.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -529.677150573041        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0721E+02  7.4813E+01  7.7089E+01  6.0884E+01  8.7470E+01  4.3021E+01 -2.9385E+01 -1.8722E+01  1.9734E+01 -6.0785E+01
            -2.9683E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1608.96440028700        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.4223E-01  1.0594E+00  1.0472E+00  9.9650E-01  1.0725E+00  8.3685E-01  9.9931E-01  9.4204E-01  8.1765E-01  8.6614E-01
             2.9839E+00
 PARAMETER:  4.0491E-02  1.5774E-01  1.4612E-01  9.6492E-02  1.6996E-01 -7.8114E-02  9.9308E-02  4.0296E-02 -1.0132E-01 -4.3713E-02
             1.1932E+00
 GRADIENT:  -2.9009E+01 -9.0982E+00 -1.6193E+01  6.3875E+00  7.4921E+00 -3.4388E+01 -3.9926E+00  4.4533E+00  4.2798E+00  1.4144E+01
             4.6671E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1616.71568060606        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.4011E-01  1.0030E+00  9.5904E-01  1.0211E+00  9.9960E-01  9.3926E-01  1.1859E+00  5.4893E-01  7.3648E-01  4.1960E-01
             2.9927E+00
 PARAMETER:  3.8246E-02  1.0298E-01  5.8177E-02  1.2088E-01  9.9599E-02  3.7335E-02  2.7051E-01 -4.9979E-01 -2.0587E-01 -7.6846E-01
             1.1962E+00
 GRADIENT:  -1.8243E+01 -9.5634E+00 -8.8932E+00  1.8005E+00  1.0360E+01  8.9246E+00  4.7419E-01  1.1695E+00  5.0880E+00  2.8102E+00
             4.0163E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1618.05486835643        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.4886E-01  1.0471E+00  8.2263E-01  9.8091E-01  9.3348E-01  9.1881E-01  1.2077E+00  4.2137E-01  6.8602E-01  3.1648E-01
             2.8836E+00
 PARAMETER:  4.7507E-02  1.4606E-01 -9.5244E-02  8.0720E-02  3.1169E-02  1.5320E-02  2.8872E-01 -7.6425E-01 -2.7684E-01 -1.0505E+00
             1.1590E+00
 GRADIENT:   5.7308E+00 -1.1632E+00 -8.4482E-01 -1.6603E+00  1.0693E+00  1.3203E+00  8.5500E-01  4.0683E-01  4.7453E-01  1.3077E+00
            -1.6933E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1618.12254102403        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.4799E-01  1.0440E+00  7.0530E-01  9.7200E-01  8.5833E-01  9.1856E-01  1.2141E+00  3.1229E-01  6.7976E-01  2.1193E-01
             2.8807E+00
 PARAMETER:  4.6588E-02  1.4308E-01 -2.4914E-01  7.1604E-02 -5.2769E-02  1.5052E-02  2.9404E-01 -1.0638E+00 -2.8602E-01 -1.4515E+00
             1.1580E+00
 GRADIENT:   4.6606E-01 -4.1745E-01 -1.3509E+00  8.8024E-01  1.9054E+00  5.9249E-01  2.4552E-01  2.0159E-01 -9.7767E-02  5.9439E-01
            -3.4768E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1618.13404096236        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  9.4824E-01  1.0414E+00  6.7073E-01  9.6913E-01  8.3351E-01  9.1797E-01  1.2140E+00  2.5658E-01  6.8317E-01  1.6413E-01
             2.8796E+00
 PARAMETER:  4.6849E-02  1.4052E-01 -2.9938E-01  6.8645E-02 -8.2105E-02  1.4411E-02  2.9395E-01 -1.2603E+00 -2.8101E-01 -1.7071E+00
             1.1577E+00
 GRADIENT:   2.1973E-01 -4.7428E-01 -6.3188E-01  2.2567E-01  7.2916E-01  2.1811E-01 -9.1891E-02  1.1204E-01  4.1270E-02  3.4185E-01
            -6.7664E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1618.26431562198        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  9.4825E-01  1.0485E+00  6.5729E-01  9.6230E-01  8.2878E-01  9.1644E-01  1.2045E+00  8.2547E-02  6.8968E-01  3.6029E-02
             2.8905E+00
 PARAMETER:  4.6858E-02  1.4736E-01 -3.1962E-01  6.1566E-02 -8.7801E-02  1.2739E-02  2.8608E-01 -2.3944E+00 -2.7152E-01 -3.2234E+00
             1.1614E+00
 GRADIENT:  -5.1097E-01 -8.2630E-01  9.6991E-01 -2.1102E+00 -1.6935E+00 -3.7467E-01 -6.7030E-01  6.6212E-03  5.3413E-01  1.5488E-02
             1.5092E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1618.47619112467        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      519
 NPARAMETR:  9.4782E-01  1.0546E+00  7.1942E-01  9.6827E-01  8.7879E-01  9.1805E-01  1.2212E+00  1.1317E-02  6.6536E-01  1.0000E-02
             2.8969E+00
 PARAMETER:  4.6414E-02  1.5315E-01 -2.2931E-01  6.7753E-02 -2.9208E-02  1.4491E-02  2.9984E-01 -4.3814E+00 -3.0742E-01 -5.7664E+00
             1.1636E+00
 GRADIENT:   8.9449E-02  5.1360E-01 -2.4697E+00  3.6511E+00  4.1592E+00  4.5325E-01  8.8010E-01  1.7645E-04 -7.8854E-01  0.0000E+00
            -1.9830E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1619.54132614710        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      698
 NPARAMETR:  9.6340E-01  1.1340E+00  1.2006E+00  9.6049E-01  1.1871E+00  9.2347E-01  1.1036E+00  1.3390E-02  7.2040E-01  1.0000E-02
             2.9759E+00
 PARAMETER:  6.2712E-02  2.2573E-01  2.8281E-01  5.9689E-02  2.7149E-01  2.0382E-02  1.9860E-01 -4.2133E+00 -2.2795E-01 -6.0692E+00
             1.1906E+00
 GRADIENT:   1.0155E+00 -8.0412E-01 -1.4105E+00  1.1814E+00  2.2624E+00  7.2744E-01  2.7108E-01  1.2607E-04  3.9917E-01  0.0000E+00
            -3.2786E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1619.62595838193        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      873
 NPARAMETR:  9.6314E-01  1.1910E+00  1.4772E+00  9.3371E-01  1.3126E+00  9.2207E-01  1.0348E+00  1.3376E-01  7.4119E-01  2.1401E-02
             2.9955E+00
 PARAMETER:  6.2439E-02  2.7480E-01  4.9012E-01  3.1410E-02  3.7200E-01  1.8868E-02  1.3416E-01 -1.9117E+00 -1.9950E-01 -3.7443E+00
             1.1971E+00
 GRADIENT:   1.2575E-01 -3.5404E-01 -5.6472E-01  2.0712E-01  1.2490E+00  3.2092E-01  1.8085E-01  9.4316E-03 -2.8030E-01  2.6444E-03
             4.4450E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1619.66078365773        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1049
 NPARAMETR:  9.6354E-01  1.2663E+00  1.8019E+00  8.9312E-01  1.4188E+00  9.2123E-01  9.4545E-01  6.5691E-02  8.0057E-01  1.0000E-02
             3.0005E+00
 PARAMETER:  6.2859E-02  3.3613E-01  6.8886E-01 -1.3033E-02  4.4982E-01  1.7958E-02  4.3905E-02 -2.6228E+00 -1.2243E-01 -5.4504E+00
             1.1988E+00
 GRADIENT:  -8.6702E-02  5.7167E-01 -1.0763E-01  7.2279E-01  2.3095E-01 -1.3574E-01 -2.3352E-02  1.4069E-03  3.4034E-03  0.0000E+00
             1.4233E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1619.66156056480        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1224
 NPARAMETR:  9.6362E-01  1.2736E+00  1.8593E+00  8.8876E-01  1.4310E+00  9.2159E-01  9.3783E-01  3.9398E-02  8.0556E-01  1.0000E-02
             3.0009E+00
 PARAMETER:  6.2937E-02  3.4181E-01  7.2018E-01 -1.7924E-02  4.5835E-01  1.8343E-02  3.5812E-02 -3.1340E+00 -1.1622E-01 -6.2454E+00
             1.1989E+00
 GRADIENT:  -1.2592E-02 -2.5982E-02  2.8453E-03 -7.4751E-02 -1.4138E-02 -7.1508E-04  7.1285E-03  4.6255E-04  4.6162E-04  0.0000E+00
            -6.2121E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1619.66415465522        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1406
 NPARAMETR:  9.6342E-01  1.2213E+00  2.0207E+00  9.2546E-01  1.4383E+00  9.2148E-01  9.5469E-01  1.0000E-02  7.9647E-01  1.0000E-02
             3.0015E+00
 PARAMETER:  6.2730E-02  2.9993E-01  8.0342E-01  2.2537E-02  4.6346E-01  1.8223E-02  5.3635E-02 -9.0166E+00 -1.2756E-01 -1.3766E+01
             1.1991E+00
 GRADIENT:  -4.3892E-02  2.8991E-01  9.4430E-03  2.4451E-01 -2.7278E-02  7.0905E-03  9.7204E-03  0.0000E+00 -1.4483E-02  0.0000E+00
            -4.9945E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1619.66801751047        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1581
 NPARAMETR:  9.6302E-01  1.1298E+00  2.1883E+00  9.8820E-01  1.4323E+00  9.2110E-01  9.9612E-01  1.0000E-02  7.7635E-01  1.0000E-02
             3.0022E+00
 PARAMETER:  6.2320E-02  2.2206E-01  8.8311E-01  8.8132E-02  4.5928E-01  1.7810E-02  9.6109E-02 -1.5436E+01 -1.5315E-01 -2.1759E+01
             1.1994E+00
 GRADIENT:  -7.5174E-02  1.0072E+00 -3.1660E-02  1.2355E+00  1.3280E-03 -3.6593E-02  5.2775E-02  0.0000E+00  4.5304E-03  0.0000E+00
             1.1848E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1619.67227031812        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1758
 NPARAMETR:  9.6265E-01  1.0335E+00  2.4203E+00  1.0541E+00  1.4332E+00  9.2066E-01  1.0372E+00  1.0000E-02  7.5971E-01  1.0000E-02
             3.0030E+00
 PARAMETER:  6.1930E-02  1.3294E-01  9.8387E-01  1.5268E-01  4.5989E-01  1.7337E-02  1.3651E-01 -2.0364E+01 -1.7482E-01 -2.7876E+01
             1.1996E+00
 GRADIENT:  -5.0473E-02  6.0460E-01 -6.0119E-02  9.2336E-01  1.6572E-01 -9.2174E-02 -4.6474E-02  0.0000E+00 -5.5199E-02  0.0000E+00
             2.1048E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1619.67625495791        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1936
 NPARAMETR:  9.6212E-01  9.2749E-01  2.6399E+00  1.1263E+00  1.4272E+00  9.2061E-01  1.0966E+00  1.0000E-02  7.4059E-01  1.0000E-02
             3.0020E+00
 PARAMETER:  6.1383E-02  2.4723E-02  1.0707E+00  2.1895E-01  4.5568E-01  1.7276E-02  1.9223E-01 -2.2366E+01 -2.0031E-01 -3.0144E+01
             1.1993E+00
 GRADIENT:  -2.1945E-01  7.0963E-01  1.5032E-02  1.1281E+00 -2.0258E-01  4.2386E-03 -1.5847E-02  0.0000E+00 -3.5328E-02  0.0000E+00
            -2.4198E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1619.67842335252        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2114
 NPARAMETR:  9.6184E-01  8.4539E-01  2.7918E+00  1.1821E+00  1.4215E+00  9.2021E-01  1.1576E+00  1.0000E-02  7.2368E-01  1.0000E-02
             3.0026E+00
 PARAMETER:  6.1095E-02 -6.7961E-02  1.1267E+00  2.6726E-01  4.5172E-01  1.6848E-02  2.4637E-01 -2.1160E+01 -2.2341E-01 -2.8325E+01
             1.1995E+00
 GRADIENT:  -7.7037E-03  7.3453E-01 -9.4820E-03  1.3065E+00 -8.5396E-02 -4.3724E-02  3.6943E-02  0.0000E+00 -9.4063E-02  0.0000E+00
            -5.2880E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1619.68055316708        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2292
 NPARAMETR:  9.6143E-01  7.6942E-01  2.9188E+00  1.2331E+00  1.4145E+00  9.2004E-01  1.2113E+00  1.0000E-02  7.1474E-01  1.0000E-02
             3.0027E+00
 PARAMETER:  6.0661E-02 -1.6212E-01  1.1712E+00  3.0955E-01  4.4679E-01  1.6660E-02  2.9166E-01 -1.7961E+01 -2.3584E-01 -2.3913E+01
             1.1995E+00
 GRADIENT:  -1.0637E-01  4.1518E-01 -7.5513E-03  7.5184E-01 -6.3449E-02 -3.1494E-03  1.0969E-01  0.0000E+00  1.8207E-01  0.0000E+00
             9.8637E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1619.68072731701        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2467
 NPARAMETR:  9.6128E-01  7.3640E-01  2.9864E+00  1.2561E+00  1.4125E+00  9.1994E-01  1.2371E+00  1.0000E-02  7.0932E-01  1.0000E-02
             3.0027E+00
 PARAMETER:  6.0506E-02 -2.0598E-01  1.1941E+00  3.2797E-01  4.4540E-01  1.6550E-02  3.1275E-01 -1.6219E+01 -2.4345E-01 -2.1552E+01
             1.1995E+00
 GRADIENT:  -1.8974E-01  7.5619E-01 -5.9914E-03  1.6635E+00 -2.3049E-01 -1.8533E-02  4.2481E-02  0.0000E+00  6.0901E-02  0.0000E+00
            -2.0157E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1619.68088104449        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2642
 NPARAMETR:  9.6119E-01  7.1110E-01  3.0510E+00  1.2736E+00  1.4124E+00  9.1987E-01  1.2564E+00  1.0000E-02  7.0537E-01  1.0000E-02
             3.0026E+00
 PARAMETER:  6.0413E-02 -2.4095E-01  1.2155E+00  3.4181E-01  4.4528E-01  1.6475E-02  3.2824E-01 -1.4700E+01 -2.4904E-01 -1.9518E+01
             1.1995E+00
 GRADIENT:  -2.0693E-01  8.1614E-01 -3.2221E-03  1.9237E+00 -2.9958E-01 -2.5761E-02 -2.0890E-02  0.0000E+00 -5.5615E-02  0.0000E+00
            -1.0991E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1619.68118352861        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2817
 NPARAMETR:  9.6112E-01  6.8324E-01  3.1380E+00  1.2928E+00  1.4138E+00  9.1981E-01  1.2770E+00  1.0000E-02  7.0161E-01  1.0000E-02
             3.0027E+00
 PARAMETER:  6.0343E-02 -2.8091E-01  1.2436E+00  3.5681E-01  4.4627E-01  1.6414E-02  3.4448E-01 -1.2847E+01 -2.5437E-01 -1.7057E+01
             1.1995E+00
 GRADIENT:  -1.7799E-01  7.1467E-01  1.7614E-04  1.7746E+00 -3.0033E-01 -2.6858E-02 -7.2790E-02  0.0000E+00 -1.3805E-01  0.0000E+00
            -1.6375E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1619.68134998440        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2993
 NPARAMETR:  9.6105E-01  6.5453E-01  3.2300E+00  1.3126E+00  1.4153E+00  9.1976E-01  1.2995E+00  1.0000E-02  6.9796E-01  1.0000E-02
             3.0028E+00
 PARAMETER:  6.0275E-02 -3.2383E-01  1.2725E+00  3.7204E-01  4.4736E-01  1.6358E-02  3.6196E-01 -1.0585E+01 -2.5959E-01 -1.4064E+01
             1.1995E+00
 GRADIENT:  -1.4291E-01  5.9830E-01  2.8727E-03  1.5547E+00 -2.8182E-01 -2.5324E-02 -1.0569E-01  0.0000E+00 -1.8360E-01  0.0000E+00
            -1.8959E-01

0ITERATION NO.:  109    OBJECTIVE VALUE:  -1619.68219207893        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     3128
 NPARAMETR:  9.6109E-01  6.4556E-01  3.2545E+00  1.3179E+00  1.4160E+00  9.1978E-01  1.3119E+00  1.0000E-02  6.9752E-01  1.0000E-02
             3.0030E+00
 PARAMETER:  6.0309E-02 -3.3764E-01  1.2800E+00  3.7604E-01  4.4781E-01  1.6384E-02  3.7149E-01 -9.8236E+00 -2.6022E-01 -1.3058E+01
             1.1996E+00
 GRADIENT:   1.3858E-01 -3.7804E-02 -2.7140E-04 -3.6478E-01  6.2481E-02  2.1128E-02  2.1859E-02  0.0000E+00  1.1995E-01  0.0000E+00
             8.1138E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3128
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3218E-03  5.9397E-05  1.9015E-05 -8.8396E-03 -2.8079E-05
 SE:             2.8995E-02  1.4354E-02  2.8177E-05  2.2426E-02  1.3523E-04
 N:                     100         100         100         100         100

 P VAL.:         9.6364E-01  9.9670E-01  4.9978E-01  6.9345E-01  8.3551E-01

 ETASHRINKSD(%)  2.8637E+00  5.1911E+01  9.9906E+01  2.4872E+01  9.9547E+01
 ETASHRINKVR(%)  5.6455E+00  7.6874E+01  1.0000E+02  4.3557E+01  9.9998E+01
 EBVSHRINKSD(%)  2.8516E+00  5.2498E+01  9.9897E+01  2.4604E+01  9.9526E+01
 EBVSHRINKVR(%)  5.6219E+00  7.7436E+01  1.0000E+02  4.3154E+01  9.9998E+01
 RELATIVEINF(%)  8.9339E+01  3.9133E-02  2.7285E-06  9.3864E-02  1.2063E-04
 EPSSHRINKSD(%)  1.8081E+01
 EPSSHRINKVR(%)  3.2893E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1619.6821920789332     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -700.74365887426052     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    45.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1619.682       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.61E-01  6.46E-01  3.25E+00  1.32E+00  1.42E+00  9.20E-01  1.31E+00  1.00E-02  6.98E-01  1.00E-02  3.00E+00
 


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
+        1.34E+03
 
 TH 2
+       -4.92E+01  3.47E+02
 
 TH 3
+       -1.97E+00 -1.04E-01  1.02E+00
 
 TH 4
+       -5.75E+01  5.00E+02 -4.33E+00  7.41E+02
 
 TH 5
+        9.28E+00 -6.15E+01 -1.01E+01 -5.05E+01  1.17E+02
 
 TH 6
+        4.10E+00 -1.05E+01  4.56E-02 -1.50E+01 -7.89E-01  2.09E+02
 
 TH 7
+        1.08E+00  3.60E+00  2.55E-01 -6.19E-01  6.76E-01  3.08E-01  6.91E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.18E+00 -1.94E+00  7.12E-01 -1.32E+01  2.71E+00  9.92E-01  2.42E+01  0.00E+00  1.30E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.88E+01 -9.58E+00  1.58E-01 -1.61E+01  2.51E-01  3.27E+00  2.52E+00  0.00E+00  1.10E+01  0.00E+00  6.77E+01
 
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
 #CPUT: Total CPU Time in Seconds,       53.085
Stop Time:
Wed Sep 29 23:16:52 CDT 2021
