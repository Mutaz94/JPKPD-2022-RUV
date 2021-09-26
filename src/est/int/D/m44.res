Sat Sep 25 05:48:52 CDT 2021
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
$DATA ../../../../data/int/D/dat44.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m44.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25838.0192989121        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1927E+02  4.8657E+02 -1.2522E+00  3.8254E+02  7.3664E+00 -2.2020E+03 -1.2589E+03 -9.0520E+01 -1.6647E+03 -4.0448E+02
            -5.3116E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -726.663416726638        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.7326E-01  1.5261E+00  9.2942E-01  1.3536E+00  1.0333E+00  4.2033E+00  4.7771E+00  1.0250E+00  2.5171E+00  1.3217E+00
             1.2388E+01
 PARAMETER:  7.2899E-02  5.2273E-01  2.6800E-02  4.0280E-01  1.3276E-01  1.5359E+00  1.6638E+00  1.2473E-01  1.0231E+00  3.7891E-01
             2.6167E+00
 GRADIENT:  -2.4738E+01  1.0121E+01 -3.9543E+01  4.2936E+01 -6.3523E+00  1.3165E+02  8.7326E+01  3.0553E+00  7.1457E+01  2.5609E+01
             5.4040E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -788.152658875034        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  8.9353E-01  7.6559E-01  5.7129E+01  3.0961E+00  2.1653E+00  2.7754E+00  1.4364E+01  1.2539E+00  2.3173E+00  1.6209E+00
             1.1970E+01
 PARAMETER: -1.2575E-02 -1.6710E-01  4.1453E+00  1.2302E+00  8.7254E-01  1.1208E+00  2.7647E+00  3.2628E-01  9.4039E-01  5.8299E-01
             2.5824E+00
 GRADIENT:  -6.9558E+01  2.2250E+01  4.0771E-01  9.6312E+01 -2.5242E+01  5.5297E+01  2.9216E+01  5.2453E-02 -8.2324E-01  3.2613E+01
             4.9281E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -962.588132766294        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.3594E+00  2.2509E-01  4.5386E+01  1.7034E+00  2.4648E+00  2.5943E+00  5.4483E+00  1.7465E+00  2.5296E+00  7.5847E-01
             8.6081E+00
 PARAMETER:  4.0705E-01 -1.3913E+00  3.9152E+00  6.3261E-01  1.0021E+00  1.0533E+00  1.7953E+00  6.5762E-01  1.0280E+00 -1.7645E-01
             2.2527E+00
 GRADIENT:   8.3234E+01 -5.2934E+00 -6.7690E+00 -3.4309E+01  5.2802E+01  6.9404E+00  3.6857E+00  7.1871E-01  2.7904E+01  5.7069E+00
             1.5155E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -996.457745000200        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.7506E-01  2.0245E-01  1.1005E+02  1.8759E+00  2.4369E+00  2.4363E+00  8.1682E+00  8.9362E-01  2.1197E+00  5.6017E-01
             7.8789E+00
 PARAMETER:  7.4742E-02 -1.4973E+00  4.8009E+00  7.2908E-01  9.9072E-01  9.9049E-01  2.2003E+00 -1.2479E-02  8.5126E-01 -4.7952E-01
             2.1642E+00
 GRADIENT:  -1.7872E+01 -3.5136E+00 -3.2815E-01  4.1085E+01  2.3503E+01 -2.0380E+00 -1.1128E+00  9.3458E-03  4.4992E+01  1.5297E+00
            -1.9719E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1001.07611247655        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  1.0347E+00  1.9103E-01  8.6413E+01  1.7272E+00  2.2767E+00  2.4338E+00  7.0111E+00  9.0407E-01  1.8478E+00  5.7674E-01
             7.8818E+00
 PARAMETER:  1.3408E-01 -1.5553E+00  4.5591E+00  6.4653E-01  9.2271E-01  9.8945E-01  2.0475E+00 -8.5109E-04  7.1399E-01 -4.5037E-01
             2.1646E+00
 GRADIENT:   4.2926E+00 -4.1768E-01 -8.6273E-01 -7.8702E+00 -2.8958E+00  7.3468E-01  6.3951E+00 -5.3944E-03  1.1028E+00  5.1097E-01
            -1.9004E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1001.41963288849        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  1.0310E+00  1.8993E-01  8.9283E+01  1.7424E+00  2.2873E+00  2.4432E+00  7.0709E+00  9.1056E-01  1.8443E+00  5.7089E-01
             7.8894E+00
 PARAMETER:  1.3055E-01 -1.5611E+00  4.5918E+00  6.5524E-01  9.2738E-01  9.9332E-01  2.0560E+00  6.3001E-03  7.1209E-01 -4.6056E-01
             2.1655E+00
 GRADIENT:   1.6884E+00  9.8391E+00  5.2772E-01 -1.1334E+01  6.0404E-01 -2.3462E+00  1.3381E+01 -6.2405E-03 -1.8594E+01  9.1393E-01
            -4.2575E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1001.61529812290        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      693
 NPARAMETR:  1.0265E+00  1.7897E-01  9.6283E+01  1.7686E+00  2.3052E+00  2.4763E+00  7.4310E+00  8.4267E-01  1.8667E+00  5.4177E-01
             7.8833E+00
 PARAMETER:  1.2619E-01 -1.6205E+00  4.6673E+00  6.7021E-01  9.3515E-01  1.0068E+00  2.1057E+00 -7.1178E-02  7.2417E-01 -5.1291E-01
             2.1647E+00
 GRADIENT:   3.1698E-01  2.3658E+01  1.8596E+00 -2.3244E+01  6.6880E+00  2.1054E+00  4.1951E+01 -9.2123E-03 -4.4255E+01  1.0741E+00
            -3.4724E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1005.82957758069        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      874
 NPARAMETR:  1.0160E+00  6.4685E-02  2.2425E+02  1.7837E+00  2.4242E+00  2.4347E+00  1.0306E+01  1.9017E-01  1.6243E+00  1.9960E-01
             7.9749E+00
 PARAMETER:  1.1584E-01 -2.6382E+00  5.5127E+00  6.7869E-01  9.8551E-01  9.8981E-01  2.4328E+00 -1.5598E+00  5.8507E-01 -1.5115E+00
             2.1763E+00
 GRADIENT:  -4.5724E+00 -8.6405E+00  1.2576E-01  2.0617E+01  8.7862E+00 -2.6390E+00 -2.6910E+01  1.7333E-04  7.9438E+00 -7.2661E-02
             5.3124E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1007.54834977488        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1051
 NPARAMETR:  1.0282E+00  4.4193E-02  3.1505E+02  1.7851E+00  2.4159E+00  2.4243E+00  1.2144E+01  1.1585E-01  1.6026E+00  1.4970E-01
             7.9596E+00
 PARAMETER:  1.2778E-01 -3.0192E+00  5.8527E+00  6.7946E-01  9.8207E-01  9.8554E-01  2.5968E+00 -2.0554E+00  5.7161E-01 -1.7991E+00
             2.1744E+00
 GRADIENT:   3.0618E-01 -4.2138E+00  9.6199E-02 -7.3468E+00  9.0622E+00 -2.9818E+00 -1.4552E+01  4.6019E-05 -1.9085E+01 -6.8326E-02
            -2.3127E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1007.87240592003        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1244             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0281E+00  4.3258E-02  3.2343E+02  1.7877E+00  2.4134E+00  2.4292E+00  1.2371E+01  1.1256E-01  1.6138E+00  1.5081E-01
             7.9588E+00
 PARAMETER:  1.2768E-01 -3.0406E+00  5.8790E+00  6.8091E-01  9.8102E-01  9.8758E-01  2.6153E+00 -2.0843E+00  5.7862E-01 -1.7917E+00
             2.1743E+00
 GRADIENT:   1.3559E+00 -3.8533E+00  1.1022E-01 -1.5163E+00  8.6273E+00  2.0224E+00 -4.5499E-01  4.6919E-05 -1.5166E+01 -6.7081E-02
             1.9522E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1007.95339408551        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1423
 NPARAMETR:  1.0275E+00  4.3258E-02  3.2343E+02  1.7877E+00  2.4133E+00  2.4508E+00  1.2371E+01  1.1255E-01  1.6139E+00  3.5999E-01
             7.9588E+00
 PARAMETER:  1.2710E-01 -3.0406E+00  5.8790E+00  6.8091E-01  9.8100E-01  9.9643E-01  2.6153E+00 -2.0843E+00  5.7863E-01 -9.2167E-01
             2.1743E+00
 GRADIENT:  -4.4118E-02 -4.7143E+00  1.0368E-01 -6.4853E+00  1.1539E+01  7.9105E-01 -1.5013E+01  4.7930E-05 -1.4862E+01  1.7447E-02
             2.1577E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1007.96583327686        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1602
 NPARAMETR:  1.0279E+00  4.3258E-02  3.2340E+02  1.7879E+00  2.4117E+00  2.4334E+00  1.2374E+01  1.1236E-01  1.6142E+00  3.6652E-01
             7.9541E+00
 PARAMETER:  1.2749E-01 -3.0406E+00  5.8789E+00  6.8105E-01  9.8033E-01  9.8927E-01  2.6156E+00 -2.0861E+00  5.7884E-01 -9.0370E-01
             2.1737E+00
 GRADIENT:   7.9255E-02 -4.7229E+00  1.0657E-01 -6.1725E+00  1.1334E+01 -1.7010E+00 -1.5028E+01  4.4122E-05 -1.4698E+01  2.2667E-02
             5.1936E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1007.97292152317        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1788
 NPARAMETR:  1.0268E+00  4.3266E-02  3.2328E+02  1.7879E+00  2.4118E+00  2.4478E+00  1.2376E+01  1.1235E-01  1.6142E+00  3.5865E-01
             7.9550E+00
 PARAMETER:  1.2648E-01 -3.0404E+00  5.8785E+00  6.8104E-01  9.8036E-01  9.9520E-01  2.6157E+00 -2.0861E+00  5.7883E-01 -9.2542E-01
             2.1738E+00
 GRADIENT:  -2.2272E-01 -4.4106E+00  1.0352E-01 -6.5274E+00  1.1180E+01  3.3184E-01 -1.4396E+01  4.4820E-05 -1.5408E+01  5.6962E-03
             1.1361E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1007.99604759046        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1967
 NPARAMETR:  1.0247E+00  4.3267E-02  3.2314E+02  1.7887E+00  2.4065E+00  2.4397E+00  1.2385E+01  1.1162E-01  1.6154E+00  2.7491E-01
             7.9456E+00
 PARAMETER:  1.2437E-01 -3.0404E+00  5.8781E+00  6.8149E-01  9.7817E-01  9.9189E-01  2.6165E+00 -2.0927E+00  5.7957E-01 -1.1913E+00
             2.1726E+00
 GRADIENT:  -7.5383E-01 -4.6253E+00  1.2447E-01 -4.9586E+00  8.3408E+00 -9.7044E-01 -1.4838E+01  5.0053E-05 -1.5122E+01 -1.2309E-01
            -2.9990E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1008.09833052246        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     2110
 NPARAMETR:  1.0238E+00  4.3276E-02  3.2293E+02  1.7892E+00  2.3634E+00  2.4118E+00  1.2393E+01  1.1116E-01  1.6161E+00  3.7260E-01
             7.9439E+00
 PARAMETER:  1.2351E-01 -3.0402E+00  5.8774E+00  6.8174E-01  9.6012E-01  9.8036E-01  2.6171E+00 -2.0968E+00  5.8000E-01 -8.8725E-01
             2.1724E+00
 GRADIENT:  -1.0916E-01 -4.0917E+00  2.2094E-01  1.4441E+00  7.3632E-01 -4.5335E-01 -8.9137E-01  5.4646E-05 -1.3449E+01 -1.3346E-01
             1.6243E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1008.11235860740        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     2180
 NPARAMETR:  1.0238E+00  4.3276E-02  3.2292E+02  1.7892E+00  2.3540E+00  2.4154E+00  1.2393E+01  1.1113E-01  1.6161E+00  4.3288E-01
             7.9312E+00
 PARAMETER:  1.2349E-01 -3.0402E+00  5.8774E+00  6.8174E-01  9.5612E-01  9.8185E-01  2.6171E+00 -2.0970E+00  5.8000E-01 -7.3729E-01
             2.1708E+00
 GRADIENT:   4.6740E-03 -4.0977E+00  2.3928E-01  1.8960E+00 -4.4946E-02  3.3593E-02 -8.6501E-01  4.3417E-05 -1.3282E+01 -8.4416E-04
            -4.9640E-02

0ITERATION NO.:   82    OBJECTIVE VALUE:  -1008.11235860740        NO. OF FUNC. EVALS.:  60
 CUMULATIVE NO. OF FUNC. EVALS.:     2240
 NPARAMETR:  1.0238E+00  4.3276E-02  3.2292E+02  1.7892E+00  2.3540E+00  2.4154E+00  1.2393E+01  1.1113E-01  1.6161E+00  4.3288E-01
             7.9312E+00
 PARAMETER:  1.2349E-01 -3.0402E+00  5.8774E+00  6.8174E-01  9.5612E-01  9.8185E-01  2.6171E+00 -2.0970E+00  5.8000E-01 -7.3729E-01
             2.1708E+00
 GRADIENT:   6.3466E+03 -1.4486E+02  6.6940E+01  1.8495E+04 -4.0990E+02 -4.0350E+02 -4.6737E+03 -2.6658E-03  1.3325E+03 -5.1234E-03
            -1.8308E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2240
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.7076E-03  2.5242E-02  1.0702E-05 -9.5528E-03 -4.9155E-03
 SE:             2.9470E-02  1.1353E-02  1.5468E-05  2.8814E-02  6.8856E-03
 N:                     100         100         100         100         100

 P VAL.:         7.9367E-01  2.6192E-02  4.8902E-01  7.4025E-01  4.7530E-01

 ETASHRINKSD(%)  1.2731E+00  6.1965E+01  9.9948E+01  3.4682E+00  7.6933E+01
 ETASHRINKVR(%)  2.5300E+00  8.5533E+01  1.0000E+02  6.8161E+00  9.4679E+01
 EBVSHRINKSD(%)  1.9771E+00  7.2794E+01  9.9929E+01  8.0552E+00  7.6827E+01
 EBVSHRINKVR(%)  3.9152E+00  9.2598E+01  1.0000E+02  1.5461E+01  9.4630E+01
 RELATIVEINF(%)  9.5942E+01  5.9156E+00  1.5056E-05  6.5288E+01  1.5728E+00
 EPSSHRINKSD(%)  8.2711E+00
 EPSSHRINKVR(%)  1.5858E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1008.1123586074050     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       645.97700116100577     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    71.53
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1008.112       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  4.33E-02  3.23E+02  1.79E+00  2.35E+00  2.42E+00  1.24E+01  1.11E-01  1.62E+00  4.33E-01  7.93E+00
 


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
+        1.23E+07
 
 TH 2
+        1.32E+02  1.76E+08
 
 TH 3
+       -1.56E-02  2.41E-01  5.44E-02
 
 TH 4
+       -1.92E+07 -1.84E+07 -4.65E-02  2.18E+06
 
 TH 5
+        1.14E+01 -2.35E+03  1.07E-01 -1.08E+06  3.87E+04
 
 TH 6
+       -7.96E+01  1.84E+01 -4.92E-03 -1.02E+06  3.42E+00  3.01E+01
 
 TH 7
+       -2.89E+00  7.45E+05  3.05E-03  4.92E+03 -3.52E+01 -6.23E+00  4.10E+03
 
 TH 8
+       -3.41E+01 -4.39E+01 -8.68E-05 -6.23E+00  1.14E+00 -3.55E+05  6.98E-03  4.05E+01
 
 TH 9
+       -3.30E+01 -2.39E+07 -7.79E-03  2.66E+06  2.25E+02 -1.26E+01  7.03E+03 -5.02E+00  3.23E+05
 
 TH10
+        4.86E+06  7.42E+01  2.32E-03  6.51E+01  3.97E+00  2.59E+05 -1.55E+00 -7.24E+00  2.64E-01  1.92E+06
 
 TH11
+       -4.42E+00  1.31E+06  2.51E-02 -1.41E+05 -2.11E+01  2.53E+00 -1.77E+00  5.08E-01 -1.83E+05  1.63E+00  6.75E+02
 
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
 #CPUT: Total CPU Time in Seconds,       88.742
Stop Time:
Sat Sep 25 05:50:23 CDT 2021
