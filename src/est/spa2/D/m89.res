Thu Sep 30 10:01:46 CDT 2021
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
$DATA ../../../../data/spa2/D/dat89.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   35154.6704254300        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.9051E+02  6.9523E+02  3.0562E+01  7.2454E+02  2.0572E+02 -3.3440E+03 -1.5438E+03 -6.8200E+01 -1.9231E+03 -9.7294E+02
            -6.6563E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -397.639402754159        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.3704E+00  1.2058E+00  9.7017E-01  1.3529E+00  1.0241E+00  1.7499E+00  1.2407E+00  9.8676E-01  1.0955E+00  1.0223E+00
             1.4737E+01
 PARAMETER:  4.1509E-01  2.8714E-01  6.9718E-02  4.0222E-01  1.2378E-01  6.5955E-01  3.1567E-01  8.6676E-02  1.9124E-01  1.2208E-01
             2.7903E+00
 GRADIENT:   5.6002E+01  6.2721E-01  6.3706E-01  1.1248E+01 -7.9327E-01  1.0606E+01 -2.0463E+01  2.7775E+00 -9.3921E+00  8.9878E+00
            -3.5996E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -489.697964138463        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.3259E+00  1.1334E+00  2.0407E+00  1.4237E+00  1.2325E+00  2.2866E+00  5.6360E+00  5.6757E-01  1.1562E+00  3.0737E-01
             1.4203E+01
 PARAMETER:  3.8208E-01  2.2525E-01  8.1327E-01  4.5325E-01  3.0901E-01  9.2707E-01  1.8292E+00 -4.6639E-01  2.4515E-01 -1.0797E+00
             2.7534E+00
 GRADIENT:   3.4975E+01 -1.1423E+00  9.8448E-01 -2.3403E+01 -3.2209E+01  4.6659E+01  4.2899E+01  4.5595E-01  8.2756E+00  6.3050E-01
             1.0072E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -503.046801838613        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      256
 NPARAMETR:  1.2343E+00  1.0410E+00  6.9285E+00  1.5915E+00  2.0171E+00  2.1382E+00  5.2024E+00  4.6359E-01  1.3011E+00  9.2169E-01
             1.3405E+01
 PARAMETER:  3.1047E-01  1.4017E-01  2.0356E+00  5.6469E-01  8.0166E-01  8.5995E-01  1.7491E+00 -6.6876E-01  3.6318E-01  1.8459E-02
             2.6956E+00
 GRADIENT:   8.2663E+00 -4.1177E+00  4.2540E-01 -2.3192E+00 -1.4458E+01  2.3524E+01 -3.9034E+01  3.6566E-02  2.8905E+00  4.4806E+00
             2.0180E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -528.958941794849        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      434
 NPARAMETR:  1.0772E+00  3.6427E-01  8.4205E+01  2.2274E+00  3.1544E+00  1.5048E+00  8.3929E+00  5.8704E+00  1.3885E+00  3.2708E-01
             1.3550E+01
 PARAMETER:  1.7434E-01 -9.0987E-01  4.5333E+00  9.0083E-01  1.2488E+00  5.0866E-01  2.2274E+00  1.8699E+00  4.2825E-01 -1.0176E+00
             2.7064E+00
 GRADIENT:  -1.4582E+01 -2.6160E-02 -9.3428E-02  4.8518E+01 -2.1466E+00  7.5182E+00 -1.4778E+01  3.9224E-02  3.4140E-01  2.7813E-01
            -2.5230E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -532.037830400062        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      614
 NPARAMETR:  1.0895E+00  3.5627E-01  2.0776E+02  2.0626E+00  3.0839E+00  1.4161E+00  8.6158E+00  3.1799E+00  1.2354E+00  1.0947E-01
             1.3671E+01
 PARAMETER:  1.8568E-01 -9.3207E-01  5.4364E+00  8.2396E-01  1.2262E+00  4.4792E-01  2.2536E+00  1.2568E+00  3.1140E-01 -2.1121E+00
             2.7153E+00
 GRADIENT:  -2.7367E-01  3.5587E-01 -8.3032E-03  3.1114E+00 -3.5236E-01  4.7987E-01  4.0045E-01  1.7063E-03 -1.5514E+00  3.2721E-02
            -5.5775E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -533.009671417798        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      794
 NPARAMETR:  1.0783E+00  3.1036E-01  7.8437E+02  2.0745E+00  3.1042E+00  1.3017E+00  9.1884E+00  2.7428E+00  1.2456E+00  3.8807E-02
             1.3750E+01
 PARAMETER:  1.7535E-01 -1.0700E+00  6.7649E+00  8.2974E-01  1.2328E+00  3.6367E-01  2.3179E+00  1.1090E+00  3.1959E-01 -3.1492E+00
             2.7211E+00
 GRADIENT:   9.2450E+00  3.5455E-01 -2.5533E-03  8.8817E-01 -3.2610E-01 -3.9793E+00  5.8687E+00  9.2766E-05  2.6291E+00  4.1161E-03
            -9.3186E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -533.421188147643        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      956
 NPARAMETR:  1.0690E+00  3.0014E-01  3.4388E+03  2.0667E+00  3.1339E+00  1.3079E+00  9.7032E+00  3.0226E+00  1.2090E+00  2.5067E-02
             1.3810E+01
 PARAMETER:  1.6677E-01 -1.1035E+00  8.2429E+00  8.2597E-01  1.2423E+00  3.6838E-01  2.3725E+00  1.2061E+00  2.8979E-01 -3.5862E+00
             2.7254E+00
 GRADIENT:   3.6531E+00  2.1304E+00 -3.5380E-04 -8.5594E+00  3.9903E-01  1.5421E-01  1.8451E+01  2.5318E-07  5.3928E-01  1.6653E-03
             1.1311E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -533.780670759864        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1124
 NPARAMETR:  1.0649E+00  2.5575E-01  5.4629E+03  2.0883E+00  3.1223E+00  1.3033E+00  1.0030E+01  3.1042E+00  1.2097E+00  1.6519E-02
             1.3808E+01
 PARAMETER:  1.6286E-01 -1.2635E+00  8.7057E+00  8.3635E-01  1.2386E+00  3.6488E-01  2.4056E+00  1.2328E+00  2.9037E-01 -4.0032E+00
             2.7253E+00
 GRADIENT:   1.1338E+00  1.4772E+00 -1.4277E-04 -6.3785E+00  1.3804E-01  3.6756E-01  1.9599E+01  9.8291E-06  8.1698E-02  7.2914E-04
             1.4302E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -533.884263589107        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1316             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0669E+00  2.3234E-01  7.6354E+03  2.1063E+00  3.1385E+00  1.3029E+00  1.0520E+01  3.0909E+00  1.2183E+00  1.5128E-02
             1.3844E+01
 PARAMETER:  1.6472E-01 -1.3596E+00  9.0405E+00  8.4494E-01  1.2438E+00  3.6456E-01  2.4533E+00  1.2284E+00  2.9750E-01 -4.0912E+00
             2.7279E+00
 GRADIENT:   6.8021E+00  6.1526E+00 -5.5491E-05  6.4277E+00  9.0565E-01  2.8345E+00  1.8767E+02 -1.4992E-05  1.0988E+00  6.4088E-04
             3.5758E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -534.015180447502        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1474
 NPARAMETR:  1.0656E+00  2.1602E-01  7.6723E+03  2.1196E+00  3.1382E+00  1.3009E+00  1.0289E+01  3.0872E+00  1.2161E+00  1.4581E-02
             1.3848E+01
 PARAMETER:  1.6355E-01 -1.4324E+00  9.0454E+00  8.5122E-01  1.2437E+00  3.6305E-01  2.4311E+00  1.2273E+00  2.9568E-01 -4.1280E+00
             2.7281E+00
 GRADIENT:  -2.6251E-01  1.5582E-01 -1.0286E-04 -4.9441E-01  7.8248E-02  9.7912E-02  1.7864E+01  4.6878E-06  6.2844E-02  5.6234E-04
            -2.2363E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -534.126832826567        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1659
 NPARAMETR:  1.0667E+00  1.8593E-01  9.0459E+03  2.1324E+00  3.1371E+00  1.2943E+00  1.0678E+01  3.2470E+00  1.2184E+00  1.0000E-02
             1.3840E+01
 PARAMETER:  1.6458E-01 -1.5824E+00  9.2101E+00  8.5726E-01  1.2433E+00  3.5799E-01  2.4682E+00  1.2777E+00  2.9756E-01 -4.5408E+00
             2.7276E+00
 GRADIENT:   1.3165E+00 -4.2564E-01 -5.4229E-05 -8.4750E-01  1.1527E-01 -3.9478E-01  2.1087E+01 -1.0709E-05  1.8795E-01  2.8753E-05
            -2.3903E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -534.152076007476        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1847             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0664E+00  1.8235E-01  6.4687E+03  2.1413E+00  3.1387E+00  1.2950E+00  1.0787E+01  3.4494E+00  1.2204E+00  1.0334E-02
             1.3872E+01
 PARAMETER:  1.6433E-01 -1.6018E+00  8.8747E+00  8.6142E-01  1.2438E+00  3.5855E-01  2.4784E+00  1.3382E+00  2.9917E-01 -4.4723E+00
             2.7299E+00
 GRADIENT:   4.7071E+00  4.1732E+00 -8.9579E-05  1.6400E+01  4.5127E-01  1.7637E+00  1.9818E+02  1.5272E-05  9.0199E-01  3.0157E-04
             3.2391E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -534.166585536181        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2030
 NPARAMETR:  1.0650E+00  1.7900E-01  4.8091E+03  2.1440E+00  3.1336E+00  1.2944E+00  1.0895E+01  3.4509E+00  1.2195E+00  1.0000E-02
             1.3867E+01
 PARAMETER:  1.6296E-01 -1.6204E+00  8.5783E+00  8.6269E-01  1.2422E+00  3.5807E-01  2.4883E+00  1.3386E+00  2.9844E-01 -4.5515E+00
             2.7295E+00
 GRADIENT:  -4.1622E-01 -1.1268E-01 -8.8467E-05  2.5434E-01 -5.6054E-02  8.9604E-02  2.3923E+01 -7.6167E-08 -5.2278E-02  0.0000E+00
            -7.0402E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -534.170428332288        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2214
 NPARAMETR:  1.0654E+00  1.7870E-01  5.9803E+03  2.1454E+00  3.1351E+00  1.2939E+00  1.0935E+01  1.5502E+00  1.2198E+00  1.0000E-02
             1.3872E+01
 PARAMETER:  1.6331E-01 -1.6221E+00  8.7962E+00  8.6334E-01  1.2427E+00  3.5766E-01  2.4920E+00  5.3839E-01  2.9871E-01 -4.5734E+00
             2.7298E+00
 GRADIENT:  -1.8001E-01  9.9968E-03 -7.2695E-05  1.0510E-01 -4.0898E-02  6.4279E-02  2.4659E+01  7.1839E-07 -6.0260E-02  0.0000E+00
            -5.5511E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -534.171852254472        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2390
 NPARAMETR:  1.0653E+00  1.7687E-01  8.9168E+03  2.1466E+00  3.1353E+00  1.2930E+00  1.0966E+01  1.5415E+00  1.2199E+00  1.0235E-02
             1.3872E+01
 PARAMETER:  1.6325E-01 -1.6324E+00  9.1957E+00  8.6388E-01  1.2427E+00  3.5699E-01  2.4948E+00  5.3278E-01  2.9875E-01 -4.4819E+00
             2.7299E+00
 GRADIENT:  -9.0237E-02 -7.4236E-02  1.5117E-04  3.1963E-01 -4.3190E-02 -3.2225E-02 -6.9884E-01  2.3216E-03 -6.4085E-02  5.9685E-04
            -3.8999E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2390
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4471E-02  7.9057E-02  1.9981E-06 -8.3004E-02  4.4540E-06
 SE:             2.5933E-02  1.9286E-02  5.1783E-07  1.6207E-02  2.7193E-05
 N:                     100         100         100         100         100

 P VAL.:         3.4536E-01  4.1468E-05  1.1412E-04  3.0372E-07  8.6989E-01

 ETASHRINKSD(%)  1.3120E+01  3.5391E+01  9.9998E+01  4.5704E+01  9.9909E+01
 ETASHRINKVR(%)  2.4519E+01  5.8257E+01  1.0000E+02  7.0520E+01  1.0000E+02
 EBVSHRINKSD(%)  1.8607E+01  3.9242E+01  9.9996E+01  3.6918E+01  9.9850E+01
 EBVSHRINKVR(%)  3.3752E+01  6.3085E+01  1.0000E+02  6.0206E+01  1.0000E+02
 RELATIVEINF(%)  5.9720E+01  2.6517E+01  2.6908E-08  1.8876E+01  3.5816E-05
 EPSSHRINKSD(%)  4.8349E+00
 EPSSHRINKVR(%)  9.4360E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -534.17185225447190     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       568.55438759113520     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    64.56
 Elapsed covariance  time in seconds:    12.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -534.172       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.07E+00  1.77E-01  8.92E+03  2.15E+00  3.14E+00  1.29E+00  1.10E+01  1.54E+00  1.22E+00  1.02E-02  1.39E+01
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.72E-02  1.09E-01  3.86E+04  2.31E-01  4.17E-01  2.40E-01  1.61E+00  4.17E-01  2.50E-01  2.51E-02  2.33E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.45E-03
 
 TH 2
+       -1.15E-03  1.19E-02
 
 TH 3
+        2.94E+02 -6.89E+02  1.49E+09
 
 TH 4
+        1.87E-02 -9.70E-03 -1.44E+03  5.32E-02
 
 TH 5
+        3.15E-03 -5.94E-03 -4.69E+03  2.95E-02  1.74E-01
 
 TH 6
+        1.54E-03  1.33E-02 -2.89E+02 -6.20E-03 -3.32E-02  5.74E-02
 
 TH 7
+        3.02E-02 -1.45E-01  2.88E+02  1.89E-01  2.90E-01 -1.94E-01  2.59E+00
 
 TH 8
+        1.10E-02 -3.69E-02  9.92E+03  3.32E-02  3.55E-02 -5.04E-02  4.76E-01  1.74E-01
 
 TH 9
+        2.90E-03  1.06E-02 -6.18E+03  6.46E-03 -1.76E-02  4.53E-02 -1.37E-01 -7.36E-02  6.23E-02
 
 TH10
+        9.65E-04 -1.78E-03  7.05E+02  2.06E-03 -1.21E-03 -4.40E-05  1.92E-02  8.88E-03 -2.34E-03  6.29E-04
 
 TH11
+        1.32E-01 -1.26E-01  6.34E+03  3.81E-01  4.03E-01 -2.98E-01  2.13E+00  5.53E-01 -2.06E-01  2.12E-02  5.44E+00
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.72E-02
 
 TH 2
+       -1.09E-01  1.09E-01
 
 TH 3
+        7.84E-02 -1.64E-01  3.86E+04
 
 TH 4
+        8.34E-01 -3.86E-01 -1.62E-01  2.31E-01
 
 TH 5
+        7.77E-02 -1.31E-01 -2.91E-01  3.07E-01  4.17E-01
 
 TH 6
+        6.62E-02  5.08E-01 -3.12E-02 -1.12E-01 -3.32E-01  2.40E-01
 
 TH 7
+        1.93E-01 -8.27E-01  4.63E-03  5.10E-01  4.31E-01 -5.03E-01  1.61E+00
 
 TH 8
+        2.72E-01 -8.12E-01  6.16E-01  3.45E-01  2.04E-01 -5.05E-01  7.11E-01  4.17E-01
 
 TH 9
+        1.20E-01  3.88E-01 -6.41E-01  1.12E-01 -1.69E-01  7.57E-01 -3.41E-01 -7.08E-01  2.50E-01
 
 TH10
+        3.96E-01 -6.51E-01  7.28E-01  3.55E-01 -1.15E-01 -7.32E-03  4.75E-01  8.50E-01 -3.74E-01  2.51E-02
 
 TH11
+        5.80E-01 -4.95E-01  7.03E-02  7.07E-01  4.14E-01 -5.34E-01  5.66E-01  5.69E-01 -3.54E-01  3.62E-01  2.33E+00
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.62E+05
 
 TH 2
+       -3.54E+05  7.90E+05
 
 TH 3
+        3.43E-01 -8.27E-01  1.32E-06
 
 TH 4
+       -1.38E+04  2.73E+04 -1.47E-02  1.90E+03
 
 TH 5
+        4.91E+03 -1.14E+04  2.04E-02 -4.26E+02  5.60E+02
 
 TH 6
+        1.44E+05 -3.15E+05  2.69E-01 -1.15E+04  2.39E+03  1.38E+05
 
 TH 7
+       -6.62E+02  1.51E+03 -9.68E-04  3.10E+01  1.46E+01 -7.92E+02  8.83E+00
 
 TH 8
+       -1.81E+04  3.69E+04 -4.44E-02  2.75E+03 -1.94E+03 -8.27E+03 -9.31E+01  9.28E+03
 
 TH 9
+       -9.86E+04  2.09E+05 -1.40E-01  9.61E+03 -1.53E+03 -9.39E+04  4.96E+02  8.84E+03  6.89E+04
 
 TH10
+       -1.69E+06  3.85E+06 -4.13E+00  1.06E+05 -3.86E+04 -1.60E+06  9.34E+03  7.40E+04  9.98E+05  2.03E+07
 
 TH11
+        9.21E+02 -1.95E+03  9.84E-04 -8.99E+01 -2.68E+00  9.67E+02 -6.17E+00 -1.93E+01 -6.98E+02 -1.01E+04  8.47E+00
 
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
 #CPUT: Total CPU Time in Seconds,       77.061
Stop Time:
Thu Sep 30 10:03:05 CDT 2021
