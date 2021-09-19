Sat Sep 18 05:30:11 CDT 2021
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
$DATA ../../../../data/int/TD1/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3847.94329068900        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.8436E+01 -7.3383E+01 -5.1106E+01 -7.9547E+01  7.9808E+00  4.7282E+00 -1.3484E+01 -1.9737E+00 -1.0978E+01  3.2676E+00
            -5.1237E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3856.78796548879        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8529E-01  1.3173E+00  1.5451E+00  8.8499E-01  1.3793E+00  9.8724E-01  1.1851E+00  1.0885E+00  9.7349E-01  1.1267E+00
             1.0883E+00
 PARAMETER:  8.5178E-02  3.7558E-01  5.3511E-01 -2.2175E-02  4.2155E-01  8.7155E-02  2.6986E-01  1.8485E-01  7.3129E-02  2.1932E-01
             1.8466E-01
 GRADIENT:   5.8985E+01  3.3700E+01  2.5640E+01 -8.5392E+00  3.7531E+01 -2.0467E-01  1.5573E+01 -5.7633E+00  9.3350E+00 -6.5895E+00
             7.3367E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3858.15261421303        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.8927E-01  1.3051E+00  1.5917E+00  9.0034E-01  1.3656E+00  9.7376E-01  1.2105E+00  1.5161E+00  9.0021E-01  1.1230E+00
             1.0634E+00
 PARAMETER:  8.9211E-02  3.6632E-01  5.6478E-01 -4.9779E-03  4.1161E-01  7.3413E-02  2.9103E-01  5.1615E-01 -5.1245E-03  2.1597E-01
             1.6148E-01
 GRADIENT:   7.0981E+01  3.9514E+01  2.0265E+01 -2.5111E+00  2.5811E+01 -6.0053E+00  1.4798E+01  3.8676E+00  4.4139E+00 -8.0662E+00
             4.5469E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3861.58977603685        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.7873E-01  1.2636E+00  1.4516E+00  9.2220E-01  1.2966E+00  9.8182E-01  1.1106E+00  1.3364E+00  9.5948E-01  1.1078E+00
             1.0511E+00
 PARAMETER:  7.8500E-02  3.3393E-01  4.7265E-01  1.9004E-02  3.5974E-01  8.1658E-02  2.0493E-01  3.8995E-01  5.8631E-02  2.0241E-01
             1.4982E-01
 GRADIENT:   4.6856E+01  3.0320E+01  1.2439E+01  3.6725E+00  1.6927E+01 -1.8587E+00  2.5358E+00  4.2786E+00  3.0562E+00 -3.2618E+00
             2.5446E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3861.63913753465        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.7842E-01  1.2625E+00  1.4492E+00  9.2264E-01  1.2951E+00  9.8209E-01  1.1101E+00  1.3324E+00  9.5941E-01  1.1077E+00
             1.0507E+00
 PARAMETER:  7.8181E-02  3.3311E-01  4.7100E-01  1.9487E-02  3.5861E-01  8.1931E-02  2.0444E-01  3.8701E-01  5.8567E-02  2.0228E-01
             1.4945E-01
 GRADIENT:   4.6154E+01  3.0029E+01  1.2334E+01  3.6036E+00  1.6636E+01 -1.7353E+00  2.4532E+00  4.2278E+00  2.9086E+00 -3.1555E+00
             2.4712E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3861.87537697891        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      559
 NPARAMETR:  9.7598E-01  1.2625E+00  1.4492E+00  9.2264E-01  1.2943E+00  9.9075E-01  1.1125E+00  1.3324E+00  9.4962E-01  1.1238E+00
             1.0413E+00
 PARAMETER:  7.5685E-02  3.3310E-01  4.7101E-01  1.9485E-02  3.5794E-01  9.0711E-02  2.0659E-01  3.8701E-01  4.8311E-02  2.1674E-01
             1.4046E-01
 GRADIENT:   2.8640E+00  2.1599E+00  1.1050E+01 -3.2616E+00  1.1832E+00 -1.4247E+00  1.6209E-01  3.2589E+00  7.1172E-01 -1.2903E+00
             6.3116E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3861.87759343313        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      753             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7601E-01  1.2624E+00  1.4494E+00  9.2262E-01  1.2943E+00  9.9451E-01  1.1119E+00  1.3323E+00  9.4863E-01  1.1239E+00
             1.0413E+00
 PARAMETER:  7.5714E-02  3.3301E-01  4.7114E-01  1.9457E-02  3.5795E-01  9.4490E-02  2.0611E-01  3.8690E-01  4.7260E-02  2.1680E-01
             1.4046E-01
 GRADIENT:   4.1313E+01  3.0564E+01  1.2996E+01  3.4842E+00  1.4888E+01  3.3334E+00  1.8814E+00  3.4766E+00  1.0636E+00 -1.3478E-01
             6.5327E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3861.98063294324        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      869
 NPARAMETR:  9.7334E-01  1.2575E+00  1.4498E+00  9.2479E-01  1.2917E+00  9.9379E-01  1.1112E+00  1.2539E+00  9.4579E-01  1.1307E+00
             1.0377E+00
 PARAMETER:  7.2979E-02  3.2915E-01  4.7142E-01  2.1806E-02  3.5593E-01  9.3770E-02  2.0548E-01  3.2629E-01  4.4262E-02  2.2287E-01
             1.3701E-01
 GRADIENT:   3.5521E+01  2.8759E+01  1.7528E+01  4.0523E+00  1.5564E+01  3.1191E+00  1.4578E+00 -5.3867E-03 -3.6434E-01  1.8519E+00
            -4.0926E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3862.00870516906        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  9.7487E-01  1.2497E+00  1.4499E+00  9.3056E-01  1.2838E+00  9.9429E-01  1.1105E+00  1.2458E+00  9.5336E-01  1.1214E+00
             1.0397E+00
 PARAMETER:  7.4549E-02  3.2288E-01  4.7148E-01  2.8035E-02  3.4979E-01  9.4270E-02  2.0478E-01  3.1976E-01  5.2235E-02  2.1454E-01
             1.3894E-01
 GRADIENT:   3.5184E-01  3.2910E-01  1.6742E+01 -3.3380E-02  1.8744E-01  2.0726E-02 -5.6290E-02 -1.2664E-01  3.3051E-02 -3.2961E-02
             1.7158E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3862.01181780471        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1207
 NPARAMETR:  9.7496E-01  1.2487E+00  1.4495E+00  9.3153E-01  1.2828E+00  9.9413E-01  1.1130E+00  1.2535E+00  9.4984E-01  1.1188E+00
             1.0400E+00
 PARAMETER:  7.4641E-02  3.2214E-01  4.7120E-01  2.9077E-02  3.4901E-01  9.4117E-02  2.0704E-01  3.2595E-01  4.8543E-02  2.1228E-01
             1.3925E-01
 GRADIENT:   5.5028E-01  5.6801E-01  1.6200E+01  6.7974E-01 -1.9402E-02 -4.1058E-02  1.0749E-01  1.5363E-01 -3.9385E-01 -4.6620E-01
             1.2322E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3862.21907150074        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1388
 NPARAMETR:  9.7622E-01  1.2444E+00  1.4208E+00  9.3588E-01  1.2717E+00  9.9349E-01  1.1222E+00  1.2142E+00  9.3637E-01  1.1021E+00
             1.0444E+00
 PARAMETER:  7.5930E-02  3.1866E-01  4.5123E-01  3.3732E-02  3.4038E-01  9.3466E-02  2.1525E-01  2.9408E-01  3.4259E-02  1.9718E-01
             1.4342E-01
 GRADIENT:   3.3450E+00  2.7231E+00  1.2517E+01  4.2355E+00 -5.0219E-02 -2.8827E-01  1.1033E+00  2.5178E-01 -3.2406E+00 -2.9752E+00
             8.7670E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3862.29823108104        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:     1538
 NPARAMETR:  9.7626E-01  1.2441E+00  1.4205E+00  9.3586E-01  1.2713E+00  9.8904E-01  1.0992E+00  1.2130E+00  9.6782E-01  1.1180E+00
             1.0382E+00
 PARAMETER:  7.5977E-02  3.1843E-01  4.5098E-01  3.3710E-02  3.4002E-01  8.8978E-02  1.9460E-01  2.9307E-01  6.7286E-02  2.1158E-01
             1.3752E-01
 GRADIENT:   4.2173E+01  2.9066E+01  1.5701E+01  1.1691E+01  1.3618E+01  1.2531E+00  1.2327E+00  3.8301E-01  1.6270E+00  1.6102E+00
            -2.0408E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3862.31099321784        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1694
 NPARAMETR:  9.7626E-01  1.2441E+00  1.4205E+00  9.3586E-01  1.2713E+00  9.9424E-01  1.1055E+00  1.2130E+00  9.6011E-01  1.1156E+00
             1.0395E+00
 PARAMETER:  7.5977E-02  3.1842E-01  4.5098E-01  3.3710E-02  3.4002E-01  9.4227E-02  2.0031E-01  2.9307E-01  5.9288E-02  2.0938E-01
             1.3872E-01
 GRADIENT:   3.5281E+00  2.1206E+00  1.3505E+01  4.5067E+00  5.8987E-02  2.6677E-02  1.0956E-02  2.0399E-01  1.6761E-02 -1.7886E-02
             1.0820E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3862.37250516592        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:     1844
 NPARAMETR:  9.7483E-01  1.2424E+00  1.4145E+00  9.3537E-01  1.2699E+00  9.9386E-01  1.1049E+00  1.2120E+00  9.6032E-01  1.1153E+00
             1.0394E+00
 PARAMETER:  7.4504E-02  3.1706E-01  4.4676E-01  3.3188E-02  3.3894E-01  9.3845E-02  1.9973E-01  2.9229E-01  5.9515E-02  2.0912E-01
             1.3866E-01
 GRADIENT:   3.8761E+01  2.7207E+01  1.4243E+01  9.8649E+00  1.3833E+01  3.1798E+00  1.7352E+00  6.2112E-01  4.7170E-01  1.1289E+00
             4.9942E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3862.37274287970        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     2010
 NPARAMETR:  9.7474E-01  1.2424E+00  1.4145E+00  9.3537E-01  1.2693E+00  9.9414E-01  1.1047E+00  1.2120E+00  9.6109E-01  1.1149E+00
             1.0393E+00
 PARAMETER:  7.4410E-02  3.1706E-01  4.4676E-01  3.3188E-02  3.3844E-01  9.4121E-02  1.9957E-01  2.9229E-01  6.0313E-02  2.0878E-01
             1.3850E-01
 GRADIENT:   3.5384E-02  6.2206E-01  1.2549E+01  2.6116E+00 -1.1123E-02 -9.5865E-03 -1.2586E-02  4.6228E-01 -2.9299E-03 -5.1888E-03
            -3.2667E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2010
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.8658E-04 -2.5182E-02 -3.3235E-02  1.9143E-02 -2.8989E-02
 SE:             2.9856E-02  2.3796E-02  1.6362E-02  2.5614E-02  2.4795E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8967E-01  2.8993E-01  4.2238E-02  4.5485E-01  2.4236E-01

 ETASHRINKSD(%)  1.0000E-10  2.0281E+01  4.5184E+01  1.4189E+01  1.6932E+01
 ETASHRINKVR(%)  1.0000E-10  3.6449E+01  6.9952E+01  2.6364E+01  3.0998E+01
 EBVSHRINKSD(%)  3.1734E-01  1.9878E+01  4.8205E+01  1.7195E+01  1.5279E+01
 EBVSHRINKVR(%)  6.3367E-01  3.5805E+01  7.3173E+01  3.1434E+01  2.8224E+01
 RELATIVEINF(%)  9.9365E+01  3.3587E+01  2.0389E+01  3.8207E+01  3.6207E+01
 EPSSHRINKSD(%)  2.0728E+01
 EPSSHRINKVR(%)  3.7159E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3862.3727428796997     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2208.2833831112889     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    56.22
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3862.373       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.24E+00  1.41E+00  9.35E-01  1.27E+00  9.94E-01  1.10E+00  1.21E+00  9.61E-01  1.11E+00  1.04E+00
 


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
+        3.93E+07
 
 TH 2
+       -9.72E+06  2.41E+06
 
 TH 3
+        1.20E+02 -7.77E+02  9.28E+05
 
 TH 4
+       -4.09E+07  1.01E+07  1.32E+03  4.26E+07
 
 TH 5
+       -8.92E+06  2.21E+06 -2.32E+03  9.29E+06  3.98E+02
 
 TH 6
+        3.89E+00 -1.21E+02  7.40E+01 -5.09E+02 -5.31E-01  1.94E+02
 
 TH 7
+        4.61E-01 -4.30E+06  4.78E+01 -1.81E+07 -3.41E+00 -1.57E+00  6.67E+01
 
 TH 8
+       -2.12E+02  2.67E+06  5.08E+03  1.12E+07  3.97E+03 -1.32E+02 -9.54E+01  2.96E+06
 
 TH 9
+        4.40E+00 -9.85E+06 -1.13E+03 -4.15E+07  1.33E+01 -1.12E+00  2.69E+01  2.02E+03  9.12E+01
 
 TH10
+        1.64E+07 -4.07E+06  2.27E+02 -1.71E+07 -2.55E+01  1.46E+00  2.57E+00 -4.03E+02  6.32E+00  9.11E+01
 
 TH11
+       -1.08E+01 -6.58E+06 -3.79E+03 -2.77E+07 -1.87E+01  4.94E-01  1.33E+01  6.76E+03  1.86E+01  1.50E+01  9.64E+02
 
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
 #CPUT: Total CPU Time in Seconds,       69.900
Stop Time:
Sat Sep 18 05:31:22 CDT 2021
