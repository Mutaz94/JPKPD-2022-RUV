Wed Sep 29 08:23:54 CDT 2021
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
$DATA ../../../../data/int/D/dat27.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m27.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12544.4991320391        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0252E+02  1.6119E+02 -6.8299E+01  1.1558E+01  3.8222E+02 -1.0142E+03 -4.2382E+02 -1.0073E+02 -9.2253E+02 -4.9295E+02
            -2.9759E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1150.04552830027        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  7.5238E-01  1.4519E+00  9.3648E-01  2.6754E+00  6.9911E-01  3.4712E+00  3.1404E+00  1.0267E+00  3.9418E+00  2.3357E+00
             1.1776E+01
 PARAMETER: -1.8452E-01  4.7287E-01  3.4373E-02  1.0841E+00 -2.5794E-01  1.3445E+00  1.2444E+00  1.2632E-01  1.4717E+00  9.4831E-01
             2.5661E+00
 GRADIENT:  -1.0988E+02  2.6579E+01 -2.2864E+01  8.2818E+01 -1.8256E+01  8.6632E+01  5.0115E+01  5.7088E+00  8.3236E+01  3.7590E+01
             7.9675E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1339.16615260455        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  6.9720E-01  8.2581E-01  4.3955E+02  3.6093E+00  3.8132E+00  4.0468E+00  8.1918E+00  7.0797E-01  2.6370E+00  2.0796E+00
             1.0002E+01
 PARAMETER: -2.6068E-01 -9.1394E-02  6.1858E+00  1.3835E+00  1.4385E+00  1.4979E+00  2.2031E+00 -2.4535E-01  1.0696E+00  8.3216E-01
             2.4028E+00
 GRADIENT:  -7.1548E+01  2.2894E+01 -8.8115E-01  1.7029E+02  8.3828E+01  1.7800E+02  5.2924E+01  2.9821E-04  1.2151E+01  4.4619E+01
             6.7578E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1717.70844766004        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.1860E+00  2.5754E-01  7.1291E+02  1.3970E+00  2.1807E+00  2.0654E+00  2.0025E+00  1.9725E+01  2.4803E+00  1.3148E+00
             5.5589E+00
 PARAMETER:  2.7059E-01 -1.2566E+00  6.6694E+00  4.3432E-01  8.7964E-01  8.2531E-01  7.9442E-01  3.0819E+00  1.0084E+00  3.7367E-01
             1.8154E+00
 GRADIENT:   3.8547E+01 -1.5259E+01 -2.7965E+00 -8.3679E+01  1.8697E+01 -7.2041E-01  2.2330E+00  3.8483E+01  2.6712E+01  3.0549E+01
             9.0771E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1742.43884763323        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0705E+00  2.5269E-01  1.9419E+02  1.7719E+00  1.8598E+00  2.0809E+00  2.6949E+00  6.8530E+00  2.1584E+00  7.0140E-01
             5.2075E+00
 PARAMETER:  1.6816E-01 -1.2756E+00  5.3688E+00  6.7206E-01  7.2045E-01  8.3281E-01  1.0913E+00  2.0247E+00  8.6935E-01 -2.5468E-01
             1.7501E+00
 GRADIENT:  -1.8547E+01  7.3257E-01 -7.1588E-01  3.1664E+01 -1.0221E+02  1.1298E+01  2.9814E+00  2.5776E+01 -8.1511E+00 -1.3557E+01
            -9.9763E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1755.11182248650        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.1042E+00  3.3957E-01  1.3953E+02  1.7058E+00  1.9591E+00  2.0713E+00  1.8226E+00  5.8978E+00  2.2826E+00  8.1266E-01
             5.3690E+00
 PARAMETER:  1.9917E-01 -9.8006E-01  5.0383E+00  6.3402E-01  7.7247E-01  8.2820E-01  7.0029E-01  1.8746E+00  9.2530E-01 -1.0745E-01
             1.7806E+00
 GRADIENT:  -4.1108E+00 -2.5527E+00  1.3012E+00  1.4387E+01 -8.0429E+00  3.1704E+00  2.8734E+00 -1.7037E+00  1.0167E+00 -2.6422E-01
            -2.0593E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1760.11222383994        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.1185E+00  7.7708E-01  5.3498E+01  1.3378E+00  1.9959E+00  2.0430E+00  6.1271E-01  4.8218E+00  2.7007E+00  9.3476E-01
             5.3670E+00
 PARAMETER:  2.1199E-01 -1.5221E-01  4.0797E+00  3.9101E-01  7.9111E-01  8.1442E-01 -3.8987E-01  1.6732E+00  1.0935E+00  3.2530E-02
             1.7803E+00
 GRADIENT:   2.6794E+00 -6.6709E+00 -3.1024E+00 -1.3451E+01  2.5861E+01 -6.4310E+00  2.1326E+00 -4.0755E+00  5.2586E+00  3.7076E+00
             1.9743E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1762.09952435999        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0970E+00  7.4119E-01  8.3371E+01  1.4217E+00  1.9600E+00  2.0615E+00  5.8029E-01  5.8819E+00  2.6158E+00  9.3401E-01
             5.3521E+00
 PARAMETER:  1.9256E-01 -1.9950E-01  4.5233E+00  4.5188E-01  7.7294E-01  8.2343E-01 -4.4422E-01  1.8719E+00  1.0616E+00  3.1731E-02
             1.7775E+00
 GRADIENT:  -4.1226E+01 -2.4505E+00 -3.1059E+00 -2.7331E+01 -1.8515E+01 -6.6676E+01  1.5550E+00 -2.7410E-01 -6.0117E+01  3.2349E+00
            -3.9900E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1768.74060447954        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  1.0984E+00  7.4043E-01  8.2347E+01  1.4268E+00  1.9568E+00  2.3462E+00  6.4190E-02  5.9213E+00  2.6268E+00  9.3353E-01
             5.4222E+00
 PARAMETER:  1.9384E-01 -2.0052E-01  4.5109E+00  4.5546E-01  7.7130E-01  9.5282E-01 -2.6459E+00  1.8786E+00  1.0658E+00  3.1218E-02
             1.7905E+00
 GRADIENT:  -3.2101E+01 -1.7661E+00 -3.6081E+00 -2.5471E+01 -2.2512E+01 -7.7579E+00  2.7414E-02  1.1012E+00 -6.2292E+01  3.3646E+00
            -8.6438E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1772.61981705884        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      829
 NPARAMETR:  1.0995E+00  7.3971E-01  8.1032E+01  1.5559E+00  1.9502E+00  2.3093E+00  1.9731E-02  5.9623E+00  3.0906E+00  8.5971E-01
             5.3858E+00
 PARAMETER:  1.9483E-01 -2.0150E-01  4.4949E+00  5.4205E-01  7.6792E-01  9.3695E-01 -3.8256E+00  1.8854E+00  1.2284E+00 -5.1162E-02
             1.7838E+00
 GRADIENT:  -3.2152E+01 -1.2785E+01 -4.2106E+00  2.5368E+00 -3.6647E+01 -1.5727E+01  3.5023E-03  8.4845E+00  6.1602E+00 -2.1136E+00
            -3.2333E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1773.68846888195        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1005
 NPARAMETR:  1.1003E+00  7.3989E-01  8.3307E+01  1.5515E+00  1.9645E+00  2.3935E+00  1.0000E-02  5.9094E+00  3.0371E+00  8.6383E-01
             5.5034E+00
 PARAMETER:  1.9554E-01 -2.0125E-01  4.5225E+00  5.3920E-01  7.7523E-01  9.7275E-01 -5.9916E+00  1.8765E+00  1.2109E+00 -4.6381E-02
             1.8054E+00
 GRADIENT:  -3.1038E+01 -1.1769E+01 -3.6843E+00  1.5615E-01 -2.9086E+01 -1.0276E+00  0.0000E+00  5.8295E+00  4.1892E-01 -2.0204E-01
             1.7421E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1776.25767803543        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1167
 NPARAMETR:  1.2639E+00  8.7243E-01  8.5600E+01  1.4815E+00  2.0373E+00  2.3919E+00  1.0000E-02  5.8820E+00  3.0192E+00  8.5744E-01
             5.4573E+00
 PARAMETER:  3.3422E-01 -3.6471E-02  4.5497E+00  4.9308E-01  8.1164E-01  9.7210E-01 -7.8222E+00  1.8719E+00  1.2050E+00 -5.3809E-02
             1.7969E+00
 GRADIENT:   7.5111E+01  1.2256E+01 -2.0580E+00  3.5996E+01  1.1468E+01  7.1303E+01  0.0000E+00 -9.4165E-01  3.3192E+01 -1.5774E-01
             1.8434E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1777.60109658783        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1361             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2110E+00  8.4031E-01  8.5274E+01  1.4770E+00  2.0216E+00  2.4580E+00  1.0000E-02  5.8941E+00  3.1612E+00  8.5082E-01
             5.4507E+00
 PARAMETER:  2.9146E-01 -7.3985E-02  4.5459E+00  4.9001E-01  8.0387E-01  9.9937E-01 -7.8222E+00  1.8740E+00  1.2509E+00 -6.1557E-02
             1.7957E+00
 GRADIENT:   5.2341E+01 -4.4890E+00 -2.2522E+00  3.1826E+01  7.5888E+00  8.6575E+01  0.0000E+00  5.9004E-01  6.1641E+01 -2.7238E-02
             2.3447E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1777.91107223421        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1517            RESET HESSIAN, TYPE II
 NPARAMETR:  1.1956E+00  8.7021E-01  8.3938E+01  1.4586E+00  2.0241E+00  2.4093E+00  1.0000E-02  5.9674E+00  3.2083E+00  8.4944E-01
             5.4219E+00
 PARAMETER:  2.7861E-01 -3.9021E-02  4.5301E+00  4.7749E-01  8.0510E-01  9.7933E-01 -7.8222E+00  1.8863E+00  1.2657E+00 -6.3172E-02
             1.7904E+00
 GRADIENT:   4.6230E+01 -1.2886E+00 -2.7663E+00  3.2540E+01  7.0360E+00  7.8323E+01  0.0000E+00  1.8900E+00  6.1576E+01 -9.7613E-01
             1.0354E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1777.97440035928        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1682
 NPARAMETR:  1.1946E+00  8.7471E-01  8.3896E+01  1.4512E+00  2.0245E+00  2.4043E+00  1.0000E-02  5.9694E+00  3.2186E+00  8.5299E-01
             5.4217E+00
 PARAMETER:  2.7782E-01 -3.3861E-02  4.5296E+00  4.7242E-01  8.0533E-01  9.7728E-01 -7.8222E+00  1.8867E+00  1.2690E+00 -5.9012E-02
             1.7904E+00
 GRADIENT:   5.9349E-01 -3.5298E+00 -3.4373E+00  4.3635E+00 -1.2995E+00  2.3810E+00  0.0000E+00 -1.3638E+00 -5.7972E+00 -9.0028E-01
            -2.1945E+01

0ITERATION NO.:   71    OBJECTIVE VALUE:  -1777.97440035928        NO. OF FUNC. EVALS.:  32
 CUMULATIVE NO. OF FUNC. EVALS.:     1714
 NPARAMETR:  1.1979E+00  8.7384E-01  8.0391E+01  1.4512E+00  2.0085E+00  2.3812E+00  1.0000E-02  6.0831E+00  3.1795E+00  8.5214E-01
             5.3468E+00
 PARAMETER:  2.7782E-01 -3.3861E-02  4.5296E+00  4.7242E-01  8.0533E-01  9.7728E-01 -7.8222E+00  1.8867E+00  1.2690E+00 -5.9012E-02
             1.7904E+00
 GRADIENT:  -1.1541E+03  3.2005E+03  6.5952E+01  5.6449E+00  3.9212E+02  3.2715E+02  0.0000E+00 -8.2764E+01  2.4577E+02  3.2022E+03
             1.4975E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1714
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.5920E-03 -5.1740E-04 -1.1054E-02 -6.1830E-03 -1.5948E-02
 SE:             2.9658E-02  1.1198E-04  9.3371E-03  2.9592E-02  1.6728E-02
 N:                     100         100         100         100         100

 P VAL.:         8.7695E-01  3.8349E-06  2.3647E-01  8.3449E-01  3.4040E-01

 ETASHRINKSD(%)  6.4233E-01  9.9625E+01  6.8719E+01  8.6396E-01  4.3958E+01
 ETASHRINKVR(%)  1.2805E+00  9.9999E+01  9.0215E+01  1.7205E+00  6.8593E+01
 EBVSHRINKSD(%)  1.0486E+00  9.9737E+01  7.6004E+01  2.0346E+00  4.4483E+01
 EBVSHRINKVR(%)  2.0861E+00  9.9999E+01  9.4242E+01  4.0278E+00  6.9179E+01
 RELATIVEINF(%)  9.7851E+01  3.2382E-04  4.1428E+00  4.5620E+01  2.1803E+01
 EPSSHRINKSD(%)  1.0599E+01
 EPSSHRINKVR(%)  2.0075E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1777.9744003592807     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -123.88504059086995     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    50.85
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1777.974       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.19E+00  8.75E-01  8.39E+01  1.45E+00  2.02E+00  2.40E+00  1.00E-02  5.97E+00  3.22E+00  8.53E-01  5.42E+00
 


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
+        7.29E+04
 
 TH 2
+        4.89E+02  1.05E+06
 
 TH 3
+        5.35E-03 -3.92E-01  6.15E-02
 
 TH 4
+       -6.38E+01  2.88E+02  5.30E-02  1.71E+04
 
 TH 5
+        3.59E+02 -1.16E+02 -3.10E-01  4.89E+00  3.06E+03
 
 TH 6
+        5.78E+00 -2.04E+01 -5.83E-03  8.85E+00 -1.65E+00  1.48E+03
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.35E-01  1.32E+01 -1.61E-01 -1.68E+00  6.80E+00  2.27E-01  0.00E+00  6.96E+01
 
 TH 9
+        2.92E+01 -1.43E+02 -2.66E-02  5.13E+00 -4.34E+00 -9.18E-01  0.00E+00  9.71E-01  4.96E+02
 
 TH10
+        6.65E+02 -2.53E+03 -5.30E-01  2.43E+02 -1.26E+02 -2.17E+01  0.00E+00  1.77E+01 -1.11E+02  1.10E+06
 
 TH11
+       -4.37E+00 -2.35E+01 -1.68E-01 -1.78E+00 -1.26E+01  4.34E-01  0.00E+00  5.57E+00  4.01E-01 -1.40E+01  1.17E+02
 
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
 #CPUT: Total CPU Time in Seconds,       67.355
Stop Time:
Wed Sep 29 08:25:03 CDT 2021
