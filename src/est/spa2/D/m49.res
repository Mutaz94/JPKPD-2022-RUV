Thu Sep 30 09:10:40 CDT 2021
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
$DATA ../../../../data/spa2/D/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17497.7434000701        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7850E+02  2.6377E+02  1.5939E+01  3.3619E+02  2.6343E+02 -1.8835E+03 -7.7688E+02 -5.2393E+01 -9.1993E+02 -6.6177E+02
            -3.5203E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -664.757663817588        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3077E+00  1.5029E+00  9.4781E-01  1.5972E+00  1.0763E+00  3.0289E+00  3.4995E+00  1.0137E+00  2.5417E+00  1.4807E+00
             1.2457E+01
 PARAMETER:  3.6824E-01  5.0743E-01  4.6395E-02  5.6826E-01  1.7353E-01  1.2082E+00  1.3526E+00  1.1356E-01  1.0329E+00  4.9253E-01
             2.6223E+00
 GRADIENT:  -1.3268E+00 -1.8273E+01 -3.5069E+01  6.3238E+01  1.2152E+01  1.0783E+02 -1.1303E+01  2.5352E+00  4.4440E+01  2.6270E+01
             3.2306E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -715.383567590692        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.3440E+00  2.3828E+00  2.2409E+00  1.3982E+00  1.7419E+00  2.4656E+00  5.8345E+00  1.0406E+00  3.9437E+00  1.6548E+00
             1.2177E+01
 PARAMETER:  3.9564E-01  9.6829E-01  9.0687E-01  4.3515E-01  6.5499E-01  1.0024E+00  1.8638E+00  1.3978E-01  1.4721E+00  6.0366E-01
             2.5995E+00
 GRADIENT:   4.0604E+01  2.4465E+01 -2.7273E+01  2.3419E+01 -1.1579E+01  1.0530E+02  1.1275E+02  5.9142E-01  5.8821E+01  2.3352E+01
             3.0918E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -799.949866897512        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  9.9498E-01  1.0808E+00  8.0145E+00  1.2901E+00  2.4200E+00  1.6818E+00  4.4318E+00  6.2608E+00  1.5454E+00  1.2300E+00
             9.6677E+00
 PARAMETER:  9.4966E-02  1.7772E-01  2.1812E+00  3.5469E-01  9.8379E-01  6.1985E-01  1.5888E+00  1.9343E+00  5.3525E-01  3.0700E-01
             2.3688E+00
 GRADIENT:  -4.8250E+01 -1.9216E+01 -6.4571E+00 -8.3274E+00  7.7397E+00 -2.8406E+00 -4.9485E+00  4.0799E+00  4.6584E+00  1.0836E+01
             1.8511E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -814.076906932521        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0208E+00  1.3739E+00  7.6166E+00  1.0130E+00  1.9814E+00  1.6709E+00  4.0189E+00  1.5529E+00  1.3702E+00  7.1383E-01
             8.3580E+00
 PARAMETER:  1.2063E-01  4.1763E-01  2.1303E+00  1.1288E-01  7.8382E-01  6.1334E-01  1.4910E+00  5.4011E-01  4.1493E-01 -2.3711E-01
             2.2232E+00
 GRADIENT:  -1.9147E+00 -1.1525E+01 -2.9641E-01 -8.1219E+00 -6.9072E+00 -7.3089E+00 -6.9196E+00  1.2872E-01  1.3006E+00  5.6313E+00
             1.6113E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -836.725102108642        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  1.0554E+00  1.0748E+00  1.4927E+01  1.4475E+00  2.2826E+00  1.7925E+00  5.9432E+00  9.0718E-01  1.3963E+00  4.3445E-01
             8.6965E+00
 PARAMETER:  1.5391E-01  1.7209E-01  2.8032E+00  4.6986E-01  9.2532E-01  6.8360E-01  1.8822E+00  2.5837E-03  4.3380E-01 -7.3367E-01
             2.2629E+00
 GRADIENT:   7.9402E+00  6.4288E+00 -1.3647E+00  8.0886E+00  4.9421E+00  9.1328E+00 -1.3564E+00  2.3813E-02 -8.1282E+00  1.4620E+00
             1.3241E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -843.259719203320        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  1.0400E+00  4.6327E-01  4.9121E+01  1.8231E+00  2.1996E+00  1.7304E+00  7.1852E+00  9.7845E-01  1.7575E+00  2.0566E-01
             8.5392E+00
 PARAMETER:  1.3920E-01 -6.6944E-01  3.9943E+00  7.0056E-01  8.8830E-01  6.4838E-01  2.0720E+00  7.8218E-02  6.6392E-01 -1.4815E+00
             2.2447E+00
 GRADIENT:   3.1997E+00 -1.4131E-01 -4.4380E-01  2.2106E+00 -2.8983E+00 -5.0917E-01 -1.0967E-01  5.5157E-03 -3.9613E+00  3.5132E-01
            -5.8067E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -844.510901653883        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      788
 NPARAMETR:  1.0336E+00  4.0975E-01  2.9047E+02  1.8499E+00  2.2581E+00  1.7266E+00  7.8713E+00  9.7275E-01  1.8751E+00  6.1497E-02
             8.5702E+00
 PARAMETER:  1.3309E-01 -7.9222E-01  5.7715E+00  7.1513E-01  9.1453E-01  6.4618E-01  2.1632E+00  7.2374E-02  7.2867E-01 -2.6888E+00
             2.2483E+00
 GRADIENT:   7.3079E+00  4.9447E+00 -7.5498E-02  1.3575E+01  1.1519E+00  9.2137E+00  1.4535E+02  1.5316E-04  1.1655E+01  3.3175E-02
             3.1346E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -845.121869260417        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      967
 NPARAMETR:  1.0369E+00  3.3984E-01  1.2312E+04  1.9588E+00  2.2737E+00  1.7378E+00  8.0536E+00  9.6690E-01  1.8480E+00  1.0000E-02
             8.5623E+00
 PARAMETER:  1.3627E-01 -9.7927E-01  9.5183E+00  7.7234E-01  9.2139E-01  6.5262E-01  2.1861E+00  6.6343E-02  7.1409E-01 -4.6472E+00
             2.2474E+00
 GRADIENT:  -6.3368E-01  5.3378E-01 -1.8921E-03  5.5876E+00 -1.7992E+00 -5.0652E-02  9.2399E+00  6.2535E-06 -1.4979E+00  0.0000E+00
            -2.6560E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -845.326445687223        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1132             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0364E+00  2.9667E-01  3.7522E+04  1.9852E+00  2.2720E+00  1.7381E+00  8.3606E+00  1.0033E+00  1.8327E+00  1.0000E-02
             8.5566E+00
 PARAMETER:  1.3574E-01 -1.1151E+00  1.0633E+01  7.8574E-01  9.2066E-01  6.5279E-01  2.2235E+00  1.0329E-01  7.0581E-01 -4.6292E+00
             2.2467E+00
 GRADIENT:   4.3650E+00  4.2470E+00 -6.0865E-04  3.7073E+01 -8.0574E-01  1.2312E+01  1.5587E+02  5.6895E-05  3.6116E+00  0.0000E+00
             2.3434E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -845.454298697218        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1296
 NPARAMETR:  1.0364E+00  2.9814E-01  5.5732E+04  1.9574E+00  2.2931E+00  1.7384E+00  8.5472E+00  1.0023E+00  1.8597E+00  1.0000E-02
             8.5736E+00
 PARAMETER:  1.3579E-01 -1.1102E+00  1.1028E+01  7.7161E-01  9.2992E-01  6.5295E-01  2.2456E+00  1.0225E-01  7.2043E-01 -4.6292E+00
             2.2487E+00
 GRADIENT:  -1.6362E+01  9.1976E-01 -4.9023E-04  2.1841E+00  4.3681E-01  8.9086E+00  1.8094E+01 -4.0168E-05 -5.2068E+00  0.0000E+00
             4.1483E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -845.536708463070        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1465
 NPARAMETR:  1.0366E+00  2.7916E-01  7.2159E+04  1.9768E+00  2.2911E+00  1.7384E+00  8.6647E+00  1.0005E+00  1.8604E+00  1.0000E-02
             8.5732E+00
 PARAMETER:  1.3592E-01 -1.1760E+00  1.1287E+01  7.8149E-01  9.2902E-01  6.5298E-01  2.2593E+00  1.0050E-01  7.2082E-01 -4.6292E+00
             2.2486E+00
 GRADIENT:  -6.6363E+00  5.1454E+00 -3.6353E-04  3.3698E+01  1.2552E+00  1.4243E+01  1.6886E+02  9.5850E-05  2.9188E+00  0.0000E+00
             2.9425E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -845.605013237503        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1658
 NPARAMETR:  1.0370E+00  2.5943E-01  1.1326E+05  1.9927E+00  2.2894E+00  1.7398E+00  8.9078E+00  9.9886E-01  1.8608E+00  1.0000E-02
             8.5727E+00
 PARAMETER:  1.3631E-01 -1.2493E+00  1.1737E+01  7.8948E-01  9.2830E-01  6.5375E-01  2.2869E+00  9.8862E-02  7.2098E-01 -4.6292E+00
             2.2486E+00
 GRADIENT:  -8.5674E+01  2.4658E+00 -2.7031E-04  3.5418E+01 -1.3555E+00  2.6449E+01  1.9199E+01 -7.1263E-05 -2.5782E+01  0.0000E+00
             1.3313E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -845.631283680427        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:     1863             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0375E+00  2.4166E-01  1.5999E+05  2.0046E+00  2.2891E+00  1.7394E+00  8.9956E+00  1.0010E+00  1.8625E+00  1.0000E-02
             8.5720E+00
 PARAMETER:  1.3679E-01 -1.3202E+00  1.2083E+01  7.9542E-01  9.2817E-01  6.5353E-01  2.2967E+00  1.0104E-01  7.2191E-01 -4.6292E+00
             2.2485E+00
 GRADIENT:  -7.7843E+01  7.0074E+00 -1.6174E-04  6.4137E+01 -5.6033E-01  1.4340E+01  1.7969E+02  1.2211E-04 -2.5472E+00  0.0000E+00
             3.8647E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -845.639727189651        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2061
 NPARAMETR:  1.0377E+00  2.3858E-01  1.5980E+05  2.0073E+00  2.2899E+00  1.7395E+00  9.0556E+00  1.0000E+00  1.8622E+00  1.0000E-02
             8.5730E+00
 PARAMETER:  1.3702E-01 -1.3331E+00  1.2082E+01  7.9678E-01  9.2849E-01  6.5361E-01  2.3034E+00  1.0002E-01  7.2174E-01 -4.6292E+00
             2.2486E+00
 GRADIENT:  -2.1313E+02  4.9227E+00 -1.8909E-04  9.8520E+01 -3.3039E+00  8.1446E+01  1.6110E+01  2.2007E-05 -6.0050E+01  0.0000E+00
             3.0804E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -845.646208452329        NO. OF FUNC. EVALS.: 207
 CUMULATIVE NO. OF FUNC. EVALS.:     2268             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0377E+00  2.3547E-01  2.3554E+05  2.0113E+00  2.2899E+00  1.7396E+00  9.1316E+00  1.0046E+00  1.8619E+00  1.0000E-02
             8.5735E+00
 PARAMETER:  1.3696E-01 -1.3462E+00  1.2470E+01  7.9879E-01  9.2853E-01  6.5363E-01  2.3117E+00  1.0459E-01  7.2157E-01 -4.6292E+00
             2.2487E+00
 GRADIENT:  -7.8343E+01  5.8332E+00 -1.1488E-04  5.7432E+01 -3.2198E-02  3.6839E+01  1.8463E+02  1.0394E-04 -5.3147E+00  0.0000E+00
             4.2464E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -845.750301583968        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:     2385
 NPARAMETR:  1.0307E+00  2.2527E-01  1.8275E+08  2.0194E+00  2.2918E+00  1.7395E+00  9.1149E+00  1.0046E+00  1.8675E+00  1.0000E-02
             8.5786E+00
 PARAMETER:  1.3023E-01 -1.3905E+00  1.9124E+01  8.0281E-01  9.2934E-01  6.5362E-01  2.3099E+00  1.0454E-01  7.2462E-01 -4.6292E+00
             2.2493E+00
 GRADIENT:  -6.7322E+00  6.7530E+00 -2.1326E-04  4.4297E+01  5.4379E-01  6.6105E+00  1.8799E+02  6.6351E-02  1.0251E+01  0.0000E+00
             2.8805E+01

0ITERATION NO.:   83    OBJECTIVE VALUE:  -845.798117457451        NO. OF FUNC. EVALS.:  82
 CUMULATIVE NO. OF FUNC. EVALS.:     2467
 NPARAMETR:  1.0288E+00  2.2471E-01  2.6841E+08  2.0179E+00  2.3050E+00  1.7393E+00  9.2813E+00  1.0056E+00  1.8674E+00  1.0000E-02
             8.5838E+00
 PARAMETER:  1.2929E-01 -1.3943E+00  1.9315E+01  8.0311E-01  9.2934E-01  6.5381E-01  2.3096E+00  1.0454E-01  7.2455E-01 -4.6292E+00
             2.2493E+00
 GRADIENT:   1.9816E+02 -1.1654E+00 -1.3855E+02  2.5217E+00 -1.9024E+01  5.7310E-01 -1.3885E+01 -2.5578E+04 -1.5079E-02  0.0000E+00
            -3.8765E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2467
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3008E-02  6.6843E-02 -2.1570E-08 -8.3202E-02 -2.6804E-05
 SE:             2.7889E-02  1.9819E-02  2.6416E-08  1.8816E-02  5.7832E-05
 N:                     100         100         100         100         100

 P VAL.:         6.4093E-01  7.4431E-04  4.1419E-01  9.7929E-06  6.4301E-01

 ETASHRINKSD(%)  6.5675E+00  3.3605E+01  1.0000E+02  3.6964E+01  9.9806E+01
 ETASHRINKVR(%)  1.2704E+01  5.5917E+01  1.0000E+02  6.0265E+01  1.0000E+02
 EBVSHRINKSD(%)  7.7737E+00  4.4172E+01  1.0000E+02  2.2696E+01  9.9719E+01
 EBVSHRINKVR(%)  1.4943E+01  6.8833E+01  1.0000E+02  4.0241E+01  9.9999E+01
 RELATIVEINF(%)  8.2800E+01  2.0046E+01  2.9464E-10  3.7383E+01  7.5236E-04
 EPSSHRINKSD(%)  1.0865E+01
 EPSSHRINKVR(%)  2.0550E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -845.79811745745064     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       256.92812238815645     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    67.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    22.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -845.798       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.24E-01  2.21E+08  2.02E+00  2.29E+00  1.74E+00  9.11E+00  1.00E+00  1.87E+00  1.00E-02  8.58E+00
 


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
+        7.48E+03
 
 TH 2
+        1.79E+04  2.38E+03
 
 TH 3
+        1.66E-07  5.01E-05  4.55E-17
 
 TH 4
+       -2.93E+05 -1.46E+05 -2.22E-08  1.32E+02
 
 TH 5
+       -2.35E+05  2.55E+03 -1.52E-07 -2.24E+02  2.70E+02
 
 TH 6
+        2.84E+03  6.86E+02  9.83E-09 -8.77E+02 -2.77E+02  1.63E+02
 
 TH 7
+        1.10E+04  1.24E+02 -2.02E-09 -3.82E+00 -1.15E+00 -1.88E+01  2.22E+00
 
 TH 8
+        2.97E+04  2.04E+06 -1.93E-06  3.64E+05 -1.74E+04 -4.57E+05 -1.84E+02  6.08E+06
 
 TH 9
+        3.18E+01  5.70E+01  6.57E-08  1.35E+02  1.52E+02  8.62E+01 -4.18E+00 -6.82E+02  2.21E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.99E+02 -4.28E+01 -1.68E-09  2.21E+03  2.71E+00  4.16E+01 -1.20E-01 -1.65E+01  2.26E+01  0.00E+00  1.27E+01
 
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
 #CPUT: Total CPU Time in Seconds,       90.409
Stop Time:
Thu Sep 30 09:12:12 CDT 2021
