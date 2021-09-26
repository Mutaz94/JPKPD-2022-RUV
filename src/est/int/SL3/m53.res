Sat Sep 25 02:23:21 CDT 2021
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
$DATA ../../../../data/int/SL3/dat53.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      986
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

 TOT. NO. OF OBS RECS:      886
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -69.3423364097031        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.0521E+01 -6.9591E+00  1.9622E+02 -5.3658E+01  1.5556E+02 -3.0795E+01 -1.2906E+02 -6.6140E+02 -1.6111E+02 -7.8916E+01
            -6.4281E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2701.29850399625        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7658E-01  1.3189E+00  1.1637E+00  8.8529E-01  1.1268E+00  1.0811E+00  9.5986E-01  9.4591E-01  7.1458E-01  1.0424E+00
             2.7818E+00
 PARAMETER:  7.6298E-02  3.7676E-01  2.5164E-01 -2.1835E-02  2.1939E-01  1.7801E-01  5.9028E-02  4.4388E-02 -2.3605E-01  1.4153E-01
             1.1231E+00
 GRADIENT:  -2.7753E+01  2.7155E+01  2.1029E+01 -4.7715E+01 -8.3033E+01  4.6100E+00  8.7197E+00  6.5626E-01 -1.4087E+01 -2.4704E+01
             5.1853E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2707.90058135695        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.9345E-01  1.5167E+00  3.0734E+00  8.4598E-01  1.7170E+00  1.1259E+00  8.7406E-01  6.5974E-01  8.1709E-01  1.4831E+00
             2.8106E+00
 PARAMETER:  9.3431E-02  5.1655E-01  1.2228E+00 -6.7263E-02  6.4056E-01  2.1859E-01 -3.4607E-02 -3.1590E-01 -1.0201E-01  4.9414E-01
             1.1334E+00
 GRADIENT:   2.0057E+00  7.8694E+01 -5.5275E+00  7.5609E+01  1.3324E+01  1.8850E+01  1.0927E+01 -7.7009E-01  1.4591E+00 -1.5229E+01
             8.5012E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2716.78099995219        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.8850E-01  1.4545E+00  4.5243E+00  8.2882E-01  1.8270E+00  1.0822E+00  8.3141E-01  2.0948E+00  7.6537E-01  1.6348E+00
             2.6631E+00
 PARAMETER:  8.8436E-02  4.7464E-01  1.6095E+00 -8.7750E-02  7.0266E-01  1.7903E-01 -8.4632E-02  8.3945E-01 -1.6739E-01  5.9150E-01
             1.0795E+00
 GRADIENT:  -1.2636E+00 -1.9118E+01 -5.7800E+00 -8.6172E+00  1.3839E+01  4.8921E+00  1.2682E+00 -9.2623E-01 -6.8241E-01 -1.4052E+00
            -1.2054E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2718.75519144328        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.8992E-01  1.4067E+00  9.0645E+00  8.8092E-01  1.9362E+00  1.0683E+00  7.8704E-01  2.9825E+00  8.1244E-01  1.7330E+00
             2.6781E+00
 PARAMETER:  8.9868E-02  4.4123E-01  2.3044E+00 -2.6793E-02  7.6071E-01  1.6608E-01 -1.3948E-01  1.1928E+00 -1.0772E-01  6.4986E-01
             1.0851E+00
 GRADIENT:   5.2200E-02  2.4719E+00 -3.8840E-01 -1.2297E-01  2.2896E+00 -3.1318E-01 -7.1986E-01 -5.7476E-01 -5.1273E-01  4.5706E-01
             4.5769E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2718.87067619960        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  9.8981E-01  1.2207E+00  1.6737E+01  1.0157E+00  1.9705E+00  1.0696E+00  8.6464E-01  3.9435E+00  7.5123E-01  1.7386E+00
             2.6702E+00
 PARAMETER:  8.9754E-02  2.9946E-01  2.9176E+00  1.1558E-01  7.7831E-01  1.6727E-01 -4.5441E-02  1.4721E+00 -1.8604E-01  6.5310E-01
             1.0822E+00
 GRADIENT:  -1.7731E-01  3.9465E+00  2.3876E-01  2.4968E+00 -1.3785E+00 -2.2188E-02  2.0025E-01  1.6310E-01 -2.1057E-01 -1.6546E-01
            -1.5995E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2718.91868346741        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      557
 NPARAMETR:  9.9409E-01  1.2776E+00  1.5646E+01  9.7718E-01  1.9883E+00  1.0753E+00  8.3269E-01  3.9238E+00  7.7767E-01  1.7456E+00
             2.6758E+00
 PARAMETER:  9.4068E-02  3.4502E-01  2.8502E+00  7.6916E-02  7.8726E-01  1.7257E-01 -8.3097E-02  1.4671E+00 -1.5145E-01  6.5710E-01
             1.0843E+00
 GRADIENT:   1.5673E+00  1.2426E+00 -1.7068E-01  2.0171E+00  5.4136E-01  1.8106E-01  8.1408E-02  5.5495E-02  1.8409E-01 -3.3736E-01
             1.3093E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2718.93710482631        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  9.9364E-01  1.3795E+00  1.3146E+01  9.0571E-01  1.9858E+00  1.0749E+00  7.8973E-01  3.7594E+00  8.1595E-01  1.7523E+00
             2.6737E+00
 PARAMETER:  9.3621E-02  4.2174E-01  2.6761E+00  9.6692E-04  7.8602E-01  1.7222E-01 -1.3606E-01  1.4243E+00 -1.0340E-01  6.6092E-01
             1.0834E+00
 GRADIENT:   7.7804E-01  1.4291E+00 -2.5806E-01  2.3215E+00  3.0474E-01  4.6217E-02 -3.4497E-01  1.9331E-01  1.2736E-01  2.1163E-02
            -6.3399E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2718.97728384218        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      910
 NPARAMETR:  9.9215E-01  1.5449E+00  1.0607E+01  7.9100E-01  1.9889E+00  1.0740E+00  7.5478E-01  3.6771E+00  8.5301E-01  1.7630E+00
             2.6748E+00
 PARAMETER:  9.2117E-02  5.3498E-01  2.4615E+00 -1.3446E-01  7.8760E-01  1.7142E-01 -1.8133E-01  1.4021E+00 -5.8979E-02  6.6701E-01
             1.0839E+00
 GRADIENT:  -2.1751E+00  3.8678E+00  2.2364E-01  2.6865E+00 -1.0910E+00 -2.9828E-01 -5.5275E-01 -8.1146E-02 -6.6477E-01  1.1297E-01
             2.4481E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2719.07843695221        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1087
 NPARAMETR:  9.9326E-01  1.7324E+00  6.6430E+00  6.5665E-01  1.9801E+00  1.0745E+00  7.1289E-01  3.2310E+00  9.5687E-01  1.7651E+00
             2.6728E+00
 PARAMETER:  9.3232E-02  6.4948E-01  1.9936E+00 -3.2060E-01  7.8314E-01  1.7183E-01 -2.3842E-01  1.2728E+00  5.5914E-02  6.6819E-01
             1.0831E+00
 GRADIENT:   7.3309E-02 -1.7190E+00 -4.2820E-01  5.1154E-01  3.8142E-01 -1.1247E-01  1.0731E+00  5.0582E-01  1.1080E+00  6.8486E-02
            -1.0191E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2719.17773153440        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1262
 NPARAMETR:  9.9327E-01  1.9174E+00  4.4032E+00  5.3166E-01  1.9792E+00  1.0748E+00  6.7927E-01  2.8705E+00  1.0249E+00  1.7721E+00
             2.6747E+00
 PARAMETER:  9.3245E-02  7.5097E-01  1.5823E+00 -5.3176E-01  7.8268E-01  1.7214E-01 -2.8674E-01  1.1545E+00  1.2457E-01  6.7215E-01
             1.0838E+00
 GRADIENT:  -2.6739E-01  9.4553E+00  2.1074E-01  3.7780E+00 -3.9398E-01 -9.3942E-02 -1.8942E-01 -5.3650E-01 -5.9658E-01  2.1639E-01
            -2.6346E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2719.28836631744        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1438
 NPARAMETR:  9.9338E-01  2.0738E+00  3.1481E+00  4.2069E-01  1.9936E+00  1.0755E+00  6.4864E-01  2.8166E+00  1.1847E+00  1.7806E+00
             2.6755E+00
 PARAMETER:  9.3360E-02  8.2939E-01  1.2468E+00 -7.6586E-01  7.8997E-01  1.7278E-01 -3.3288E-01  1.1355E+00  2.6945E-01  6.7698E-01
             1.0842E+00
 GRADIENT:  -1.0891E-01  3.9548E+00  9.0150E-02  1.4381E+00 -3.5286E-01  1.1522E-01  4.0789E-01 -4.3949E-02  3.6165E-01 -2.2424E-01
             7.2076E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2719.34702302084        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1613
 NPARAMETR:  9.9365E-01  2.2306E+00  1.9291E+00  3.1211E-01  2.0053E+00  1.0751E+00  6.2552E-01  2.7711E+00  1.3228E+00  1.7906E+00
             2.6735E+00
 PARAMETER:  9.3629E-02  9.0227E-01  7.5706E-01 -1.0644E+00  7.9580E-01  1.7239E-01 -3.6918E-01  1.1193E+00  3.7971E-01  6.8256E-01
             1.0834E+00
 GRADIENT:   2.4019E-01  5.5304E+00 -5.5034E-02  1.2131E+00 -2.2907E-01 -1.0845E-01 -5.9003E-01  1.2415E-01 -3.2096E-01 -7.2004E-02
            -1.4250E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2719.37084484293        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1789
 NPARAMETR:  9.9324E-01  2.2994E+00  1.3419E+00  2.6204E-01  2.0008E+00  1.0754E+00  6.1563E-01  2.5433E+00  1.4701E+00  1.7877E+00
             2.6763E+00
 PARAMETER:  9.3219E-02  9.3263E-01  3.9409E-01 -1.2393E+00  7.9357E-01  1.7270E-01 -3.8510E-01  1.0335E+00  4.8532E-01  6.8090E-01
             1.0844E+00
 GRADIENT:  -5.8496E-01 -5.4955E-01 -2.6088E-01  5.4333E-01 -5.6824E-01  1.3711E-02  1.8122E-01  5.2718E-01  5.4682E-01 -1.2685E-01
             7.8874E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2719.47200076580        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1968
 NPARAMETR:  9.9465E-01  2.4025E+00  6.1660E-01  1.8979E-01  1.9880E+00  1.0763E+00  6.0259E-01  1.7272E+00  1.6867E+00  1.7788E+00
             2.6733E+00
 PARAMETER:  9.4638E-02  9.7651E-01 -3.8354E-01 -1.5618E+00  7.8714E-01  1.7354E-01 -4.0652E-01  6.4648E-01  6.2276E-01  6.7594E-01
             1.0833E+00
 GRADIENT:   1.8920E+00 -5.6495E-01 -2.0839E-02 -2.4426E-01 -6.9569E-01  2.5967E-01  4.3078E-01  2.9166E-01  1.4331E-01 -5.5565E-01
            -1.8239E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2719.55960184560        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2146
 NPARAMETR:  9.8913E-01  2.3203E+00  8.0790E-01  2.5258E-01  1.9537E+00  1.0719E+00  6.0933E-01  1.0951E+00  1.5611E+00  1.7665E+00
             2.6841E+00
 PARAMETER:  8.9075E-02  9.4168E-01 -1.1332E-01 -1.2760E+00  7.6973E-01  1.6947E-01 -3.9540E-01  1.9088E-01  5.4536E-01  6.6900E-01
             1.0873E+00
 GRADIENT:  -9.2152E+00  1.9496E+01  9.3854E-02  2.4308E+00  1.5861E+00 -1.4121E+00  7.5404E-01  1.0243E-01 -1.6803E+00  2.5052E+00
             8.0523E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2720.10085353369        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2324
 NPARAMETR:  9.9624E-01  2.2600E+00  6.5662E-01  2.8198E-01  1.8938E+00  1.0768E+00  6.1038E-01  4.6886E-01  1.6017E+00  1.7123E+00
             2.6703E+00
 PARAMETER:  9.6229E-02  9.1538E-01 -3.2065E-01 -1.1659E+00  7.3858E-01  1.7397E-01 -3.9367E-01 -6.5746E-01  5.7107E-01  6.3782E-01
             1.0822E+00
 GRADIENT:   4.4748E+00 -6.9055E+00 -1.0993E+00 -1.2823E+00  8.6397E+00  3.4921E-01  3.5721E+00  1.5700E-02 -5.3476E+00  2.7098E+00
             1.4769E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2721.43804888247        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2501
 NPARAMETR:  9.9132E-01  2.3378E+00  5.7061E-01  2.2841E-01  1.9061E+00  1.0750E+00  5.6739E-01  1.3646E-01  2.2540E+00  1.7017E+00
             2.6668E+00
 PARAMETER:  9.1285E-02  9.4921E-01 -4.6106E-01 -1.3766E+00  7.4508E-01  1.7233E-01 -4.6671E-01 -1.8918E+00  9.1271E-01  6.3160E-01
             1.0809E+00
 GRADIENT:  -4.8264E+00 -1.6678E+00 -4.6378E-01  3.6158E+00 -1.2166E+00 -3.8453E-01 -2.1403E+00  8.3091E-03  2.3111E+00  6.8897E-01
             5.2490E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2722.33182465714        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2677
 NPARAMETR:  9.9320E-01  2.5042E+00  3.7444E-01  1.1467E-01  2.0254E+00  1.0766E+00  5.6056E-01  1.4132E-01  3.3374E+00  1.7590E+00
             2.6553E+00
 PARAMETER:  9.3174E-02  1.0180E+00 -8.8233E-01 -2.0657E+00  8.0579E-01  1.7385E-01 -4.7882E-01 -1.8567E+00  1.3052E+00  6.6476E-01
             1.0766E+00
 GRADIENT:  -1.2181E+00 -3.0669E+00 -6.7685E-01  2.2436E+00  1.6615E+00  1.3670E-02 -2.6810E+00  1.4207E-02  3.2381E+00  9.8882E-02
             1.5064E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2722.42959488696        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2852
 NPARAMETR:  9.9396E-01  2.5413E+00  3.3664E-01  8.8327E-02  2.0419E+00  1.0767E+00  5.6433E-01  1.5115E-01  3.7838E+00  1.7694E+00
             2.6516E+00
 PARAMETER:  9.3943E-02  1.0327E+00 -9.8875E-01 -2.3267E+00  8.1390E-01  1.7394E-01 -4.7212E-01 -1.7895E+00  1.4307E+00  6.7064E-01
             1.0752E+00
 GRADIENT:   3.0018E-01  2.4776E+00  1.9141E-01 -1.5078E-01 -2.0874E-01 -1.0088E-03  1.1862E-01  1.5362E-02 -3.1599E-01 -3.2106E-02
            -3.3080E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2722.43287552967        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     3031
 NPARAMETR:  9.9390E-01  2.5387E+00  3.3589E-01  8.9217E-02  2.0420E+00  1.0768E+00  5.6428E-01  1.1912E-01  3.7541E+00  1.7701E+00
             2.6522E+00
 PARAMETER:  9.3878E-02  1.0316E+00 -9.9096E-01 -2.3167E+00  8.1393E-01  1.7401E-01 -4.7221E-01 -2.0276E+00  1.4229E+00  6.7104E-01
             1.0754E+00
 GRADIENT:   2.0125E-01 -1.4009E+00 -4.2033E-03 -1.2483E-01  1.1179E-01  4.0474E-02  2.5405E-02  9.6253E-03 -1.6388E-01  6.9073E-02
             1.4109E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2722.43815489307        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3206
 NPARAMETR:  9.9379E-01  2.5413E+00  3.3392E-01  8.7895E-02  2.0429E+00  1.0767E+00  5.6430E-01  1.0895E-02  3.7914E+00  1.7701E+00
             2.6520E+00
 PARAMETER:  9.3776E-02  1.0327E+00 -9.9685E-01 -2.3316E+00  8.1439E-01  1.7393E-01 -4.7216E-01 -4.4194E+00  1.4327E+00  6.7103E-01
             1.0753E+00
 GRADIENT:  -5.5506E-03  3.8811E-02 -7.0152E-03  1.2737E-02 -8.0135E-03 -2.0315E-03 -1.0506E-03  8.1035E-05  1.4769E-02 -2.9075E-03
            -7.0986E-03

0ITERATION NO.:  107    OBJECTIVE VALUE:  -2722.43816111082        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     3263
 NPARAMETR:  9.9379E-01  2.5413E+00  3.3396E-01  8.7877E-02  2.0430E+00  1.0767E+00  5.6430E-01  1.0000E-02  3.7920E+00  1.7701E+00
             2.6520E+00
 PARAMETER:  9.3775E-02  1.0327E+00 -9.9675E-01 -2.3318E+00  8.1440E-01  1.7393E-01 -4.7217E-01 -4.6442E+00  1.4329E+00  6.7102E-01
             1.0753E+00
 GRADIENT:  -5.8019E-03  4.8670E-02 -2.5143E-03  8.6372E-03 -7.5880E-03 -1.9981E-03 -1.7296E-03  0.0000E+00  1.1378E-02 -3.2529E-03
            -8.0094E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3263
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4805E-03 -2.7613E-02 -4.1565E-06  3.5134E-02 -2.0627E-02
 SE:             2.9459E-02  2.5102E-02  1.1609E-05  1.4756E-02  2.6531E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5992E-01  2.7131E-01  7.2030E-01  1.7266E-02  4.3688E-01

 ETASHRINKSD(%)  1.3092E+00  1.5906E+01  9.9961E+01  5.0565E+01  1.1118E+01
 ETASHRINKVR(%)  2.6012E+00  2.9281E+01  1.0000E+02  7.5562E+01  2.1000E+01
 EBVSHRINKSD(%)  1.3999E+00  1.3158E+01  9.9935E+01  5.9996E+01  8.3813E+00
 EBVSHRINKVR(%)  2.7802E+00  2.4585E+01  1.0000E+02  8.3997E+01  1.6060E+01
 RELATIVEINF(%)  9.7164E+01  2.2512E+01  3.2852E-05  4.6706E+00  6.8251E+01
 EPSSHRINKSD(%)  1.5998E+01
 EPSSHRINKVR(%)  2.9437E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          886
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1628.3590808386800     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2722.4381611108156     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1094.0790802721356     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    79.97
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2722.438       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  2.54E+00  3.34E-01  8.79E-02  2.04E+00  1.08E+00  5.64E-01  1.00E-02  3.79E+00  1.77E+00  2.65E+00
 


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
+        9.53E+02
 
 TH 2
+       -1.37E+01  5.23E+02
 
 TH 3
+       -7.66E-01  1.41E+02  2.05E+02
 
 TH 4
+       -1.69E+01 -3.10E+02 -1.10E+03  8.85E+03
 
 TH 5
+       -1.74E+00 -5.95E+00  1.03E+01 -6.85E+01  7.33E+01
 
 TH 6
+        4.42E+00 -4.21E+00  1.23E-01 -1.02E+01 -6.89E-01  1.63E+02
 
 TH 7
+        3.03E+00  7.48E+01  9.11E+01 -6.98E+02  9.38E+00 -1.17E+00  4.78E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.80E-01 -3.24E+01 -3.79E+01  2.98E+02 -4.11E+00 -2.20E-01 -2.61E+01  0.00E+00  1.15E+01
 
 TH10
+       -4.86E-01  2.63E+00  1.03E+01 -3.54E+01 -4.70E+00  7.66E-02  4.25E+00  0.00E+00 -1.29E+00  4.21E+01
 
 TH11
+       -1.24E+01 -3.59E+00  2.11E+01 -1.47E+02  2.80E+00  2.10E+00  2.46E+01  0.00E+00 -3.71E+00  5.65E+00  1.68E+02
 
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
 #CPUT: Total CPU Time in Seconds,       93.473
Stop Time:
Sat Sep 25 02:24:56 CDT 2021
