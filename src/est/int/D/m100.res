Wed Sep 29 10:23:47 CDT 2021
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
$DATA ../../../../data/int/D/dat100.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   65847.9925542920        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4395E+03  1.0795E+03  2.2646E+01  1.1408E+03 -7.4912E+01 -3.5402E+03 -2.1632E+03 -7.2575E+01 -2.9241E+03 -7.3054E+02
            -1.2804E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -275.328361226817        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.3504E-01  2.3068E+00  9.0834E-01  1.5774E+00  9.6784E-01  4.7435E+00  4.6457E+00  9.7929E-01  1.6782E+00  1.1782E+00
             1.3397E+01
 PARAMETER:  3.2838E-02  9.3586E-01  3.8581E-03  5.5579E-01  6.7312E-02  1.6568E+00  1.6359E+00  7.9072E-02  6.1773E-01  2.6399E-01
             2.6950E+00
 GRADIENT:  -2.1863E+01  6.6491E+01 -3.0233E+01  1.5244E+02 -2.8820E+01  2.0549E+02  2.1622E+01  4.4435E+00 -8.8365E-01  2.3194E+01
            -1.1895E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -345.440673301505        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  8.8442E-01  4.2135E+00  9.3615E+00  1.9237E+00  2.1276E+00  2.9165E+00  9.6637E+00  7.4417E-01  1.3068E+00  1.2528E+00
             1.3804E+01
 PARAMETER: -2.2824E-02  1.5383E+00  2.3366E+00  7.5426E-01  8.5502E-01  1.1704E+00  2.3684E+00 -1.9548E-01  3.6760E-01  3.2535E-01
             2.7249E+00
 GRADIENT:  -3.6837E+01  5.3758E+01 -9.2474E+00  1.4229E+02 -7.9637E+00  1.1064E+02  1.9994E+01  2.8329E-01 -2.2412E+01  1.9855E+01
            -7.8024E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -450.595195007463        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0695E+00  1.9796E+00  5.9652E+00  8.6901E-01  2.3066E+00  1.7781E+00  3.8293E+00  7.2474E-01  6.9432E-01  5.2235E-01
             1.4630E+01
 PARAMETER:  1.6718E-01  7.8288E-01  1.8859E+00 -4.0395E-02  9.3576E-01  6.7556E-01  1.4427E+00 -2.2195E-01 -2.6483E-01 -5.4943E-01
             2.7831E+00
 GRADIENT:   1.3606E+01 -1.0083E+00 -2.5746E+00 -7.2276E+00  1.2552E+01  2.4723E+01 -1.8411E+01  1.3030E-01  4.4339E+00  4.1261E+00
             8.0420E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -459.387307268291        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0214E+00  1.6275E+00  8.1936E+00  1.0198E+00  2.2358E+00  1.5958E+00  4.3896E+00  6.1050E-01  7.1325E-01  3.4464E-01
             1.4412E+01
 PARAMETER:  1.2113E-01  5.8706E-01  2.2034E+00  1.1961E-01  9.0462E-01  5.6737E-01  1.5792E+00 -3.9348E-01 -2.3792E-01 -9.6525E-01
             2.7681E+00
 GRADIENT:  -3.8030E+00 -7.9584E-01 -1.2028E+00 -2.8250E+00  2.1512E+00  2.7065E-01  1.7508E+00  7.5680E-02  4.3834E+00  1.7589E+00
             5.9323E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -460.531459857003        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      395
 NPARAMETR:  1.0047E+00  1.4623E+00  9.6761E+00  1.0047E+00  2.2630E+00  1.6239E+00  4.4714E+00  4.3640E-01  5.0731E-01  1.5753E-01
             1.3963E+01
 PARAMETER:  1.0472E-01  4.8001E-01  2.3697E+00  1.0471E-01  9.1667E-01  5.8481E-01  1.5977E+00 -7.2919E-01 -5.7863E-01 -1.7481E+00
             2.7364E+00
 GRADIENT:  -7.5115E+00 -1.1507E+01 -4.1521E-01 -4.9961E+00 -5.9148E-01 -5.2930E-01 -3.9928E+01  2.8308E-02  1.3308E+00  3.5263E-01
            -2.9566E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -466.099955715152        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      572
 NPARAMETR:  1.0280E+00  1.2518E+00  1.4858E+01  1.1787E+00  2.3753E+00  1.6085E+00  5.7261E+00  4.9716E-01  5.0942E-01  1.5792E-01
             1.4459E+01
 PARAMETER:  1.2762E-01  3.2462E-01  2.7986E+00  2.6440E-01  9.6514E-01  5.7532E-01  1.8450E+00 -5.9885E-01 -5.7448E-01 -1.7457E+00
             2.7713E+00
 GRADIENT:  -1.8242E+00 -2.3831E+00 -8.6039E-01 -6.4270E-01  2.9906E+00  4.1228E-01  1.0717E+00  2.2873E-02 -8.8481E-01  3.1621E-01
             6.8104E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -466.250375071424        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      750
 NPARAMETR:  1.0353E+00  1.3855E+00  1.8036E+01  1.1566E+00  2.4075E+00  1.6058E+00  5.5475E+00  4.2089E-01  4.8610E-01  8.3882E-02
             1.4558E+01
 PARAMETER:  1.3471E-01  4.2608E-01  2.9924E+00  2.4548E-01  9.7860E-01  5.7364E-01  1.8133E+00 -7.6538E-01 -6.2134E-01 -2.3783E+00
             2.7781E+00
 GRADIENT:   2.2358E-01  5.4126E-01 -3.5506E-01  1.9644E+00  7.8811E-01  2.9108E-02 -1.0952E+00  9.8624E-03 -7.0975E-01  8.8283E-02
             5.2757E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -466.789315671163        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      928
 NPARAMETR:  1.0315E+00  1.2035E+00  5.8345E+01  1.2551E+00  2.5055E+00  1.6031E+00  5.8638E+00  4.9733E-01  7.2779E-01  2.2263E-02
             1.4477E+01
 PARAMETER:  1.3104E-01  2.8519E-01  4.1664E+00  3.2718E-01  1.0185E+00  5.7195E-01  1.8688E+00 -5.9850E-01 -2.1774E-01 -3.7048E+00
             2.7725E+00
 GRADIENT:   1.8631E-01 -1.0719E-01 -5.1618E-02 -1.2507E+00  1.9765E-01 -2.4731E-01  2.1015E-01  1.3317E-03  4.2657E-01  6.1431E-03
            -8.2611E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -466.825954050047        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1103
 NPARAMETR:  1.0315E+00  1.1958E+00  1.4083E+02  1.2639E+00  2.5429E+00  1.6030E+00  5.8889E+00  3.8524E-01  7.3246E-01  1.0000E-02
             1.4488E+01
 PARAMETER:  1.3098E-01  2.7880E-01  5.0475E+00  3.3418E-01  1.0333E+00  5.7190E-01  1.8731E+00 -8.5388E-01 -2.1135E-01 -5.4228E+00
             2.7733E+00
 GRADIENT:   2.9202E-02  1.6109E-02 -5.7740E-03 -3.7462E-01  1.8096E-01 -5.8268E-02  4.7677E-02  1.3490E-04  1.3672E-01  0.0000E+00
             1.7024E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -466.844263420710        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1281
 NPARAMETR:  1.0313E+00  1.1773E+00  3.2828E+02  1.2769E+00  2.5719E+00  1.6025E+00  5.9274E+00  3.0869E-01  7.5463E-01  1.0000E-02
             1.4489E+01
 PARAMETER:  1.3084E-01  2.6319E-01  5.8939E+00  3.4446E-01  1.0446E+00  5.7153E-01  1.8796E+00 -1.0754E+00 -1.8153E-01 -7.0289E+00
             2.7734E+00
 GRADIENT:   3.9223E-03 -6.3010E-02 -7.5324E-03 -6.4864E-01  1.6967E+00 -4.8595E-02  1.8843E-01  1.6437E-05  2.5755E-01  0.0000E+00
             2.8432E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -466.895776263217        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1463
 NPARAMETR:  1.0317E+00  1.1536E+00  7.5630E+02  1.2854E+00  2.5652E+00  1.6031E+00  6.0366E+00  3.1197E-01  7.4652E-01  1.0000E-02
             1.4514E+01
 PARAMETER:  1.3122E-01  2.4286E-01  6.7284E+00  3.5111E-01  1.0421E+00  5.7197E-01  1.8978E+00 -1.0648E+00 -1.9233E-01 -7.0289E+00
             2.7751E+00
 GRADIENT:  -9.3318E-02 -1.2099E-02 -3.1418E-04 -1.0258E+00  3.0201E-01  2.7766E-01  2.6781E+00  3.9748E-06 -1.1885E-01  0.0000E+00
             2.4850E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -466.910559424627        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1643
 NPARAMETR:  1.0320E+00  1.1498E+00  2.1741E+03  1.2910E+00  2.5630E+00  1.6018E+00  6.0801E+00  2.9923E-01  7.6282E-01  1.0000E-02
             1.4490E+01
 PARAMETER:  1.3148E-01  2.3961E-01  7.7844E+00  3.5539E-01  1.0412E+00  5.7110E-01  1.9050E+00 -1.1066E+00 -1.7073E-01 -7.0289E+00
             2.7735E+00
 GRADIENT:   6.2674E-01  4.4597E-01  1.7267E-04 -1.4060E+00 -2.1585E-01 -5.3521E-02  3.9122E+00  1.5808E-06  8.3497E-02  0.0000E+00
            -2.5353E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -466.912176035487        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1818
 NPARAMETR:  1.0321E+00  1.1364E+00  2.2061E+03  1.2995E+00  2.5675E+00  1.6028E+00  6.0624E+00  2.9925E-01  7.7288E-01  1.0000E-02
             1.4517E+01
 PARAMETER:  1.3162E-01  2.2789E-01  7.7990E+00  3.6199E-01  1.0429E+00  5.7173E-01  1.9021E+00 -1.1065E+00 -1.5763E-01 -7.0289E+00
             2.7753E+00
 GRADIENT:  -3.5956E-02  6.8911E-02 -1.7218E-04 -3.3887E-01  1.2851E-01  2.1798E-01  2.2365E+00  1.3049E-06 -4.3867E-02  0.0000E+00
             2.4003E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -466.931059873642        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1997
 NPARAMETR:  1.0311E+00  1.1214E+00  3.5303E+03  1.3066E+00  2.5665E+00  1.6036E+00  6.1585E+00  2.9900E-01  7.7749E-01  1.0000E-02
             1.4507E+01
 PARAMETER:  1.3064E-01  2.1462E-01  8.2691E+00  3.6740E-01  1.0426E+00  5.7228E-01  1.9178E+00 -1.1073E+00 -1.5168E-01 -7.0289E+00
             2.7746E+00
 GRADIENT:  -2.3795E-01  3.7366E-01 -1.3870E-04 -1.2153E+00  2.2269E-02  3.9780E-01  4.7411E+00  1.9738E-06 -5.3892E-02  0.0000E+00
             1.7104E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -466.932998831985        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2172
 NPARAMETR:  1.0326E+00  1.1118E+00  4.0933E+03  1.3145E+00  2.5685E+00  1.6022E+00  6.1306E+00  2.9899E-01  7.8920E-01  1.0000E-02
             1.4525E+01
 PARAMETER:  1.3212E-01  2.0597E-01  8.4171E+00  3.7342E-01  1.0433E+00  5.7135E-01  1.9133E+00 -1.1073E+00 -1.3674E-01 -7.0289E+00
             2.7759E+00
 GRADIENT:   1.2613E-01  1.0816E-01 -1.9662E-04  8.7033E-02  5.8831E-02  2.0929E-01  2.8974E+00  5.2919E-07 -1.6505E-01  0.0000E+00
             2.6010E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -466.949006426014        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2358
 NPARAMETR:  1.0309E+00  1.0922E+00  5.9647E+03  1.3187E+00  2.5653E+00  1.6012E+00  6.1815E+00  2.9324E-01  7.9848E-01  1.0000E-02
             1.4477E+01
 PARAMETER:  1.3043E-01  1.8823E-01  8.7936E+00  3.7665E-01  1.0421E+00  5.7075E-01  1.9216E+00 -1.1268E+00 -1.2504E-01 -7.0289E+00
             2.7725E+00
 GRADIENT:   2.0894E-01  6.6931E-03 -1.2119E-04 -6.2308E-01 -1.0841E-01 -4.9207E-02  3.8523E+00  3.3397E-08 -1.2288E-01  0.0000E+00
            -1.2871E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -466.953941796722        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2543             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0314E+00  1.0806E+00  1.4499E+04  1.3244E+00  2.5672E+00  1.6005E+00  6.2062E+00  3.0685E-01  8.0796E-01  1.0000E-02
             1.4482E+01
 PARAMETER:  1.3089E-01  1.7752E-01  9.6818E+00  3.8094E-01  1.0428E+00  5.7032E-01  1.9255E+00 -1.0814E+00 -1.1324E-01 -7.0289E+00
             2.7729E+00
 GRADIENT:   4.3177E+00  9.3458E-01 -8.1912E-05  5.6398E+00  2.2602E+00  3.8051E+00  6.7974E+01 -9.1442E-07  8.4886E-02  0.0000E+00
             4.9422E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -466.956910252556        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2723
 NPARAMETR:  1.0308E+00  1.0824E+00  1.4124E+04  1.3292E+00  2.5665E+00  1.6023E+00  6.2211E+00  3.0324E-01  8.1454E-01  1.0000E-02
             1.4490E+01
 PARAMETER:  1.3037E-01  1.7921E-01  9.6556E+00  3.8457E-01  1.0425E+00  5.7142E-01  1.9280E+00 -1.0932E+00 -1.0513E-01 -7.0289E+00
             2.7734E+00
 GRADIENT:  -2.1321E-01  2.2926E-01 -8.6356E-05 -3.7057E-01 -9.5078E-02  2.0121E-01  4.3184E+00  4.4051E-07 -5.5225E-02  0.0000E+00
            -1.6188E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -466.959046350040        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     2918             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0316E+00  1.0709E+00  1.7171E+04  1.3329E+00  2.5665E+00  1.6002E+00  6.2287E+00  2.9689E-01  8.1685E-01  1.0000E-02
             1.4477E+01
 PARAMETER:  1.3107E-01  1.6847E-01  9.8510E+00  3.8739E-01  1.0425E+00  5.7010E-01  1.9292E+00 -1.1144E+00 -1.0230E-01 -7.0289E+00
             2.7726E+00
 GRADIENT:   4.3293E+00  1.0345E+00 -7.0477E-05  7.0664E+00  1.9814E+00  3.7793E+00  6.8140E+01  2.7245E-06 -1.2811E-01  0.0000E+00
             4.8168E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -466.960542692331        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     3086             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0310E+00  1.0703E+00  3.3079E+04  1.3348E+00  2.5661E+00  1.6013E+00  6.2358E+00  2.9505E-01  8.1969E-01  1.0000E-02
             1.4483E+01
 PARAMETER:  1.3049E-01  1.6791E-01  1.0507E+01  3.8879E-01  1.0424E+00  5.7079E-01  1.9303E+00 -1.1206E+00 -9.8831E-02 -7.0289E+00
             2.7730E+00
 GRADIENT:   3.7329E+00  1.0930E+00 -3.4708E-05  7.2063E+00  1.9000E+00  4.0242E+00  6.8286E+01  1.7184E-06 -1.2151E-01  0.0000E+00
             4.8855E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -466.961880020756        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     3277             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0310E+00  1.0679E+00  3.7648E+04  1.3361E+00  2.5668E+00  1.6014E+00  6.2483E+00  2.9199E-01  8.2258E-01  1.0000E-02
             1.4486E+01
 PARAMETER:  1.3053E-01  1.6570E-01  1.0636E+01  3.8975E-01  1.0426E+00  5.7086E-01  1.9323E+00 -1.1310E+00 -9.5308E-02 -7.0289E+00
             2.7732E+00
 GRADIENT:   3.2844E+00  1.0744E+00 -3.8514E-05  6.9640E+00  2.0341E+00  4.3263E+00  6.8805E+01  4.3419E-06 -8.2145E-02  0.0000E+00
             4.9475E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -466.963192502890        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     3445
 NPARAMETR:  1.0308E+00  1.0653E+00  3.7830E+04  1.3383E+00  2.5667E+00  1.6013E+00  6.2544E+00  2.9157E-01  8.3112E-01  1.0000E-02
             1.4484E+01
 PARAMETER:  1.3038E-01  1.6323E-01  1.0641E+01  3.9137E-01  1.0426E+00  5.7079E-01  1.9333E+00 -1.1325E+00 -8.4975E-02 -7.0289E+00
             2.7731E+00
 GRADIENT:   3.6106E+00  1.0824E+00 -4.2635E-05  6.4024E+00  2.1217E+00  3.9801E+00  6.9064E+01  1.3445E-06  1.7327E-01  0.0000E+00
             4.9591E+01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -466.963976864355        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     3617             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0309E+00  1.0586E+00  3.8556E+04  1.3397E+00  2.5670E+00  1.6015E+00  6.2568E+00  2.9344E-01  8.3043E-01  1.0000E-02
             1.4486E+01
 PARAMETER:  1.3045E-01  1.5695E-01  1.0660E+01  3.9248E-01  1.0428E+00  5.7092E-01  1.9337E+00 -1.1261E+00 -8.5806E-02 -7.0289E+00
             2.7732E+00
 GRADIENT:   3.1780E+00  8.4180E-01 -4.3609E-05  6.6515E+00  2.1650E+00  4.1814E+00  6.8988E+01  2.1118E-06  4.8165E-02  0.0000E+00
             4.9911E+01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -466.964516345209        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3797            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0310E+00  1.0603E+00  3.9541E+04  1.3410E+00  2.5669E+00  1.6014E+00  6.2637E+00  2.9016E-01  8.3465E-01  1.0000E-02
             1.4486E+01
 PARAMETER:  1.3057E-01  1.5857E-01  1.0685E+01  3.9344E-01  1.0427E+00  5.7088E-01  1.9348E+00 -1.1373E+00 -8.0745E-02 -7.0289E+00
             2.7732E+00
 GRADIENT:   3.0699E+00  1.0238E+00 -4.4330E-05  6.6977E+00  2.1374E+00  4.2633E+00  6.9221E+01  2.2934E-06  9.5283E-02  0.0000E+00
             4.9924E+01

0ITERATION NO.:  123    OBJECTIVE VALUE:  -466.964543323866        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:     3875
 NPARAMETR:  1.0307E+00  1.0586E+00  4.0192E+04  1.3446E+00  2.5659E+00  1.6017E+00  6.2621E+00  2.8867E-01  8.3480E-01  1.0000E-02
             1.4483E+01
 PARAMETER:  1.3057E-01  1.5853E-01  1.0685E+01  3.9345E-01  1.0427E+00  5.7088E-01  1.9348E+00 -1.1373E+00 -7.9745E-02 -7.0289E+00
             2.7732E+00
 GRADIENT:   9.6058E-02  8.5329E-02 -1.1696E-04 -6.7821E-01  5.2721E-02 -1.8894E-02  2.6779E-02  2.6025E-03  5.4369E-02  0.0000E+00
             1.7809E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3875
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3679E-02  3.0764E-02 -1.1256E-08 -6.3660E-02  1.5784E-05
 SE:             2.6784E-02  2.4279E-02  8.1209E-08  1.1556E-02  8.1602E-05
 N:                     100         100         100         100         100

 P VAL.:         3.7666E-01  2.0513E-01  8.8976E-01  3.6250E-08  8.4663E-01

 ETASHRINKSD(%)  1.0270E+01  1.8661E+01  1.0000E+02  6.1285E+01  9.9727E+01
 ETASHRINKVR(%)  1.9485E+01  3.3839E+01  1.0000E+02  8.5011E+01  9.9999E+01
 EBVSHRINKSD(%)  1.3022E+01  1.2968E+01  1.0000E+02  6.6044E+01  9.9637E+01
 EBVSHRINKVR(%)  2.4349E+01  2.4255E+01  1.0000E+02  8.8470E+01  9.9999E+01
 RELATIVEINF(%)  7.5226E+01  4.0604E+01  3.1502E-10  6.0770E+00  2.7551E-04
 EPSSHRINKSD(%)  2.5712E+00
 EPSSHRINKVR(%)  5.0764E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -466.96454332386628     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       1187.1248164445444     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   138.06
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    17.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -466.965       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.06E+00  3.95E+04  1.34E+00  2.57E+00  1.60E+00  6.26E+00  2.90E-01  8.35E-01  1.00E-02  1.45E+01
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.21E+02
 
 TH 2
+       -1.30E+01  5.72E+00
 
 TH 3
+        4.42E-07 -2.71E-09  1.94E-15
 
 TH 4
+       -1.20E+02  3.10E+01 -3.85E-08  1.79E+02
 
 TH 5
+        1.47E+01 -3.72E+00  1.48E-08 -2.12E+01  2.58E+00
 
 TH 6
+       -4.26E+01  4.15E+00 -2.53E-07  2.07E+01 -4.06E+00  3.87E+01
 
 TH 7
+        7.97E+00 -1.80E+00 -2.11E-10 -1.07E+01  1.24E+00 -7.12E-01  6.59E-01
 
 TH 8
+       -1.35E+00  2.74E-01 -1.20E-09  1.63E+00 -1.98E-01  3.03E-01 -9.79E-02  1.62E-02
 
 TH 9
+        4.05E+01 -1.06E+01  1.12E-08 -6.12E+01  7.24E+00 -6.87E+00  3.67E+00 -5.58E-01  2.09E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.49E+00 -1.18E+00 -2.30E-08 -5.97E+00  6.24E-01  1.69E+00  3.68E-01 -5.29E-02  2.07E+00  0.00E+00  7.04E-01
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.79E+02
 
 TH 2
+        5.36E+00  2.14E+01
 
 TH 3
+        4.34E-07 -5.07E-07  6.75E-13
 
 TH 4
+       -1.27E+01  2.60E+01  1.87E-07  1.44E+02
 
 TH 5
+        1.11E+00 -3.28E+00 -4.89E-08 -1.39E+01  2.10E+01
 
 TH 6
+       -1.65E+01 -9.90E-02 -2.15E-07  1.54E+00  6.43E-01  3.99E+01
 
 TH 7
+        1.51E+00  2.41E+00 -1.42E-08 -9.10E+00  4.80E-01  1.51E-01  2.84E+00
 
 TH 8
+        2.22E-01 -2.80E+00 -1.31E-06  3.93E-01 -4.80E-02  3.67E-01 -1.64E-02  1.77E+00
 
 TH 9
+        2.23E+00 -5.49E+00 -1.72E-07 -4.41E+01  3.92E+00 -5.36E-01  2.08E+00 -5.52E+00  4.01E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.98E+00 -2.20E+00 -4.18E-09 -8.52E+00  6.04E-01  1.65E+00  2.49E-01  5.62E-03  2.63E+00  0.00E+00  3.96E+00
 
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
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.13E+02
 
 TH 2
+        3.96E+01  1.96E+01
 
 TH 3
+        1.36E-08  7.75E-09  8.66E-17
 
 TH 4
+        8.62E+01  2.62E+01  4.71E-08  1.54E+02
 
 TH 5
+       -1.11E+01 -2.24E+00 -2.95E-08 -2.15E+01  1.21E+01
 
 TH 6
+        2.26E+00  7.91E-01  9.36E-09 -1.06E+01 -3.35E+00  4.39E+01
 
 TH 7
+       -4.47E+00  2.69E+00 -3.48E-09 -1.03E+01  2.41E+00  9.05E-01  3.09E+00
 
 TH 8
+        1.62E-03 -7.81E-04 -1.39E-12  6.61E-03 -1.77E-04 -1.67E-03 -3.66E-04  8.49E-06
 
 TH 9
+       -8.63E+00 -2.00E+00 -9.73E-09 -3.29E+01  4.31E+00  2.21E+00  1.74E+00 -4.84E-03  1.64E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.92E+01 -1.90E+00 -8.60E-09 -2.24E+01  4.99E+00 -6.93E+00  3.17E+00 -4.02E-03  4.28E+00  0.00E+00  6.45E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,      155.600
Stop Time:
Wed Sep 29 10:26:24 CDT 2021
