Sat Sep 18 09:36:22 CDT 2021
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
$DATA ../../../../data/spa/A2/dat4.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      500
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

 TOT. NO. OF OBS RECS:      400
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
 RAW OUTPUT FILE (FILE): m4.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -872.142225072536        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1364E+02  2.0095E+01  8.0244E+01 -5.3473E+01  1.4036E+02 -2.4928E+01 -2.5130E+01 -6.0009E+01 -7.5115E+01 -9.2422E+01
            -1.2508E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1322.60968266619        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0207E+00  7.8821E-01  6.8362E-01  1.1718E+00  7.0275E-01  9.7541E-01  8.8849E-01  1.2180E+00  1.1162E+00  7.9594E-01
             3.2740E+00
 PARAMETER:  1.2045E-01 -1.3799E-01 -2.8036E-01  2.5853E-01 -2.5276E-01  7.5105E-02 -1.8230E-02  2.9721E-01  2.0996E-01 -1.2824E-01
             1.2860E+00
 GRADIENT:   6.0644E+01 -2.2308E+01 -3.2836E+01  1.0187E+01  3.4149E+01 -2.0983E+01  4.5304E+00  1.4785E+01  2.0188E+01  1.7862E+01
             1.1117E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1333.22850650011        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0287E+00  6.5904E-01  2.7768E-01  1.1501E+00  3.4378E-01  1.1335E+00  8.5502E-01  1.0345E+00  1.0090E+00  4.2490E-01
             2.8524E+00
 PARAMETER:  1.2826E-01 -3.1697E-01 -1.1813E+00  2.3981E-01 -9.6775E-01  2.2532E-01 -5.6629E-02  1.3390E-01  1.0894E-01 -7.5591E-01
             1.1482E+00
 GRADIENT:   5.5841E+01  1.1712E+02  6.9398E+01  1.1693E+02 -1.8482E+02  2.0676E+01 -8.1045E+00  2.1940E+00 -8.0937E+00  1.2074E+00
             7.5433E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1364.90651370364        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.5675E-01  6.6986E-01  1.9472E-01  9.7668E-01  3.2635E-01  1.1021E+00  9.7881E-01  1.0806E+00  1.0260E+00  3.0182E-01
             2.0464E+00
 PARAMETER:  5.5783E-02 -3.0069E-01 -1.5362E+00  7.6400E-02 -1.0198E+00  1.9717E-01  7.8585E-02  1.7747E-01  1.2564E-01 -1.0979E+00
             8.1607E-01
 GRADIENT:  -1.7583E+01  5.9327E+00  1.5310E+01  4.9913E+00 -2.3613E+00  1.4147E+01 -2.3603E+00 -2.2955E+01 -1.7456E+01  3.0067E+00
            -1.5012E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1365.67534051839        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:      331
 NPARAMETR:  9.5573E-01  6.6852E-01  1.9167E-01  9.7281E-01  3.2399E-01  1.1009E+00  9.7552E-01  1.1103E+00  1.0321E+00  2.9500E-01
             2.0300E+00
 PARAMETER:  5.4720E-02 -3.0268E-01 -1.5520E+00  7.2435E-02 -1.0270E+00  1.9614E-01  7.5219E-02  2.0462E-01  1.3157E-01 -1.1208E+00
             8.0803E-01
 GRADIENT:  -1.8055E+01  6.5471E+00  1.4556E+01  3.6015E+00 -3.8379E+00  1.3947E+01 -1.9991E+00 -2.3394E+01 -1.7731E+01  2.9227E+00
            -1.4117E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1366.48035910773        NO. OF FUNC. EVALS.: 213
 CUMULATIVE NO. OF FUNC. EVALS.:      544             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8365E-01  6.6330E-01  1.9179E-01  9.7273E-01  3.2713E-01  1.1009E+00  1.0170E+00  1.1105E+00  1.0322E+00  1.7483E-01
             2.0306E+00
 PARAMETER:  8.3516E-02 -3.1053E-01 -1.5513E+00  7.2349E-02 -1.0174E+00  1.9616E-01  1.1690E-01  2.0481E-01  1.3171E-01 -1.6439E+00
             8.0834E-01
 GRADIENT:   3.3874E+01 -1.6921E+01  1.6818E+00  2.4325E+00  3.0076E+01  1.3866E+01 -2.9980E+00 -2.5198E+01 -1.5188E+01  7.0163E-01
            -1.6134E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1366.89905951366        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      615
 NPARAMETR:  9.6421E-01  6.6108E-01  1.9180E-01  9.7273E-01  3.2374E-01  1.1009E+00  1.0572E+00  1.1105E+00  1.0322E+00  1.0120E-01
             2.0306E+00
 PARAMETER:  6.3550E-02 -3.1388E-01 -1.5513E+00  7.2348E-02 -1.0278E+00  1.9616E-01  1.5560E-01  2.0481E-01  1.3171E-01 -2.1906E+00
             8.0835E-01
 GRADIENT:  -2.6241E+00  4.7991E-01  1.1846E+01  1.3902E+00 -6.7108E-01  1.3897E+01  1.4636E-02 -2.6086E+01 -1.5024E+01  1.8340E-01
            -1.4510E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1367.10408698393        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      691
 NPARAMETR:  9.6380E-01  6.7432E-01  1.9179E-01  9.7273E-01  3.2782E-01  1.1009E+00  1.0500E+00  1.1105E+00  1.0322E+00  1.0000E-02
             2.0307E+00
 PARAMETER:  6.3127E-02 -2.9404E-01 -1.5513E+00  7.2347E-02 -1.0153E+00  1.9617E-01  1.4874E-01  2.0481E-01  1.3171E-01 -7.7201E+00
             8.0839E-01
 GRADIENT:  -2.7101E+00  2.9375E+00  1.0090E+01  9.7258E+00 -2.7475E-02  1.4005E+01 -6.4283E-01 -2.6867E+01 -1.5075E+01  0.0000E+00
            -1.5926E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1367.21487191434        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      879            RESET HESSIAN, TYPE II
 NPARAMETR:  9.6897E-01  6.8867E-01  1.9179E-01  9.7272E-01  3.3299E-01  1.1009E+00  1.0428E+00  1.1112E+00  1.0323E+00  1.0000E-02
             2.0308E+00
 PARAMETER:  6.8483E-02 -2.7299E-01 -1.5514E+00  7.2337E-02 -9.9965E-01  1.9614E-01  1.4195E-01  2.0541E-01  1.3176E-01 -1.0989E+01
             8.0841E-01
 GRADIENT:   7.7533E+00  1.4029E+00  4.9126E+00  1.8406E+01  1.0358E+01  1.4293E+01  1.0646E-01 -2.6611E+01 -1.4817E+01  0.0000E+00
            -1.7529E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1367.22551198081        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  9.6907E-01  6.8848E-01  1.9209E-01  9.7261E-01  3.3309E-01  1.1011E+00  1.0442E+00  1.1112E+00  1.0322E+00  1.0000E-02
             2.0324E+00
 PARAMETER:  6.8581E-02 -2.7326E-01 -1.5498E+00  7.2233E-02 -9.9934E-01  1.9632E-01  1.4327E-01  2.0548E-01  1.3165E-01 -1.0989E+01
             8.0922E-01
 GRADIENT:   6.7897E-02 -1.7591E-01  2.4093E+00  1.6082E+01  7.2475E-02  1.2596E+01  1.8956E-01 -2.6547E+01 -1.4985E+01  0.0000E+00
            -1.7908E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1367.30096724338        NO. OF FUNC. EVALS.: 221
 CUMULATIVE NO. OF FUNC. EVALS.:     1282
 NPARAMETR:  9.6907E-01  6.8848E-01  1.9208E-01  9.7258E-01  3.3310E-01  1.1009E+00  1.0428E+00  1.1139E+00  1.0325E+00  1.0000E-02
             2.0327E+00
 PARAMETER:  6.8581E-02 -2.7326E-01 -1.5498E+00  7.2195E-02 -9.9932E-01  1.9611E-01  1.4196E-01  2.0789E-01  1.3201E-01 -1.0989E+01
             8.0938E-01
 GRADIENT:   8.2721E-02 -3.5879E-01  2.4013E+00  1.6049E+01  2.1306E-01  1.2546E+01  8.2645E-02 -2.6489E+01 -1.4959E+01  0.0000E+00
            -1.7588E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1367.37454241276        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1476            RESET HESSIAN, TYPE II
 NPARAMETR:  9.6917E-01  6.8830E-01  1.9238E-01  9.7245E-01  3.3320E-01  1.1010E+00  1.0427E+00  1.1164E+00  1.0326E+00  1.0000E-02
             2.0347E+00
 PARAMETER:  6.8681E-02 -2.7354E-01 -1.5483E+00  7.2061E-02 -9.9901E-01  1.9618E-01  1.4177E-01  2.1012E-01  1.3209E-01 -1.0989E+01
             8.1033E-01
 GRADIENT:   7.8807E+00  8.7305E-01  5.6499E+00  1.6686E+01  1.0280E+01  1.4371E+01  1.6182E-01 -2.6263E+01 -1.4498E+01  0.0000E+00
            -1.6368E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1367.37670883291        NO. OF FUNC. EVALS.: 224
 CUMULATIVE NO. OF FUNC. EVALS.:     1700
 NPARAMETR:  9.6921E-01  6.8820E-01  1.9253E-01  9.7240E-01  3.3318E-01  1.1011E+00  1.0415E+00  1.1163E+00  1.0327E+00  1.0000E-02
             2.0355E+00
 PARAMETER:  6.8731E-02 -2.7367E-01 -1.5475E+00  7.2011E-02 -9.9907E-01  1.9627E-01  1.4068E-01  2.1001E-01  1.3215E-01 -1.0989E+01
             8.1074E-01
 GRADIENT:   1.7852E-01 -3.8415E-01  3.2506E+00  1.4773E+01 -5.9050E-01  1.2636E+01 -1.1189E-01 -2.6337E+01 -1.4795E+01  0.0000E+00
            -1.6913E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1367.41348458587        NO. OF FUNC. EVALS.: 211
 CUMULATIVE NO. OF FUNC. EVALS.:     1911             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6913E-01  6.8855E-01  1.9260E-01  9.7236E-01  3.3324E-01  1.1010E+00  1.0419E+00  1.1176E+00  1.0328E+00  1.0000E-02
             2.0361E+00
 PARAMETER:  6.8646E-02 -2.7317E-01 -1.5471E+00  7.1968E-02 -9.9889E-01  1.9625E-01  1.4108E-01  2.1119E-01  1.3224E-01 -1.0989E+01
             8.1102E-01
 GRADIENT:   7.7173E+00  1.5873E+00  6.3652E+00  1.6213E+01  8.8388E+00  1.4410E+01  9.4366E-02 -2.6198E+01 -1.4424E+01  0.0000E+00
            -1.6031E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1367.41477286940        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2088
 NPARAMETR:  9.6904E-01  6.9084E-01  1.9260E-01  9.7236E-01  3.3404E-01  1.1010E+00  1.0404E+00  1.1176E+00  1.0328E+00  1.0000E-02
             2.0361E+00
 PARAMETER:  6.8550E-02 -2.6985E-01 -1.5471E+00  7.1968E-02 -9.9648E-01  1.9625E-01  1.3961E-01  2.1119E-01  1.3224E-01 -1.0989E+01
             8.1102E-01
 GRADIENT:  -4.1862E-02  1.0563E-02  2.8561E+00  1.6173E+01  1.0597E-01  1.2676E+01  6.9281E-02 -2.6252E+01 -1.4756E+01  0.0000E+00
            -1.6922E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1367.42377175615        NO. OF FUNC. EVALS.: 207
 CUMULATIVE NO. OF FUNC. EVALS.:     2295
 NPARAMETR:  9.6915E-01  6.9053E-01  1.9287E-01  9.7226E-01  3.3399E-01  1.1012E+00  1.0399E+00  1.1177E+00  1.0327E+00  1.0000E-02
             2.0376E+00
 PARAMETER:  6.8668E-02 -2.7029E-01 -1.5457E+00  7.1871E-02 -9.9664E-01  1.9642E-01  1.3912E-01  2.1130E-01  1.3214E-01 -1.0989E+01
             8.1177E-01
 GRADIENT:   3.7970E-02  3.6113E-01  3.5972E+00  1.5293E+01 -1.0568E+00  1.2736E+01 -7.3699E-02 -2.6212E+01 -1.4698E+01  0.0000E+00
            -1.6655E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1367.47890600502        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     2467
 NPARAMETR:  9.6915E-01  6.9053E-01  1.9287E-01  9.7224E-01  3.3410E-01  1.1011E+00  1.0405E+00  1.1198E+00  1.0328E+00  1.0000E-02
             2.0378E+00
 PARAMETER:  6.8666E-02 -2.7029E-01 -1.5457E+00  7.1842E-02 -9.9631E-01  1.9631E-01  1.3968E-01  2.1316E-01  1.3231E-01 -1.0989E+01
             8.1189E-01
 GRADIENT:   7.7746E+00  1.2276E+00  6.0653E+00  1.6790E+01  1.0057E+01  1.4493E+01  1.6714E-01 -2.6037E+01 -1.4298E+01  0.0000E+00
            -1.5747E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1367.48774643716        NO. OF FUNC. EVALS.: 212
 CUMULATIVE NO. OF FUNC. EVALS.:     2679             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6915E-01  6.9056E-01  1.9289E-01  9.7223E-01  3.3410E-01  1.1010E+00  1.0403E+00  1.1200E+00  1.0330E+00  1.0000E-02
             2.0379E+00
 PARAMETER:  6.8663E-02 -2.7026E-01 -1.5456E+00  7.1833E-02 -9.9632E-01  1.9623E-01  1.3950E-01  2.1336E-01  1.3246E-01 -1.0989E+01
             8.1192E-01
 GRADIENT:   7.7627E+00  1.3026E+00  6.1495E+00  1.6747E+01  9.8986E+00  1.4465E+01  1.4980E-01 -2.6037E+01 -1.4267E+01  0.0000E+00
            -1.5711E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1367.48797654563        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2856
 NPARAMETR:  9.6906E-01  6.9153E-01  1.9289E-01  9.7223E-01  3.3438E-01  1.1010E+00  1.0387E+00  1.1200E+00  1.0330E+00  1.0000E-02
             2.0379E+00
 PARAMETER:  6.8576E-02 -2.6885E-01 -1.5456E+00  7.1833E-02 -9.9547E-01  1.9623E-01  1.3801E-01  2.1336E-01  1.3246E-01 -1.0989E+01
             8.1192E-01
 GRADIENT:  -5.8361E-02  5.4939E-02  3.2362E+00  1.5861E+01 -9.7750E-02  1.2706E+01 -4.2223E-02 -2.6133E+01 -1.4641E+01  0.0000E+00
            -1.6472E+01

0ITERATION NO.:   92    OBJECTIVE VALUE:  -1367.48799041010        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:     2958
 NPARAMETR:  9.6909E-01  6.9159E-01  1.9289E-01  9.7223E-01  3.3441E-01  1.1009E+00  1.0389E+00  1.1200E+00  1.0331E+00  1.0000E-02
             2.0379E+00
 PARAMETER:  6.8601E-02 -2.6876E-01 -1.5456E+00  7.1833E-02 -9.9538E-01  1.9623E-01  1.3819E-01  2.1336E-01  1.3246E-01 -1.0989E+01
             8.1192E-01
 GRADIENT:  -1.2040E-02  1.7955E-02 -3.5836E+04  2.7711E+05  3.0102E-02  1.2690E+01 -6.5978E-03  2.5970E+05 -1.4623E+01  0.0000E+00
            -3.4189E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2958
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3626E-04  4.0327E-03 -1.3920E-02 -1.6496E-02  1.3895E-04
 SE:             2.8509E-02  2.3379E-02  2.0063E-02  2.7180E-02  3.4640E-04
 N:                     100         100         100         100         100

 P VAL.:         9.9339E-01  8.6305E-01  4.8779E-01  5.4392E-01  6.8832E-01

 ETASHRINKSD(%)  4.4928E+00  2.1676E+01  3.2787E+01  8.9431E+00  9.8840E+01
 ETASHRINKVR(%)  8.7838E+00  3.8654E+01  5.4824E+01  1.7086E+01  9.9987E+01
 EBVSHRINKSD(%)  1.3654E+00  2.1355E+01  4.4092E+01  1.2669E+01  9.8859E+01
 EBVSHRINKVR(%)  2.7122E+00  3.8150E+01  6.8743E+01  2.3733E+01  9.9987E+01
 RELATIVEINF(%)  9.4923E+01  3.6474E+00  9.1222E+00  2.8961E+01  6.5773E-04
 EPSSHRINKSD(%)  3.8848E+01
 EPSSHRINKVR(%)  6.2605E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1367.4879904100976     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -632.33716384635943     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    49.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1367.488       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.69E-01  6.92E-01  1.93E-01  9.72E-01  3.34E-01  1.10E+00  1.04E+00  1.12E+00  1.03E+00  1.00E-02  2.04E+00
 


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
+        9.52E+02
 
 TH 2
+       -9.24E+00  1.35E+03
 
 TH 3
+       -5.83E+03  1.61E+04  1.56E+08
 
 TH 4
+        3.39E+04 -4.86E+04  4.98E+04  1.47E+09
 
 TH 5
+        1.04E+02 -3.84E+03 -2.28E+04  4.65E+04  1.35E+04
 
 TH 6
+        6.62E+08 -3.45E+08 -7.27E+03  5.26E+04 -1.56E+08  2.97E+08
 
 TH 7
+       -9.26E-01  4.17E+01 -1.20E+03  1.99E+03  4.46E+01  3.05E-01  7.18E+01
 
 TH 8
+        2.13E+04 -2.21E+04  3.39E+04  5.96E+08  1.87E+04 -2.68E+08  2.79E+02  2.43E+08
 
 TH 9
+       -1.05E+09 -2.83E+01 -1.91E+03 -2.42E+04  1.07E+02 -4.69E+08  1.68E+01 -3.65E+04  7.40E+08
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.37E+02  2.45E+03 -2.57E+04  6.55E+03 -2.79E+03 -9.29E+02 -2.25E+02  5.37E+03 -8.00E+02  0.00E+00  5.07E+06
 
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
 #CPUT: Total CPU Time in Seconds,       57.641
Stop Time:
Sat Sep 18 09:37:21 CDT 2021
