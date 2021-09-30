Wed Sep 29 18:20:57 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat67.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1652.56499016208        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1549E+02 -7.1303E+01  4.8142E+00 -1.0572E+02 -6.3600E+00  4.1812E+01 -1.4667E+01  4.9868E+00 -1.2352E+01  5.9135E+00
             6.4671E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1663.82701778513        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9460E-01  1.0757E+00  1.0435E+00  1.0703E+00  1.0534E+00  1.0490E+00  1.0941E+00  9.7327E-01  1.0203E+00  9.8640E-01
             9.7825E-01
 PARAMETER:  9.4582E-02  1.7302E-01  1.4258E-01  1.6791E-01  1.5206E-01  1.4788E-01  1.8994E-01  7.2906E-02  1.2006E-01  8.6304E-02
             7.8012E-02
 GRADIENT:   2.2789E-01 -3.0164E+00  2.3210E+00 -2.0802E+00  1.0594E+00  9.6692E+00 -7.3940E+00 -1.7082E-01  3.3788E+00 -5.8817E+00
            -6.0963E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1664.66851390447        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  9.9850E-01  1.0207E+00  1.0973E+00  1.1093E+00  1.0727E+00  1.0233E+00  1.2802E+00  8.9156E-01  9.3654E-01  1.0608E+00
             9.8944E-01
 PARAMETER:  9.8503E-02  1.2044E-01  1.9288E-01  2.0376E-01  1.7017E-01  1.2307E-01  3.4705E-01 -1.4786E-02  3.4437E-02  1.5902E-01
             8.9382E-02
 GRADIENT:   1.0071E+01  1.3444E+00  1.6145E+00  8.0633E-02  7.0547E+00  3.1252E-01 -3.3968E-01 -1.6295E+00 -1.6862E+00 -1.2721E+00
            -2.6796E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.94204157143        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.9251E-01  8.9033E-01  1.1159E+00  1.1949E+00  1.0080E+00  1.0212E+00  1.3792E+00  9.1270E-01  9.0832E-01  1.0229E+00
             9.9109E-01
 PARAMETER:  9.2486E-02 -1.6167E-02  2.0968E-01  2.7807E-01  1.0799E-01  1.2096E-01  4.2152E-01  8.6565E-03  3.8370E-03  1.2262E-01
             9.1048E-02
 GRADIENT:  -8.2975E-01  5.0976E+00  1.6918E+00  5.2470E+00 -4.4313E+00 -4.6890E-02 -4.1072E-01  7.4113E-02 -2.7853E-01  4.7077E-01
            -7.2507E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1665.02474394381        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.8935E-01  7.0882E-01  1.2427E+00  1.3127E+00  9.9621E-01  1.0169E+00  1.5661E+00  1.0064E+00  8.6914E-01  1.0345E+00
             9.9116E-01
 PARAMETER:  8.9296E-02 -2.4415E-01  3.1731E-01  3.7209E-01  9.6206E-02  1.1678E-01  5.4861E-01  1.0636E-01 -4.0246E-02  1.3395E-01
             9.1120E-02
 GRADIENT:  -2.4144E+00  4.1920E+00  2.2162E+00  4.3486E+00 -4.5533E+00 -7.2523E-01  1.8898E-03  1.3018E-01 -4.9837E-01 -7.4233E-02
            -1.6105E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1665.04102082167        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.8714E-01  5.5501E-01  1.3760E+00  1.4190E+00  9.9920E-01  1.0145E+00  1.7592E+00  1.1135E+00  8.4480E-01  1.0587E+00
             9.9188E-01
 PARAMETER:  8.7060E-02 -4.8876E-01  4.1916E-01  4.4998E-01  9.9198E-02  1.1444E-01  6.6486E-01  2.0749E-01 -6.8652E-02  1.5705E-01
             9.1843E-02
 GRADIENT:  -1.9240E+00  6.6364E+00  3.7836E+00  1.1048E+01 -6.2596E+00 -7.3557E-01  8.1094E-01 -2.2620E-01 -6.6575E-01 -2.1210E-01
            -3.9365E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1665.05307273876        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.8578E-01  4.5645E-01  1.4751E+00  1.4894E+00  1.0070E+00  1.0130E+00  1.9341E+00  1.1968E+00  8.2891E-01  1.0767E+00
             9.9251E-01
 PARAMETER:  8.5681E-02 -6.8428E-01  4.8875E-01  4.9836E-01  1.0694E-01  1.1296E-01  7.5963E-01  2.7962E-01 -8.7638E-02  1.7390E-01
             9.2479E-02
 GRADIENT:  -1.3135E+00  8.1409E+00  5.0664E+00  1.8839E+01 -6.9829E+00 -7.3058E-01  1.3433E+00 -7.0275E-01 -6.8818E-01 -3.2837E-01
            -6.3277E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1665.06619894157        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  9.8507E-01  3.9281E-01  1.5396E+00  1.5341E+00  1.0126E+00  1.0123E+00  2.0950E+00  1.2474E+00  8.1769E-01  1.0872E+00
             9.9302E-01
 PARAMETER:  8.4954E-02 -8.3444E-01  5.3151E-01  5.2794E-01  1.1254E-01  1.1221E-01  8.3956E-01  3.2102E-01 -1.0127E-01  1.8361E-01
             9.2993E-02
 GRADIENT:  -5.7940E-01  8.2592E+00  5.7206E+00  2.3246E+01 -6.3825E+00 -6.4884E-01  1.6355E+00 -1.3193E+00 -6.4898E-01 -5.6513E-01
            -8.6076E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1665.08446609337        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1409
 NPARAMETR:  9.8456E-01  3.3205E-01  1.6073E+00  1.5764E+00  1.0203E+00  1.0120E+00  2.2955E+00  1.3089E+00  8.0705E-01  1.0979E+00
             9.9339E-01
 PARAMETER:  8.4443E-02 -1.0025E+00  5.7453E-01  5.5512E-01  1.2011E-01  1.1190E-01  9.3094E-01  3.6922E-01 -1.1437E-01  1.9338E-01
             9.3367E-02
 GRADIENT:   4.9920E-01  7.8246E+00  5.6033E+00  2.6252E+01 -5.1894E+00 -4.1842E-01  1.9094E+00 -1.5748E+00 -5.1730E-01 -6.0526E-01
            -9.1935E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1665.15698572857        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1584
 NPARAMETR:  9.8384E-01  2.4620E-01  1.7010E+00  1.6343E+00  1.0298E+00  1.0116E+00  2.6897E+00  1.3991E+00  7.9245E-01  1.1098E+00
             9.9366E-01
 PARAMETER:  8.3708E-02 -1.3016E+00  6.3121E-01  5.9120E-01  1.2933E-01  1.1156E-01  1.0894E+00  4.3581E-01 -1.3262E-01  2.0421E-01
             9.3643E-02
 GRADIENT:   1.9378E+00  6.5080E+00  4.8653E+00  2.7323E+01 -3.0850E+00 -5.3284E-02  2.1997E+00 -1.7519E+00 -4.5267E-01 -6.6876E-01
            -9.3266E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1665.32742065164        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1759
 NPARAMETR:  9.8279E-01  1.4455E-01  1.8088E+00  1.7009E+00  1.0386E+00  1.0111E+00  3.5028E+00  1.5057E+00  7.7663E-01  1.1213E+00
             9.9384E-01
 PARAMETER:  8.2637E-02 -1.8341E+00  6.9268E-01  6.3117E-01  1.3783E-01  1.1106E-01  1.3536E+00  5.0926E-01 -1.5279E-01  2.1451E-01
             9.3819E-02
 GRADIENT:   3.1217E+00  4.4258E+00  3.5752E+00  2.4664E+01 -6.8967E-01  3.1390E-01  2.4082E+00 -1.6815E+00 -5.3197E-01 -7.1627E-01
            -8.4449E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1665.71890296819        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1937
 NPARAMETR:  9.8089E-01  5.3215E-02  1.8455E+00  1.7538E+00  1.0260E+00  1.0099E+00  5.5500E+00  1.5589E+00  7.6281E-01  1.1164E+00
             9.9327E-01
 PARAMETER:  8.0701E-02 -2.8334E+00  7.1273E-01  6.6180E-01  1.2564E-01  1.0983E-01  1.8138E+00  5.4396E-01 -1.7075E-01  2.1013E-01
             9.3252E-02
 GRADIENT:   2.1911E+00  1.2872E+00  1.6923E+00  1.5475E+01 -3.4854E-01  3.5692E-01  3.3571E-01 -6.7570E-01  2.1239E-01 -1.2254E-01
            -4.5250E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1665.96361839697        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2113
 NPARAMETR:  9.7912E-01  1.3939E-02  1.8085E+00  1.7713E+00  1.0042E+00  1.0082E+00  1.0060E+01  1.5378E+00  7.5729E-01  1.1053E+00
             9.9278E-01
 PARAMETER:  7.8904E-02 -4.1731E+00  6.9248E-01  6.7170E-01  1.0415E-01  1.0814E-01  2.4086E+00  5.3035E-01 -1.7801E-01  2.0015E-01
             9.2758E-02
 GRADIENT:  -5.3040E-02  2.8431E-01  2.4091E-01  5.2043E+00 -7.3022E-01 -4.4204E-02  1.4629E-02 -1.8373E-02  6.5920E-02  2.1869E-01
            -3.1887E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1666.14300875368        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2296             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7956E-01  1.0000E-02  1.7980E+00  1.7570E+00  9.9823E-01  1.0083E+00  1.1768E+01  1.5313E+00  7.5575E-01  1.1025E+00
             9.9061E-01
 PARAMETER:  7.9345E-02 -4.6675E+00  6.8669E-01  6.6358E-01  9.8227E-02  1.0829E-01  2.5654E+00  5.2613E-01 -1.8004E-01  1.9756E-01
             9.0569E-02
 GRADIENT:   4.0583E+02  0.0000E+00  1.1343E+01  1.2468E+03  4.4077E+00  5.4301E+01  7.4418E-01  2.1626E+00  2.8636E+01  2.1853E+00
             3.3591E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1666.14536967150        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2484             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7964E-01  1.0000E-02  1.7917E+00  1.7563E+00  9.9474E-01  1.0085E+00  1.2012E+01  1.5235E+00  7.5652E-01  1.1027E+00
             9.9187E-01
 PARAMETER:  7.9434E-02 -4.6675E+00  6.8316E-01  6.6322E-01  9.4725E-02  1.0846E-01  2.5859E+00  5.2103E-01 -1.7902E-01  1.9773E-01
             9.1835E-02
 GRADIENT:   4.0490E+02  0.0000E+00  1.2028E+01  1.2421E+03  2.8259E+00  5.4277E+01  7.8924E-01  2.1350E+00  2.8755E+01  2.5276E+00
             9.0439E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1666.14985972671        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2666
 NPARAMETR:  9.7954E-01  1.0000E-02  1.7841E+00  1.7580E+00  9.9661E-01  1.0084E+00  1.1917E+01  1.5198E+00  7.5654E-01  1.0967E+00
             9.9181E-01
 PARAMETER:  7.9331E-02 -4.6675E+00  6.7894E-01  6.6415E-01  9.6609E-02  1.0839E-01  2.5779E+00  5.1859E-01 -1.7900E-01  1.9231E-01
             9.1775E-02
 GRADIENT:   1.2836E+00  0.0000E+00 -1.5404E-01 -2.4771E+01  1.7934E+00  1.4651E-01  1.4251E-02  1.4743E-01  8.3308E-02 -3.0289E-01
            -2.3523E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1666.15600055302        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2858             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7962E-01  1.0000E-02  1.7795E+00  1.7557E+00  9.9186E-01  1.0085E+00  1.2010E+01  1.5131E+00  7.5672E-01  1.0991E+00
             9.9169E-01
 PARAMETER:  7.9405E-02 -4.6675E+00  6.7632E-01  6.6287E-01  9.1827E-02  1.0843E-01  2.5857E+00  5.1414E-01 -1.7876E-01  1.9445E-01
             9.1660E-02
 GRADIENT:   4.0484E+02  0.0000E+00  1.1422E+01  1.2409E+03  3.8330E+00  5.4242E+01  7.8969E-01  2.0488E+00  2.8737E+01  2.2141E+00
             8.5577E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1666.15972236066        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3046
 NPARAMETR:  9.7961E-01  1.0000E-02  1.7746E+00  1.7555E+00  9.9090E-01  1.0085E+00  1.2012E+01  1.5090E+00  7.5681E-01  1.0974E+00
             9.9164E-01
 PARAMETER:  7.9398E-02 -4.6675E+00  6.7359E-01  6.6273E-01  9.0860E-02  1.0843E-01  2.5859E+00  5.1145E-01 -1.7865E-01  1.9298E-01
             9.1609E-02
 GRADIENT:   1.4958E+00  0.0000E+00  1.1504E+00 -2.8478E+01 -9.2599E-01  1.8132E-01  1.7674E-02  1.7973E-01  1.4707E-01  3.3665E-01
             3.9321E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1666.16335505764        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3238             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7960E-01  1.0000E-02  1.7639E+00  1.7550E+00  9.9024E-01  1.0084E+00  1.2009E+01  1.5008E+00  7.5693E-01  1.0922E+00
             9.9149E-01
 PARAMETER:  7.9387E-02 -4.6675E+00  6.6755E-01  6.6245E-01  9.0195E-02  1.0841E-01  2.5857E+00  5.0603E-01 -1.7848E-01  1.8823E-01
             9.1456E-02
 GRADIENT:   4.0477E+02  0.0000E+00  9.5860E+00  1.2391E+03  7.3839E+00  5.4169E+01  7.9219E-01  1.8697E+00  2.8675E+01  1.2493E+00
             7.4388E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1666.16597907766        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3426
 NPARAMETR:  9.7959E-01  1.0000E-02  1.7612E+00  1.7548E+00  9.8871E-01  1.0084E+00  1.2014E+01  1.4979E+00  7.5700E-01  1.0928E+00
             9.9147E-01
 PARAMETER:  7.9382E-02 -4.6675E+00  6.6601E-01  6.6235E-01  8.8645E-02  1.0840E-01  2.5861E+00  5.0406E-01 -1.7839E-01  1.8874E-01
             9.1434E-02
 GRADIENT:   1.5186E+00  0.0000E+00 -4.3301E-02 -2.8706E+01  1.2313E+00  1.8140E-01  1.7233E-02  9.7312E-02  1.2328E-01 -1.9769E-01
            -3.2491E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1666.16815750192        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3618             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7959E-01  1.0000E-02  1.7579E+00  1.7546E+00  9.8644E-01  1.0084E+00  1.2021E+01  1.4940E+00  7.5711E-01  1.0935E+00
             9.9150E-01
 PARAMETER:  7.9375E-02 -4.6675E+00  6.6411E-01  6.6222E-01  8.6344E-02  1.0840E-01  2.5867E+00  5.0142E-01 -1.7825E-01  1.8939E-01
             9.1463E-02
 GRADIENT:   4.0467E+02  0.0000E+00  1.0463E+01  1.2383E+03  5.3044E+00  5.4144E+01  7.9346E-01  1.8975E+00  2.8693E+01  1.7880E+00
             8.2160E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1666.16969155996        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3806
 NPARAMETR:  9.7958E-01  1.0000E-02  1.7534E+00  1.7543E+00  9.8489E-01  1.0084E+00  1.2023E+01  1.4894E+00  7.5719E-01  1.0932E+00
             9.9143E-01
 PARAMETER:  7.9367E-02 -4.6675E+00  6.6158E-01  6.6208E-01  8.4775E-02  1.0839E-01  2.5868E+00  4.9841E-01 -1.7814E-01  1.8911E-01
             9.1389E-02
 GRADIENT:   1.5024E+00  0.0000E+00  5.3586E-01 -2.8594E+01 -2.9667E-01  1.8967E-01  1.7447E-02  1.0591E-01  1.3410E-01  1.8618E-01
             9.6734E-03

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1666.17088133438        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3998             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7957E-01  1.0000E-02  1.7486E+00  1.7541E+00  9.8436E-01  1.0084E+00  1.2024E+01  1.4860E+00  7.5726E-01  1.0907E+00
             9.9139E-01
 PARAMETER:  7.9363E-02 -4.6675E+00  6.5881E-01  6.6195E-01  8.4234E-02  1.0838E-01  2.5869E+00  4.9608E-01 -1.7805E-01  1.8678E-01
             9.1357E-02
 GRADIENT:   4.0467E+02  0.0000E+00  9.8845E+00  1.2372E+03  6.2787E+00  5.4128E+01  7.9477E-01  1.8259E+00  2.8669E+01  1.5104E+00
             7.9012E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1666.17127081339        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     4186
 NPARAMETR:  9.7957E-01  1.0000E-02  1.7457E+00  1.7539E+00  9.8376E-01  1.0084E+00  1.2025E+01  1.4837E+00  7.5730E-01  1.0899E+00
             9.9136E-01
 PARAMETER:  7.9360E-02 -4.6675E+00  6.5715E-01  6.6187E-01  8.3628E-02  1.0838E-01  2.5870E+00  4.9454E-01 -1.7800E-01  1.8608E-01
             9.1318E-02
 GRADIENT:   1.5897E+00  0.0000E+00 -2.9314E-01 -2.8795E+01  1.1962E+00  2.0371E-01  1.7275E-02  7.7889E-02  1.2861E-01 -2.0342E-01
            -1.6545E-02

0ITERATION NO.:  118    OBJECTIVE VALUE:  -1666.17193209698        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:     4290
 NPARAMETR:  9.7957E-01  1.0000E-02  1.7406E+00  1.7537E+00  9.8296E-01  1.0084E+00  1.2025E+01  1.4792E+00  7.5739E-01  1.0891E+00
             9.9130E-01
 PARAMETER:  7.9358E-02 -4.6675E+00  6.5747E-01  6.6186E-01  8.2617E-02  1.0838E-01  2.5872E+00  4.9418E-01 -1.7796E-01  1.8719E-01
             9.1337E-02
 GRADIENT:   1.8605E-03  0.0000E+00  3.7016E-01  1.6015E-01 -5.4614E-02  1.4999E-03  2.6829E-05  7.5073E-02 -8.2897E-03  7.1702E-02
             9.6353E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4290
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.1847E-04  3.8356E-04 -3.5817E-02 -7.2820E-03 -4.9090E-02
 SE:             2.9844E-02  1.8058E-03  1.8992E-02  2.9234E-02  1.9732E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9416E-01  8.3180E-01  5.9310E-02  8.0329E-01  1.2854E-02

 ETASHRINKSD(%)  2.0350E-02  9.3950E+01  3.6375E+01  2.0611E+00  3.3895E+01
 ETASHRINKVR(%)  4.0696E-02  9.9634E+01  5.9518E+01  4.0798E+00  5.6301E+01
 EBVSHRINKSD(%)  4.0253E-01  9.4147E+01  3.9686E+01  2.3645E+00  3.0455E+01
 EBVSHRINKVR(%)  8.0344E-01  9.9657E+01  6.3622E+01  4.6732E+00  5.1635E+01
 RELATIVEINF(%)  9.3509E+01  8.9684E-03  1.0043E+01  2.9349E+00  7.0854E+00
 EPSSHRINKSD(%)  4.5295E+01
 EPSSHRINKVR(%)  7.0074E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1666.1719320969758     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -931.02110553323757     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    65.25
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     7.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1666.172       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.00E-02  1.75E+00  1.75E+00  9.83E-01  1.01E+00  1.20E+01  1.48E+00  7.57E-01  1.09E+00  9.91E-01
 


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
+        1.15E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.60E+00  0.00E+00  4.06E+01
 
 TH 4
+        4.69E+00  0.00E+00  1.39E+00  6.44E+02
 
 TH 5
+       -1.89E+01  0.00E+00 -1.36E+02 -8.51E+01  4.86E+02
 
 TH 6
+       -1.65E+01  0.00E+00  1.41E+01  1.98E+01 -3.85E+01  2.31E+02
 
 TH 7
+        2.25E-03  0.00E+00  7.02E-04  1.24E-02 -1.33E-03  1.62E-03  5.79E-06
 
 TH 8
+        5.03E+00  0.00E+00 -3.57E-01 -6.20E+00 -3.01E+00 -2.56E+00 -8.41E-05  1.24E+00
 
 TH 9
+        4.88E+00  0.00E+00  1.91E+00  1.07E+02  1.62E+00  9.00E+00  4.04E-02 -1.64E+00  2.83E+02
 
 TH10
+        9.01E+00  0.00E+00  1.76E+01 -1.27E+00 -6.65E+01 -3.51E+00 -6.25E-04  1.61E+00 -6.88E+00  1.06E+01
 
 TH11
+        5.71E+01  0.00E+00 -2.53E+01 -1.73E+01  2.69E+01  1.56E+01  4.14E-03  1.51E+01  1.82E+01  8.62E+00  2.16E+02
 
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
+        1.13E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.72E+00  0.00E+00  6.39E+01
 
 TH 4
+       -7.55E+00  0.00E+00 -1.50E+01  6.01E+02
 
 TH 5
+        1.89E+00  0.00E+00 -1.32E+02 -4.35E+01  4.91E+02
 
 TH 6
+        2.33E-01  0.00E+00 -4.56E-02 -2.24E+00 -8.08E-01  1.93E+02
 
 TH 7
+        1.48E-03  0.00E+00  9.22E-04 -2.93E-03 -2.02E-03  3.44E-04  2.47E-03
 
 TH 8
+       -1.49E-01  0.00E+00 -1.80E+01 -3.04E+00 -9.88E+00 -2.80E-02  1.82E-03  2.14E+01
 
 TH 9
+        2.92E+00  0.00E+00  5.66E+00 -5.71E-01 -1.05E+00 -5.15E-01  4.60E-02  7.13E-01  3.17E+02
 
 TH10
+        1.44E+00  0.00E+00  8.97E-01  1.31E-01 -6.72E+01  1.73E-01  1.59E-03  1.14E+01 -5.35E-01  5.48E+01
 
 TH11
+       -8.62E+00  0.00E+00 -9.69E+00 -1.02E+01 -1.11E+01  1.76E+00  2.67E-03  1.15E+01  9.86E+00  7.62E+00  2.08E+02
 
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
+        1.13E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.63E+00  0.00E+00  6.28E+01
 
 TH 4
+       -2.37E+01  0.00E+00 -2.35E+00  6.01E+02
 
 TH 5
+        3.69E+01  0.00E+00 -1.33E+02 -3.13E+00  4.99E+02
 
 TH 6
+        1.82E+01  0.00E+00 -2.48E+00 -2.16E+01  2.60E+01  1.58E+02
 
 TH 7
+        3.60E-03  0.00E+00  2.22E-03 -2.25E-02  4.61E-04 -2.73E-03  1.97E-05
 
 TH 8
+       -6.89E+00  0.00E+00 -1.78E+01 -2.20E+01 -8.96E+00  3.65E-01  2.53E-03  2.01E+01
 
 TH 9
+        4.68E+00  0.00E+00  1.89E+01 -1.27E+02 -1.88E+01 -1.20E+01  7.63E-02  3.84E+00  3.78E+02
 
 TH10
+       -1.67E+01  0.00E+00  3.05E+00 -1.57E+01 -1.07E+02 -1.63E+01 -4.38E-03  1.38E+01 -1.38E+01  9.20E+01
 
 TH11
+       -7.44E+01  0.00E+00  1.58E+00  4.19E+00 -4.67E+01 -8.43E+00  1.63E-03  8.70E+00  4.57E+00  2.29E+01  2.06E+02
 
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
 #CPUT: Total CPU Time in Seconds,       72.637
Stop Time:
Wed Sep 29 18:22:12 CDT 2021
