Wed Sep 29 08:38:54 CDT 2021
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
$DATA ../../../../data/int/D/dat37.csv ignore=@
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   28524.9187585340        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.7610E+02  4.2080E+02  7.5655E+01  2.7151E+02  1.5270E+02 -1.5473E+03 -8.6761E+02 -1.1056E+02 -1.3159E+03 -5.7118E+02
            -5.9780E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -892.736876968847        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1975E+00  1.9218E+00  7.7869E-01  2.0604E+00  9.1570E-01  4.6533E+00  6.1893E+00  1.0176E+00  2.6572E+00  1.9627E+00
             1.2858E+01
 PARAMETER:  2.8028E-01  7.5328E-01 -1.5014E-01  8.2289E-01  1.1934E-02  1.6376E+00  1.9228E+00  1.1745E-01  1.0773E+00  7.7433E-01
             2.6540E+00
 GRADIENT:  -9.6215E+00  2.5147E+01 -5.4814E+01  1.0694E+02  6.7292E+00  1.8585E+02  1.2680E+02  4.6904E+00  6.5010E+01  4.1363E+01
             5.0680E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1012.07027743575        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.1498E+00  1.6790E+00  4.4520E+01  2.7045E+00  2.5726E+00  3.0905E+00  1.1809E+01  8.4798E-01  2.5294E+00  2.4592E+00
             1.2736E+01
 PARAMETER:  2.3956E-01  6.1817E-01  3.8959E+00  1.0949E+00  1.0449E+00  1.2283E+00  2.5689E+00 -6.4899E-02  1.0280E+00  9.9984E-01
             2.6444E+00
 GRADIENT:  -3.1707E+01  2.9468E+01 -2.7974E+00  5.5101E+01  4.3199E+00  1.1546E+02  9.6935E+01  2.6315E-02  5.6291E+01  6.5371E+01
             5.4229E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1198.98380564250        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1227E+00  3.4379E-01  1.2421E+02  1.8646E+00  2.6637E+00  1.9527E+00  1.4284E+00  3.1423E-01  2.4214E+00  8.6333E-01
             9.5352E+00
 PARAMETER:  2.1577E-01 -9.6774E-01  4.9220E+00  7.2305E-01  1.0797E+00  7.6922E-01  4.5658E-01 -1.0576E+00  9.8433E-01 -4.6962E-02
             2.3550E+00
 GRADIENT:  -5.2243E+01 -9.9955E+00 -2.9430E+00 -2.7163E+01  4.7198E+01 -4.6027E+01  1.4872E+00  1.2109E-02 -1.1825E+01  8.5974E+00
             3.2956E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1229.14193587818        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      322
 NPARAMETR:  1.1672E+00  1.9064E-01  2.3099E+02  1.9173E+00  2.5300E+00  2.2374E+00  1.0384E+00  1.5091E-01  2.3705E+00  5.9128E-01
             8.1550E+00
 PARAMETER:  2.5461E-01 -1.5574E+00  5.5424E+00  7.5091E-01  1.0282E+00  9.0534E-01  1.3764E-01 -1.7911E+00  9.6311E-01 -4.2547E-01
             2.1986E+00
 GRADIENT:  -2.0510E+01 -6.4562E+00 -1.3125E+00 -4.5635E+01  2.8170E+01 -1.9874E+01  3.0097E-01  1.1030E-03 -1.9345E+01  2.2965E+00
             2.1357E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1241.26171144949        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      499
 NPARAMETR:  1.2070E+00  5.9491E-01  2.0548E+02  1.9693E+00  2.5400E+00  2.1752E+00  1.8912E+00  1.7776E-01  2.7053E+00  3.2473E-01
             8.1520E+00
 PARAMETER:  2.8816E-01 -4.1935E-01  5.4254E+00  7.7767E-01  1.0322E+00  8.7713E-01  7.3721E-01 -1.6273E+00  1.0952E+00 -1.0247E+00
             2.1983E+00
 GRADIENT:  -1.8761E+00  8.3367E-01 -1.4513E+00  1.8153E+01  1.3225E+01 -2.5154E+01  4.3750E+00  2.2518E-03 -2.0670E+01  5.2004E-02
            -1.7699E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1256.62340361976        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      677
 NPARAMETR:  1.2108E+00  1.4614E+00  9.6416E+01  1.0954E+00  2.5253E+00  2.2836E+00  8.5997E-01  6.2709E-01  4.3794E+00  2.5538E-01
             8.2150E+00
 PARAMETER:  2.9126E-01  4.7940E-01  4.6687E+00  1.9109E-01  1.0264E+00  9.2577E-01 -5.0860E-02 -3.6667E-01  1.5769E+00 -1.2650E+00
             2.2060E+00
 GRADIENT:  -3.3582E-01  5.1485E+00 -3.4353E+00  8.6854E+00  9.5858E+00 -6.8302E+00  1.7128E+00  9.7109E-02  2.2883E+00 -3.0078E-01
            -2.6902E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1259.51262071057        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      853
 NPARAMETR:  1.2089E+00  1.3970E+00  5.7728E+02  1.0172E+00  2.5074E+00  2.3262E+00  6.7984E-01  1.3112E-01  4.3640E+00  2.3527E-01
             8.1830E+00
 PARAMETER:  2.8973E-01  4.3430E-01  6.4583E+00  1.1702E-01  1.0192E+00  9.4425E-01 -2.8590E-01 -1.9316E+00  1.5734E+00 -1.3470E+00
             2.2021E+00
 GRADIENT:   1.3032E-01 -1.9838E+00 -2.2124E-01 -1.1917E+00 -8.4585E-01  3.2858E-01  7.7434E-01  2.2705E-04  8.8490E-01 -2.3667E-01
            -8.9785E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1259.74703858373        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  1.2129E+00  1.3110E+00  1.4007E+05  1.0984E+00  2.5271E+00  2.3205E+00  5.2743E-01  1.0000E-02  4.2095E+00  1.9846E-01
             8.1905E+00
 PARAMETER:  2.9305E-01  3.7081E-01  1.1950E+01  1.9381E-01  1.0271E+00  9.4178E-01 -5.3973E-01 -6.8628E+00  1.5374E+00 -1.5172E+00
             2.2030E+00
 GRADIENT:  -8.0085E+00  1.7650E+00 -5.7480E-04  5.6588E-01 -8.4213E-02  4.1556E+00 -1.7493E-01  0.0000E+00 -4.3975E-01 -1.4187E-01
             1.9986E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1259.84556028985        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1196
 NPARAMETR:  1.2096E+00  1.3101E+00  1.8248E+05  1.0964E+00  2.5279E+00  2.3351E+00  5.0595E-01  1.0000E-02  4.2262E+00  3.6737E-01
             8.1906E+00
 PARAMETER:  2.9027E-01  3.7008E-01  1.2214E+01  1.9206E-01  1.0274E+00  9.4805E-01 -5.8133E-01 -6.8628E+00  1.5413E+00 -9.0139E-01
             2.2030E+00
 GRADIENT:  -8.3772E+01  4.9641E+01 -3.7696E-04  1.3757E+01  6.9960E+00  9.9213E+01 -4.3501E+00  0.0000E+00  3.8794E+01 -7.3178E-03
             5.6806E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1259.86862030803        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1353
 NPARAMETR:  1.2133E+00  1.3184E+00  1.8604E+05  1.0944E+00  2.5186E+00  2.3372E+00  5.0683E-01  1.0000E-02  4.2234E+00  4.0182E-01
             8.1853E+00
 PARAMETER:  2.9332E-01  3.7640E-01  1.2234E+01  1.9022E-01  1.0237E+00  9.4894E-01 -5.7959E-01 -6.8628E+00  1.5406E+00 -8.1174E-01
             2.2023E+00
 GRADIENT:  -9.4811E+01  1.7540E+01 -3.5114E-04  4.7623E+00  1.4500E+00  7.7577E+01 -1.5733E+00  0.0000E+00 -5.8720E+00  1.6350E-02
             2.2661E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1259.87485440418        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1528
 NPARAMETR:  1.2134E+00  1.3214E+00  1.9407E+05  1.0922E+00  2.5077E+00  2.3357E+00  5.0761E-01  1.0000E-02  4.2300E+00  4.0014E-01
             8.1795E+00
 PARAMETER:  2.9344E-01  3.7872E-01  1.2276E+01  1.8824E-01  1.0194E+00  9.4830E-01 -5.7805E-01 -6.8628E+00  1.5422E+00 -8.1594E-01
             2.2016E+00
 GRADIENT:  -5.5707E+01  1.8554E+01 -3.0214E-04  4.7665E+00 -5.8208E-01  4.0869E+01 -1.8153E+00  0.0000E+00 -6.3469E+00 -6.1968E-02
             9.7287E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1259.88017181982        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1712             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2095E+00  1.3249E+00  2.1386E+05  1.0861E+00  2.5043E+00  2.3333E+00  5.1431E-01  1.0000E-02  4.2350E+00  4.1571E-01
             8.1739E+00
 PARAMETER:  2.9021E-01  3.8133E-01  1.2373E+01  1.8255E-01  1.0180E+00  9.4728E-01 -5.6493E-01 -6.8628E+00  1.5434E+00 -7.7776E-01
             2.2009E+00
 GRADIENT:  -3.7474E+01  3.1086E+01 -2.4468E-04  9.2715E+00  3.3899E+00  7.5492E+01 -2.1833E+00  0.0000E+00  4.7451E+01  2.8363E-02
             4.8324E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1259.88124011005        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     1907
 NPARAMETR:  1.2100E+00  1.3249E+00  1.7181E+05  1.0851E+00  2.5034E+00  2.3332E+00  5.2053E-01  1.0000E-02  4.2370E+00  4.2225E-01
             8.1734E+00
 PARAMETER:  2.9064E-01  3.8137E-01  1.2154E+01  1.8168E-01  1.0176E+00  9.4725E-01 -5.5291E-01 -6.8628E+00  1.5439E+00 -7.6217E-01
             2.2009E+00
 GRADIENT:  -1.1322E+01  2.9342E+00 -3.1898E-04  1.2422E+00 -2.6586E-01  6.2400E+00 -3.8441E-01  0.0000E+00 -3.5227E-01  1.4199E-02
             1.6937E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1259.88162816413        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2104             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2100E+00  1.3260E+00  3.0936E+05  1.0837E+00  2.5040E+00  2.3332E+00  5.2288E-01  1.0000E-02  4.2386E+00  4.1808E-01
             8.1738E+00
 PARAMETER:  2.9062E-01  3.8217E-01  1.2742E+01  1.8034E-01  1.0179E+00  9.4725E-01 -5.4840E-01 -6.8628E+00  1.5442E+00 -7.7207E-01
             2.2009E+00
 GRADIENT:  -3.2564E+02  1.0905E+02 -1.7158E-04  2.6071E+01  3.5424E+00  2.6017E+02 -9.5088E+00  0.0000E+00  1.7089E+01  1.8465E-01
             9.9826E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1259.88178938462        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2301
 NPARAMETR:  1.2101E+00  1.3268E+00  3.8926E+05  1.0835E+00  2.5040E+00  2.3333E+00  5.2290E-01  1.0000E-02  4.2394E+00  4.1812E-01
             8.1738E+00
 PARAMETER:  2.9066E-01  3.8279E-01  1.2972E+01  1.8023E-01  1.0179E+00  9.4726E-01 -5.4836E-01 -6.8628E+00  1.5444E+00 -7.7200E-01
             2.2009E+00
 GRADIENT:  -2.2017E+02 -1.5299E+01 -1.3739E-04 -1.4379E+01  5.5241E-02  1.8346E+02 -6.8736E-01  0.0000E+00  7.8634E+00 -1.3433E-01
             7.7847E+01

0ITERATION NO.:   78    OBJECTIVE VALUE:  -1259.88187246905        NO. OF FUNC. EVALS.: 108
 CUMULATIVE NO. OF FUNC. EVALS.:     2409
 NPARAMETR:  1.2099E+00  1.3276E+00  5.3118E+05  1.0831E+00  2.5047E+00  2.3332E+00  5.2498E-01  1.0000E-02  4.2419E+00  4.2197E-01
             8.1744E+00
 PARAMETER:  2.9068E-01  3.8274E-01  1.3244E+01  1.8012E-01  1.0180E+00  9.4727E-01 -5.4780E-01 -6.8628E+00  1.5447E+00 -7.7051E-01
             2.2010E+00
 GRADIENT:   1.5874E-02 -1.0821E-01 -2.4490E-04  1.2440E-02 -2.5685E-02  2.7574E-03 -1.6136E-02  0.0000E+00 -3.4031E-02 -1.5189E-02
            -2.7056E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2409
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.0288E-03 -3.5978E-02  1.9473E-09  1.1521E-02 -4.8169E-04
 SE:             2.8914E-02  7.8153E-03  6.5293E-10  2.7034E-02  6.5824E-03
 N:                     100         100         100         100         100

 P VAL.:         7.5484E-01  4.1591E-06  2.8601E-03  6.6998E-01  9.4166E-01

 ETASHRINKSD(%)  3.1344E+00  7.3818E+01  1.0000E+02  9.4326E+00  7.7948E+01
 ETASHRINKVR(%)  6.1705E+00  9.3145E+01  1.0000E+02  1.7975E+01  9.5137E+01
 EBVSHRINKSD(%)  2.8622E+00  8.0998E+01  1.0000E+02  5.6547E+00  7.8708E+01
 EBVSHRINKVR(%)  5.6425E+00  9.6389E+01  1.0000E+02  1.0990E+01  9.5466E+01
 RELATIVEINF(%)  9.4168E+01  1.7753E+00  0.0000E+00  4.3997E+01  3.8422E+00
 EPSSHRINKSD(%)  6.8367E+00
 EPSSHRINKVR(%)  1.3206E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1259.8818724690543     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       394.20748729935644     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    79.39
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1259.882       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.21E+00  1.33E+00  5.11E+05  1.08E+00  2.50E+00  2.33E+00  5.23E-01  1.00E-02  4.24E+00  4.19E-01  8.17E+00
 


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
+        1.31E+02
 
 TH 2
+       -2.84E+00  1.26E+02
 
 TH 3
+        8.80E-08 -2.93E-08  5.73E-15
 
 TH 4
+       -1.77E+00  4.17E+01  2.85E-09  5.37E+01
 
 TH 5
+       -2.07E+00 -1.56E+01  7.11E-09 -6.55E+00  4.17E+01
 
 TH 6
+        6.69E-01  1.41E+00 -1.81E-08 -1.25E+00 -4.92E-01  3.25E+01
 
 TH 7
+        2.04E+00 -3.18E+01 -7.35E-08 -9.78E+00  3.29E+00 -8.67E-01  1.71E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.04E-01 -1.74E+01 -2.32E-10  2.79E+00  1.62E+00 -4.80E-01  3.88E+00  0.00E+00  7.70E+00
 
 TH10
+       -1.95E+00 -7.01E+00 -1.55E-07  2.26E+00  5.42E+00  2.13E-01  2.27E+00  0.00E+00  9.69E-01  3.99E+00
 
 TH11
+       -5.54E+00 -9.30E+00  8.03E-10 -4.30E+00  3.12E-01  6.65E-01  1.80E+00  0.00E+00  9.64E-01  1.45E+00  1.59E+01
 
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
 #CPUT: Total CPU Time in Seconds,       96.091
Stop Time:
Wed Sep 29 08:40:32 CDT 2021
