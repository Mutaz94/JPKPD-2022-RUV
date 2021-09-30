Thu Sep 30 02:29:00 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat94.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2098.24363322063        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2673E+02 -2.5245E+01 -4.0720E+01  1.1455E+01  9.9818E+01  3.9141E+01 -1.2187E+00  1.3399E+01 -1.3089E+01 -4.9713E+00
            -8.6434E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2105.53458300220        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0011E+00  9.4896E-01  1.0337E+00  1.0818E+00  9.2214E-01  1.0357E+00  1.0033E+00  9.2057E-01  1.1112E+00  9.7857E-01
             1.0283E+00
 PARAMETER:  1.0115E-01  4.7612E-02  1.3310E-01  1.7862E-01  1.8945E-02  1.3510E-01  1.0325E-01  1.7235E-02  2.0547E-01  7.8332E-02
             1.2794E-01
 GRADIENT:   3.5873E+00  7.9713E+00 -2.9164E+00  1.8538E+00  5.2697E+00  7.0391E+00  5.0986E+00  1.0230E+01  1.2160E+01  2.4239E+00
             1.4965E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2109.44797292469        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  1.0012E+00  8.0549E-01  9.2907E-01  1.1815E+00  8.1985E-01  9.9471E-01  1.0027E+00  4.8863E-01  1.0556E+00  1.0530E+00
             1.0088E+00
 PARAMETER:  1.0117E-01 -1.1630E-01  2.6427E-02  2.6679E-01 -9.8631E-02  9.4692E-02  1.0267E-01 -6.1616E-01  1.5415E-01  1.5163E-01
             1.0877E-01
 GRADIENT:   5.7784E+00  1.8065E+01 -1.9653E+01  3.7147E+01 -1.9498E+00 -8.4704E+00  2.5593E+00  3.9323E+00  1.0817E+01  1.8863E+01
             7.2016E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2112.60988372117        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  9.9798E-01  6.8074E-01  9.8097E-01  1.2362E+00  7.9672E-01  1.0143E+00  1.1481E+00  3.0523E-01  9.6801E-01  1.0032E+00
             1.0020E+00
 PARAMETER:  9.7980E-02 -2.8457E-01  8.0785E-02  3.1201E-01 -1.2725E-01  1.1423E-01  2.3811E-01 -1.0867E+00  6.7492E-02  1.0318E-01
             1.0197E-01
 GRADIENT:   2.7845E+00  1.0408E+01  1.1213E+01  1.7359E+00 -1.6765E+01  2.2623E-01  1.0363E+00  1.6880E-01  1.1919E+00  1.8368E+00
            -2.6156E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2113.15377054527        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      712             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9596E-01  5.9778E-01  9.9012E-01  1.2823E+00  7.8495E-01  1.0124E+00  1.0422E+00  2.5349E-01  9.4595E-01  1.0144E+00
             1.0064E+00
 PARAMETER:  9.5952E-02 -4.1453E-01  9.0067E-02  3.4869E-01 -1.4214E-01  1.1229E-01  1.4135E-01 -1.2724E+00  4.4433E-02  1.1433E-01
             1.0641E-01
 GRADIENT:   4.1311E+02  6.2429E+01 -1.0415E+00  4.1589E+02  1.8270E+01  4.9117E+01  1.6442E+00  2.3767E-01  8.6631E+00  9.9458E-01
             1.3851E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2113.17770703541        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  9.9591E-01  5.9779E-01  9.9606E-01  1.2821E+00  7.8570E-01  1.0124E+00  9.9925E-01  2.7924E-01  9.5448E-01  1.0166E+00
             1.0061E+00
 PARAMETER:  9.5905E-02 -4.1452E-01  9.6056E-02  3.4849E-01 -1.4118E-01  1.1231E-01  9.9250E-02 -1.1757E+00  5.3408E-02  1.1643E-01
             1.0610E-01
 GRADIENT:   8.2147E-01  3.1909E+00 -1.1064E+00  6.0161E-01 -3.2831E-01  1.0207E-01  1.0180E-03  2.3573E-02 -3.4832E-01 -2.1005E-02
             1.7625E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2113.20061131135        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  9.9592E-01  5.9172E-01  1.0074E+00  1.2813E+00  7.9171E-01  1.0122E+00  9.4556E-01  2.8811E-01  9.6016E-01  1.0259E+00
             1.0061E+00
 PARAMETER:  9.5907E-02 -4.2473E-01  1.0733E-01  3.4785E-01 -1.3356E-01  1.1215E-01  4.4019E-02 -1.1444E+00  5.9348E-02  1.2561E-01
             1.0610E-01
 GRADIENT:   1.2169E+00  3.5695E-02 -2.2419E+00 -7.5817E+00  3.6695E+00  1.2093E-01  8.0831E-02  1.9586E-02  7.2747E-01 -4.7156E-02
             2.6430E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2113.41457963238        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1249
 NPARAMETR:  9.9191E-01  4.8022E-01  1.0411E+00  1.3578E+00  7.7008E-01  1.0090E+00  7.6749E-01  3.3608E-01  9.2889E-01  1.0368E+00
             1.0060E+00
 PARAMETER:  9.1881E-02 -6.3352E-01  1.4023E-01  4.0585E-01 -1.6126E-01  1.0895E-01 -1.6464E-01 -9.9041E-01  2.6232E-02  1.3614E-01
             1.0601E-01
 GRADIENT:  -3.3566E+00  4.6950E+00  2.0611E+00  9.7724E+00 -7.8316E+00 -3.9308E-01 -1.5648E-01 -1.6159E-01 -3.5369E-01  6.2991E-01
            -1.1604E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2113.47504950637        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  9.9127E-01  3.9737E-01  1.0721E+00  1.4072E+00  7.6151E-01  1.0081E+00  6.3799E-01  4.4350E-01  9.0328E-01  1.0289E+00
             1.0069E+00
 PARAMETER:  9.1229E-02 -8.2290E-01  1.6958E-01  4.4158E-01 -1.7246E-01  1.0808E-01 -3.4943E-01 -7.1305E-01 -1.7240E-03  1.2846E-01
             1.0685E-01
 GRADIENT:  -1.3028E+00  2.7206E+00 -6.3037E-02  6.8344E+00 -6.4156E+00 -1.9601E-01 -8.4605E-02  2.9427E-02 -3.8199E-01  6.3695E-01
            -3.4596E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2113.51578217626        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1609
 NPARAMETR:  9.9172E-01  3.8307E-01  1.0786E+00  1.4116E+00  7.6289E-01  1.0084E+00  6.7628E-01  4.5170E-01  8.9776E-01  1.0263E+00
             1.0071E+00
 PARAMETER:  9.1687E-02 -8.5955E-01  1.7564E-01  4.4472E-01 -1.7064E-01  1.0837E-01 -2.9114E-01 -6.9473E-01 -7.8510E-03  1.2599E-01
             1.0711E-01
 GRADIENT:   3.9524E-01  4.1285E-01 -1.9741E+00 -2.1414E+00 -1.1755E+00  2.4542E-02 -3.6390E-02 -3.7277E-02 -1.1371E-01 -2.4773E-01
            -2.0400E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2113.52699893545        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1791
 NPARAMETR:  9.9172E-01  3.7996E-01  1.0847E+00  1.4126E+00  7.6447E-01  1.0085E+00  7.2118E-01  4.5947E-01  8.9575E-01  1.0277E+00
             1.0073E+00
 PARAMETER:  9.1687E-02 -8.6768E-01  1.8132E-01  4.4545E-01 -1.6858E-01  1.0851E-01 -2.2686E-01 -6.7768E-01 -1.0098E-02  1.2732E-01
             1.0725E-01
 GRADIENT:   5.7665E-01  6.8970E-02 -6.2707E-01 -5.1358E+00 -1.8377E+00  1.0142E-01 -5.3161E-05 -5.0270E-02  3.9737E-03 -1.2797E-01
            -1.8323E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2113.53264399415        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1978
 NPARAMETR:  9.9214E-01  3.8135E-01  1.0882E+00  1.4117E+00  7.6584E-01  1.0087E+00  7.1733E-01  4.6352E-01  8.9594E-01  1.0281E+00
             1.0071E+00
 PARAMETER:  9.2108E-02 -8.6403E-01  1.8455E-01  4.4482E-01 -1.6678E-01  1.0863E-01 -2.3222E-01 -6.6891E-01 -9.8852E-03  1.2771E-01
             1.0707E-01
 GRADIENT:   1.4476E+00  1.7578E-01  7.4268E-01 -5.9498E+00 -2.9531E+00  1.4138E-01 -8.3409E-04 -8.0005E-02 -1.1524E-01 -2.5308E-01
            -1.5346E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2113.53672046129        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2165
 NPARAMETR:  9.9215E-01  3.8226E-01  1.0901E+00  1.4114E+00  7.6712E-01  1.0087E+00  7.1291E-01  4.6635E-01  8.9631E-01  1.0289E+00
             1.0070E+00
 PARAMETER:  9.2121E-02 -8.6166E-01  1.8624E-01  4.4457E-01 -1.6512E-01  1.0866E-01 -2.3840E-01 -6.6282E-01 -9.4651E-03  1.2849E-01
             1.0702E-01
 GRADIENT:   1.4531E+00  1.5817E-01  6.5198E-01 -5.9531E+00 -2.5013E+00  1.4278E-01 -3.1573E-04 -6.4950E-02 -1.1567E-01 -2.3348E-01
            -1.3975E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2113.53976908969        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2352
 NPARAMETR:  9.9216E-01  3.8305E-01  1.0913E+00  1.4110E+00  7.6818E-01  1.0087E+00  7.0755E-01  4.6835E-01  8.9669E-01  1.0297E+00
             1.0070E+00
 PARAMETER:  9.2133E-02 -8.5959E-01  1.8736E-01  4.4433E-01 -1.6373E-01  1.0867E-01 -2.4594E-01 -6.5853E-01 -9.0437E-03  1.2927E-01
             1.0699E-01
 GRADIENT:   1.4566E+00  1.0815E-01  3.5051E-01 -5.9396E+00 -1.8770E+00  1.4325E-01 -6.2905E-04 -4.4614E-02 -1.1671E-01 -1.9036E-01
            -1.0670E-01

0ITERATION NO.:   69    OBJECTIVE VALUE:  -2113.54093352651        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     2487
 NPARAMETR:  9.9204E-01  3.8342E-01  1.0913E+00  1.4114E+00  7.6915E-01  1.0086E+00  7.0647E-01  4.7084E-01  8.9702E-01  1.0305E+00
             1.0071E+00
 PARAMETER:  9.2009E-02 -8.5863E-01  1.8738E-01  4.4459E-01 -1.6247E-01  1.0861E-01 -2.4748E-01 -6.5324E-01 -8.6797E-03  1.3007E-01
             1.0706E-01
 GRADIENT:  -3.0886E-01 -1.0271E-01 -6.5300E-01  1.6406E+00 -7.1246E-02 -3.4749E-02 -4.0861E-04  1.9762E-02 -7.4681E-02 -2.4775E-02
             7.3187E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2487
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.6070E-04 -8.4338E-03 -1.4498E-02 -2.3700E-03 -1.8347E-02
 SE:             2.9861E-02  4.7937E-03  9.5255E-03  2.9264E-02  2.5354E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8235E-01  7.8520E-02  1.2801E-01  9.3545E-01  4.6931E-01

 ETASHRINKSD(%)  1.0000E-10  8.3940E+01  6.8088E+01  1.9613E+00  1.5059E+01
 ETASHRINKVR(%)  1.0000E-10  9.7421E+01  8.9816E+01  3.8842E+00  2.7851E+01
 EBVSHRINKSD(%)  3.4478E-01  8.4837E+01  7.0434E+01  2.1151E+00  1.2684E+01
 EBVSHRINKVR(%)  6.8838E-01  9.7701E+01  9.1258E+01  4.1855E+00  2.3759E+01
 RELATIVEINF(%)  9.7824E+01  1.9693E-01  1.7494E+00  1.0749E+01  1.3187E+01
 EPSSHRINKSD(%)  3.3264E+01
 EPSSHRINKVR(%)  5.5463E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2113.5409335265149     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1194.6024003218422     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.25
 Elapsed covariance  time in seconds:     6.42
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2113.541       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.92E-01  3.83E-01  1.09E+00  1.41E+00  7.69E-01  1.01E+00  7.06E-01  4.71E-01  8.97E-01  1.03E+00  1.01E+00
 


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
 
         3.09E-02  2.79E-01  1.51E-01  1.79E-01  7.43E-02  7.31E-02  1.18E+00  3.09E-01  1.23E-01  1.17E-01  5.10E-02
 


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
+        9.53E-04
 
 TH 2
+        1.88E-03  7.77E-02
 
 TH 3
+       -4.63E-04 -2.30E-02  2.28E-02
 
 TH 4
+       -1.22E-03 -4.85E-02  1.63E-02  3.19E-02
 
 TH 5
+        2.00E-04  1.05E-02  4.10E-03 -5.60E-03  5.52E-03
 
 TH 6
+       -7.74E-05  4.20E-05 -4.11E-04  3.11E-04 -5.83E-05  5.34E-03
 
 TH 7
+       -2.54E-03 -6.41E-02  2.88E-02  4.13E-02 -1.03E-03  6.25E-03  1.39E+00
 
 TH 8
+       -8.66E-04 -3.71E-02  3.14E-02  2.52E-02  2.13E-03 -8.94E-04  9.31E-02  9.56E-02
 
 TH 9
+        5.00E-04  2.76E-02 -7.64E-03 -1.74E-02  3.88E-03 -2.42E-05 -5.34E-02 -1.60E-02  1.50E-02
 
 TH10
+       -1.62E-04 -1.23E-03  3.01E-03  1.17E-03  2.54E-03  4.39E-04 -1.48E-02 -1.22E-02  1.15E-03  1.38E-02
 
 TH11
+        8.36E-05  2.05E-03 -6.39E-04 -1.32E-03 -9.37E-05 -2.34E-04 -1.21E-03 -4.71E-04  2.20E-04 -7.05E-04  2.60E-03
 
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
+        3.09E-02
 
 TH 2
+        2.18E-01  2.79E-01
 
 TH 3
+       -9.92E-02 -5.46E-01  1.51E-01
 
 TH 4
+       -2.20E-01 -9.73E-01  6.04E-01  1.79E-01
 
 TH 5
+        8.73E-02  5.05E-01  3.66E-01 -4.22E-01  7.43E-02
 
 TH 6
+       -3.43E-02  2.06E-03 -3.72E-02  2.38E-02 -1.07E-02  7.31E-02
 
 TH 7
+       -6.97E-02 -1.95E-01  1.62E-01  1.96E-01 -1.18E-02  7.25E-02  1.18E+00
 
 TH 8
+       -9.07E-02 -4.31E-01  6.72E-01  4.57E-01  9.27E-02 -3.95E-02  2.55E-01  3.09E-01
 
 TH 9
+        1.32E-01  8.08E-01 -4.12E-01 -7.94E-01  4.26E-01 -2.69E-03 -3.70E-01 -4.23E-01  1.23E-01
 
 TH10
+       -4.46E-02 -3.74E-02  1.70E-01  5.57E-02  2.91E-01  5.12E-02 -1.07E-01 -3.37E-01  8.01E-02  1.17E-01
 
 TH11
+        5.31E-02  1.44E-01 -8.30E-02 -1.44E-01 -2.47E-02 -6.29E-02 -2.02E-02 -2.99E-02  3.52E-02 -1.18E-01  5.10E-02
 
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
+        1.14E+03
 
 TH 2
+       -5.95E+01  4.13E+02
 
 TH 3
+       -9.67E+01  1.76E+02  4.69E+02
 
 TH 4
+        4.79E+01  4.05E+02 -9.59E+01  7.09E+02
 
 TH 5
+        1.43E+02 -4.67E+02 -7.12E+02 -1.71E+01  1.59E+03
 
 TH 6
+        9.03E+00 -2.81E+01  1.46E+01 -5.10E+01  3.00E-02  1.94E+02
 
 TH 7
+        1.48E+00  5.12E-01  1.86E+00  1.01E+00 -4.75E+00 -1.12E+00  9.11E-01
 
 TH 8
+        1.26E+01 -2.61E+00 -5.81E+01  9.36E+00  2.10E+01 -1.07E+00 -5.52E-01  3.25E+01
 
 TH 9
+        5.89E+01 -8.20E+01 -6.26E+01  4.42E+01  7.82E+01 -4.59E+00  4.95E+00  1.02E+01  2.43E+02
 
 TH10
+        7.14E+00  5.16E+01  1.37E+00  8.37E+00 -1.58E+02 -8.48E+00  5.69E-01  3.58E+01 -1.43E+01  1.39E+02
 
 TH11
+        1.67E+01 -7.39E+01 -9.69E+01  1.06E+01  1.89E+02  1.45E+01  3.05E-01  7.27E+00  5.18E+01  3.02E+00  4.30E+02
 
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
 #CPUT: Total CPU Time in Seconds,       43.749
Stop Time:
Thu Sep 30 02:29:45 CDT 2021
