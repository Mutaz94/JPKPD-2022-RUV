Wed Sep 29 12:38:45 CDT 2021
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
$DATA ../../../../data/spa/A2/dat20.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -887.310015969539        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5716E+02  6.4656E+00  3.4739E+01  5.8103E+00  8.4779E+01  3.1453E+01 -6.4855E+00 -4.0745E+01 -2.2025E+01 -1.6915E+01
            -1.3887E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1354.57911482581        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0876E+00  1.0178E+00  9.0478E-01  1.0398E+00  8.9218E-01  1.1631E+00  9.5670E-01  1.1746E+00  1.0394E+00  8.0013E-01
             2.2064E+00
 PARAMETER:  1.8400E-01  1.1769E-01 -6.6428E-05  1.3906E-01 -1.4087E-02  2.5106E-01  5.5735E-02  2.6092E-01  1.3860E-01 -1.2298E-01
             8.9134E-01
 GRADIENT:   3.1333E+02  1.4103E+01  1.6389E+01  1.2198E+01 -1.3572E+01  5.5742E+01  3.9994E+00 -8.5121E+00  4.0713E+00  1.3153E+01
            -9.8649E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1369.75275291690        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0480E+00  8.5884E-01  5.2382E-01  1.1102E+00  6.1055E-01  9.6362E-01  1.3501E+00  6.2868E-01  7.5363E-01  2.1502E-01
             2.2636E+00
 PARAMETER:  1.4692E-01 -5.2170E-02 -5.4661E-01  2.0453E-01 -3.9340E-01  6.2939E-02  4.0020E-01 -3.6413E-01 -1.8285E-01 -1.4370E+00
             9.1696E-01
 GRADIENT:   2.4733E+02  4.0967E+01  1.3899E+01  7.5557E+01  1.5812E+01 -2.1747E+00  1.0785E+01 -1.3776E+01 -2.3808E+01  1.9694E-01
            -9.4708E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1385.12051353897        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.8940E-01  9.1004E-01  5.1376E-01  1.0403E+00  6.2970E-01  9.4616E-01  1.1085E+00  9.8694E-01  8.3921E-01  1.8814E-01
             2.4054E+00
 PARAMETER:  8.9342E-02  5.7354E-03 -5.6600E-01  1.3947E-01 -3.6251E-01  4.4656E-02  2.0302E-01  8.6850E-02 -7.5293E-02 -1.5705E+00
             9.7773E-01
 GRADIENT:   1.2430E+01 -9.0995E+00  1.0671E+00 -5.9528E+00  1.1194E+01 -2.4279E+00 -8.3233E-01 -3.0012E+00 -6.4048E+00  1.3630E+00
            -3.7111E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1391.83179915864        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  9.6876E-01  1.1127E+00  2.3934E-01  8.5358E-01  4.9559E-01  9.0300E-01  8.4536E-01  1.1512E+00  8.4034E-01  8.7337E-02
             2.5342E+00
 PARAMETER:  6.8262E-02  2.0679E-01 -1.3299E+00 -5.8313E-02 -6.0201E-01 -2.0338E-03 -6.7993E-02  2.4081E-01 -7.3949E-02 -2.3380E+00
             1.0299E+00
 GRADIENT:  -2.4870E+01  1.0886E+02  4.6509E+01  1.8136E+01 -1.2281E+02 -1.7951E+01 -7.1268E+00 -1.1851E+01 -1.3560E+01  4.6266E-01
             3.6823E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1421.62132882730        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      651
 NPARAMETR:  9.6981E-01  1.4822E+00  1.9388E-01  6.4646E-01  6.5975E-01  9.6451E-01  7.8415E-01  2.3004E+00  1.1245E+00  1.0213E-01
             2.0268E+00
 PARAMETER:  6.9342E-02  4.9356E-01 -1.5405E+00 -3.3625E-01 -3.1589E-01  6.3861E-02 -1.4316E-01  9.3306E-01  2.1731E-01 -2.1815E+00
             8.0648E-01
 GRADIENT:   1.3802E+01  5.9881E+01  1.4702E+01  2.7050E+01 -2.9981E+01  1.1504E+01  3.7717E+00  2.4372E+00 -1.2497E+01  4.4461E-01
             8.7502E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1429.69308182497        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      826
 NPARAMETR:  9.5630E-01  1.8170E+00  1.0091E-01  4.3371E-01  7.9794E-01  9.2479E-01  6.8822E-01  1.9908E+00  1.8188E+00  2.3807E-02
             2.0940E+00
 PARAMETER:  5.5311E-02  6.9716E-01 -2.1935E+00 -7.3539E-01 -1.2572E-01  2.1806E-02 -2.7365E-01  7.8854E-01  6.9818E-01 -3.6378E+00
             8.3908E-01
 GRADIENT:  -1.1964E+01  4.0332E+01 -2.0172E+00  2.4078E+01 -9.1399E+00 -4.1735E+00  1.0114E+00  1.1578E+00  1.8313E+00  6.7893E-03
             6.7037E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1432.04422824309        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      991
 NPARAMETR:  9.5819E-01  2.0109E+00  6.2505E-02  3.1850E-01  9.0987E-01  9.2820E-01  6.6254E-01  1.7523E+00  2.3415E+00  1.0000E-02
             2.1525E+00
 PARAMETER:  5.7288E-02  7.9858E-01 -2.6725E+00 -1.0441E+00  5.5423E-03  2.5492E-02 -3.1168E-01  6.6096E-01  9.5080E-01 -4.8376E+00
             8.6662E-01
 GRADIENT:   7.2728E+01  2.3725E+02  7.9923E+00  3.6311E+01  7.4305E+00  2.9770E+00  4.6718E+00  5.5200E+00  5.1359E+00  0.0000E+00
             9.1921E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1432.66955527600        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1170
 NPARAMETR:  9.6451E-01  1.9795E+00  6.4657E-02  3.1261E-01  8.9437E-01  9.3140E-01  6.5463E-01  1.7485E+00  2.3411E+00  1.0000E-02
             2.1528E+00
 PARAMETER:  6.3868E-02  7.8285E-01 -2.6387E+00 -1.0628E+00 -1.1636E-02  2.8937E-02 -3.2368E-01  6.5873E-01  9.5062E-01 -4.8376E+00
             8.6675E-01
 GRADIENT:   1.1118E+01 -5.7792E+00 -3.5924E+00  8.9449E+00 -8.0081E+00 -8.4393E-01 -2.1373E+00  4.8485E+00  1.1189E+00  0.0000E+00
             5.7292E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1433.34048685044        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1347
 NPARAMETR:  9.6240E-01  1.9878E+00  6.4356E-02  3.0537E-01  9.0018E-01  9.3473E-01  6.6391E-01  1.5279E+00  2.3545E+00  1.0000E-02
             2.1311E+00
 PARAMETER:  6.1673E-02  7.8703E-01 -2.6433E+00 -1.0862E+00 -5.1586E-03  3.2499E-02 -3.0961E-01  5.2389E-01  9.5633E-01 -4.8376E+00
             8.5665E-01
 GRADIENT:   5.3374E+00 -2.5425E+00  8.6259E-01  5.8706E+00 -2.9334E+00  1.9886E-01  5.3215E-01  1.1863E+00 -8.3564E-01  0.0000E+00
             1.1731E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1434.74512646423        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1516
 NPARAMETR:  9.5603E-01  2.0469E+00  4.8684E-02  2.5368E-01  9.3706E-01  9.3463E-01  6.6040E-01  1.0528E+00  2.6849E+00  1.0000E-02
             2.1501E+00
 PARAMETER:  5.5030E-02  8.1632E-01 -2.9224E+00 -1.2717E+00  3.4995E-02  3.2393E-02 -3.1492E-01  1.5144E-01  1.0876E+00 -4.8376E+00
             8.6552E-01
 GRADIENT:   7.4563E+01  2.0773E+02  1.4730E+01  1.0621E+01 -9.2064E+00  6.2220E+00  4.9406E+00 -1.5464E+00  7.4370E+00  0.0000E+00
             2.7048E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1434.76319424026        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1703
 NPARAMETR:  9.5651E-01  2.0554E+00  4.9394E-02  2.5531E-01  9.3660E-01  9.3336E-01  6.6143E-01  1.0520E+00  2.6996E+00  1.0000E-02
             2.1410E+00
 PARAMETER:  5.5532E-02  8.2045E-01 -2.9079E+00 -1.2653E+00  3.4498E-02  3.1033E-02 -3.1334E-01  1.5068E-01  1.0931E+00 -4.8376E+00
             8.6125E-01
 GRADIENT:  -6.5591E+00 -3.3927E+00  1.0525E+00 -4.7859E-01 -1.7777E+01 -2.4621E-01  1.2979E+00 -1.9111E+00 -7.2678E-01  0.0000E+00
            -6.3407E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1434.98708924998        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     1850
 NPARAMETR:  9.5803E-01  2.0580E+00  4.9233E-02  2.5617E-01  9.4127E-01  9.3269E-01  6.5597E-01  1.1320E+00  2.7131E+00  1.0000E-02
             2.1500E+00
 PARAMETER:  5.7120E-02  8.2173E-01 -2.9112E+00 -1.2619E+00  3.9475E-02  3.0318E-02 -3.2164E-01  2.2402E-01  1.0981E+00 -4.8376E+00
             8.6548E-01
 GRADIENT:   7.8170E+01  2.2501E+02  1.5534E+01  1.3876E+01 -1.0306E+01  5.3430E+00  3.9495E+00 -7.6835E-01  9.2104E+00  0.0000E+00
             2.9006E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1435.07587438579        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2028
 NPARAMETR:  9.5871E-01  2.0594E+00  4.9589E-02  2.5486E-01  9.4511E-01  9.3402E-01  6.5611E-01  1.1703E+00  2.7202E+00  1.0000E-02
             2.1581E+00
 PARAMETER:  5.7837E-02  8.2243E-01 -2.9040E+00 -1.2671E+00  4.3548E-02  3.1743E-02 -3.2143E-01  2.5729E-01  1.1007E+00 -4.8376E+00
             8.6923E-01
 GRADIENT:  -2.4268E+00 -1.2229E+01  7.2111E-01  1.3643E+00 -8.9493E+00  2.4739E-01  4.2291E-01 -5.4686E-01 -3.1660E-01  0.0000E+00
            -1.6176E+00

0ITERATION NO.:   67    OBJECTIVE VALUE:  -1435.08563775147        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     2093
 NPARAMETR:  9.5947E-01  2.0735E+00  5.0754E-02  2.5718E-01  9.4484E-01  9.3298E-01  6.5431E-01  1.1709E+00  2.7453E+00  1.0000E-02
             2.1434E+00
 PARAMETER:  5.7819E-02  8.2262E-01 -2.9041E+00 -1.2682E+00  4.4068E-02  3.1628E-02 -3.2170E-01  2.5987E-01  1.1010E+00 -4.8376E+00
             8.6940E-01
 GRADIENT:  -5.8883E+03 -1.4409E+03 -3.9510E+02 -9.3271E+02  1.1767E+04  1.1242E-01  3.6625E-01  4.5213E+03 -1.0810E+03  0.0000E+00
             1.3483E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2093
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3959E-03 -1.9290E-02  1.5953E-02  2.0523E-02 -4.0591E-04
 SE:             2.9291E-02  2.6355E-02  7.5036E-03  2.1897E-02  2.6335E-04
 N:                     100         100         100         100         100

 P VAL.:         9.3481E-01  4.6422E-01  3.3504E-02  3.4863E-01  1.2324E-01

 ETASHRINKSD(%)  1.8718E+00  1.1707E+01  7.4862E+01  2.6641E+01  9.9118E+01
 ETASHRINKVR(%)  3.7085E+00  2.2043E+01  9.3681E+01  4.6185E+01  9.9992E+01
 EBVSHRINKSD(%)  2.0300E+00  1.2521E+01  7.4966E+01  2.3087E+01  9.9121E+01
 EBVSHRINKVR(%)  4.0189E+00  2.3474E+01  9.3733E+01  4.0843E+01  9.9992E+01
 RELATIVEINF(%)  9.3743E+01  2.1986E+01  4.3044E+00  2.2310E+01  2.0217E-03
 EPSSHRINKSD(%)  3.4175E+01
 EPSSHRINKVR(%)  5.6670E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1435.0856377514726     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -699.93481118773445     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1435.086       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.59E-01  2.06E+00  4.96E-02  2.55E-01  9.46E-01  9.34E-01  6.56E-01  1.17E+00  2.72E+00  1.00E-02  2.16E+00
 


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
+        3.20E+06
 
 TH 2
+       -7.84E+01  1.07E+04
 
 TH 3
+       -8.35E+02  7.47E+02  1.37E+06
 
 TH 4
+        1.15E+04  7.71E+02  5.95E+03  2.87E+05
 
 TH 5
+       -1.30E+03 -4.70E+02 -1.12E+03  6.86E+01  3.29E+06
 
 TH 6
+       -4.86E+02 -3.46E+01 -3.31E+02 -1.60E+02  4.93E+02  2.08E+02
 
 TH 7
+       -1.56E+02 -7.79E+00 -1.30E+02 -3.12E+01  1.82E+02  2.98E+00  6.61E+05
 
 TH 8
+        3.92E+03  2.19E+02  2.44E+03  1.17E+03 -3.98E+03  1.54E+02 -4.57E+05  3.15E+05
 
 TH 9
+       -1.63E+01  5.85E+03  1.52E+03  3.98E+02 -5.43E+01 -1.48E+01 -2.20E+00  1.28E+02  3.36E+03
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.52E+01  4.78E+01  8.58E+02 -5.93E+02  4.57E+01  2.73E+01  2.63E+01 -1.98E+02  3.90E+01  0.00E+00  8.36E+03
 
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
 #CPUT: Total CPU Time in Seconds,       36.327
Stop Time:
Wed Sep 29 12:39:24 CDT 2021
