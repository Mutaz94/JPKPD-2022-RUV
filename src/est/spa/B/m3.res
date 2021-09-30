Wed Sep 29 10:54:03 CDT 2021
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
$DATA ../../../../data/spa/B/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1718.52591078141        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0483E+02 -4.2664E+01 -4.9017E+01  2.0417E+01  5.0633E+01  7.0798E+01  7.2858E+00  1.1584E+01  1.8811E+01  3.1068E+01
             3.2591E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1728.22698797019        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0217E+00  1.0958E+00  1.1255E+00  1.0014E+00  1.0346E+00  1.0001E+00  9.5510E-01  9.5906E-01  9.6231E-01  8.2771E-01
             1.0298E+00
 PARAMETER:  1.2144E-01  1.9147E-01  2.1819E-01  1.0141E-01  1.3398E-01  1.0012E-01  5.4062E-02  5.8197E-02  6.1579E-02 -8.9088E-02
             1.2941E-01
 GRADIENT:  -8.4896E-01  6.6299E+00  5.8418E+00  7.4589E+00 -8.2755E-01 -1.3434E-01  1.8299E+00 -2.5841E+00 -2.4572E+00  1.5078E+00
             4.4599E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1728.70314521588        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0278E+00  1.1863E+00  1.0234E+00  9.4500E-01  1.0195E+00  1.0056E+00  9.1784E-01  1.0131E+00  1.0002E+00  7.5963E-01
             1.0160E+00
 PARAMETER:  1.2746E-01  2.7081E-01  1.2317E-01  4.3435E-02  1.1928E-01  1.0561E-01  1.4265E-02  1.1305E-01  1.0022E-01 -1.7492E-01
             1.1589E-01
 GRADIENT:   1.0813E+01  1.6987E+01  8.4553E+00  1.3271E+01 -1.3970E+01  1.3950E+00  1.1330E+00 -5.4923E-01 -2.0944E+00 -1.6773E+00
            -1.3472E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1730.72833884511        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  1.0225E+00  1.5523E+00  6.3741E-01  6.9307E-01  1.0364E+00  1.0077E+00  7.5779E-01  6.3820E-01  1.2520E+00  7.7260E-01
             1.0143E+00
 PARAMETER:  1.2227E-01  5.3974E-01 -3.5034E-01 -2.6662E-01  1.3573E-01  1.0769E-01 -1.7735E-01 -3.4910E-01  3.2476E-01 -1.5800E-01
             1.1423E-01
 GRADIENT:  -5.9085E+00  6.5784E+00 -2.9162E+00  1.4005E+01 -1.0819E+01  8.4446E-01 -1.7652E+00  1.3443E+00  2.3901E+00  4.2323E+00
             1.5088E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1731.78078739064        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  1.0296E+00  1.8184E+00  4.9837E-01  5.2432E-01  1.1201E+00  1.0058E+00  7.0626E-01  4.2573E-01  1.4768E+00  7.7012E-01
             1.0175E+00
 PARAMETER:  1.2919E-01  6.9797E-01 -5.9641E-01 -5.4565E-01  2.1344E-01  1.0573E-01 -2.4777E-01 -7.5396E-01  4.8989E-01 -1.6121E-01
             1.1730E-01
 GRADIENT:   8.7178E+00  3.4132E+01  4.4428E+00  1.0957E+01 -1.2506E+01 -1.3367E-01 -9.2226E-02  3.5176E-01 -3.1530E+00 -1.7930E+00
            -1.4057E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1732.27785733897        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      913
 NPARAMETR:  1.0262E+00  1.9114E+00  4.2764E-01  4.3254E-01  1.1746E+00  1.0080E+00  6.7483E-01  9.5366E-02  1.6920E+00  8.2240E-01
             1.0238E+00
 PARAMETER:  1.2588E-01  7.4783E-01 -7.4947E-01 -7.3807E-01  2.6089E-01  1.0795E-01 -2.9330E-01 -2.2500E+00  6.2590E-01 -9.5529E-02
             1.2357E-01
 GRADIENT:   2.0112E+00 -4.2998E+01 -2.6801E+00 -6.2036E+00  4.9849E+00  8.7383E-01  7.7147E-01  3.0224E-02  1.6361E+00  2.4503E+00
             2.7196E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1732.41627192875        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  1.0266E+00  1.9222E+00  4.3043E-01  4.3554E-01  1.1741E+00  1.0066E+00  6.7317E-01  3.2576E-02  1.6822E+00  8.0483E-01
             1.0187E+00
 PARAMETER:  1.2629E-01  7.5347E-01 -7.4298E-01 -7.3117E-01  2.6050E-01  1.0658E-01 -2.9576E-01 -3.3242E+00  6.2008E-01 -1.1712E-01
             1.1850E-01
 GRADIENT:   2.7541E+00 -1.5452E+01  3.2849E-01 -1.4218E+00  1.8539E+00  2.7491E-01 -4.0470E-01  2.9940E-03 -1.9237E-01 -9.5424E-01
            -1.0342E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1732.42262200269        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1257             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0268E+00  1.9213E+00  4.2979E-01  4.3593E-01  1.1735E+00  1.0068E+00  6.7376E-01  1.0000E-02  1.6812E+00  8.0679E-01
             1.0196E+00
 PARAMETER:  1.2642E-01  7.5301E-01 -7.4445E-01 -7.3027E-01  2.5999E-01  1.0678E-01 -2.9487E-01 -5.2309E+00  6.1952E-01 -1.1469E-01
             1.1942E-01
 GRADIENT:   5.6922E+02  9.7318E+02  4.3230E+00  1.2000E+02  1.7258E+01  7.4952E+01  2.0329E+01  0.0000E+00  2.4397E+01 -3.6062E-02
             5.1221E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1732.42497020432        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1442             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0268E+00  1.9202E+00  4.2991E-01  4.3691E-01  1.1727E+00  1.0068E+00  6.7404E-01  1.0000E-02  1.6794E+00  8.0570E-01
             1.0196E+00
 PARAMETER:  1.2641E-01  7.5242E-01 -7.4418E-01 -7.2802E-01  2.5929E-01  1.0678E-01 -2.9447E-01 -5.2309E+00  6.1842E-01 -1.1605E-01
             1.1945E-01
 GRADIENT:   5.6905E+02  9.7136E+02  4.2201E+00  1.2032E+02  1.7466E+01  7.4930E+01  2.0272E+01  0.0000E+00  2.4423E+01 -7.9465E-02
             5.3938E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1732.42643232100        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1631             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0268E+00  1.9188E+00  4.2997E-01  4.3768E-01  1.1721E+00  1.0068E+00  6.7430E-01  1.0000E-02  1.6774E+00  8.0483E-01
             1.0196E+00
 PARAMETER:  1.2641E-01  7.5173E-01 -7.4405E-01 -7.2627E-01  2.5883E-01  1.0678E-01 -2.9408E-01 -5.2309E+00  6.1727E-01 -1.1713E-01
             1.1946E-01
 GRADIENT:   5.6900E+02  9.6868E+02  4.0367E+00  1.2048E+02  1.8056E+01  7.4919E+01  2.0213E+01  0.0000E+00  2.4411E+01 -1.1902E-01
             5.6799E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1732.42792961525        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1815
 NPARAMETR:  1.0266E+00  1.9196E+00  4.3041E-01  4.3786E-01  1.1709E+00  1.0067E+00  6.7510E-01  1.0000E-02  1.6774E+00  8.1053E-01
             1.0205E+00
 PARAMETER:  1.2627E-01  7.5211E-01 -7.4301E-01 -7.2586E-01  2.5775E-01  1.0670E-01 -2.9289E-01 -5.2309E+00  6.1723E-01 -1.1007E-01
             1.2026E-01
 GRADIENT:   2.6625E+00 -1.4174E+01  6.7041E-02 -1.0040E+00 -9.0305E-01  3.2571E-01  2.4774E-01  0.0000E+00  2.6013E-01  4.2092E-01
             2.3609E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1732.43087576538        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2008
 NPARAMETR:  1.0268E+00  1.9180E+00  4.3067E-01  4.3836E-01  1.1706E+00  1.0068E+00  6.7507E-01  1.0000E-02  1.6759E+00  8.0818E-01
             1.0202E+00
 PARAMETER:  1.2641E-01  7.5130E-01 -7.4242E-01 -7.2470E-01  2.5753E-01  1.0679E-01 -2.9294E-01 -5.2309E+00  6.1636E-01 -1.1297E-01
             1.1999E-01
 GRADIENT:   2.9836E+00 -1.5533E+01 -5.1256E-02 -1.1766E+00  1.3874E-01  3.6054E-01  1.0291E-01  0.0000E+00  2.3111E-01  1.0945E-01
             6.0716E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1732.43449365197        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2202             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0268E+00  1.9153E+00  4.3122E-01  4.4031E-01  1.1698E+00  1.0068E+00  6.7544E-01  1.0000E-02  1.6721E+00  8.0557E-01
             1.0199E+00
 PARAMETER:  1.2640E-01  7.4986E-01 -7.4113E-01 -7.2028E-01  2.5687E-01  1.0678E-01 -2.9239E-01 -5.2309E+00  6.1407E-01 -1.1620E-01
             1.1973E-01
 GRADIENT:   5.6850E+02  9.6255E+02  3.9291E+00  1.2083E+02  1.7195E+01  7.4854E+01  2.0220E+01  0.0000E+00  2.4491E+01  2.7084E-01
             8.8572E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1732.43555997173        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2393
 NPARAMETR:  1.0268E+00  1.9139E+00  4.3147E-01  4.4123E-01  1.1694E+00  1.0068E+00  6.7574E-01  1.0000E-02  1.6705E+00  8.0520E-01
             1.0199E+00
 PARAMETER:  1.2640E-01  7.4914E-01 -7.4056E-01 -7.1820E-01  2.5653E-01  1.0678E-01 -2.9195E-01 -5.2309E+00  6.1312E-01 -1.1666E-01
             1.1972E-01
 GRADIENT:   2.9373E+00 -1.6694E+01 -7.4910E-01 -3.1328E-01  2.6327E+00  3.4767E-01 -1.3552E-02  0.0000E+00  4.3251E-01 -1.8038E-01
            -3.1515E-02

0ITERATION NO.:   68    OBJECTIVE VALUE:  -1732.43772812524        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     2498
 NPARAMETR:  1.0267E+00  1.9092E+00  4.3117E-01  4.4426E-01  1.1691E+00  1.0068E+00  6.7676E-01  1.0000E-02  1.6666E+00  8.0549E-01
             1.0201E+00
 PARAMETER:  1.2641E-01  7.4917E-01 -7.3795E-01 -7.1852E-01  2.5563E-01  1.0679E-01 -2.9165E-01 -5.2309E+00  6.1223E-01 -1.1541E-01
             1.1976E-01
 GRADIENT:   7.6146E-03  1.5214E+00  9.9667E-02 -6.0781E-01 -2.0102E-01  1.7472E-03 -4.2655E-02  0.0000E+00  4.5467E-02  1.9001E-02
            -9.6580E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2498
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.4014E-04 -3.3793E-02 -2.2148E-04  3.2168E-02 -3.9322E-02
 SE:             2.9871E-02  2.4295E-02  8.4800E-05  2.2809E-02  2.0960E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9091E-01  1.6425E-01  9.0054E-03  1.5844E-01  6.0651E-02

 ETASHRINKSD(%)  1.0000E-10  1.8608E+01  9.9716E+01  2.3589E+01  2.9780E+01
 ETASHRINKVR(%)  1.0000E-10  3.3754E+01  9.9999E+01  4.1613E+01  5.0692E+01
 EBVSHRINKSD(%)  4.3022E-01  1.8420E+01  9.9737E+01  2.5046E+01  2.8967E+01
 EBVSHRINKVR(%)  8.5858E-01  3.3448E+01  9.9999E+01  4.3819E+01  4.9543E+01
 RELATIVEINF(%)  9.9115E+01  5.6943E+00  8.0233E-05  4.8293E+00  1.2036E+01
 EPSSHRINKSD(%)  4.3821E+01
 EPSSHRINKVR(%)  6.8439E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1732.4377281252353     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -997.28690156149707     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.41
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1732.438       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.91E+00  4.33E-01  4.41E-01  1.17E+00  1.01E+00  6.76E-01  1.00E-02  1.67E+00  8.06E-01  1.02E+00
 


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
+        1.03E+03
 
 TH 2
+       -6.93E+00  4.63E+02
 
 TH 3
+        5.55E+00  1.94E+02  4.42E+02
 
 TH 4
+       -1.62E+01  3.40E+02 -3.79E+02  1.19E+03
 
 TH 5
+       -5.32E+00 -2.29E+02 -3.81E+02  3.23E+02  6.39E+02
 
 TH 6
+       -3.46E-01 -1.05E+00  1.96E+00 -3.54E+00 -5.66E-01  1.93E+02
 
 TH 7
+        4.75E-01  6.98E+00 -9.17E+00 -2.08E+01 -1.25E+01  7.20E-02  2.11E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.09E+00 -1.78E+01 -3.12E+01  6.57E+01 -2.89E+00 -9.70E-02  1.30E+01  0.00E+00  3.07E+01
 
 TH10
+        3.14E-01 -1.80E+01 -3.92E+01 -5.25E+00 -7.66E+01  1.04E-01  2.45E+01  0.00E+00  7.05E+00  9.00E+01
 
 TH11
+       -6.54E+00 -2.26E+01 -2.72E+01  4.00E+00 -1.43E+01  2.10E+00  1.31E+01  0.00E+00  4.92E+00  2.39E+01  2.04E+02
 
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
 #CPUT: Total CPU Time in Seconds,       42.618
Stop Time:
Wed Sep 29 10:54:48 CDT 2021
