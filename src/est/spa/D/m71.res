Sat Sep 18 15:35:00 CDT 2021
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
$DATA ../../../../data/spa/D/dat71.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13399.3104542946        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2233E-01  2.3493E+02 -8.2055E+01 -1.8122E+02  3.2860E+02 -2.5638E+03 -7.7968E+02 -8.0002E+01 -1.8826E+03 -7.3699E+02
            -2.4074E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -545.204705301075        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4185E+00  9.0552E-01  8.0249E-01  2.3960E+00  1.9722E+00  3.2795E+00  1.5259E+00  9.3804E-01  2.1597E+00  1.1853E+00
             1.2052E+01
 PARAMETER:  4.4961E-01  7.5393E-04 -1.2003E-01  9.7381E-01  7.7913E-01  1.2877E+00  5.2255E-01  3.6034E-02  8.6995E-01  2.7001E-01
             2.5892E+00
 GRADIENT:   1.0102E+01  4.2609E+01 -3.9016E+01  7.3063E+01 -1.9175E+01  5.4145E+01 -2.5205E-01  8.2362E+00 -2.4722E+01  3.3865E+00
             1.8401E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -578.687746558809        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.3998E+00  1.0081E+00  1.0769E+00  2.3075E+00  4.3010E+00  3.1880E+00  2.6124E+00  3.4414E-01  2.9018E+00  1.1643E+00
             1.1412E+01
 PARAMETER:  4.3632E-01  1.0810E-01  1.7407E-01  9.3617E-01  1.5589E+00  1.2594E+00  1.0603E+00 -9.6672E-01  1.1653E+00  2.5215E-01
             2.5347E+00
 GRADIENT:   9.9045E+00  3.7315E+01 -1.0792E+01  4.9365E+01 -1.4954E+01  3.5328E+01  7.1996E+00  2.9404E-01  8.4258E+00  9.1289E-01
             6.5838E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -607.914112968552        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  1.2499E+00  6.4931E-01  7.8661E-01  1.7623E+00  5.2370E+00  2.8066E+00  4.9722E-01  1.7785E-02  2.4197E+00  9.0427E+00
             1.1446E+01
 PARAMETER:  3.2303E-01 -3.3185E-01 -1.4002E-01  6.6662E-01  1.7557E+00  1.1320E+00 -5.9872E-01 -3.9294E+00  9.8366E-01  2.3020E+00
             2.5376E+00
 GRADIENT:  -3.1481E+00  9.8467E+00  1.0639E+01 -7.2424E+00 -1.3738E+01  1.5789E+01  9.1063E-01  2.7252E-04  1.3945E+01  6.2753E+00
             1.0250E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -661.631192968483        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  9.5706E-01  1.2391E-01  1.2609E-01  1.2366E+00  9.2805E+00  3.0146E+00  7.2160E-01  1.0000E-02  7.0939E-01  4.2019E+00
             1.1488E+01
 PARAMETER:  5.6114E-02 -1.9882E+00 -1.9707E+00  3.1234E-01  2.3279E+00  1.2035E+00 -2.2628E-01 -1.1876E+01 -2.4335E-01  1.5355E+00
             2.5413E+00
 GRADIENT:   2.6087E+01  3.8165E+01 -1.0191E+02  1.4105E+02 -3.7324E+01  6.0974E+01  8.5192E+00  0.0000E+00  5.4313E+00  3.4130E+01
             1.6151E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -746.693851649950        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  4.0125E-01  1.0000E-02  1.6875E-02  2.3833E-01  6.6383E+00  1.7869E+00  8.4506E-01  1.0000E-02  1.4272E-01  1.0448E+00
             1.0378E+01
 PARAMETER: -8.1318E-01 -4.7436E+00 -3.9819E+00 -1.3341E+00  1.9929E+00  6.8050E-01 -6.8352E-02 -1.7607E+01 -1.8468E+00  1.4381E-01
             2.4397E+00
 GRADIENT:  -9.4626E-01  0.0000E+00 -9.9517E-01  1.3768E+01 -1.7170E+01 -3.3318E+01  5.6859E-01  0.0000E+00  3.8373E-01  2.8366E+00
             4.0216E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -749.713048977365        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      494
 NPARAMETR:  3.7441E-01  1.0000E-02  1.4352E-02  2.0787E-01  7.1251E+00  1.9296E+00  7.3153E-01  1.0000E-02  1.3590E-01  9.9159E-01
             9.8965E+00
 PARAMETER: -8.8240E-01 -4.9399E+00 -4.1438E+00 -1.4709E+00  2.0636E+00  7.5733E-01 -2.1262E-01 -1.8402E+01 -1.8959E+00  9.1558E-02
             2.3922E+00
 GRADIENT:  -1.0230E+00  0.0000E+00 -6.9824E+00  4.1930E+00 -1.8630E+00 -4.7015E+00  8.5953E-01  0.0000E+00  1.1225E-01 -5.2606E-01
             4.1899E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -750.143765490353        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      669
 NPARAMETR:  3.8705E-01  1.0000E-02  1.5582E-02  2.2102E-01  7.2020E+00  1.9652E+00  5.4946E-01  1.0000E-02  1.5295E-01  1.1476E+00
             9.7941E+00
 PARAMETER: -8.4920E-01 -4.8038E+00 -4.0617E+00 -1.4095E+00  2.0744E+00  7.7560E-01 -4.9882E-01 -1.8108E+01 -1.7776E+00  2.3771E-01
             2.3818E+00
 GRADIENT:   2.9928E-01  0.0000E+00 -7.1976E-01 -1.6142E-02 -9.3030E-01 -1.0186E-01  2.9036E-01  0.0000E+00 -3.5112E-02  4.0151E-01
            -2.2673E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -750.253357909892        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      846
 NPARAMETR:  3.8597E-01  1.0000E-02  1.5445E-02  2.1923E-01  7.2510E+00  1.9733E+00  2.3653E-01  1.0000E-02  1.7184E-01  1.1669E+00
             9.8264E+00
 PARAMETER: -8.5200E-01 -4.7811E+00 -4.0704E+00 -1.4176E+00  2.0811E+00  7.7969E-01 -1.3417E+00 -1.7981E+01 -1.6612E+00  2.5439E-01
             2.3851E+00
 GRADIENT:   1.5407E-01  0.0000E+00  3.3286E-01 -9.2751E-01 -4.5971E-01  1.5023E+00  5.0675E-02  0.0000E+00 -3.7463E-02  5.0237E-01
             1.5974E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -752.414806528857        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1025
 NPARAMETR:  3.7516E-01  1.0226E-02  1.4493E-02  2.0404E-01  7.7410E+00  1.9474E+00  1.0000E-02  1.0000E-02  7.4421E-01  1.2521E+00
             9.5529E+00
 PARAMETER: -8.8039E-01 -4.4828E+00 -4.1341E+00 -1.4895E+00  2.1465E+00  7.6647E-01 -1.0996E+01 -1.5669E+01 -1.9543E-01  3.2482E-01
             2.3568E+00
 GRADIENT:  -8.9588E+00  8.4772E+00 -1.7092E+01  1.5933E+01 -1.1528E+01 -5.3593E+00  0.0000E+00  0.0000E+00  3.9614E+00 -1.3205E+00
             3.7897E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -755.636025233450        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1201
 NPARAMETR:  3.7304E-01  1.0642E-02  1.4368E-02  1.9860E-01  9.9571E+00  1.9809E+00  1.0000E-02  1.0000E-02  6.9467E-01  2.2649E+00
             9.0026E+00
 PARAMETER: -8.8607E-01 -4.4429E+00 -4.1428E+00 -1.5165E+00  2.3983E+00  7.8353E-01 -1.0330E+01 -1.6711E+01 -2.6432E-01  9.1754E-01
             2.2975E+00
 GRADIENT:  -1.3886E+00  7.5741E+00 -9.4549E-01 -5.5967E+00 -1.1916E+00  3.8326E-02  0.0000E+00  0.0000E+00 -2.2254E+00 -1.5845E+00
            -4.2644E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -756.093341458730        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1381
 NPARAMETR:  3.7615E-01  1.0568E-02  1.4553E-02  2.0121E-01  1.0856E+01  1.9849E+00  1.0000E-02  1.0000E-02  6.9428E-01  2.8077E+00
             9.0594E+00
 PARAMETER: -8.7776E-01 -4.4499E+00 -4.1299E+00 -1.5034E+00  2.4847E+00  7.8556E-01 -1.0398E+01 -1.7005E+01 -2.6488E-01  1.1324E+00
             2.3038E+00
 GRADIENT:   1.6101E-02  5.6373E+00 -4.9784E+00  5.0416E-01 -2.8635E-01  3.1769E-01  0.0000E+00  0.0000E+00 -1.7014E+00 -8.9516E-01
            -1.8607E-01

0ITERATION NO.:   58    OBJECTIVE VALUE:  -756.094006603349        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:     1506
 NPARAMETR:  3.7615E-01  1.0581E-02  1.4538E-02  2.0129E-01  1.0848E+01  1.9833E+00  1.0000E-02  1.0000E-02  6.9458E-01  2.8086E+00
             9.0537E+00
 PARAMETER: -8.7776E-01 -4.4499E+00 -4.1299E+00 -1.5034E+00  2.4847E+00  7.8555E-01 -1.0398E+01 -1.7005E+01 -2.6438E-01  1.1324E+00
             2.3038E+00
 GRADIENT:   3.4635E-03 -1.3430E+02  1.3991E+02 -4.0701E+02  2.4915E+02  2.4338E-01  0.0000E+00  0.0000E+00  2.3506E+03 -2.7573E+02
             2.7069E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1506
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4449E-03  2.5344E-06  1.7147E-04 -1.9731E-02 -4.6321E-03
 SE:             2.9298E-02  1.4444E-05  2.5518E-04  2.1032E-02  4.2367E-03
 N:                     100         100         100         100         100

 P VAL.:         9.6067E-01  8.6071E-01  5.0163E-01  3.4818E-01  2.7425E-01

 ETASHRINKSD(%)  1.8497E+00  9.9952E+01  9.9145E+01  2.9540E+01  8.5806E+01
 ETASHRINKVR(%)  3.6651E+00  1.0000E+02  9.9993E+01  5.0353E+01  9.7985E+01
 EBVSHRINKSD(%)  1.7414E+00  9.9911E+01  9.9176E+01  3.1788E+01  8.8446E+01
 EBVSHRINKVR(%)  3.4525E+00  1.0000E+02  9.9993E+01  5.3471E+01  9.8665E+01
 RELATIVEINF(%)  8.2701E+00  1.9635E-05  4.0629E-05  2.7358E-01  8.8395E-01
 EPSSHRINKSD(%)  1.3249E+01
 EPSSHRINKVR(%)  2.4742E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -756.09400660334916     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -20.943180039610979     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -756.094       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.76E-01  1.06E-02  1.46E-02  2.01E-01  1.09E+01  1.98E+00  1.00E-02  1.00E-02  6.95E-01  2.81E+00  9.06E+00
 


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
+        1.43E+06
 
 TH 2
+       -1.22E+04  7.05E+07
 
 TH 3
+       -9.46E+03  1.03E+05  4.29E+07
 
 TH 4
+       -1.85E+03 -3.86E+05  3.44E+04  1.67E+06
 
 TH 5
+        1.96E+01  1.60E+02  9.18E+04  6.57E+02  2.12E+02
 
 TH 6
+        3.03E+05  6.83E+01  1.80E+02 -5.48E+01  1.47E-01  4.28E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.81E+03  8.73E+04 -6.81E+04  1.36E+04 -1.52E+02  2.07E+01  0.00E+00  0.00E+00  4.60E+06
 
 TH10
+       -1.63E+02 -3.10E+03  2.40E+03 -5.01E+02  5.39E+00 -1.24E+00  0.00E+00  0.00E+00  1.29E+03  1.54E+04
 
 TH11
+        9.68E+00 -3.25E+02  8.10E+02 -8.42E+01  5.67E-01  7.88E-01  0.00E+00  0.00E+00 -1.89E+02 -4.51E+00  3.63E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.695
Stop Time:
Sat Sep 18 15:35:34 CDT 2021
