Sat Sep 25 09:40:02 CDT 2021
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
$DATA ../../../../data/spa/S1/dat4.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1609.54356451505        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.3992E+01 -2.9215E+01 -1.6041E+01  1.7834E+00  8.3940E+01 -2.6725E+01 -3.0578E+00 -9.2908E+00  1.8396E+00 -2.9110E+01
            -2.2560E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1615.60060501627        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  9.7884E-01  1.0079E+00  9.6171E-01  1.0071E+00  9.3280E-01  1.0661E+00  1.0046E+00  1.0361E+00  9.9341E-01  1.0650E+00
             1.0491E+00
 PARAMETER:  7.8613E-02  1.0782E-01  6.0959E-02  1.0711E-01  3.0438E-02  1.6397E-01  1.0456E-01  1.3545E-01  9.3386E-02  1.6293E-01
             1.4793E-01
 GRADIENT:   1.4045E+00  1.1247E+00  2.2303E+00  4.0551E+00  2.6365E+00 -4.1556E+00 -2.2515E-01 -3.0063E+00  1.8524E+00 -2.1042E+00
             2.8415E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1615.74307506330        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      315
 NPARAMETR:  9.7750E-01  1.0460E+00  9.3234E-01  9.8184E-01  9.3380E-01  1.0740E+00  1.0154E+00  1.1062E+00  9.8957E-01  1.0480E+00
             1.0254E+00
 PARAMETER:  7.7242E-02  1.4498E-01  2.9947E-02  8.1673E-02  3.1506E-02  1.7143E-01  1.1527E-01  2.0092E-01  8.9513E-02  1.4685E-01
             1.2513E-01
 GRADIENT:  -1.6808E+00  1.5721E+00 -1.8623E-01  3.6644E+00  1.9776E+00 -1.5459E+00  6.3653E-01 -3.1515E-01 -3.4209E-01 -1.8287E+00
            -5.9258E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1615.90870686363        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      492
 NPARAMETR:  9.7936E-01  1.1531E+00  9.0469E-01  9.1441E-01  9.6745E-01  1.0792E+00  9.1939E-01  1.1466E+00  1.0550E+00  1.0807E+00
             1.0435E+00
 PARAMETER:  7.9144E-02  2.4246E-01 -1.6830E-04  1.0525E-02  6.6912E-02  1.7618E-01  1.5959E-02  2.3680E-01  1.5355E-01  1.7760E-01
             1.4257E-01
 GRADIENT:   8.1318E-01  3.7300E+00  1.7113E+00  3.2429E+00 -2.9488E+00  1.7417E-02 -7.1706E-02 -3.8890E-01 -1.3125E-01  1.8751E-01
             8.3559E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1616.10342008999        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  9.8166E-01  1.4197E+00  6.4720E-01  7.3586E-01  9.7577E-01  1.0809E+00  8.3321E-01  9.4635E-01  1.1931E+00  1.0498E+00
             1.0464E+00
 PARAMETER:  8.1487E-02  4.5045E-01 -3.3510E-01 -2.0672E-01  7.5476E-02  1.7776E-01 -8.2465E-02  4.4854E-02  2.7652E-01  1.4864E-01
             1.4534E-01
 GRADIENT:   1.7465E+00  9.7731E+00  1.4181E+00  5.7772E+00 -3.7542E+00 -4.0884E-01 -5.6840E-01 -2.5351E-01 -5.8241E-01  3.0366E-01
             1.3811E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1616.19633851740        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      848
 NPARAMETR:  9.8087E-01  1.5337E+00  5.3509E-01  6.5209E-01  9.8324E-01  1.0828E+00  8.0696E-01  7.9726E-01  1.2718E+00  1.0343E+00
             1.0425E+00
 PARAMETER:  8.0689E-02  5.2768E-01 -5.2532E-01 -3.2757E-01  8.3097E-02  1.7952E-01 -1.1448E-01 -1.2658E-01  3.4045E-01  1.3376E-01
             1.4160E-01
 GRADIENT:  -3.3149E-01  2.4310E+00  1.5876E-01  1.3297E+00 -4.0848E-01  1.9784E-02  1.4739E-01  1.3229E-03 -2.5206E-01 -2.0297E-01
            -3.2876E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1616.20877309364        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1023
 NPARAMETR:  9.8104E-01  1.5954E+00  4.8938E-01  6.0855E-01  9.9590E-01  1.0829E+00  7.8820E-01  7.2595E-01  1.3261E+00  1.0390E+00
             1.0439E+00
 PARAMETER:  8.0858E-02  5.6710E-01 -6.1461E-01 -3.9668E-01  9.5894E-02  1.7965E-01 -1.3801E-01 -2.2027E-01  3.8224E-01  1.3830E-01
             1.4299E-01
 GRADIENT:  -5.7062E-02 -2.3451E-01 -9.2887E-02 -8.1023E-02  6.5075E-02  1.5999E-02 -9.0391E-03  7.6850E-02  1.3498E-02  2.8945E-02
            -3.4395E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1616.21369225034        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1200
 NPARAMETR:  9.8120E-01  1.6274E+00  4.4718E-01  5.8487E-01  9.8966E-01  1.0829E+00  7.8245E-01  5.5680E-01  1.3480E+00  1.0292E+00
             1.0447E+00
 PARAMETER:  8.1021E-02  5.8697E-01 -7.0479E-01 -4.3637E-01  8.9602E-02  1.7962E-01 -1.4532E-01 -4.8554E-01  3.9863E-01  1.2876E-01
             1.4371E-01
 GRADIENT:   1.7075E-01  9.1535E-02 -8.5275E-02  1.0430E-01 -1.1121E-01 -1.9784E-02 -3.7593E-02  2.3984E-02 -3.2969E-03  7.2779E-02
             1.6053E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1616.21548980134        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1376
 NPARAMETR:  9.8107E-01  1.6644E+00  4.0959E-01  5.5781E-01  9.9030E-01  1.0829E+00  7.7527E-01  3.4057E-01  1.3789E+00  1.0245E+00
             1.0449E+00
 PARAMETER:  8.0887E-02  6.0944E-01 -7.9260E-01 -4.8374E-01  9.0248E-02  1.7967E-01 -1.5455E-01 -9.7713E-01  4.2129E-01  1.2425E-01
             1.4389E-01
 GRADIENT:   8.1528E-03 -7.6817E-01 -2.3888E-01 -1.0032E-01  4.4570E-01  5.4049E-03  5.7154E-02  1.5815E-02  1.1967E-02  2.2311E-02
             7.3136E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1616.21894932654        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1554
 NPARAMETR:  9.8107E-01  1.6744E+00  3.9643E-01  5.5047E-01  9.8763E-01  1.0829E+00  7.7377E-01  1.6102E-01  1.3861E+00  1.0211E+00
             1.0449E+00
 PARAMETER:  8.0890E-02  6.1546E-01 -8.2525E-01 -4.9699E-01  8.7550E-02  1.7962E-01 -1.5648E-01 -1.7262E+00  4.2651E-01  1.2089E-01
             1.4388E-01
 GRADIENT:   5.7318E-02  2.1764E-01 -2.7617E-02  1.0452E-01 -3.6322E-02 -6.5981E-03  3.7094E-02  3.2012E-03 -1.3909E-02  2.9540E-03
             1.7356E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1616.22031249129        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1729
 NPARAMETR:  9.8104E-01  1.6795E+00  3.9228E-01  5.4665E-01  9.8859E-01  1.0829E+00  7.7248E-01  3.6671E-02  1.3914E+00  1.0215E+00
             1.0449E+00
 PARAMETER:  8.0856E-02  6.1853E-01 -8.3577E-01 -5.0394E-01  8.8529E-02  1.7962E-01 -1.5815E-01 -3.2058E+00  4.3029E-01  1.2123E-01
             1.4394E-01
 GRADIENT:   1.9256E-02 -1.8521E-01 -2.4694E-02 -5.6211E-02  9.7732E-02 -1.4709E-03  3.0116E-04  1.8706E-04 -1.2575E-02  1.2906E-04
             1.2197E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1616.22043482014        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1905
 NPARAMETR:  9.8103E-01  1.6786E+00  3.9248E-01  5.4730E-01  9.8794E-01  1.0829E+00  7.7276E-01  1.1557E-02  1.3904E+00  1.0210E+00
             1.0449E+00
 PARAMETER:  8.0845E-02  6.1798E-01 -8.3528E-01 -5.0275E-01  8.7869E-02  1.7963E-01 -1.5778E-01 -4.3605E+00  4.2957E-01  1.2076E-01
             1.4389E-01
 GRADIENT:  -3.7450E-03  1.9174E-02  9.1460E-03  2.0517E-03 -1.1902E-02  1.2374E-03  4.5270E-04  1.8294E-05  4.5476E-04 -2.8532E-03
            -3.2637E-03

0ITERATION NO.:   58    OBJECTIVE VALUE:  -1616.22043531758        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:     2002
 NPARAMETR:  9.8103E-01  1.6787E+00  3.9244E-01  5.4727E-01  9.8795E-01  1.0829E+00  7.7275E-01  1.0526E-02  1.3904E+00  1.0210E+00
             1.0449E+00
 PARAMETER:  8.0846E-02  6.1801E-01 -8.3537E-01 -5.0281E-01  8.7879E-02  1.7962E-01 -1.5780E-01 -4.4540E+00  4.2960E-01  1.2077E-01
             1.4389E-01
 GRADIENT:  -1.7211E-02  7.9968E-02  6.3652E-03  9.4502E-03 -1.3845E-02 -5.4684E-03  3.0838E-04  1.0397E-05  1.7526E-05 -1.4222E-03
            -2.6675E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2002
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8690E-04 -3.2466E-02 -2.9085E-04  2.7105E-02 -3.5794E-02
 SE:             2.9855E-02  2.3600E-02  1.0966E-04  2.3408E-02  2.2596E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9501E-01  1.6892E-01  7.9966E-03  2.4690E-01  1.1318E-01

 ETASHRINKSD(%)  1.0000E-10  2.0937E+01  9.9633E+01  2.1580E+01  2.4300E+01
 ETASHRINKVR(%)  1.0000E-10  3.7490E+01  9.9999E+01  3.8503E+01  4.2696E+01
 EBVSHRINKSD(%)  4.0675E-01  2.0671E+01  9.9696E+01  2.2887E+01  2.2603E+01
 EBVSHRINKVR(%)  8.1184E-01  3.7070E+01  9.9999E+01  4.0536E+01  4.0097E+01
 RELATIVEINF(%)  9.9135E+01  4.3006E+00  9.2088E-05  3.8243E+00  1.0591E+01
 EPSSHRINKSD(%)  4.5253E+01
 EPSSHRINKVR(%)  7.0028E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1616.2204353175753     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -881.06960875383709     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.57
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.27
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1616.220       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.68E+00  3.92E-01  5.47E-01  9.88E-01  1.08E+00  7.73E-01  1.05E-02  1.39E+00  1.02E+00  1.04E+00
 


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
+        9.78E+02
 
 TH 2
+       -6.41E+00  4.36E+02
 
 TH 3
+        1.32E+00  2.09E+02  6.87E+02
 
 TH 4
+       -1.25E+01  3.44E+02 -4.77E+02  1.15E+03
 
 TH 5
+       -3.29E+00 -2.42E+02 -5.01E+02  3.41E+02  6.56E+02
 
 TH 6
+       -4.35E-01 -9.78E-01  1.00E+00 -3.28E+00 -3.08E+00  1.68E+02
 
 TH 7
+       -1.37E+00  3.21E+00 -3.43E+01 -4.73E+00 -1.27E+01 -1.50E-01  1.44E+02
 
 TH 8
+        6.98E-02  7.81E-01 -2.36E+00  2.07E+00  2.43E+00 -1.24E+00  2.49E+00  2.55E+00
 
 TH 9
+        2.25E+00 -1.89E+01 -4.14E+01  6.04E+01 -2.82E+00 -2.96E-01  2.03E+01 -1.72E+00  4.52E+01
 
 TH10
+       -1.31E+00 -1.77E+01 -4.76E+01 -2.99E+00 -6.05E+01 -9.92E-01  1.37E+01 -6.91E+00  6.50E+00  7.95E+01
 
 TH11
+       -5.89E+00 -1.62E+01 -2.54E+01  1.35E+00 -2.93E+00  1.08E+00  1.23E+01  3.54E+00  4.57E+00  1.54E+01  1.90E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.904
Stop Time:
Sat Sep 25 09:40:35 CDT 2021
