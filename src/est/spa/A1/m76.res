Sat Sep 25 08:14:31 CDT 2021
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
$DATA ../../../../data/spa/A1/dat76.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m76.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -956.948498452628        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2541E+02  4.1016E+01  3.9961E+01  4.5843E+00  1.8088E+01 -4.1471E+01 -1.6202E+01 -9.2547E+00 -3.0235E+01 -5.2814E+01
            -1.2419E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1368.14212388498        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0411E+00  1.0035E+00  1.1280E+00  1.0803E+00  1.0579E+00  1.0820E+00  9.2111E-01  9.1946E-01  9.4351E-01  7.8064E-01
             3.1989E+00
 PARAMETER:  1.4031E-01  1.0353E-01  2.2044E-01  1.7726E-01  1.5631E-01  1.7883E-01  1.7824E-02  1.6026E-02  4.1850E-02 -1.4764E-01
             1.2628E+00
 GRADIENT:   9.1461E+01  3.7009E+01 -4.3801E+00  6.2969E+01 -6.4270E+00  1.9494E+00  2.2289E+00  3.0623E+00  6.0923E-01  9.0859E+00
             9.6119E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1381.66851774347        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.9723E-01  9.8089E-01  1.1088E+00  1.0299E+00  1.0260E+00  1.0730E+00  8.0054E-01  2.8605E-01  1.0617E+00  8.0689E-01
             2.8314E+00
 PARAMETER:  9.7224E-02  8.0708E-02  2.0331E-01  1.2944E-01  1.2565E-01  1.7048E-01 -1.2247E-01 -1.1516E+00  1.5989E-01 -1.1457E-01
             1.1408E+00
 GRADIENT:   3.3072E+01  6.3953E+00  4.6721E+00  5.1210E+00 -6.7268E+00  3.2925E+00  9.7638E-01  3.3014E-01  6.1646E+00  5.5384E+00
             3.7516E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1384.70972219540        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.7371E-01  7.5080E-01  8.7632E-01  1.1519E+00  8.1546E-01  1.0657E+00  1.0343E+00  2.0511E-01  9.4189E-01  6.5354E-01
             2.6075E+00
 PARAMETER:  7.3356E-02 -1.8662E-01 -3.2028E-02  2.4145E-01 -1.0401E-01  1.6363E-01  1.3372E-01 -1.4842E+00  4.0133E-02 -3.2536E-01
             1.0584E+00
 GRADIENT:  -5.2837E+00  5.2258E+00 -4.9636E-01  8.5639E+00  2.9330E-01 -6.0067E-01 -8.2483E-01  9.5665E-02 -1.5787E+00 -7.1680E-01
            -5.5782E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1385.40440096503        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  9.7413E-01  5.5388E-01  8.7094E-01  1.2634E+00  7.3865E-01  1.0655E+00  1.2714E+00  5.1371E-02  8.8481E-01  6.6748E-01
             2.6000E+00
 PARAMETER:  7.3786E-02 -4.9081E-01 -3.8177E-02  3.3379E-01 -2.0292E-01  1.6347E-01  3.4012E-01 -2.8687E+00 -2.2384E-02 -3.0424E-01
             1.0555E+00
 GRADIENT:  -5.7190E-01  4.9249E+00  3.8549E+00  7.8701E+00 -7.4107E+00 -2.0260E-01 -1.0496E-01  9.4941E-03 -5.3303E-01 -3.5009E-01
            -1.8753E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1385.89578877783        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.7063E-01  3.3269E-01  8.4287E-01  1.3794E+00  6.6341E-01  1.0637E+00  1.7428E+00  1.0000E-02  8.3753E-01  6.7183E-01
             2.5860E+00
 PARAMETER:  7.0185E-02 -1.0005E+00 -7.0938E-02  4.2163E-01 -3.1036E-01  1.6176E-01  6.5548E-01 -6.1542E+00 -7.7297E-02 -2.9775E-01
             1.0501E+00
 GRADIENT:  -1.5116E-01  1.0413E+00  8.1607E-01  5.3798E+00 -2.3771E+00 -2.4893E-01 -1.2287E-01  0.0000E+00  1.1424E+00 -4.7869E-01
            -1.2709E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1385.96082423331        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  9.6970E-01  2.9853E-01  8.5027E-01  1.3962E+00  6.5961E-01  1.0639E+00  1.8860E+00  1.0000E-02  8.2368E-01  6.8846E-01
             2.5784E+00
 PARAMETER:  6.9228E-02 -1.1089E+00 -6.2196E-02  4.3376E-01 -3.1611E-01  1.6192E-01  7.3445E-01 -6.9706E+00 -9.3974E-02 -2.7330E-01
             1.0472E+00
 GRADIENT:  -1.2976E-01  3.4814E-01  4.2829E-01  1.0115E+00 -9.0234E-01  1.2965E-02  4.6899E-02  0.0000E+00 -2.4675E-01  6.6831E-02
            -2.1054E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1386.07636877280        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      560
 NPARAMETR:  9.7140E-01  2.3449E-01  8.6297E-01  1.4395E+00  6.5089E-01  1.0671E+00  2.1756E+00  1.0000E-02  8.0972E-01  6.9884E-01
             2.5789E+00
 PARAMETER:  7.0979E-02 -1.3503E+00 -4.7374E-02  4.6429E-01 -3.2942E-01  1.6491E-01  8.7731E-01 -8.7398E+00 -1.1106E-01 -2.5833E-01
             1.0474E+00
 GRADIENT:  -5.8054E-01  7.9257E-01  1.3511E+00  3.2249E+00 -3.5561E+00  2.6137E-01 -1.9876E-02  0.0000E+00 -3.5727E-01  1.4335E-01
            -4.8908E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1386.13650822619        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      736
 NPARAMETR:  9.6954E-01  1.5566E-01  8.8909E-01  1.4847E+00  6.4804E-01  1.0643E+00  2.6643E+00  1.0000E-02  7.9486E-01  7.0989E-01
             2.5805E+00
 PARAMETER:  6.9071E-02 -1.7601E+00 -1.7559E-02  4.9519E-01 -3.3380E-01  1.6235E-01  1.0799E+00 -1.1986E+01 -1.2959E-01 -2.4265E-01
             1.0480E+00
 GRADIENT:   9.0445E-02  1.0000E-01  6.2005E-02  1.2011E+00 -4.1492E-01 -1.4000E-01 -6.4270E-02  0.0000E+00  3.8984E-02  4.1225E-02
            -1.1766E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1386.14494048930        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      911
 NPARAMETR:  9.6824E-01  1.0657E-01  8.9586E-01  1.5110E+00  6.4099E-01  1.0639E+00  3.2722E+00  1.0000E-02  7.8562E-01  7.1486E-01
             2.5785E+00
 PARAMETER:  6.7722E-02 -2.1390E+00 -9.9698E-03  5.1279E-01 -3.4474E-01  1.6199E-01  1.2855E+00 -1.5162E+01 -1.4128E-01 -2.3567E-01
             1.0472E+00
 GRADIENT:   1.7762E-01 -1.2170E-02 -5.7069E-02 -4.0037E-01  1.4354E-01 -2.2189E-03  6.2935E-03  0.0000E+00 -1.9541E-02 -2.4138E-02
            -1.3940E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1386.14510111642        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1087
 NPARAMETR:  9.6796E-01  9.9393E-02  8.9692E-01  1.5152E+00  6.3994E-01  1.0639E+00  3.3930E+00  1.0000E-02  7.8445E-01  7.1595E-01
             2.5787E+00
 PARAMETER:  6.7432E-02 -2.2087E+00 -8.7865E-03  5.1557E-01 -3.4638E-01  1.6190E-01  1.3217E+00 -1.5755E+01 -1.4278E-01 -2.3414E-01
             1.0473E+00
 GRADIENT:  -2.3814E-02  8.3130E-03  1.1744E-02  4.1977E-02 -2.3165E-02  2.2307E-03  1.0140E-02  0.0000E+00  1.1450E-02  3.7603E-03
             2.1935E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1386.14510350394        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1269
 NPARAMETR:  9.6794E-01  9.8302E-02  8.9704E-01  1.5158E+00  6.3980E-01  1.0638E+00  3.4008E+00  1.0000E-02  7.8424E-01  7.1605E-01
             2.5786E+00
 PARAMETER:  6.7416E-02 -2.2197E+00 -8.6495E-03  5.1592E-01 -3.4660E-01  1.6189E-01  1.3240E+00 -1.5843E+01 -1.4304E-01 -2.3401E-01
             1.0472E+00
 GRADIENT:   1.2286E-02 -6.9671E-03 -3.2191E-02 -9.0211E-02  5.3589E-02  7.5270E-04 -5.4177E-03  0.0000E+00 -1.0645E-02 -3.8970E-03
            -1.2060E-02

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1386.14510747524        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1326
 NPARAMETR:  9.6794E-01  9.8361E-02  8.9707E-01  1.5158E+00  6.3978E-01  1.0638E+00  3.4053E+00  1.0000E-02  7.8426E-01  7.1609E-01
             2.5786E+00
 PARAMETER:  6.7413E-02 -2.2191E+00 -8.6226E-03  5.1594E-01 -3.4662E-01  1.6188E-01  1.3253E+00 -1.5843E+01 -1.4301E-01 -2.3396E-01
             1.0472E+00
 GRADIENT:  -4.7216E-03  3.7658E-03  2.5542E-02  2.2834E-02 -4.2168E-02 -2.7390E-04  2.0528E-03  0.0000E+00  6.2913E-04  8.7339E-04
             3.4263E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1326
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -7.5425E-05  8.4198E-04  8.4199E-05 -1.2258E-02 -1.7369E-02
 SE:             2.9155E-02  5.4598E-03  1.6163E-04  2.6326E-02  1.7774E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9794E-01  8.7744E-01  6.0242E-01  6.4147E-01  3.2846E-01

 ETASHRINKSD(%)  2.3270E+00  8.1709E+01  9.9459E+01  1.1806E+01  4.0455E+01
 ETASHRINKVR(%)  4.5998E+00  9.6654E+01  9.9997E+01  2.2219E+01  6.4544E+01
 EBVSHRINKSD(%)  2.1855E+00  8.2512E+01  9.9431E+01  1.1369E+01  4.0914E+01
 EBVSHRINKVR(%)  4.3233E+00  9.6942E+01  9.9997E+01  2.1445E+01  6.5089E+01
 RELATIVEINF(%)  8.6561E+01  6.6914E-02  1.7367E-04  2.5875E+00  1.3308E+00
 EPSSHRINKSD(%)  2.9872E+01
 EPSSHRINKVR(%)  5.0820E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1386.1451074752399     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -650.99428091150173     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.84
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1386.145       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  9.84E-02  8.97E-01  1.52E+00  6.40E-01  1.06E+00  3.41E+00  1.00E-02  7.84E-01  7.16E-01  2.58E+00
 


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
+        9.97E+02
 
 TH 2
+       -5.82E+01  3.46E+02
 
 TH 3
+        3.41E+00  1.62E+02  4.28E+02
 
 TH 4
+       -4.22E+01  3.54E+02 -3.53E+01  6.17E+02
 
 TH 5
+        2.46E+01 -4.38E+02 -8.37E+02 -1.28E+02  1.75E+03
 
 TH 6
+        6.82E-01 -7.84E+00  1.04E+01 -1.21E+01 -2.30E+00  1.55E+02
 
 TH 7
+       -1.04E-02  3.24E+00  5.02E-01 -2.28E-01  1.19E-01  7.94E-02  1.95E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.89E+00 -2.53E+01  1.91E+01 -1.33E+01 -7.98E+00  2.60E+00  9.49E-01  0.00E+00  1.99E+02
 
 TH10
+       -6.30E+00  5.77E+00 -2.89E+00 -8.07E+00 -1.70E+01 -6.26E-01  1.30E-01  0.00E+00  5.79E+00  4.30E+01
 
 TH11
+       -1.30E+01 -1.51E+00 -4.99E+00 -8.13E+00 -1.06E+01  2.86E+00  1.56E-01  0.00E+00  1.05E+01  2.44E+01  4.70E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.195
Stop Time:
Sat Sep 25 08:14:54 CDT 2021
