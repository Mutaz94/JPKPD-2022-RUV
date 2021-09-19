Sat Sep 18 10:08:36 CDT 2021
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
$DATA ../../../../data/spa/A2/dat92.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -926.706332254307        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.0926E+01 -3.5385E+01 -1.0902E+00 -5.1211E+01  1.2111E+02  1.9312E+01 -4.8649E+01 -1.9020E+01 -8.1503E+01 -7.3143E+01
            -1.2270E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1365.25055990249        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0011E+00  1.0948E+00  1.1452E+00  1.0514E+00  1.0811E+00  8.3816E-01  9.9509E-01  9.4994E-01  9.4163E-01  8.3914E-01
             3.2035E+00
 PARAMETER:  1.0112E-01  1.9061E-01  2.3562E-01  1.5014E-01  1.7798E-01 -7.6542E-02  9.5075E-02  4.8640E-02  3.9855E-02 -7.5381E-02
             1.2642E+00
 GRADIENT:   4.0117E+00 -1.5393E+00 -1.6608E+01  1.4059E+01  4.5982E+00 -2.8546E+01  2.6488E+00  3.5913E+00  5.2704E+00  1.3103E+01
             5.8212E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1373.53327987703        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.9785E-01  8.8701E-01  1.5930E+00  1.2056E+00  1.1141E+00  9.3566E-01  9.7924E-01  8.6992E-01  9.1456E-01  4.7014E-01
             3.0697E+00
 PARAMETER:  9.7850E-02 -1.9903E-02  5.6564E-01  2.8698E-01  2.0806E-01  3.3499E-02  7.9024E-02 -3.9357E-02  1.0688E-02 -6.5472E-01
             1.2216E+00
 GRADIENT:   9.1515E+00  1.7790E+01  4.1712E-01  3.4969E+01 -1.3772E+01  1.0105E+01  8.6243E-01  1.4606E+00  3.1424E+00  2.2075E+00
             1.8028E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1375.94195610315        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.9210E-01  9.8940E-01  2.0793E+00  1.1111E+00  1.2854E+00  9.0208E-01  5.5757E-01  2.2168E-01  1.0504E+00  2.1138E-01
             3.0549E+00
 PARAMETER:  9.2073E-02  8.9345E-02  8.3202E-01  2.0538E-01  3.5105E-01 -3.0478E-03 -4.8417E-01 -1.4065E+00  1.4916E-01 -1.4541E+00
             1.2167E+00
 GRADIENT:   2.1190E+00 -7.9239E+00 -1.2634E+00 -6.1014E+00  7.8345E+00 -3.3666E-01  9.6568E-01  5.6442E-02  6.0551E+00  2.5452E-01
             3.8836E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1376.24098846938        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.9178E-01  9.5793E-01  2.2114E+00  1.1350E+00  1.2792E+00  9.0361E-01  4.2074E-01  5.7476E-02  1.0245E+00  8.2863E-02
             3.0543E+00
 PARAMETER:  9.1744E-02  5.7021E-02  8.9364E-01  2.2663E-01  3.4624E-01 -1.3619E-03 -7.6573E-01 -2.7564E+00  1.2421E-01 -2.3906E+00
             1.2166E+00
 GRADIENT:   1.4603E+00  1.3711E+00  5.5468E-01  1.8350E+00 -1.3383E+00  5.4703E-01  2.1331E-02  3.3707E-03  4.0367E-01  2.9140E-02
            -2.2349E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1376.32256501078        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  9.9009E-01  8.5811E-01  2.2429E+00  1.1974E+00  1.2596E+00  9.0186E-01  2.8724E-01  1.0000E-02  9.8613E-01  2.1014E-02
             3.0568E+00
 PARAMETER:  9.0043E-02 -5.3017E-02  9.0776E-01  2.8013E-01  3.3083E-01 -3.2965E-03 -1.1474E+00 -4.6559E+00  8.6034E-02 -3.7626E+00
             1.2174E+00
 GRADIENT:  -9.1124E-01 -6.2858E-01  2.9172E-02 -1.2722E+00  5.3520E-01  6.6583E-02  1.0627E-02  0.0000E+00 -2.4659E-01  1.7590E-03
             3.6099E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1376.45445861693        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      556
 NPARAMETR:  9.9188E-01  6.6139E-01  2.3234E+00  1.3367E+00  1.2174E+00  9.0301E-01  1.0344E-01  1.0000E-02  8.9818E-01  1.0000E-02
             3.0582E+00
 PARAMETER:  9.1843E-02 -3.1342E-01  9.4303E-01  3.9022E-01  2.9674E-01 -2.0185E-03 -2.1687E+00 -9.5397E+00 -7.3828E-03 -7.3756E+00
             1.2178E+00
 GRADIENT:   1.3990E+00  3.4652E+00  1.9294E-01  8.5361E+00 -1.7351E+00  4.5279E-02  1.0311E-03  0.0000E+00 -2.0545E-01  0.0000E+00
            -7.5902E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1376.54776237397        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  9.8970E-01  4.8285E-01  2.3573E+00  1.4502E+00  1.1748E+00  9.0189E-01  1.4144E-02  1.0000E-02  8.2391E-01  1.0000E-02
             3.0586E+00
 PARAMETER:  8.9651E-02 -6.2805E-01  9.5751E-01  4.7168E-01  2.6112E-01 -3.2604E-03 -4.1585E+00 -1.8150E+01 -9.3690E-02 -1.4004E+01
             1.2180E+00
 GRADIENT:   7.4894E-02  7.4062E-01 -1.7418E-03  2.4029E+00 -2.6985E-01 -2.4783E-02  5.3861E-05  0.0000E+00 -2.0006E-01  0.0000E+00
            -4.8028E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1376.57474179848        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  9.8855E-01  3.7439E-01  2.3561E+00  1.5198E+00  1.1445E+00  9.0140E-01  1.0000E-02  1.0000E-02  7.8479E-01  1.0000E-02
             3.0579E+00
 PARAMETER:  8.8482E-02 -8.8245E-01  9.5699E-01  5.1856E-01  2.3496E-01 -3.8047E-03 -6.1135E+00 -2.6388E+01 -1.4233E-01 -2.0371E+01
             1.2177E+00
 GRADIENT:  -2.2375E-01  4.4506E-01  1.5806E-01  1.6924E+00 -7.4061E-01 -3.2093E-02  0.0000E+00  0.0000E+00 -2.2156E-01  0.0000E+00
            -9.9864E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1376.59235021258        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  9.8740E-01  2.5483E-01  2.3731E+00  1.5978E+00  1.1171E+00  9.0092E-01  1.0000E-02  1.0000E-02  7.4525E-01  1.0000E-02
             3.0578E+00
 PARAMETER:  8.7323E-02 -1.2671E+00  9.6421E-01  5.6865E-01  2.1076E-01 -4.3367E-03 -9.5039E+00 -4.0394E+01 -1.9404E-01 -3.1247E+01
             1.2177E+00
 GRADIENT:  -3.7654E-01  3.3826E-01  4.6209E-02  2.5311E+00 -4.0711E-01 -4.8356E-02  0.0000E+00  0.0000E+00 -2.3401E-01  0.0000E+00
            -1.1935E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1376.59827613034        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1260
 NPARAMETR:  9.8679E-01  2.0263E-01  2.3607E+00  1.6300E+00  1.1013E+00  9.0074E-01  1.0000E-02  1.0000E-02  7.3063E-01  1.0000E-02
             3.0569E+00
 PARAMETER:  8.6705E-02 -1.4964E+00  9.5896E-01  5.8858E-01  1.9647E-01 -4.5363E-03 -1.1698E+01 -4.9386E+01 -2.1385E-01 -3.8236E+01
             1.2174E+00
 GRADIENT:  -2.9173E-01  1.1782E-01  8.9794E-02  1.0377E+00 -3.8641E-01 -6.1692E-03  0.0000E+00  0.0000E+00  8.7571E-03  0.0000E+00
            -1.6727E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1376.60168065897        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1435
 NPARAMETR:  9.8635E-01  1.5175E-01  2.3731E+00  1.6629E+00  1.0911E+00  9.0053E-01  1.0000E-02  1.0000E-02  7.1559E-01  1.0000E-02
             3.0568E+00
 PARAMETER:  8.6254E-02 -1.7855E+00  9.6418E-01  6.0858E-01  1.8721E-01 -4.7668E-03 -1.4593E+01 -6.1160E+01 -2.3465E-01 -4.7412E+01
             1.2174E+00
 GRADIENT:  -1.3325E-01  6.5647E-02  1.0578E-02  9.4835E-01 -8.0139E-02  1.8547E-03  0.0000E+00  0.0000E+00  2.6307E-02  0.0000E+00
            -1.7762E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1376.60328300109        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1611
 NPARAMETR:  9.8593E-01  1.1274E-01  2.3581E+00  1.6873E+00  1.0783E+00  9.0031E-01  1.0000E-02  1.0000E-02  7.0482E-01  1.0000E-02
             3.0571E+00
 PARAMETER:  8.5834E-02 -2.0827E+00  9.5787E-01  6.2310E-01  1.7538E-01 -5.0161E-03 -1.7663E+01 -7.3627E+01 -2.4981E-01 -5.7127E+01
             1.2175E+00
 GRADIENT:   6.8164E-03  6.2291E-02  7.7995E-02  1.0846E+00 -3.8356E-01 -1.5393E-03  0.0000E+00  0.0000E+00 -3.7911E-02  0.0000E+00
            -1.1790E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1376.60408186278        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1787
 NPARAMETR:  9.8583E-01  9.2339E-02  2.3710E+00  1.7002E+00  1.0759E+00  9.0019E-01  1.0000E-02  1.0000E-02  6.9899E-01  1.0000E-02
             3.0574E+00
 PARAMETER:  8.5732E-02 -2.2823E+00  9.6331E-01  6.3074E-01  1.7316E-01 -5.1484E-03 -1.9785E+01 -8.2201E+01 -2.5812E-01 -6.3823E+01
             1.2176E+00
 GRADIENT:   3.0447E-01  5.0258E-03  1.8394E-02  5.5355E-02 -2.4404E-02  6.5679E-03  0.0000E+00  0.0000E+00 -1.3079E-05  0.0000E+00
            -1.9469E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1376.60442664295        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1962
 NPARAMETR:  9.8543E-01  7.3429E-02  2.3639E+00  1.7122E+00  1.0699E+00  9.0005E-01  1.0000E-02  1.0000E-02  6.9396E-01  1.0000E-02
             3.0574E+00
 PARAMETER:  8.5322E-02 -2.5114E+00  9.6030E-01  6.3776E-01  1.6759E-01 -5.3097E-03 -2.2241E+01 -9.2128E+01 -2.6535E-01 -7.1569E+01
             1.2176E+00
 GRADIENT:  -1.6921E-01  1.9945E-02  2.0236E-02  5.3982E-01 -1.2912E-01 -1.5470E-02  0.0000E+00  0.0000E+00 -4.3627E-02  0.0000E+00
            -8.3826E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1376.60462338833        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2140
 NPARAMETR:  9.8534E-01  5.9527E-02  2.3642E+00  1.7209E+00  1.0666E+00  9.0001E-01  1.0000E-02  1.0000E-02  6.9034E-01  1.0000E-02
             3.0573E+00
 PARAMETER:  8.5234E-02 -2.7213E+00  9.6046E-01  6.4284E-01  1.6448E-01 -5.3499E-03 -2.4524E+01 -1.0134E+02 -2.7057E-01 -7.8758E+01
             1.2175E+00
 GRADIENT:   3.5939E-02  9.8208E-03  1.5908E-02  3.2609E-01 -8.8742E-02  2.2568E-04  0.0000E+00  0.0000E+00 -2.0340E-02  0.0000E+00
            -2.8361E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1376.60469776143        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2317
 NPARAMETR:  9.8517E-01  4.8698E-02  2.3667E+00  1.7278E+00  1.0645E+00  8.9992E-01  1.0000E-02  1.0000E-02  6.8743E-01  1.0000E-02
             3.0574E+00
 PARAMETER:  8.5059E-02 -2.9221E+00  9.6149E-01  6.4684E-01  1.6248E-01 -5.4444E-03 -2.6726E+01 -1.1021E+02 -2.7479E-01 -8.5688E+01
             1.2176E+00
 GRADIENT:  -1.0167E-01  5.4013E-03  2.8187E-03  2.2784E-01 -2.6150E-02 -8.3786E-03  0.0000E+00  0.0000E+00 -2.0060E-02  0.0000E+00
            -4.0052E-04

0ITERATION NO.:   84    OBJECTIVE VALUE:  -1376.60472232250        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     2459
 NPARAMETR:  9.8513E-01  4.2768E-02  2.3638E+00  1.7312E+00  1.0626E+00  8.9992E-01  1.0000E-02  1.0000E-02  6.8605E-01  1.0000E-02
             3.0574E+00
 PARAMETER:  8.5032E-02 -3.0489E+00  9.6033E-01  6.4888E-01  1.6065E-01 -5.4543E-03 -2.8124E+01 -1.1584E+02 -2.7683E-01 -9.0087E+01
             1.2175E+00
 GRADIENT:   3.7476E-02  1.5839E-03  7.0982E-03  1.2503E-01 -3.9443E-02 -1.4837E-04  0.0000E+00  0.0000E+00 -4.9526E-03  0.0000E+00
            -1.2487E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2459
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7623E-04 -1.3084E-05  6.0287E-05 -1.0103E-02 -1.7127E-05
 SE:             2.8781E-02  6.9223E-06  5.8142E-05  2.5647E-02  1.5628E-04
 N:                     100         100         100         100         100

 P VAL.:         9.9234E-01  5.8734E-02  2.9978E-01  6.9363E-01  9.1273E-01

 ETASHRINKSD(%)  3.5812E+00  9.9977E+01  9.9805E+01  1.4078E+01  9.9476E+01
 ETASHRINKVR(%)  7.0342E+00  1.0000E+02  1.0000E+02  2.6174E+01  9.9997E+01
 EBVSHRINKSD(%)  3.5566E+00  9.9978E+01  9.9784E+01  1.3901E+01  9.9458E+01
 EBVSHRINKVR(%)  6.9866E+00  1.0000E+02  1.0000E+02  2.5869E+01  9.9997E+01
 RELATIVEINF(%)  8.6450E+01  1.2317E-07  2.8086E-05  2.3734E+00  1.2536E-04
 EPSSHRINKSD(%)  2.2116E+01
 EPSSHRINKVR(%)  3.9341E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1376.6047223224998     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -641.45389575876163     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.55
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1376.605       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  4.29E-02  2.36E+00  1.73E+00  1.06E+00  9.00E-01  1.00E-02  1.00E-02  6.86E-01  1.00E-02  3.06E+00
 


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
+        1.31E+03
 
 TH 2
+       -6.46E+01  3.05E+02
 
 TH 3
+        4.69E-02  1.11E+01  8.33E+00
 
 TH 4
+       -8.23E+01  4.02E+02  1.90E+00  5.78E+02
 
 TH 5
+       -5.63E+00 -1.54E+02 -4.45E+01 -1.22E+02  2.83E+02
 
 TH 6
+       -7.38E+00 -1.46E-01  1.79E+00 -1.75E+01 -4.41E+00  2.19E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.63E+00 -7.67E+01  2.64E+00 -1.73E+01  1.74E+01  1.99E+00  0.00E+00  0.00E+00  2.23E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.79E+01 -1.15E+01  3.94E-01 -1.06E+01  1.74E+00  5.14E+00  0.00E+00  0.00E+00  1.85E+01  0.00E+00  4.50E+01
 
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
 #CPUT: Total CPU Time in Seconds,       34.428
Stop Time:
Sat Sep 18 10:09:12 CDT 2021
