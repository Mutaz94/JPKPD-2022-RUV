Wed Sep 29 03:49:06 CDT 2021
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
$DATA ../../../../data/int/SL3/dat2.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      986
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      886
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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -495.033264901265        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6484E+02  5.3718E+01  9.5842E+01  1.3110E+02  9.3157E+01  3.5192E+01 -9.0497E+01 -1.6426E+02 -6.1606E+01 -2.7139E+01
            -6.2201E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2696.17911685139        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0584E+00  1.3027E+00  1.0385E+00  9.1449E-01  1.1598E+00  9.5770E-01  9.8918E-01  8.8930E-01  7.2997E-01  1.0939E+00
             2.7186E+00
 PARAMETER:  1.5678E-01  3.6446E-01  1.3775E-01  1.0616E-02  2.4825E-01  5.6779E-02  8.9120E-02 -1.7317E-02 -2.1475E-01  1.8979E-01
             1.1001E+00
 GRADIENT:   2.2567E+02  1.0976E+02 -6.6733E+00  7.0491E+01 -1.3545E+01 -1.1934E+01  1.5797E+01  3.6364E+00 -1.3759E+01 -2.1896E+01
             3.4451E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2702.61738962523        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0242E+00  1.2878E+00  1.1028E+00  9.1992E-01  1.1940E+00  9.5397E-01  1.0045E+00  2.2278E-01  9.3675E-01  1.2964E+00
             2.7327E+00
 PARAMETER:  1.2392E-01  3.5294E-01  1.9789E-01  1.6530E-02  2.7730E-01  5.2879E-02  1.0444E-01 -1.4016E+00  3.4664E-02  3.5960E-01
             1.1053E+00
 GRADIENT:   1.3106E+02  7.3227E+01 -1.2539E+00  7.0261E+01  3.8957E+00 -7.0882E+00  2.6589E+01  2.2189E-01  1.3826E+01  1.5525E+01
             6.3763E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2715.40682913017        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      317
 NPARAMETR:  9.9385E-01  1.5270E+00  1.2063E+00  7.3751E-01  1.4567E+00  9.7468E-01  6.1196E-01  3.9285E-01  1.0070E+00  1.4384E+00
             2.6631E+00
 PARAMETER:  9.3832E-02  5.2331E-01  2.8754E-01 -2.0448E-01  4.7617E-01  7.4359E-02 -3.9109E-01 -8.3434E-01  1.0699E-01  4.6352E-01
             1.0795E+00
 GRADIENT:  -7.7357E+00 -1.8100E+01 -8.4603E-01  1.7239E+01  7.4239E+00 -1.9043E+00 -8.6263E+00  4.5230E-02 -6.1614E+00 -1.7980E+00
            -1.8392E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2722.22649382703        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      492
 NPARAMETR:  9.9924E-01  2.0376E+00  1.2576E+00  4.2621E-01  1.9016E+00  9.7359E-01  5.9584E-01  5.9422E-01  1.2781E+00  1.7244E+00
             2.6544E+00
 PARAMETER:  9.9241E-02  8.1176E-01  3.2918E-01 -7.5282E-01  7.4267E-01  7.3237E-02 -4.1779E-01 -4.2051E-01  3.4537E-01  6.4488E-01
             1.0762E+00
 GRADIENT:   3.1222E+00  5.2523E+01 -5.3286E-01  2.2796E+01  6.1053E+00 -3.1202E+00 -9.4234E+00 -2.1877E-01 -3.4349E+00  2.0105E+00
             3.6845E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2724.73062913599        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      672
 NPARAMETR:  9.9737E-01  2.3212E+00  1.0587E+00  2.1501E-01  2.0833E+00  9.7992E-01  5.8164E-01  3.3201E+00  1.5497E+00  1.8418E+00
             2.6315E+00
 PARAMETER:  9.7371E-02  9.4207E-01  1.5707E-01 -1.4371E+00  8.3396E-01  7.9718E-02 -4.4190E-01  1.3000E+00  5.3808E-01  7.1073E-01
             1.0675E+00
 GRADIENT:   6.6992E-01 -2.7150E+00 -7.8951E-01  1.1718E+00 -1.5160E+00 -7.1668E-01 -9.3186E+00  2.9154E-01  1.6051E-01 -5.1836E-02
            -7.3930E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2724.96311167673        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      830
 NPARAMETR:  9.9632E-01  2.3378E+00  1.0922E+00  2.0220E-01  2.1050E+00  9.8244E-01  6.0508E-01  3.4710E+00  1.5288E+00  1.8529E+00
             2.6383E+00
 PARAMETER:  9.6310E-02  9.4922E-01  1.8822E-01 -1.4985E+00  8.4432E-01  8.2289E-02 -4.0239E-01  1.3444E+00  5.2451E-01  7.1676E-01
             1.0701E+00
 GRADIENT:   6.0535E+01  3.5171E+02 -3.3006E-01  8.8986E+00  3.0209E+01  5.9927E+00  5.7533E+00 -9.4974E-02  9.9417E-01  8.9166E+00
             2.0405E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2725.16074235751        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  9.9857E-01  2.3497E+00  1.5338E+00  2.0179E-01  2.1231E+00  9.8217E-01  6.0547E-01  5.2867E+00  1.3139E+00  1.8666E+00
             2.6349E+00
 PARAMETER:  9.8567E-02  9.5427E-01  5.2776E-01 -1.5005E+00  8.5288E-01  8.2010E-02 -4.0174E-01  1.7652E+00  3.7298E-01  7.2412E-01
             1.0689E+00
 GRADIENT:   2.9149E+00  1.2732E+01 -4.1631E-01  2.9162E+00 -3.4996E+00 -2.5807E-01 -2.4249E+00  1.1606E+00  9.8799E-01 -8.7246E-01
            -2.0352E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2725.28674433010        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1188
 NPARAMETR:  9.9713E-01  2.3470E+00  1.9520E+00  1.9965E-01  2.1533E+00  9.8290E-01  6.1705E-01  6.2976E+00  1.0968E+00  1.8857E+00
             2.6336E+00
 PARAMETER:  9.7128E-02  9.5314E-01  7.6884E-01 -1.5112E+00  8.6698E-01  8.2755E-02 -3.8281E-01  1.9402E+00  1.9236E-01  7.3430E-01
             1.0684E+00
 GRADIENT:  -1.6325E-01 -3.6948E-01 -8.2189E-02  2.5918E-02  6.8433E-02 -3.5886E-02  2.3853E-01  1.5514E-01  1.2192E-01 -1.1073E-01
             3.1584E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2725.29706981463        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1365
 NPARAMETR:  9.9719E-01  2.3429E+00  1.9562E+00  2.0244E-01  2.1513E+00  9.8297E-01  6.1867E-01  6.2914E+00  1.0233E+00  1.8847E+00
             2.6334E+00
 PARAMETER:  9.7182E-02  9.5140E-01  7.7100E-01 -1.4973E+00  8.6605E-01  8.2822E-02 -3.8019E-01  1.9392E+00  1.2305E-01  7.3377E-01
             1.0683E+00
 GRADIENT:  -3.8581E-02  1.6358E-01 -8.4883E-02  9.1765E-02 -6.8677E-02 -6.2390E-03 -1.8276E-01  1.4590E-01 -4.2868E-03 -1.1581E-01
            -1.0737E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2725.29807208742        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1544
 NPARAMETR:  9.9720E-01  2.3426E+00  1.9646E+00  2.0266E-01  2.1515E+00  9.8299E-01  6.1885E-01  6.2754E+00  1.0350E+00  1.8852E+00
             2.6332E+00
 PARAMETER:  9.7197E-02  9.5127E-01  7.7529E-01 -1.4962E+00  8.6616E-01  8.2843E-02 -3.7990E-01  1.9366E+00  1.3440E-01  7.3402E-01
             1.0682E+00
 GRADIENT:  -1.6738E-04  1.4145E-01 -1.7429E-02 -5.8680E-02 -1.3794E-02  1.7736E-04  4.8475E-03  1.1570E-02 -3.8521E-03  4.0217E-04
            -1.5197E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2725.30520259613        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1725
 NPARAMETR:  9.9719E-01  2.3386E+00  2.0394E+00  2.0589E-01  2.1527E+00  9.8300E-01  6.1868E-01  6.2710E+00  1.0666E+00  1.8849E+00
             2.6334E+00
 PARAMETER:  9.7189E-02  9.4956E-01  8.1264E-01 -1.4804E+00  8.6673E-01  8.2852E-02 -3.8017E-01  1.9359E+00  1.6445E-01  7.3389E-01
             1.0683E+00
 GRADIENT:  -4.9934E-02  1.0096E+00 -1.0929E-01  2.8959E-01  9.2408E-02 -6.5923E-03 -4.5329E-02  2.2504E-01  8.5858E-02 -1.3941E-01
             7.5192E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2725.38344507529        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1895
 NPARAMETR:  9.9680E-01  2.2790E+00  2.7461E+00  2.4438E-01  2.1506E+00  9.8291E-01  6.3597E-01  6.2556E+00  9.2536E-01  1.8826E+00
             2.6327E+00
 PARAMETER:  9.6798E-02  9.2375E-01  1.1102E+00 -1.3090E+00  8.6575E-01  8.2766E-02 -3.5260E-01  1.9335E+00  2.2429E-02  7.3263E-01
             1.0680E+00
 GRADIENT:  -7.5286E-01 -6.7419E+00 -3.0400E-01 -4.2454E-01  2.8934E-02 -3.0032E-02  1.3278E+00  8.7073E-01  1.0896E-01 -3.7478E-01
             6.4671E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2725.40077208585        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2073
 NPARAMETR:  9.9702E-01  2.2810E+00  2.5020E+00  2.4526E-01  2.1456E+00  9.8280E-01  6.3034E-01  5.8150E+00  9.7810E-01  1.8775E+00
             2.6335E+00
 PARAMETER:  9.7014E-02  9.2463E-01  1.0171E+00 -1.3055E+00  8.6344E-01  8.2647E-02 -3.6149E-01  1.8604E+00  7.7861E-02  7.2996E-01
             1.0683E+00
 GRADIENT:  -4.2304E-01  2.0529E-01  1.6436E-01 -2.5586E-01  6.0101E-01 -5.6439E-02 -3.3327E-01 -3.2905E-01 -9.4199E-02  9.0905E-02
             8.0761E-01

0ITERATION NO.:   69    OBJECTIVE VALUE:  -2725.40346729151        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     2204
 NPARAMETR:  9.9726E-01  2.2801E+00  2.3720E+00  2.4479E-01  2.1388E+00  9.8295E-01  6.3080E-01  5.6697E+00  9.8716E-01  1.8759E+00
             2.6328E+00
 PARAMETER:  9.7183E-02  9.2470E-01  9.6355E-01 -1.3071E+00  8.6063E-01  8.2784E-02 -3.6077E-01  1.8389E+00  8.8079E-02  7.2869E-01
             1.0680E+00
 GRADIENT:  -1.9150E-01  2.4757E+00 -4.0336E-03  5.4929E-02  2.3733E-01 -8.5544E-03 -1.0223E-04  4.5433E-01  3.3429E-03 -1.2288E-01
            -3.4054E-02
 NUMSIGDIG:         3.1         3.2         3.7         3.7         3.3         3.6         6.0         2.6         1.9         3.2
                    4.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2204
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0553E-03 -1.0342E-02 -2.0592E-02 -2.0976E-04 -1.8320E-02
 SE:             2.9351E-02  2.7010E-02  7.5774E-03  7.6520E-03  2.6667E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7132E-01  7.0180E-01  6.5774E-03  9.7813E-01  4.9209E-01

 ETASHRINKSD(%)  1.6702E+00  9.5130E+00  7.4615E+01  7.4365E+01  1.0662E+01
 ETASHRINKVR(%)  3.3124E+00  1.8121E+01  9.3556E+01  9.3428E+01  2.0187E+01
 EBVSHRINKSD(%)  1.6266E+00  9.6618E+00  8.2394E+01  7.6903E+01  8.0036E+00
 EBVSHRINKVR(%)  3.2268E+00  1.8390E+01  9.6900E+01  9.4665E+01  1.5367E+01
 RELATIVEINF(%)  9.6695E+01  2.0334E+00  7.8298E-01  1.2386E-01  5.7722E+01
 EPSSHRINKSD(%)  1.5934E+01
 EPSSHRINKVR(%)  2.9329E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          886
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1628.3590808386800     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2725.4034672915068     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1097.0443864528268     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    57.16
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2725.403       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.97E-01  2.28E+00  2.37E+00  2.45E-01  2.14E+00  9.83E-01  6.31E-01  5.69E+00  9.88E-01  1.88E+00  2.63E+00
 


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
+        1.12E+03
 
 TH 2
+       -1.79E+01  4.40E+02
 
 TH 3
+       -9.28E-02  3.33E-01  1.83E+00
 
 TH 4
+       -2.38E+01  6.14E+02 -5.30E+01  3.53E+03
 
 TH 5
+       -2.25E+00 -1.34E+01 -5.87E-01 -2.14E+01  5.92E+01
 
 TH 6
+        6.05E+00 -6.16E+00  5.80E-02 -9.81E+00 -8.74E-01  1.92E+02
 
 TH 7
+        4.67E+00  2.24E+00  5.31E-01 -9.30E+01 -4.44E-01 -2.41E+00  3.35E+02
 
 TH 8
+        3.74E-01 -1.36E+00 -9.65E+00  1.01E+02 -1.05E+01 -1.25E-01 -9.92E-01  3.69E+00
 
 TH 9
+        2.71E-01 -2.62E+00 -1.54E+00  3.68E+01 -9.83E-01 -1.21E-01  1.19E+01  3.49E+00  2.81E+00
 
 TH10
+       -2.99E-02 -1.28E-01  1.14E+00 -4.04E+01 -3.16E+00  4.36E-01 -8.28E-01 -1.47E+01 -1.53E+00  4.10E+01
 
 TH11
+       -1.60E+01 -1.59E+01  1.47E+00 -8.37E+01  1.80E+00  2.72E+00  1.43E+01 -1.03E+01 -5.67E-01  6.46E+00  1.69E+02
 
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
 #CPUT: Total CPU Time in Seconds,       70.068
Stop Time:
Wed Sep 29 03:50:18 CDT 2021
