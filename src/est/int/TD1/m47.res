Sat Sep 18 05:13:15 CDT 2021
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
$DATA ../../../../data/int/TD1/dat47.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3439.40227239724        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.5981E+01 -1.1151E+02  1.7525E+02 -1.1827E+02  5.2888E+01  2.5842E+01 -8.7250E+00 -5.1803E+02 -1.3687E+02 -7.8546E+00
            -2.1348E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3624.56283854567        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.2544E-01  1.1278E+00  1.0486E+00  9.8959E-01  1.0514E+00  9.6100E-01  9.9918E-01  1.6717E+00  8.8163E-01  1.0451E+00
             1.1194E+00
 PARAMETER:  2.2515E-02  2.2028E-01  1.4742E-01  8.9538E-02  1.5014E-01  6.0215E-02  9.9180E-02  6.1384E-01 -2.5986E-02  1.4416E-01
             2.1279E-01
 GRADIENT:  -1.1729E+02 -4.6543E-01  5.3559E+01 -4.4150E+01 -8.2147E+00  2.1443E+00 -1.5601E+00 -2.2680E+02 -5.5579E+01 -1.0181E+01
             7.7541E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3627.66884270717        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  9.4396E-01  1.0882E+00  1.0267E+00  9.9338E-01  1.1040E+00  8.7677E-01  1.3754E+00  1.7784E+00  7.7314E-01  1.0797E+00
             1.1204E+00
 PARAMETER:  4.2329E-02  1.8449E-01  1.2638E-01  9.3358E-02  1.9894E-01 -3.1516E-02  4.1877E-01  6.7573E-01 -1.5730E-01  1.7668E-01
             2.1369E-01
 GRADIENT:  -9.0516E+01 -4.3575E+01  1.9205E+01 -1.3685E+01  4.5507E+01 -3.2179E+01  3.2696E+01 -2.0741E+02 -3.8761E+01 -1.5118E+01
             8.5457E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3669.87801331796        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      240
 NPARAMETR:  9.6076E-01  1.1968E+00  1.2165E+00  9.4001E-01  1.2216E+00  9.2746E-01  1.3330E+00  2.2527E+00  8.4131E-01  1.1217E+00
             1.0901E+00
 PARAMETER:  5.9971E-02  2.7964E-01  2.9598E-01  3.8131E-02  3.0018E-01  2.4696E-02  3.8742E-01  9.1214E-01 -7.2791E-02  2.1489E-01
             1.8627E-01
 GRADIENT:  -2.9854E+01 -1.1174E+01  3.6392E+01 -2.5116E+01  5.8164E+01 -4.5642E+00  3.1730E+01 -1.2657E+02 -8.1215E+00 -1.8927E+01
             3.5874E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3678.52834847605        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      316
 NPARAMETR:  9.6282E-01  1.2054E+00  1.2100E+00  9.3829E-01  1.2252E+00  9.2842E-01  1.3186E+00  2.3978E+00  8.0600E-01  1.1692E+00
             1.0839E+00
 PARAMETER:  6.2115E-02  2.8681E-01  2.9065E-01  3.6301E-02  3.0309E-01  2.5732E-02  3.7656E-01  9.7454E-01 -1.1567E-01  2.5634E-01
             1.8054E-01
 GRADIENT:  -2.3685E+01 -8.8785E+00  2.9820E+01 -1.8901E+01  4.7031E+01 -3.8777E+00  2.5888E+01 -1.0290E+02 -5.7435E+00 -1.5192E+01
             2.8900E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3679.46353984206        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      460
 NPARAMETR:  9.6584E-01  1.2111E+00  1.2100E+00  9.4032E-01  1.2252E+00  9.3075E-01  1.2909E+00  2.3983E+00  8.1645E-01  1.1695E+00
             1.0839E+00
 PARAMETER:  6.5239E-02  2.9153E-01  2.9064E-01  3.8469E-02  3.0313E-01  2.8237E-02  3.5537E-01  9.7476E-01 -1.0279E-01  2.5655E-01
             1.8053E-01
 GRADIENT:  -1.5477E+01 -4.3690E+00  2.9768E+01 -1.2571E+01  4.7946E+01 -2.5584E+00  2.3029E+01 -1.0200E+02 -6.2800E+00 -1.3795E+01
             2.8193E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3679.80735065202        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      609             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6588E-01  1.2111E+00  1.2101E+00  9.4035E-01  1.2253E+00  9.4645E-01  1.2902E+00  2.3989E+00  9.0872E-01  1.1695E+00
             1.0838E+00
 PARAMETER:  6.5284E-02  2.9152E-01  2.9071E-01  3.8499E-02  3.0321E-01  4.4959E-02  3.5476E-01  9.7501E-01  4.2794E-03  2.5661E-01
             1.8048E-01
 GRADIENT:  -1.3604E+01 -4.0010E+00  2.9827E+01 -1.1657E+01  5.3741E+01  4.1130E+00  2.9815E+01 -9.7446E+01  2.7711E+00 -6.4210E+00
             2.8724E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3679.83383807668        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  9.6588E-01  1.2111E+00  1.2101E+00  9.4035E-01  1.2253E+00  9.4271E-01  1.2902E+00  2.3989E+00  8.9120E-01  1.1695E+00
             1.0838E+00
 PARAMETER:  6.5284E-02  2.9152E-01  2.9071E-01  3.8499E-02  3.0321E-01  4.1005E-02  3.5476E-01  9.7501E-01 -1.5190E-02  2.5661E-01
             1.8048E-01
 GRADIENT:  -4.8763E+01 -2.9166E+01  2.8535E+01 -1.8120E+01  4.1033E+01 -1.1287E+00  2.5466E+01 -1.0166E+02  4.8307E-01 -9.0628E+00
             2.8271E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3679.87420902775        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      864
 NPARAMETR:  9.6593E-01  1.2109E+00  1.2103E+00  9.4040E-01  1.2255E+00  9.4527E-01  1.2904E+00  2.4001E+00  8.8953E-01  1.1697E+00
             1.0837E+00
 PARAMETER:  6.5334E-02  2.9137E-01  2.9085E-01  3.8549E-02  3.0336E-01  4.3714E-02  3.5494E-01  9.7551E-01 -1.7060E-02  2.5674E-01
             1.8039E-01
 GRADIENT:  -1.3611E+01 -4.2386E+00  2.9717E+01 -1.1979E+01  5.2816E+01  3.6310E+00  2.8536E+01 -9.8161E+01  7.4652E-01 -7.7349E+00
             2.8438E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3688.84117432913        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      982
 NPARAMETR:  9.6956E-01  1.2174E+00  1.1633E+00  9.4415E-01  1.2058E+00  9.4322E-01  1.2173E+00  2.6201E+00  8.9408E-01  1.1851E+00
             1.0774E+00
 PARAMETER:  6.9089E-02  2.9674E-01  2.5123E-01  4.2534E-02  2.8718E-01  4.1544E-02  2.9667E-01  1.0632E+00 -1.1961E-02  2.6983E-01
             1.7452E-01
 GRADIENT:  -3.9943E+00  2.2835E+00  1.8292E+01 -3.6893E+00  3.5779E+01  3.0130E+00  1.8753E+01 -5.7157E+01 -1.6165E+00 -4.1101E+00
             1.9137E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3689.88553397066        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1112
 NPARAMETR:  9.6965E-01  1.2151E+00  1.1604E+00  9.4582E-01  1.2018E+00  9.4595E-01  1.2038E+00  2.6541E+00  9.0058E-01  1.1814E+00
             1.0760E+00
 PARAMETER:  6.9177E-02  2.9483E-01  2.4879E-01  4.4301E-02  2.8381E-01  4.4431E-02  2.8544E-01  1.0761E+00 -4.7135E-03  2.6674E-01
             1.7326E-01
 GRADIENT:  -3.4233E+00  2.1282E+00  1.7113E+01 -4.2945E+00  3.4263E+01  4.1589E+00  1.7536E+01 -5.2108E+01 -1.7699E+00 -3.7826E+00
             1.7376E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3689.88954369876        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1272
 NPARAMETR:  9.6965E-01  1.2151E+00  1.1604E+00  9.4582E-01  1.2018E+00  9.4887E-01  1.2038E+00  2.6544E+00  9.0058E-01  1.1814E+00
             1.0760E+00
 PARAMETER:  6.9179E-02  2.9483E-01  2.4879E-01  4.4302E-02  2.8381E-01  4.7521E-02  2.8544E-01  1.0762E+00 -4.7128E-03  2.6674E-01
             1.7326E-01
 GRADIENT:  -3.8467E+01 -2.3732E+01  1.6043E+01 -1.0491E+01  2.3050E+01  1.6740E+00  1.5173E+01 -5.5656E+01 -2.1866E+00 -5.2975E+00
             1.7078E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3690.46969241750        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1429
 NPARAMETR:  9.7041E-01  1.2163E+00  1.1576E+00  9.4618E-01  1.2003E+00  9.6983E-01  1.1962E+00  2.6780E+00  9.0121E-01  1.1830E+00
             1.0753E+00
 PARAMETER:  6.9961E-02  2.9578E-01  2.4631E-01  4.4674E-02  2.8261E-01  6.9369E-02  2.7916E-01  1.0851E+00 -4.0119E-03  2.6803E-01
             1.7261E-01
 GRADIENT:   3.9139E-01  3.3571E+00  1.5840E+01 -3.7416E+00  3.2277E+01  1.3790E+01  1.6632E+01 -4.9014E+01 -1.9310E+00 -3.2516E+00
             1.6471E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3690.91118318031        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1603
 NPARAMETR:  9.8354E-01  1.2162E+00  1.1576E+00  9.4619E-01  1.2004E+00  9.4591E-01  1.1962E+00  2.6783E+00  9.1528E-01  1.1876E+00
             1.0753E+00
 PARAMETER:  8.3406E-02  2.9575E-01  2.4633E-01  4.4685E-02  2.8264E-01  4.4395E-02  2.7919E-01  1.0852E+00  1.1478E-02  2.7196E-01
             1.7260E-01
 GRADIENT:  -3.5432E+00 -2.3125E+01  1.0231E+01 -1.0471E+01  1.6941E+01  1.0813E+00  1.5274E+01 -5.4209E+01 -3.3980E-01 -4.3303E-02
             1.5368E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -3690.91633639056        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1779
 NPARAMETR:  9.8487E-01  1.2162E+00  1.1576E+00  9.4619E-01  1.2004E+00  9.4298E-01  1.1962E+00  2.6783E+00  9.1809E-01  1.1877E+00
             1.0753E+00
 PARAMETER:  8.4757E-02  2.9575E-01  2.4633E-01  4.4685E-02  2.8264E-01  4.1285E-02  2.7919E-01  1.0852E+00  1.4540E-02  2.7200E-01
             1.7260E-01
 GRADIENT:  -2.1122E-01 -2.3135E+01  9.9017E+00 -1.0490E+01  1.6752E+01 -1.3430E-01  1.5490E+01 -5.4246E+01  6.7836E-03  3.6922E-01
             1.5286E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -3691.78615043988        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1972
 NPARAMETR:  9.7483E-01  1.2172E+00  1.1571E+00  9.4611E-01  1.1996E+00  9.1152E-01  1.1955E+00  2.7459E+00  9.0238E-01  1.1943E+00
             1.0753E+00
 PARAMETER:  7.4507E-02  2.9656E-01  2.4588E-01  4.4599E-02  2.8199E-01  7.3617E-03  2.7858E-01  1.1101E+00 -2.7179E-03  2.7753E-01
             1.7263E-01
 GRADIENT:  -2.7581E+01 -2.2343E+01  1.5255E+01 -8.4981E+00  1.9562E+01 -1.4310E+01  1.4155E+01 -4.3173E+01 -1.0964E+00 -4.2944E+00
             1.7254E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -3692.03292389268        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2154
 NPARAMETR:  9.7492E-01  1.2173E+00  1.1570E+00  9.4608E-01  1.1995E+00  9.4402E-01  1.1954E+00  2.7451E+00  9.0258E-01  1.1957E+00
             1.0754E+00
 PARAMETER:  7.4596E-02  2.9664E-01  2.4582E-01  4.4572E-02  2.8191E-01  4.2394E-02  2.7850E-01  1.1098E+00 -2.4943E-03  2.7871E-01
             1.7268E-01
 GRADIENT:  -2.5437E+01 -2.2392E+01  1.5213E+01 -8.5512E+00  1.9234E+01 -6.2975E-03  1.4157E+01 -4.3353E+01 -1.0237E+00 -3.9818E+00
             1.7512E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -3692.34876753436        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     2283
 NPARAMETR:  9.7565E-01  1.2194E+00  1.1518E+00  9.4656E-01  1.1974E+00  9.4384E-01  1.1828E+00  2.7451E+00  9.0346E-01  1.1969E+00
             1.0744E+00
 PARAMETER:  7.5351E-02  2.9832E-01  2.4131E-01  4.5079E-02  2.8015E-01  4.2197E-02  2.6787E-01  1.1098E+00 -1.5236E-03  2.7971E-01
             1.7177E-01
 GRADIENT:  -2.3534E+01 -2.0367E+01  8.4529E+00 -9.5297E+00  1.1976E+01 -2.8276E-02  1.2906E+01 -4.6633E+01 -1.6457E+00  1.5153E-01
             1.5432E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -3692.48264302527        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2468
 NPARAMETR:  9.8487E-01  1.2193E+00  1.1518E+00  9.4658E-01  1.1975E+00  9.4358E-01  1.1828E+00  2.7458E+00  9.1813E-01  1.1959E+00
             1.0744E+00
 PARAMETER:  8.4757E-02  2.9826E-01  2.4136E-01  4.5101E-02  2.8022E-01  4.1924E-02  2.6792E-01  1.1101E+00  1.4583E-02  2.7887E-01
             1.7173E-01
 GRADIENT:  -2.1632E-01 -2.0267E+01  9.9692E+00 -9.0811E+00  1.4354E+01  9.5039E-02  1.4287E+01 -4.5252E+01  3.4833E-02 -1.1888E-02
             1.5376E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -3692.54128886424        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2648
 NPARAMETR:  9.8498E-01  1.2194E+00  1.1518E+00  9.4657E-01  1.1974E+00  9.4412E-01  1.1828E+00  2.7567E+00  9.2572E-01  1.1678E+00
             1.0744E+00
 PARAMETER:  8.4870E-02  2.9837E-01  2.4131E-01  4.5093E-02  2.8013E-01  4.2494E-02  2.6785E-01  1.1140E+00  2.2813E-02  2.5515E-01
             1.7172E-01
 GRADIENT:   1.7293E-02 -1.8551E+01  1.4390E+01 -8.5824E+00  2.1711E+01  2.8807E-01  1.5027E+01 -4.0390E+01 -3.9553E-01 -7.1779E+00
             1.4474E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -3692.84139097385        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2844             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8527E-01  1.2199E+00  1.1515E+00  9.4655E-01  1.1970E+00  9.4662E-01  1.1824E+00  2.7972E+00  9.3268E-01  1.1353E+00
             1.0743E+00
 PARAMETER:  8.5164E-02  2.9877E-01  2.4110E-01  4.5066E-02  2.7984E-01  4.5145E-02  2.6757E-01  1.1286E+00  3.0309E-02  2.2693E-01
             1.7169E-01
 GRADIENT:   3.6160E+01  1.0182E+01  1.4035E+01 -3.7218E+00  3.4834E+01  4.8537E+00  1.7619E+01 -3.1524E+01 -1.6296E-01 -1.0894E+01
             1.3519E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -3692.84376372110        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3021
 NPARAMETR:  9.8527E-01  1.2199E+00  1.1515E+00  9.4655E-01  1.1970E+00  9.4301E-01  1.1824E+00  2.7972E+00  9.3315E-01  1.1353E+00
             1.0743E+00
 PARAMETER:  8.5165E-02  2.9876E-01  2.4110E-01  4.5065E-02  2.7984E-01  4.1320E-02  2.6757E-01  1.1286E+00  3.0807E-02  2.2693E-01
             1.7169E-01
 GRADIENT:   7.2980E-01 -1.6226E+01  1.2982E+01 -9.9062E+00  2.3748E+01 -1.9702E-01  1.5528E+01 -3.5185E+01 -5.4973E-01 -1.2102E+01
             1.3205E+01

0ITERATION NO.:  106    OBJECTIVE VALUE:  -3692.84376372110        NO. OF FUNC. EVALS.:  31
 CUMULATIVE NO. OF FUNC. EVALS.:     3052
 NPARAMETR:  9.8528E-01  1.2199E+00  1.1515E+00  9.4654E-01  1.1970E+00  9.4310E-01  1.1824E+00  2.7972E+00  9.3314E-01  1.1353E+00
             1.0743E+00
 PARAMETER:  8.5165E-02  2.9876E-01  2.4110E-01  4.5065E-02  2.7984E-01  4.1320E-02  2.6757E-01  1.1286E+00  3.0807E-02  2.2693E-01
             1.7169E-01
 GRADIENT:  -2.5181E+05  8.4261E+04  5.2237E+04  2.5179E+05  2.1000E+01 -1.8761E-01  1.2976E+01 -1.1886E+02  2.5181E+05  1.1095E+05
             3.8378E+00
 NUMSIGDIG:         3.3         3.3         3.3         3.3         7.2         2.1         7.5         5.9         3.3         3.3
                    8.2

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3052
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0495E-05 -2.6115E-02 -4.4485E-02  3.4356E-02 -6.2964E-02
 SE:             2.9906E-02  2.2404E-02  2.7774E-02  2.5141E-02  2.4662E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9839E-01  2.4375E-01  1.0923E-01  1.7177E-01  1.0678E-02

 ETASHRINKSD(%)  1.0000E-10  2.4944E+01  6.9535E+00  1.5775E+01  1.7379E+01
 ETASHRINKVR(%)  1.0000E-10  4.3666E+01  1.3424E+01  2.9061E+01  3.1737E+01
 EBVSHRINKSD(%)  3.2652E-01  1.8971E+01  1.6633E+01  1.9315E+01  2.0214E+01
 EBVSHRINKVR(%)  6.5197E-01  3.4342E+01  3.0500E+01  3.4900E+01  3.6342E+01
 RELATIVEINF(%)  9.9346E+01  3.2436E+01  6.4482E+01  3.3924E+01  3.6155E+01
 EPSSHRINKSD(%)  2.3660E+01
 EPSSHRINKVR(%)  4.1722E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3692.8437637211009     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2038.7544039526902     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    93.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3692.844       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  1.22E+00  1.15E+00  9.47E-01  1.20E+00  9.43E-01  1.18E+00  2.80E+00  9.33E-01  1.14E+00  1.07E+00
 


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
+        6.48E+08
 
 TH 2
+        1.75E+08  9.48E+07
 
 TH 3
+       -2.30E+08  1.24E+01  1.63E+08
 
 TH 4
+        1.42E+03  2.87E+02  2.40E+08  1.41E+09
 
 TH 5
+        9.27E-01 -5.16E+07 -4.33E+01  1.94E+02  1.12E+08
 
 TH 6
+       -1.72E+03  4.60E+02  6.07E+02  1.78E+03 -2.22E+00  7.08E+08
 
 TH 7
+        2.11E+00  2.51E+01  7.17E+07  2.10E+08 -1.23E+01 -1.09E-01  1.26E+08
 
 TH 8
+       -8.13E-02  5.48E+06 -1.41E+01  2.83E+00 -2.23E+04 -3.18E-01 -6.29E+06  6.29E+05
 
 TH 9
+        2.63E+03 -7.05E+02 -9.24E+02 -2.73E+03  2.30E+01  1.80E+03  3.10E+01  5.88E+00  7.23E+08
 
 TH10
+        1.80E+03 -5.12E+02 -6.29E+02 -1.85E+03 -3.67E+01  6.56E+02  3.68E-01 -7.37E-01 -9.82E+02  9.48E+07
 
 TH11
+       -9.60E+00 -3.35E+01  7.84E+03 -3.61E+08 -3.15E+03  9.63E-01 -1.08E+08  1.08E+07  1.03E+01  2.52E+01  3.70E+08
 
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
 #CPUT: Total CPU Time in Seconds,      108.627
Stop Time:
Sat Sep 18 05:15:05 CDT 2021
