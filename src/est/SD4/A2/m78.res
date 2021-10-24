Sun Oct 24 02:19:26 CDT 2021
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
$DATA ../../../../data/SD4/A2/dat78.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1044.88594681702        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7121E+02 -1.6806E+01 -5.7846E+00 -4.1821E+00  1.0292E+02  4.4657E+01 -2.6904E+01  4.6820E+00 -2.9034E+01 -4.8846E+01
            -1.1214E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1385.76091467705        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1954E+00  1.1248E+00  1.2531E+00  1.0552E+00  1.1032E+00  1.0646E+00  1.0498E+00  8.5395E-01  1.0371E+00  7.7628E-01
             2.7471E+00
 PARAMETER:  2.7851E-01  2.1761E-01  3.2562E-01  1.5369E-01  1.9820E-01  1.6259E-01  1.4858E-01 -5.7883E-02  1.3638E-01 -1.5324E-01
             1.1106E+00
 GRADIENT:   4.9793E+02  2.6099E+01 -6.2989E+00  4.5313E+01  1.1040E+01 -1.5876E+01  5.2322E+00  2.2571E+00  5.7434E+00  3.7564E+00
             3.8741E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1402.53214072813        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.1248E+00  9.1306E-01  3.6689E+00  1.1442E+00  1.3788E+00  1.0242E+00  5.1679E-01  5.1487E-01  9.8922E-01  9.5584E-01
             2.7408E+00
 PARAMETER:  2.1760E-01  9.0508E-03  1.3999E+00  2.3467E-01  4.2123E-01  1.2391E-01 -5.6011E-01 -5.6384E-01  8.9166E-02  5.4839E-02
             1.1082E+00
 GRADIENT:   3.7923E+02 -2.4574E+01  3.3215E-01 -2.9670E+01  1.7407E+01  1.2644E+00  8.3848E-01  3.7754E-02 -5.5122E+00  3.3447E+00
             3.1988E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1422.97963359739        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.6970E-01  1.0404E+00  1.5201E+00  1.0297E+00  1.1016E+00  9.0082E-01  4.1290E-01  2.1895E-01  1.2175E+00  9.1572E-01
             2.4417E+00
 PARAMETER:  6.9231E-02  1.3957E-01  5.1875E-01  1.2922E-01  1.9676E-01 -4.4469E-03 -7.8454E-01 -1.4189E+00  2.9676E-01  1.1952E-02
             9.9270E-01
 GRADIENT:   2.5412E+01  1.1799E+01  5.5058E+00  3.9156E+00 -8.4822E+00 -1.2115E+01  7.7001E-01  6.3944E-02  1.3489E+00 -2.4079E+00
            -9.4967E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1424.30214096467        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.8525E-01  9.9658E-01  1.5103E+00  1.0747E+00  1.1013E+00  9.5058E-01  4.7021E-01  1.9425E-01  1.1722E+00  9.1947E-01
             2.4991E+00
 PARAMETER:  8.5144E-02  9.6579E-02  5.1233E-01  1.7204E-01  1.9646E-01  4.9316E-02 -6.5458E-01 -1.5386E+00  2.5891E-01  1.6044E-02
             1.0159E+00
 GRADIENT:  -6.4626E+00  3.4846E+00  9.6191E-01  5.5334E-01 -1.9037E+00  2.2101E+00  2.6348E-03  6.0127E-02 -1.4778E+00 -5.2159E-01
            -2.9324E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1424.70564955371        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  9.8564E-01  7.4767E-01  1.5660E+00  1.2429E+00  1.0312E+00  9.4164E-01  3.0104E-01  5.2953E-02  1.0592E+00  9.1995E-01
             2.5016E+00
 PARAMETER:  8.5537E-02 -1.9080E-01  5.4855E-01  3.1748E-01  1.3075E-01  3.9863E-02 -1.1005E+00 -2.8383E+00  1.5750E-01  1.6562E-02
             1.0169E+00
 GRADIENT:  -2.1360E+00  6.0263E+00  2.3074E+00  1.0825E+01 -5.6106E+00 -6.7238E-01  1.2491E-01  4.7388E-03  6.7979E-01 -4.8926E-02
             1.3352E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1424.92888408683        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  9.8564E-01  5.5838E-01  1.4950E+00  1.3540E+00  9.5584E-01  9.4045E-01  9.8495E-02  1.0000E-02  9.7071E-01  9.2485E-01
             2.4725E+00
 PARAMETER:  8.5538E-02 -4.8271E-01  5.0211E-01  4.0304E-01  5.4839E-02  3.8608E-02 -2.2177E+00 -4.8750E+00  7.0269E-02  2.1880E-02
             1.0052E+00
 GRADIENT:   2.6831E+00  1.3863E+00  1.6184E-01  1.3700E+00 -1.9970E+00 -6.1911E-01  1.0553E-02  0.0000E+00 -1.2351E+00  1.4765E-02
             5.2449E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1424.96709307932        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  9.8296E-01  4.5039E-01  1.5374E+00  1.4255E+00  9.3972E-01  9.4123E-01  3.2587E-02  1.0000E-02  9.2562E-01  9.3640E-01
             2.4601E+00
 PARAMETER:  8.2811E-02 -6.9764E-01  5.3012E-01  4.5449E-01  3.7828E-02  3.9437E-02 -3.3238E+00 -6.5130E+00  2.2709E-02  3.4289E-02
             1.0002E+00
 GRADIENT:  -8.6280E-01  9.7198E-01  1.8381E-01  3.2051E+00 -5.2290E-02 -1.3487E-01  1.0332E-03  0.0000E+00 -3.6313E-01  1.2108E-02
            -1.5255E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1424.98281863321        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1075
 NPARAMETR:  9.8196E-01  3.5061E-01  1.5053E+00  1.4867E+00  9.0213E-01  9.4145E-01  1.0000E-02  1.0000E-02  8.8807E-01  9.3521E-01
             2.4468E+00
 PARAMETER:  8.1791E-02 -9.4808E-01  5.0897E-01  4.9655E-01 -3.0012E-03  3.9669E-02 -4.8183E+00 -8.6479E+00 -1.8706E-02  3.3012E-02
             9.9478E-01
 GRADIENT:  -3.5412E-01  9.3095E-01  1.0223E-01  3.9448E+00 -5.1688E-02  1.8597E-01  0.0000E+00  0.0000E+00 -5.6131E-01 -3.3101E-02
            -5.2758E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1425.03041832575        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1253
 NPARAMETR:  9.8033E-01  2.1412E-01  1.3998E+00  1.5628E+00  8.3128E-01  9.4129E-01  1.0000E-02  1.0000E-02  8.4678E-01  9.2662E-01
             2.4229E+00
 PARAMETER:  8.0132E-02 -1.4412E+00  4.3634E-01  5.4648E-01 -8.4785E-02  3.9495E-02 -8.2610E+00 -1.3291E+01 -6.6316E-02  2.3790E-02
             9.8495E-01
 GRADIENT:   6.1048E-01  6.0900E-01  3.5812E-02  2.2816E+00 -7.5919E-01  5.4404E-01  0.0000E+00  0.0000E+00 -4.8651E-01 -6.7233E-02
            -6.7996E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1425.04057552028        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1428
 NPARAMETR:  9.7921E-01  1.6204E-01  1.3549E+00  1.5933E+00  8.0077E-01  9.3967E-01  1.0000E-02  1.0000E-02  8.3195E-01  9.2275E-01
             2.4138E+00
 PARAMETER:  7.8987E-02 -1.7199E+00  4.0374E-01  5.6579E-01 -1.2218E-01  3.7776E-02 -1.0366E+01 -1.6027E+01 -8.3981E-02  1.9598E-02
             9.8120E-01
 GRADIENT:  -3.2136E-01  1.0031E+00  1.3710E+00  7.0736E+00 -4.0157E+00  5.7402E-02  0.0000E+00  0.0000E+00 -7.9221E-01  5.8397E-02
            -5.2183E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1425.08000508887        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1607
 NPARAMETR:  9.7770E-01  1.0381E-01  1.3030E+00  1.6241E+00  7.6836E-01  9.3647E-01  1.0000E-02  1.0000E-02  8.1841E-01  9.1713E-01
             2.4060E+00
 PARAMETER:  7.7449E-02 -2.1652E+00  3.6470E-01  5.8494E-01 -1.6349E-01  3.4360E-02 -1.3953E+01 -2.0560E+01 -1.0039E-01  1.3494E-02
             9.7795E-01
 GRADIENT:  -1.5876E+00  7.3583E-01  1.4822E+00  8.1026E+00 -4.5698E+00 -9.5957E-01  0.0000E+00  0.0000E+00 -4.1342E-01  1.2792E-01
             3.9620E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1425.18981896771        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1785
 NPARAMETR:  9.7664E-01  3.3985E-02  1.2801E+00  1.6613E+00  7.4566E-01  9.3780E-01  1.0000E-02  1.0000E-02  7.9999E-01  9.2045E-01
             2.3919E+00
 PARAMETER:  7.6366E-02 -3.2818E+00  3.4691E-01  6.0757E-01 -1.9349E-01  3.5780E-02 -2.3469E+01 -3.2280E+01 -1.2315E-01  1.7110E-02
             9.7207E-01
 GRADIENT:  -6.3615E-01  1.8721E-01  7.6333E-01  4.7107E+00 -2.2679E+00 -1.7803E-01  0.0000E+00  0.0000E+00 -3.8722E-01  1.0050E-01
            -1.5180E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1425.21910589894        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1962
 NPARAMETR:  9.7641E-01  1.1888E-02  1.3054E+00  1.6762E+00  7.5075E-01  9.3820E-01  1.0000E-02  1.0000E-02  7.9233E-01  9.2399E-01
             2.3926E+00
 PARAMETER:  7.6128E-02 -4.3322E+00  3.6650E-01  6.1655E-01 -1.8669E-01  3.6205E-02 -3.2694E+01 -4.3465E+01 -1.3278E-01  2.0949E-02
             9.7238E-01
 GRADIENT:  -2.2204E-01  5.2794E-02  2.9071E-01  3.6617E+00 -5.6189E-01  6.6725E-02  0.0000E+00  0.0000E+00 -4.2560E-01 -3.4007E-02
            -2.8861E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1425.23143581887        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2146
 NPARAMETR:  9.7639E-01  1.0000E-02  1.3000E+00  1.6730E+00  7.5010E-01  9.3770E-01  1.0000E-02  1.0000E-02  7.9283E-01  9.2379E-01
             2.3922E+00
 PARAMETER:  7.6208E-02 -4.8537E+00  3.6599E-01  6.1445E-01 -1.8830E-01  3.5883E-02 -3.7351E+01 -4.9096E+01 -1.3201E-01  2.1308E-02
             9.7203E-01
 GRADIENT:   2.1672E-01  0.0000E+00  6.2981E-01 -4.0708E-01 -6.8479E-01  6.5123E-02  0.0000E+00  0.0000E+00  3.2905E-02  3.4951E-02
            -8.9916E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2146
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.2045E-04 -2.8855E-06  6.7719E-05 -1.2241E-02 -2.5864E-02
 SE:             2.9062E-02  1.5372E-06  1.1621E-04  2.7203E-02  1.9765E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9395E-01  6.0513E-02  5.6009E-01  6.5270E-01  1.9069E-01

 ETASHRINKSD(%)  2.6402E+00  9.9995E+01  9.9611E+01  8.8679E+00  3.3784E+01
 ETASHRINKVR(%)  5.2106E+00  1.0000E+02  9.9998E+01  1.6949E+01  5.6154E+01
 EBVSHRINKSD(%)  2.4211E+00  9.9995E+01  9.9574E+01  8.4253E+00  3.4383E+01
 EBVSHRINKVR(%)  4.7837E+00  1.0000E+02  9.9998E+01  1.6141E+01  5.6944E+01
 RELATIVEINF(%)  9.0238E+01  8.6199E-09  1.2511E-04  4.2415E+00  1.8989E+00
 EPSSHRINKSD(%)  3.1816E+01
 EPSSHRINKVR(%)  5.3509E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1425.2314358188726     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -690.08060925513439     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1425.231       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  1.00E-02  1.30E+00  1.67E+00  7.50E-01  9.38E-01  1.00E-02  1.00E-02  7.93E-01  9.24E-01  2.39E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       61.003
Stop Time:
Sun Oct 24 02:19:39 CDT 2021
