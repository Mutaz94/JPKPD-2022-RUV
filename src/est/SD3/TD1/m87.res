Sun Oct 24 00:35:08 CDT 2021
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
$DATA ../../../../data/SD3/TD1/dat87.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1677.48719457934        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9383E+02 -7.2562E+01 -5.3527E+00 -1.0203E+01  9.0298E+01  4.4535E+01 -3.9021E+00 -1.6951E+02 -9.3093E+00  1.5044E+01
            -7.3350E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2083.78500052826        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:      116
 NPARAMETR:  8.7312E-01  1.0469E+00  1.0214E+00  9.9096E-01  9.4921E-01  1.0717E+00  1.0000E+00  1.1330E+00  1.0022E+00  9.6997E-01
             1.5301E+00
 PARAMETER: -3.5684E-02  1.4579E-01  1.2117E-01  9.0923E-02  4.7878E-02  1.6925E-01  1.0003E-01  2.2484E-01  1.0222E-01  6.9507E-02
             5.2534E-01
 GRADIENT:  -3.0590E+02 -4.8893E+01 -2.3498E+01 -6.0294E+01 -2.9519E+01 -2.5298E+01  2.8015E+00  1.8985E+01  4.1528E+00  2.7883E+01
             2.6366E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2120.05535025868        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      281
 NPARAMETR:  1.0333E+00  1.0721E+00  1.0892E+00  1.0219E+00  1.0012E+00  1.0036E+00  9.3749E-01  5.8041E-01  9.8000E-01  8.3129E-01
             1.5312E+00
 PARAMETER:  1.3279E-01  1.6959E-01  1.8546E-01  1.2165E-01  1.0124E-01  1.0361E-01  3.5446E-02 -4.4402E-01  7.9793E-02 -8.4780E-02
             5.2604E-01
 GRADIENT:   1.7694E+01 -9.0379E+00 -1.0037E+01 -1.6347E+00  1.1958E+01 -1.4850E+00 -3.1331E+00  1.9774E+00 -4.4125E+00 -1.5399E-01
             2.4945E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2120.72993236997        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      456
 NPARAMETR:  1.0252E+00  1.0762E+00  1.1350E+00  1.0209E+00  1.0226E+00  1.0147E+00  9.5990E-01  3.4521E-01  1.0066E+00  8.8129E-01
             1.5309E+00
 PARAMETER:  1.2489E-01  1.7341E-01  2.2662E-01  1.2069E-01  1.2231E-01  1.1455E-01  5.9078E-02 -9.6360E-01  1.0654E-01 -2.6366E-02
             5.2584E-01
 GRADIENT:   1.0055E+00 -5.1208E+00 -1.3704E+00 -3.3306E+00  5.3974E+00  2.9939E+00  4.2182E-01  3.0558E-01  1.0926E+00  1.4306E+00
             2.4940E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2122.14958115353        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      634
 NPARAMETR:  1.0313E+00  1.2602E+00  8.7261E-01  8.9814E-01  9.7718E-01  1.0314E+00  8.6362E-01  1.0000E-02  1.0763E+00  7.7416E-01
             1.5203E+00
 PARAMETER:  1.3080E-01  3.3126E-01 -3.6270E-02 -7.4307E-03  7.6918E-02  1.3089E-01 -4.6622E-02 -1.0723E+01  1.7356E-01 -1.5598E-01
             5.1890E-01
 GRADIENT:   8.8240E+00  4.5258E+00 -8.3110E-01  1.8410E+00 -8.6428E+00  8.0410E+00 -3.3926E+00  0.0000E+00 -2.1706E+00 -1.9789E+00
             2.4558E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2126.89805532629        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      812
 NPARAMETR:  1.0185E+00  1.8183E+00  7.8554E-01  5.6374E-01  1.2407E+00  9.8498E-01  7.3118E-01  1.0000E-02  1.5675E+00  9.7952E-01
             1.4756E+00
 PARAMETER:  1.1829E-01  6.9792E-01 -1.4139E-01 -4.7317E-01  3.1568E-01  8.4863E-02 -2.1309E-01 -5.3873E+01  5.4950E-01  7.9306E-02
             4.8907E-01
 GRADIENT:  -2.1733E+01  4.7393E+01  8.0960E+00  1.3442E+01 -1.0760E+01 -1.1118E+01  3.4798E-01  0.0000E+00  3.6138E+00 -3.7877E-01
             2.3478E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2137.28146954588        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  1.0331E+00  2.3612E+00  3.7277E-01  2.1390E-01  1.4383E+00  1.0170E+00  6.8971E-01  1.0000E-02  2.5865E+00  1.1186E+00
             1.3695E+00
 PARAMETER:  1.3260E-01  9.5916E-01 -8.8678E-01 -1.4422E+00  4.6348E-01  1.1689E-01 -2.7149E-01 -1.6103E+02  1.0503E+00  2.1212E-01
             4.1442E-01
 GRADIENT:   1.0074E+01  7.1661E+01  3.6687E+00  3.2903E+00 -2.4031E-01  9.7178E-01 -2.4851E+00  0.0000E+00 -5.5398E-01 -3.3226E+00
             2.0296E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2144.81255594242        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:     1142
 NPARAMETR:  1.0320E+00  2.3748E+00  3.5473E-01  2.0117E-01  1.4478E+00  1.0183E+00  6.8703E-01  1.0000E-02  2.6620E+00  1.1373E+00
             1.3163E+00
 PARAMETER:  1.3151E-01  9.6492E-01 -9.3640E-01 -1.5036E+00  4.7004E-01  1.1813E-01 -2.7538E-01 -1.6833E+02  1.0791E+00  2.2862E-01
             3.7481E-01
 GRADIENT:   8.9633E+00  6.3818E+01  2.9879E+00  2.1671E+00  1.3021E+00  1.3198E+00 -3.6943E+00  0.0000E+00 -5.2734E-01 -3.5593E+00
             1.8582E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2144.90368608523        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     1287
 NPARAMETR:  1.0327E+00  2.3631E+00  3.5640E-01  1.9966E-01  1.4512E+00  1.0067E+00  6.9496E-01  1.0000E-02  2.6477E+00  1.1651E+00
             1.3187E+00
 PARAMETER:  1.3217E-01  9.5997E-01 -9.3172E-01 -1.5111E+00  4.7240E-01  1.0669E-01 -2.6390E-01 -1.6833E+02  1.0737E+00  2.5285E-01
             3.7667E-01
 GRADIENT:   1.0812E+01  3.9872E+01  2.6853E+00 -5.6702E-01  3.3915E-01 -3.1199E+00 -2.2646E-01  0.0000E+00 -6.2844E-01  9.3707E-02
             1.8837E+02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2156.90951311402        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1424
 NPARAMETR:  1.0132E+00  2.1982E+00  3.6239E-01  2.0020E-01  1.4307E+00  9.9041E-01  6.8136E-01  1.0000E-02  2.5112E+00  1.1651E+00
             1.1961E+00
 PARAMETER:  1.1311E-01  8.8764E-01 -9.1503E-01 -1.5084E+00  4.5815E-01  9.0362E-02 -2.8366E-01 -1.6833E+02  1.0208E+00  2.5284E-01
             2.7907E-01
 GRADIENT:   3.5123E+02  1.0606E+03 -7.5814E-01  2.4386E+01  2.6203E+01  4.1615E+01  8.8349E+00  0.0000E+00  2.0892E+01  3.7786E+00
             1.4759E+02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2178.99922615824        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1601
 NPARAMETR:  1.0233E+00  2.2913E+00  3.8326E-01  2.2953E-01  1.4293E+00  1.0165E+00  6.9501E-01  1.0000E-02  2.4786E+00  1.2147E+00
             9.6300E-01
 PARAMETER:  1.2307E-01  9.2911E-01 -8.5904E-01 -1.3717E+00  4.5719E-01  1.1638E-01 -2.6383E-01 -1.6833E+02  1.0077E+00  2.9453E-01
             6.2299E-02
 GRADIENT:  -9.2343E-01 -3.8197E+00 -1.9907E-01 -1.6979E+00 -1.6134E+00  2.3069E-02 -3.6364E+00  0.0000E+00  1.2763E-01 -8.2563E-02
            -3.2647E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2179.69331933879        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1778            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0253E+00  2.2507E+00  3.9463E-01  2.3990E-01  1.4282E+00  1.0178E+00  7.0690E-01  1.0000E-02  2.4330E+00  1.2160E+00
             9.6625E-01
 PARAMETER:  1.2503E-01  9.1124E-01 -8.2982E-01 -1.3275E+00  4.5640E-01  1.1760E-01 -2.4686E-01 -1.6833E+02  9.8912E-01  2.9557E-01
             6.5665E-02
 GRADIENT:   6.4308E+02  1.9047E+03 -6.2779E-01  8.6680E+01  2.9661E+01  9.0433E+01  3.1037E+01  0.0000E+00  4.0945E+01  4.9814E+00
             1.8790E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2182.03215950429        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     1922
 NPARAMETR:  1.0214E+00  2.1005E+00  5.1130E-01  3.2929E-01  1.3718E+00  1.0131E+00  6.9737E-01  1.0000E-02  2.1525E+00  1.1440E+00
             9.6180E-01
 PARAMETER:  1.2118E-01  8.4217E-01 -5.7080E-01 -1.0108E+00  4.1611E-01  1.1306E-01 -2.6044E-01 -1.6833E+02  8.6664E-01  2.3449E-01
             6.1056E-02
 GRADIENT:   6.1956E+02  1.5767E+03  2.4649E+00  1.0949E+02  3.6387E+01  8.6389E+01  2.2164E+01  0.0000E+00  4.8017E+01  7.3558E-01
            -5.9905E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2182.53715099898        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2099
 NPARAMETR:  1.0220E+00  2.1169E+00  4.8225E-01  3.4052E-01  1.3459E+00  1.0148E+00  7.1735E-01  1.0000E-02  2.0061E+00  1.1447E+00
             9.6306E-01
 PARAMETER:  1.2172E-01  8.4996E-01 -6.2930E-01 -9.7729E-01  3.9709E-01  1.1469E-01 -2.3219E-01 -1.6833E+02  7.9619E-01  2.3514E-01
             6.2358E-02
 GRADIENT:  -3.6200E+00 -3.9380E+00 -2.2692E-01  4.3908E-02  7.6718E-01 -5.7508E-01 -4.4953E-01  0.0000E+00  4.6092E-01  6.6135E-01
             1.7303E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2182.87448440766        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2286
 NPARAMETR:  1.0244E+00  2.0942E+00  4.9504E-01  3.5343E-01  1.3356E+00  1.0168E+00  7.1770E-01  1.0000E-02  1.9860E+00  1.1351E+00
             9.6236E-01
 PARAMETER:  1.2410E-01  8.3916E-01 -6.0312E-01 -9.4007E-01  3.8939E-01  1.1671E-01 -2.3170E-01 -1.6833E+02  7.8611E-01  2.2676E-01
             6.1631E-02
 GRADIENT:   1.6956E+00 -7.6819E+00 -4.2065E-01  1.1059E-01  7.8792E-01  2.6815E-01 -4.5530E-01  0.0000E+00  1.9344E+00  7.7359E-01
            -1.4767E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2183.28703882557        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     2434
 NPARAMETR:  1.0226E+00  2.0525E+00  5.0523E-01  3.6948E-01  1.3213E+00  1.0128E+00  7.1800E-01  1.0000E-02  1.9378E+00  1.1258E+00
             9.6125E-01
 PARAMETER:  1.2235E-01  8.1907E-01 -5.8274E-01 -8.9566E-01  3.7863E-01  1.1272E-01 -2.3128E-01 -1.6833E+02  7.6158E-01  2.1845E-01
             6.0476E-02
 GRADIENT:   6.2702E+02  1.4667E+03 -4.8021E-01  1.1943E+02  2.6975E+01  8.5368E+01  2.3053E+01  0.0000E+00  3.8233E+01  4.3694E+00
             8.0334E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2183.47575741037        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     2595             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0285E+00  2.0446E+00  5.1773E-01  3.7731E-01  1.3163E+00  1.0213E+00  7.2515E-01  1.0000E-02  1.8914E+00  1.1141E+00
             9.6197E-01
 PARAMETER:  1.2808E-01  8.1521E-01 -5.5831E-01 -8.7468E-01  3.7486E-01  1.2107E-01 -2.2137E-01 -1.6833E+02  7.3734E-01  2.0807E-01
             6.1232E-02
 GRADIENT:   6.6788E+02  1.4487E+03  1.8165E+00  1.2022E+02  2.4780E+01  9.3583E+01  2.3370E+01  0.0000E+00  3.4664E+01  2.0672E+00
             1.0255E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2183.67713406241        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     2737
 NPARAMETR:  1.0242E+00  2.0248E+00  5.2014E-01  3.9038E-01  1.3095E+00  1.0158E+00  7.2379E-01  1.0000E-02  1.8784E+00  1.1141E+00
             9.6190E-01
 PARAMETER:  1.2395E-01  8.0546E-01 -5.5365E-01 -8.4064E-01  3.6964E-01  1.1570E-01 -2.2325E-01 -1.6833E+02  7.3042E-01  2.0808E-01
             6.1158E-02
 GRADIENT:   6.3712E+02  1.4047E+03 -1.3112E+00  1.2457E+02  2.9941E+01  8.8238E+01  2.2411E+01  0.0000E+00  3.7442E+01  3.9031E+00
             1.6113E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2183.75753183558        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     2900
 NPARAMETR:  1.0246E+00  2.0197E+00  5.3601E-01  3.9234E-01  1.3012E+00  1.0172E+00  7.2473E-01  1.0000E-02  1.8308E+00  1.0983E+00
             9.6146E-01
 PARAMETER:  1.2428E-01  8.0295E-01 -5.2359E-01 -8.3564E-01  3.6331E-01  1.1703E-01 -2.2195E-01 -1.6833E+02  7.0477E-01  1.9372E-01
             6.0700E-02
 GRADIENT:   2.3412E+00 -1.7741E+01  2.4016E+00 -8.1165E+00 -2.0151E+00  4.8219E-01 -6.5404E-01  0.0000E+00 -2.0477E+00 -1.4791E+00
            -5.9578E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2183.96243493392        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     3062
 NPARAMETR:  1.0233E+00  2.0048E+00  5.3324E-01  4.0582E-01  1.2973E+00  1.0162E+00  7.2665E-01  1.0000E-02  1.8349E+00  1.1039E+00
             9.6181E-01
 PARAMETER:  1.2299E-01  7.9556E-01 -5.2879E-01 -8.0186E-01  3.6027E-01  1.1608E-01 -2.1931E-01 -1.6833E+02  7.0700E-01  1.9886E-01
             6.1064E-02
 GRADIENT:  -5.3291E-01 -1.9818E+01 -2.7707E+00 -6.2868E-01  6.1427E+00  6.0669E-02 -5.6892E-02  0.0000E+00  2.6479E+00  1.3269E+00
             5.7418E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2184.10254325917        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     3229             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0268E+00  1.9935E+00  5.4042E-01  4.1252E-01  1.2918E+00  1.0184E+00  7.2956E-01  1.0000E-02  1.8114E+00  1.0971E+00
             9.6163E-01
 PARAMETER:  1.2649E-01  7.8987E-01 -5.1540E-01 -7.8548E-01  3.5606E-01  1.1824E-01 -2.1531E-01 -1.6833E+02  6.9408E-01  1.9265E-01
             6.0872E-02
 GRADIENT:   6.5572E+02  1.3410E+03 -7.6731E-02  1.2811E+02  2.7374E+01  9.0457E+01  2.1566E+01  0.0000E+00  3.5520E+01  2.7381E+00
             1.3637E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2184.16608604683        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     3391
 NPARAMETR:  1.0255E+00  1.9895E+00  5.5505E-01  4.1695E-01  1.2824E+00  1.0175E+00  7.2888E-01  1.0000E-02  1.7821E+00  1.0903E+00
             9.6117E-01
 PARAMETER:  1.2516E-01  7.8787E-01 -4.8870E-01 -7.7478E-01  3.4875E-01  1.1731E-01 -2.1624E-01 -1.6833E+02  6.7782E-01  1.8642E-01
             6.0395E-02
 GRADIENT:   4.3160E+00 -6.5883E+00  2.6424E+00 -4.9678E+00 -5.1610E+00  6.0614E-01 -4.7645E-01  0.0000E+00 -9.3963E-01 -3.7762E-01
            -5.9702E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -2184.33380660187        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3578             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0250E+00  1.9727E+00  5.5835E-01  4.2422E-01  1.2786E+00  1.0173E+00  7.3100E-01  1.0000E-02  1.7632E+00  1.0876E+00
             9.6131E-01
 PARAMETER:  1.2474E-01  7.7941E-01 -4.8277E-01 -7.5750E-01  3.4576E-01  1.1717E-01 -2.1334E-01 -1.6833E+02  6.6713E-01  1.8400E-01
             6.0538E-02
 GRADIENT:   6.4325E+02  1.3002E+03  3.1318E+00  1.2696E+02  1.9991E+01  8.9636E+01  2.0463E+01  0.0000E+00  3.2447E+01  1.9858E+00
             9.0677E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -2184.44974461979        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3765             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0254E+00  1.9589E+00  5.5797E-01  4.3273E-01  1.2769E+00  1.0175E+00  7.3418E-01  1.0000E-02  1.7583E+00  1.0863E+00
             9.6146E-01
 PARAMETER:  1.2507E-01  7.7240E-01 -4.8345E-01 -7.3765E-01  3.4445E-01  1.1731E-01 -2.0900E-01 -1.6833E+02  6.6434E-01  1.8277E-01
             6.0702E-02
 GRADIENT:   6.4509E+02  1.2682E+03 -1.3767E-01  1.3002E+02  2.7306E+01  8.9595E+01  2.0525E+01  0.0000E+00  3.4784E+01  2.7240E+00
             1.5285E+00

0ITERATION NO.:  120    OBJECTIVE VALUE:  -2184.57817360076        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3950
 NPARAMETR:  1.0248E+00  1.9528E+00  5.6761E-01  4.4195E-01  1.2693E+00  1.0171E+00  7.3309E-01  1.0000E-02  1.7332E+00  1.0789E+00
             9.6109E-01
 PARAMETER:  1.2448E-01  7.6928E-01 -4.6631E-01 -7.1655E-01  3.3844E-01  1.1697E-01 -2.1049E-01 -1.6833E+02  6.4997E-01  1.7594E-01
             6.0311E-02
 GRADIENT:   2.8939E+00 -8.7460E+00 -3.4561E-01 -1.3467E+00  1.9963E+00  4.7253E-01 -4.0213E-01  0.0000E+00  6.0756E-01  1.6426E-02
            -1.7314E-01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -2184.76276552578        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     4117
 NPARAMETR:  1.0227E+00  1.9255E+00  5.7378E-01  4.5674E-01  1.2622E+00  1.0155E+00  7.3081E-01  1.0000E-02  1.7001E+00  1.0646E+00
             9.6063E-01
 PARAMETER:  1.2246E-01  7.5519E-01 -4.5550E-01 -6.8365E-01  3.3282E-01  1.1543E-01 -2.1361E-01 -1.6833E+02  6.3070E-01  1.6262E-01
             5.9835E-02
 GRADIENT:   6.2630E+02  1.2041E+03 -8.1092E-01  1.3398E+02  3.1852E+01  8.7943E+01  1.7999E+01  0.0000E+00  3.2953E+01  6.0185E-01
             6.0404E-01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -2184.84561043190        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     4298
 NPARAMETR:  1.0246E+00  1.9247E+00  5.8109E-01  4.5959E-01  1.2545E+00  1.0167E+00  7.3872E-01  1.0000E-02  1.6930E+00  1.0690E+00
             9.6106E-01
 PARAMETER:  1.2427E-01  7.5479E-01 -4.4285E-01 -6.7742E-01  3.2676E-01  1.1658E-01 -2.0283E-01 -1.6833E+02  6.2650E-01  1.6669E-01
             6.0285E-02
 GRADIENT:   2.5025E+00 -8.3680E+00 -1.4555E-01 -1.9193E+00  3.5668E-01  3.4054E-01  2.3412E-01  0.0000E+00  6.8046E-01  2.0171E-01
            -7.9349E-04

0ITERATION NO.:  135    OBJECTIVE VALUE:  -2184.97390608368        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     4465
 NPARAMETR:  1.0234E+00  1.9019E+00  5.9610E-01  4.6846E-01  1.2428E+00  1.0160E+00  7.4057E-01  1.0000E-02  1.6626E+00  1.0638E+00
             9.6091E-01
 PARAMETER:  1.2313E-01  7.4285E-01 -4.1735E-01 -6.5831E-01  3.1739E-01  1.1589E-01 -2.0034E-01 -1.6833E+02  6.0841E-01  1.6184E-01
             6.0120E-02
 GRADIENT:   6.3102E+02  1.1594E+03  4.5899E+00  1.3016E+02  1.4214E+01  8.8239E+01  1.8348E+01  0.0000E+00  3.1246E+01  2.0185E+00
             8.9170E-01

0ITERATION NO.:  140    OBJECTIVE VALUE:  -2185.05643530057        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     4644
 NPARAMETR:  1.0250E+00  1.9010E+00  5.9262E-01  4.7630E-01  1.2448E+00  1.0170E+00  7.4068E-01  1.0000E-02  1.6611E+00  1.0606E+00
             9.6090E-01
 PARAMETER:  1.2470E-01  7.4238E-01 -4.2321E-01 -6.4171E-01  3.1900E-01  1.1684E-01 -2.0019E-01 -1.6833E+02  6.0750E-01  1.5885E-01
             6.0118E-02
 GRADIENT:   3.5194E+00 -6.7826E+00 -9.0608E-01 -1.8687E-01  2.6059E+00  4.5771E-01 -1.2488E-01  0.0000E+00  9.6224E-01  7.2749E-02
            -7.3695E-02

0ITERATION NO.:  145    OBJECTIVE VALUE:  -2185.14202369206        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4829
 NPARAMETR:  1.0246E+00  1.8912E+00  6.0113E-01  4.8428E-01  1.2400E+00  1.0169E+00  7.4221E-01  1.0000E-02  1.6473E+00  1.0572E+00
             9.6087E-01
 PARAMETER:  1.2435E-01  7.3722E-01 -4.0894E-01 -6.2509E-01  3.1513E-01  1.1672E-01 -1.9812E-01 -1.6833E+02  5.9911E-01  1.5558E-01
             6.0081E-02
 GRADIENT:   2.7930E+00 -2.6225E+00 -1.8503E-01  5.4627E-01  1.3882E+00  4.2370E-01 -1.7589E-01  0.0000E+00  9.0258E-01 -1.2447E-01
            -1.9953E-01

0ITERATION NO.:  150    OBJECTIVE VALUE:  -2185.27659298575        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     4975
 NPARAMETR:  1.0221E+00  1.8689E+00  6.0883E-01  4.9137E-01  1.2316E+00  1.0139E+00  7.4251E-01  1.0000E-02  1.6287E+00  1.0514E+00
             9.6029E-01
 PARAMETER:  1.2190E-01  7.2537E-01 -3.9622E-01 -6.1057E-01  3.0829E-01  1.1385E-01 -1.9772E-01 -1.6833E+02  5.8779E-01  1.5011E-01
             5.9483E-02
 GRADIENT:   6.2213E+02  1.0964E+03  2.1202E+00  1.3411E+02  2.0805E+01  8.6040E+01  1.7026E+01  0.0000E+00  3.2446E+01  1.5892E+00
             5.6191E-01

0ITERATION NO.:  155    OBJECTIVE VALUE:  -2185.32836874636        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     5138
 NPARAMETR:  1.0252E+00  1.8630E+00  6.1068E-01  4.9491E-01  1.2300E+00  1.0193E+00  7.4629E-01  1.0000E-02  1.6157E+00  1.0513E+00
             9.6089E-01
 PARAMETER:  1.2493E-01  7.2216E-01 -3.9318E-01 -6.0337E-01  3.0701E-01  1.1915E-01 -1.9264E-01 -1.6833E+02  5.7976E-01  1.4998E-01
             6.0110E-02
 GRADIENT:   4.2218E+00 -1.5988E+01 -4.7012E-01 -4.9161E+00  1.8887E+00  1.4421E+00  1.0417E-01  0.0000E+00  4.0798E-01  3.0066E-01
             1.6739E-01

0ITERATION NO.:  160    OBJECTIVE VALUE:  -2185.39043665163        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     5280
 NPARAMETR:  1.0227E+00  1.8497E+00  6.1204E-01  5.0383E-01  1.2270E+00  1.0135E+00  7.4577E-01  1.0000E-02  1.6089E+00  1.0485E+00
             9.6073E-01
 PARAMETER:  1.2249E-01  7.1500E-01 -3.9096E-01 -5.8551E-01  3.0461E-01  1.1339E-01 -1.9333E-01 -1.6833E+02  5.7557E-01  1.4733E-01
             5.9937E-02
 GRADIENT:  -1.1490E+00 -1.9786E+01 -3.2754E+00 -2.1755E+00  7.9163E+00 -8.9549E-01 -2.9292E-01  0.0000E+00  1.8082E+00  6.5458E-01
             2.9546E-01

0ITERATION NO.:  165    OBJECTIVE VALUE:  -2185.45308709556        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     5442
 NPARAMETR:  1.0254E+00  1.8469E+00  6.2053E-01  5.0548E-01  1.2215E+00  1.0165E+00  7.4847E-01  1.0000E-02  1.5942E+00  1.0457E+00
             9.6078E-01
 PARAMETER:  1.2509E-01  7.1351E-01 -3.7719E-01 -5.8225E-01  3.0009E-01  1.1635E-01 -1.8972E-01 -1.6833E+02  5.6639E-01  1.4469E-01
             5.9991E-02
 GRADIENT:   4.6445E+00 -1.3713E+01  4.3694E-01 -5.5371E+00 -3.0586E-01  3.3606E-01  2.0556E-02  0.0000E+00  1.0438E-01  2.1798E-01
            -2.8455E-03

0ITERATION NO.:  170    OBJECTIVE VALUE:  -2185.51902409846        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     5627
 NPARAMETR:  1.0248E+00  1.8337E+00  6.2369E-01  5.1207E-01  1.2192E+00  1.0170E+00  7.5003E-01  1.0000E-02  1.5837E+00  1.0441E+00
             9.6101E-01
 PARAMETER:  1.2451E-01  7.0633E-01 -3.7210E-01 -5.6929E-01  2.9816E-01  1.1684E-01 -1.8764E-01 -1.6833E+02  5.5979E-01  1.4315E-01
             6.0234E-02
 GRADIENT:   3.4289E+00 -1.9793E+01 -1.2765E+00 -5.6046E+00  3.8520E+00  5.4509E-01  2.0066E-02  0.0000E+00  6.6165E-01  5.1609E-01
             4.3925E-01

0ITERATION NO.:  175    OBJECTIVE VALUE:  -2185.60084303225        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     5767
 NPARAMETR:  1.0221E+00  1.8283E+00  6.2864E-01  5.2173E-01  1.2135E+00  1.0132E+00  7.5000E-01  1.0000E-02  1.5702E+00  1.0392E+00
             9.6080E-01
 PARAMETER:  1.2187E-01  7.0337E-01 -3.6419E-01 -5.5061E-01  2.9355E-01  1.1311E-01 -1.8768E-01 -1.6833E+02  5.5122E-01  1.3849E-01
             6.0012E-02
 GRADIENT:  -2.4369E+00 -8.2416E+00 -1.0280E+00 -1.1124E+00  2.9098E+00 -9.8525E-01 -4.4142E-01  0.0000E+00  6.6897E-01  1.7646E-01
             4.6543E-02

0ITERATION NO.:  180    OBJECTIVE VALUE:  -2185.65918272820        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     5953
 NPARAMETR:  1.0245E+00  1.8195E+00  6.3097E-01  5.2688E-01  1.2110E+00  1.0157E+00  7.5265E-01  1.0000E-02  1.5633E+00  1.0373E+00
             9.6087E-01
 PARAMETER:  1.2416E-01  6.9858E-01 -3.6049E-01 -5.4077E-01  2.9141E-01  1.1560E-01 -1.8415E-01 -1.6833E+02  5.4678E-01  1.3665E-01
             6.0087E-02
 GRADIENT:   2.6520E+00 -1.0553E+01 -1.9915E+00 -6.9969E-01  4.9032E+00  4.0241E-02 -4.3184E-02  0.0000E+00  1.1865E+00  3.2892E-01
             2.2232E-01

0ITERATION NO.:  185    OBJECTIVE VALUE:  -2185.73209548963        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     6136             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0250E+00  1.8022E+00  6.3665E-01  5.3375E-01  1.2057E+00  1.0168E+00  7.5518E-01  1.0000E-02  1.5486E+00  1.0336E+00
             9.6101E-01
 PARAMETER:  1.2470E-01  6.8902E-01 -3.5153E-01 -5.2783E-01  2.8702E-01  1.1667E-01 -1.8080E-01 -1.6833E+02  5.3734E-01  1.3308E-01
             6.0225E-02
 GRADIENT:   6.4046E+02  9.6681E+02 -4.1772E-01  1.3561E+02  2.5646E+01  8.8144E+01  1.5197E+01  0.0000E+00  3.1638E+01  1.7962E+00
             1.5194E+00

0ITERATION NO.:  190    OBJECTIVE VALUE:  -2185.79421109814        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     6325             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0247E+00  1.8011E+00  6.4389E-01  5.4245E-01  1.2004E+00  1.0168E+00  7.5591E-01  1.0000E-02  1.5384E+00  1.0294E+00
             9.6063E-01
 PARAMETER:  1.2440E-01  6.8840E-01 -3.4022E-01 -5.1166E-01  2.8263E-01  1.1668E-01 -1.7983E-01 -1.6833E+02  5.3076E-01  1.2895E-01
             5.9835E-02
 GRADIENT:   6.3861E+02  9.7430E+02  1.3829E+00  1.3905E+02  2.1194E+01  8.8245E+01  1.4924E+01  0.0000E+00  3.1271E+01  1.1678E+00
             7.7157E-01

0ITERATION NO.:  195    OBJECTIVE VALUE:  -2185.82925036087        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     6492
 NPARAMETR:  1.0243E+00  1.7959E+00  6.4701E-01  5.4352E-01  1.1979E+00  1.0165E+00  7.5735E-01  1.0000E-02  1.5280E+00  1.0295E+00
             9.6086E-01
 PARAMETER:  1.2406E-01  6.8551E-01 -3.3539E-01 -5.0969E-01  2.8060E-01  1.1640E-01 -1.7794E-01 -1.6833E+02  5.2399E-01  1.2911E-01
             6.0069E-02
 GRADIENT:   2.5468E+00 -4.3041E+00  1.4588E-01 -1.4220E+00  1.3810E-01  4.0073E-01 -1.1264E-01  0.0000E+00  1.6960E-01  1.0716E-01
            -1.7247E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     6492
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1028E-04 -3.7326E-02 -2.6095E-04  3.0585E-02 -3.6472E-02
 SE:             2.9899E-02  2.2752E-02  9.9759E-05  2.3461E-02  2.3157E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8105E-01  1.0088E-01  8.9030E-03  1.9235E-01  1.1526E-01

 ETASHRINKSD(%)  1.0000E-10  2.3779E+01  9.9666E+01  2.1403E+01  2.2422E+01
 ETASHRINKVR(%)  1.0000E-10  4.1904E+01  9.9999E+01  3.8226E+01  3.9817E+01
 EBVSHRINKSD(%)  3.0171E-01  2.2504E+01  9.9707E+01  2.3811E+01  1.9974E+01
 EBVSHRINKVR(%)  6.0251E-01  3.9944E+01  9.9999E+01  4.1952E+01  3.5959E+01
 RELATIVEINF(%)  9.9319E+01  4.9113E+00  1.4739E-04  4.9307E+00  1.7373E+01
 EPSSHRINKSD(%)  3.3668E+01
 EPSSHRINKVR(%)  5.6000E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2185.8292503608714     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1266.8907171561987     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.60
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2185.829       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.80E+00  6.47E-01  5.44E-01  1.20E+00  1.02E+00  7.57E-01  1.00E-02  1.53E+00  1.03E+00  9.61E-01
 


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
 #CPUT: Total CPU Time in Seconds,      242.122
Stop Time:
Sun Oct 24 00:35:46 CDT 2021
