Sun Oct 24 01:15:51 CDT 2021
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
$DATA ../../../../data/SD3/D2/dat46.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2008.88318185431        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2691E+02 -9.6460E+01 -5.8389E+01 -5.7770E+01  7.2477E+01 -2.8630E+01 -5.0728E+01 -3.9648E+00 -4.6627E+01  1.4749E+01
            -1.1438E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2024.75697411648        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  1.0083E+00  1.2585E+00  1.3366E+00  1.0145E+00  1.1914E+00  1.6055E+00  1.7092E+00  1.1069E+00  1.4092E+00  8.3707E-01
             9.8965E-01
 PARAMETER:  1.0829E-01  3.2990E-01  3.9013E-01  1.1440E-01  2.7511E-01  5.7346E-01  6.3605E-01  2.0153E-01  4.4299E-01 -7.7843E-02
             8.9596E-02
 GRADIENT:   2.3065E+01  2.5786E+01 -7.9348E+00  5.3246E+01  5.5455E+01  9.1675E+01  2.5637E+01 -1.5929E+01  4.1009E+01 -2.0757E+01
            -1.6549E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2035.68177841597        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0146E+00  9.5866E-01  2.6420E+00  1.1956E+00  1.3637E+00  1.5086E+00  1.8504E+00  2.2929E+00  1.2788E+00  9.7868E-01
             1.0229E+00
 PARAMETER:  1.1454E-01  5.7783E-02  1.0715E+00  2.7867E-01  4.1023E-01  5.1121E-01  7.1542E-01  9.2980E-01  3.4595E-01  7.8446E-02
             1.2265E-01
 GRADIENT:   3.3466E+01  8.2968E+00 -9.8578E+00  1.7803E+01  5.0014E+01  7.4944E+01  2.1409E+01  9.6482E+00  3.5693E+01 -2.6852E+01
             5.3971E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2053.45552138131        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.8202E-01  7.7730E-01  3.2303E+00  1.3057E+00  1.3486E+00  1.1820E+00  1.7589E+00  2.2908E+00  1.0167E+00  1.1865E+00
             1.0050E+00
 PARAMETER:  8.1861E-02 -1.5193E-01  1.2726E+00  3.6677E-01  3.9909E-01  2.6725E-01  6.6467E-01  9.2890E-01  1.1653E-01  2.7100E-01
             1.0498E-01
 GRADIENT:   2.4583E+00  7.9711E+00  2.4088E+00  3.1875E+00 -6.4346E+00 -7.6798E+00  1.7505E+00  1.5928E+00 -2.4379E+00  1.0931E+00
            -2.9937E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2055.44896550956        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.7615E-01  4.9063E-01  3.5828E+00  1.5066E+00  1.3235E+00  1.1800E+00  1.2705E+00  2.3698E+00  1.0512E+00  1.1852E+00
             1.0038E+00
 PARAMETER:  7.5857E-02 -6.1207E-01  1.3761E+00  5.0984E-01  3.8027E-01  2.6549E-01  3.3938E-01  9.6282E-01  1.4990E-01  2.6987E-01
             1.0379E-01
 GRADIENT:  -4.1543E+00  6.8677E+00 -7.3052E-01  2.0950E+01 -1.7048E+00 -7.1260E+00  1.6440E+00  1.1303E+00  5.5630E+00  7.2442E-01
            -2.6725E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2056.27354545626        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  9.7664E-01  3.0451E-01  3.6326E+00  1.6211E+00  1.2887E+00  1.1978E+00  7.2787E-01  2.3547E+00  9.8514E-01  1.1588E+00
             1.0059E+00
 PARAMETER:  7.6367E-02 -1.0891E+00  1.3900E+00  5.8311E-01  3.5363E-01  2.8047E-01 -2.1763E-01  9.5642E-01  8.5026E-02  2.4736E-01
             1.0591E-01
 GRADIENT:  -7.9231E-01  2.1576E+00 -1.3152E+00  1.1295E+01  8.6947E-01 -3.4026E-01  3.3748E-01  9.8173E-03  5.3704E-01  3.1448E-01
             2.4798E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2056.50010755630        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  9.7774E-01  2.8582E-01  3.7038E+00  1.6183E+00  1.2880E+00  1.2033E+00  2.8390E-01  2.3745E+00  9.7805E-01  1.1566E+00
             1.0056E+00
 PARAMETER:  7.7489E-02 -1.1524E+00  1.4094E+00  5.8140E-01  3.5309E-01  2.8508E-01 -1.1591E+00  9.6477E-01  7.7810E-02  2.4547E-01
             1.0557E-01
 GRADIENT:   1.1395E+00 -6.3616E-01 -1.2319E-02 -9.8522E+00 -1.7848E-01  1.6175E+00  5.9174E-02 -2.8804E-01 -5.7007E-01  9.2778E-02
             4.9961E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2056.57842401665        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1230
 NPARAMETR:  9.7602E-01  3.0286E-01  3.6993E+00  1.6093E+00  1.2899E+00  1.1971E+00  2.6648E-02  2.3811E+00  9.8737E-01  1.1552E+00
             1.0057E+00
 PARAMETER:  7.5726E-02 -1.0945E+00  1.4081E+00  5.7577E-01  3.5458E-01  2.7994E-01 -3.5250E+00  9.6757E-01  8.7290E-02  2.4429E-01
             1.0570E-01
 GRADIENT:  -1.8616E+00 -1.5054E-01 -2.4086E-02 -5.6215E+00 -1.1812E+00 -5.2257E-01  8.4034E-04  1.8846E-01  7.7517E-02 -2.4352E-01
            -1.3468E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2056.61329042759        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1369
 NPARAMETR:  9.7591E-01  3.0926E-01  3.6953E+00  1.6017E+00  1.2911E+00  1.2038E+00  1.0000E-02  2.3788E+00  9.8679E-01  1.1579E+00
             1.0058E+00
 PARAMETER:  7.5611E-02 -1.0736E+00  1.4071E+00  5.7110E-01  3.5547E-01  2.8551E-01 -5.3816E+00  9.6660E-01  8.6701E-02  2.4658E-01
             1.0577E-01
 GRADIENT:  -2.0677E+00 -5.1963E-01  7.6665E-02 -9.6891E+00 -1.2147E+00  1.7814E+00  0.0000E+00  1.6489E-01 -1.1885E+00 -2.9825E-02
             2.6092E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2056.64686413328        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1555
 NPARAMETR:  9.7496E-01  3.2415E-01  3.6986E+00  1.5952E+00  1.2943E+00  1.2042E+00  1.0000E-02  2.3786E+00  9.9161E-01  1.1602E+00
             1.0059E+00
 PARAMETER:  7.4637E-02 -1.0266E+00  1.4079E+00  5.6703E-01  3.5795E-01  2.8585E-01 -5.3598E+00  9.6651E-01  9.1570E-02  2.4860E-01
             1.0591E-01
 GRADIENT:  -3.8010E+00  3.5504E-01  2.0458E-01 -4.3695E+00 -1.6057E+00  1.8968E+00  0.0000E+00  1.4086E-02 -1.7613E+00 -5.4819E-02
            -6.1764E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2056.67814420160        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1743             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7787E-01  3.3497E-01  3.6993E+00  1.5889E+00  1.2963E+00  1.2080E+00  1.0000E-02  2.3812E+00  9.9464E-01  1.1624E+00
             1.0059E+00
 PARAMETER:  7.7620E-02 -9.9370E-01  1.4081E+00  5.6304E-01  3.5955E-01  2.8899E-01 -5.3598E+00  9.6763E-01  9.4628E-02  2.5052E-01
             1.0586E-01
 GRADIENT:   3.9332E+02  4.1492E+01  7.4549E+00  8.6402E+02  1.3663E+01  2.0068E+02  0.0000E+00  7.9698E+00  1.1635E+01  1.8081E+00
             8.3582E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2056.69580565042        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1903
 NPARAMETR:  9.7948E-01  3.3425E-01  3.6889E+00  1.5849E+00  1.2982E+00  1.2072E+00  1.0000E-02  2.3810E+00  1.0029E+00  1.1625E+00
             1.0060E+00
 PARAMETER:  7.9262E-02 -9.9588E-01  1.4053E+00  5.6049E-01  3.6096E-01  2.8833E-01 -5.3598E+00  9.6752E-01  1.0294E-01  2.5061E-01
             1.0602E-01
 GRADIENT:   3.1481E+00 -7.7443E-01 -2.6131E-01 -8.9396E+00  1.4849E-01  2.8698E+00  0.0000E+00  2.6717E-01  9.5745E-01 -1.3713E-01
             7.1549E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2056.71966919892        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2090             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7852E-01  3.4217E-01  3.6884E+00  1.5790E+00  1.3001E+00  1.2048E+00  1.0000E-02  2.3814E+00  1.0059E+00  1.1642E+00
             1.0061E+00
 PARAMETER:  7.8281E-02 -9.7245E-01  1.4052E+00  5.5678E-01  3.6247E-01  2.8628E-01 -5.3598E+00  9.6768E-01  1.0584E-01  2.5207E-01
             1.0607E-01
 GRADIENT:   3.9417E+02  4.0797E+01  6.9648E+00  8.4177E+02  1.6095E+01  1.9739E+02  0.0000E+00  8.0647E+00  1.5550E+01  1.6656E+00
             1.1414E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2056.73235538372        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2280
 NPARAMETR:  9.7852E-01  3.4695E-01  3.6881E+00  1.5756E+00  1.3009E+00  1.2048E+00  1.0000E-02  2.3804E+00  1.0071E+00  1.1650E+00
             1.0061E+00
 PARAMETER:  7.8290E-02 -9.5858E-01  1.4051E+00  5.5463E-01  3.6309E-01  2.8631E-01 -5.3598E+00  9.6729E-01  1.0709E-01  2.5270E-01
             1.0612E-01
 GRADIENT:   1.4950E+00 -7.4247E-01 -1.4079E-01 -9.3820E+00  1.1399E-01  2.0432E+00  0.0000E+00  1.9577E-01  4.5603E-01 -1.0713E-01
             8.6699E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2056.74031172275        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2465
 NPARAMETR:  9.7859E-01  3.5112E-01  3.6862E+00  1.5706E+00  1.3023E+00  1.2048E+00  1.0000E-02  2.3805E+00  1.0094E+00  1.1660E+00
             1.0062E+00
 PARAMETER:  7.8353E-02 -9.4663E-01  1.4046E+00  5.5148E-01  3.6413E-01  2.8634E-01 -5.3598E+00  9.6730E-01  1.0934E-01  2.5360E-01
             1.0620E-01
 GRADIENT:   1.5481E+00 -1.2104E+00 -1.7034E-01 -1.2116E+01  5.6021E-01  2.0530E+00  0.0000E+00  2.1037E-01  6.2139E-01 -1.1279E-01
             1.7680E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2056.75622135386        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2652
 NPARAMETR:  9.7841E-01  3.5781E-01  3.6840E+00  1.5697E+00  1.3031E+00  1.2035E+00  1.0000E-02  2.3793E+00  1.0117E+00  1.1667E+00
             1.0062E+00
 PARAMETER:  7.8177E-02 -9.2774E-01  1.4040E+00  5.5086E-01  3.6472E-01  2.8523E-01 -5.3598E+00  9.6679E-01  1.1165E-01  2.5419E-01
             1.0618E-01
 GRADIENT:   1.1652E+00 -3.4793E-01 -1.2890E-01 -7.0744E+00  4.9748E-02  1.5864E+00  0.0000E+00  1.4543E-01  4.0047E-01 -9.6252E-02
             4.8292E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2056.76666238014        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2846             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7863E-01  3.6837E-01  3.6834E+00  1.5654E+00  1.3035E+00  1.2049E+00  1.0000E-02  2.3786E+00  1.0117E+00  1.1681E+00
             1.0061E+00
 PARAMETER:  7.8403E-02 -8.9867E-01  1.4038E+00  5.4814E-01  3.6505E-01  2.8636E-01 -5.3598E+00  9.6651E-01  1.1162E-01  2.5537E-01
             1.0612E-01
 GRADIENT:   3.9403E+02  4.5321E+01  7.2352E+00  8.1502E+02  1.4475E+01  1.9760E+02  0.0000E+00  7.8315E+00  1.4173E+01  1.8681E+00
             8.8884E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2056.77152119768        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3026
 NPARAMETR:  9.7849E-01  3.6715E-01  3.6794E+00  1.5651E+00  1.3048E+00  1.2036E+00  1.0000E-02  2.3785E+00  1.0153E+00  1.1682E+00
             1.0063E+00
 PARAMETER:  7.8260E-02 -9.0198E-01  1.4028E+00  5.4794E-01  3.6604E-01  2.8529E-01 -5.3598E+00  9.6648E-01  1.1521E-01  2.5545E-01
             1.0624E-01
 GRADIENT:   1.1563E+00  1.5456E-01 -1.5867E-01 -4.3856E+00 -5.2709E-02  1.5939E+00  0.0000E+00  1.3584E-01  1.9469E-01 -7.2537E-02
             1.0622E-02

0ITERATION NO.:   86    OBJECTIVE VALUE:  -2056.77152119768        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     3054
 NPARAMETR:  9.7899E-01  3.6694E-01  3.6855E+00  1.5565E+00  1.3047E+00  1.2069E+00  1.0000E-02  2.3793E+00  1.0142E+00  1.1697E+00
             1.0062E+00
 PARAMETER:  7.8260E-02 -9.0198E-01  1.4028E+00  5.4794E-01  3.6604E-01  2.8529E-01 -5.3598E+00  9.6648E-01  1.1521E-01  2.5545E-01
             1.0624E-01
 GRADIENT:  -3.1246E-01  1.6759E-02 -7.2923E-02  4.6378E+00  2.9478E-02 -4.3440E-01  0.0000E+00 -1.7236E-02  1.6796E-01 -7.9824E-02
             1.0754E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3054
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9470E-04 -1.6165E-04 -4.4612E-02 -7.5173E-03 -5.2749E-02
 SE:             2.9943E-02  5.3544E-05  1.7308E-02  2.9597E-02  2.0775E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9481E-01  2.5364E-03  9.9525E-03  7.9950E-01  1.1116E-02

 ETASHRINKSD(%)  1.0000E-10  9.9821E+01  4.2014E+01  8.4615E-01  3.0400E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.6377E+01  1.6851E+00  5.1558E+01
 EBVSHRINKSD(%)  2.4349E-01  9.9831E+01  5.1039E+01  1.0927E+00  2.4357E+01
 EBVSHRINKVR(%)  4.8638E-01  1.0000E+02  7.6028E+01  2.1734E+00  4.2781E+01
 RELATIVEINF(%)  9.9117E+01  2.6068E-05  1.0680E+01  9.3253E+00  2.3460E+01
 EPSSHRINKSD(%)  3.3928E+01
 EPSSHRINKVR(%)  5.6345E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2056.7715211976774     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1137.8329879930047     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2056.772       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.78E-01  3.67E-01  3.68E+00  1.57E+00  1.30E+00  1.20E+00  1.00E-02  2.38E+00  1.02E+00  1.17E+00  1.01E+00
 


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
 #CPUT: Total CPU Time in Seconds,      110.773
Stop Time:
Sun Oct 24 01:16:10 CDT 2021
