Sun Oct 24 02:09:00 CDT 2021
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
$DATA ../../../../data/SD4/A2/dat6.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1240.75156442972        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1232E+02 -4.6362E+01 -5.0427E+01  2.7323E+01  1.8702E+02  4.3936E+01 -1.7628E+01 -4.6191E+00 -1.5229E+01 -4.5980E+01
            -7.3419E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1456.36683713616        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0292E+00  1.0053E+00  1.0278E+00  1.0646E+00  8.8536E-01  1.0285E+00  1.0173E+00  9.8558E-01  1.0477E+00  7.4079E-01
             1.8171E+00
 PARAMETER:  1.2876E-01  1.0530E-01  1.2741E-01  1.6256E-01 -2.1756E-02  1.2809E-01  1.1715E-01  8.5475E-02  1.4658E-01 -2.0004E-01
             6.9723E-01
 GRADIENT:   2.0637E+02  2.1676E+01  5.0238E+00  5.7220E+01  2.8608E+01  3.4422E+01 -2.3420E+00 -4.3444E+00  7.0718E+00 -1.7482E+01
            -8.8592E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1458.59774068771        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0193E+00  9.7088E-01  4.9925E-01  1.0419E+00  5.9145E-01  1.0243E+00  1.1312E+00  8.7367E-01  8.9910E-01  4.3900E-01
             1.8135E+00
 PARAMETER:  1.1912E-01  7.0450E-02 -5.9465E-01  1.4109E-01 -4.2518E-01  1.2403E-01  2.2324E-01 -3.5051E-02 -6.3633E-03 -7.2327E-01
             6.9527E-01
 GRADIENT:   1.5457E+02  4.9334E+01  1.3689E+01  9.2286E+01  5.6085E+00  3.1704E+01 -5.9395E+00  3.9185E+00  2.5207E-01 -2.0575E+00
            -7.2406E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1465.63093833371        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      335
 NPARAMETR:  1.0160E+00  8.8380E-01  6.1082E-01  1.0624E+00  6.3961E-01  9.7862E-01  1.2116E+00  6.7120E-01  8.5771E-01  5.5118E-01
             1.9675E+00
 PARAMETER:  1.1589E-01 -2.3528E-02 -3.9295E-01  1.6055E-01 -3.4690E-01  7.8391E-02  2.9195E-01 -2.9869E-01 -5.3486E-02 -4.9570E-01
             7.7679E-01
 GRADIENT:   1.0892E+01 -1.2621E+01  2.9094E+00 -1.9247E+01  4.3964E+00  4.6575E+00 -7.3061E+00  9.1386E-01 -4.9990E+00 -2.7520E+00
            -3.5548E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1467.56187875099        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      510
 NPARAMETR:  1.0129E+00  8.0878E-01  6.9001E-01  1.1347E+00  6.5359E-01  9.6045E-01  1.3443E+00  3.9645E-01  8.4797E-01  6.7150E-01
             2.1086E+00
 PARAMETER:  1.1283E-01 -1.1222E-01 -2.7105E-01  2.2637E-01 -3.2528E-01  5.9650E-02  3.9584E-01 -8.2521E-01 -6.4905E-02 -2.9824E-01
             8.4604E-01
 GRADIENT:   2.2006E+00  3.8045E+00  4.0256E+00 -5.9604E-01 -7.5933E+00 -5.3934E-02 -3.6757E-01  5.2072E-02 -1.0181E+00  7.2397E-01
            -2.5044E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1467.65894117940        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      685
 NPARAMETR:  1.0094E+00  6.8010E-01  7.6225E-01  1.2177E+00  6.5132E-01  9.6021E-01  1.5226E+00  3.3990E-01  8.2198E-01  7.1164E-01
             2.1201E+00
 PARAMETER:  1.0940E-01 -2.8552E-01 -1.7148E-01  2.9696E-01 -3.2876E-01  5.9399E-02  5.2039E-01 -9.7909E-01 -9.6044E-02 -2.4018E-01
             8.5145E-01
 GRADIENT:  -1.9606E+00  1.7989E+00  2.2818E+00 -3.8977E-01 -3.3732E+00  8.8195E-01  2.5871E-01 -4.8130E-02 -7.9566E-01 -1.7774E-01
             8.2980E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1467.73415756689        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  1.0088E+00  4.9879E-01  8.6190E-01  1.3387E+00  6.4752E-01  9.5071E-01  1.8050E+00  2.8724E-01  8.0420E-01  7.9212E-01
             2.1144E+00
 PARAMETER:  1.0875E-01 -5.9557E-01 -4.8611E-02  3.9173E-01 -3.3460E-01  4.9450E-02  6.9054E-01 -1.1474E+00 -1.1790E-01 -1.3305E-01
             8.4879E-01
 GRADIENT:   3.9951E+00  3.4994E+00  2.3473E+00  6.6562E+00 -5.4622E+00 -1.4720E+00  8.0433E-01  2.5500E-02  2.5752E-01  2.9087E-01
            -1.8092E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1468.16447626425        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1044
 NPARAMETR:  1.0013E+00  2.6011E-01  1.0409E+00  1.5040E+00  6.6723E-01  9.4712E-01  2.3558E+00  2.9316E-01  7.7696E-01  8.6336E-01
             2.1474E+00
 PARAMETER:  1.0131E-01 -1.2467E+00  1.4007E-01  5.0814E-01 -3.0462E-01  4.5673E-02  9.5690E-01 -1.1270E+00 -1.5236E-01 -4.6921E-02
             8.6427E-01
 GRADIENT:  -1.6097E+00  3.9524E+00  3.0402E+00  2.0764E+01 -6.8625E+00 -3.0390E-01  1.0997E+00  1.4397E-02  1.9531E-02 -6.1321E-01
             4.1145E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1468.65576366191        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1222
 NPARAMETR:  9.9668E-01  1.4496E-01  1.0539E+00  1.5626E+00  6.4815E-01  9.5187E-01  2.7089E+00  3.7623E-01  7.6018E-01  8.6331E-01
             2.1413E+00
 PARAMETER:  9.6674E-02 -1.8313E+00  1.5247E-01  5.4637E-01 -3.3363E-01  5.0675E-02  1.0965E+00 -8.7756E-01 -1.7420E-01 -4.6985E-02
             8.6142E-01
 GRADIENT:  -6.6782E+00  1.2088E+00  2.7730E+00  8.1678E+00 -5.1108E+00  2.0606E+00 -4.2615E-01  5.1700E-02 -1.5416E+00 -8.3377E-01
             1.8980E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1468.83847360569        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1398
 NPARAMETR:  9.9781E-01  6.4079E-02  1.0635E+00  1.6095E+00  6.3600E-01  9.4505E-01  3.2518E+00  5.5921E-01  7.5188E-01  8.6709E-01
             2.1188E+00
 PARAMETER:  9.7809E-02 -2.6476E+00  1.6158E-01  5.7594E-01 -3.5256E-01  4.3484E-02  1.2792E+00 -4.8123E-01 -1.8517E-01 -4.2608E-02
             8.5086E-01
 GRADIENT:   2.7205E-01  5.9879E-01  7.5995E-01  9.8912E+00 -4.0463E+00 -6.4247E-01 -1.2261E-01  3.6005E-01 -7.5323E-01  1.3534E+00
            -3.8621E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1468.90483399858        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1575
 NPARAMETR:  9.9559E-01  3.1993E-02  1.1356E+00  1.6321E+00  6.5607E-01  9.4822E-01  3.9319E+00  7.3705E-01  7.4529E-01  8.4017E-01
             2.1218E+00
 PARAMETER:  9.5581E-02 -3.3423E+00  2.2714E-01  5.8987E-01 -3.2148E-01  4.6828E-02  1.4691E+00 -2.0509E-01 -1.9399E-01 -7.4155E-02
             8.5224E-01
 GRADIENT:  -2.7196E+00  1.6535E-01 -2.6737E-03  5.7292E-01 -2.3666E-01  9.5550E-01 -3.4237E-02  2.8995E-01 -3.6034E-02  3.7447E-01
             8.3356E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1468.91819269645        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1752
 NPARAMETR:  9.9604E-01  1.7357E-02  1.1645E+00  1.6441E+00  6.6266E-01  9.4637E-01  4.7376E+00  8.1765E-01  7.4250E-01  8.2235E-01
             2.1171E+00
 PARAMETER:  9.6028E-02 -3.9538E+00  2.5226E-01  5.9722E-01 -3.1150E-01  4.4876E-02  1.6555E+00 -1.0133E-01 -1.9773E-01 -9.5588E-02
             8.5004E-01
 GRADIENT:  -7.3783E-01  9.4266E-02 -4.1942E-02  1.6877E+00  1.9731E-01  2.8246E-01 -1.3618E-02  1.8445E-01  3.1482E-01 -3.2370E-01
            -9.9936E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1469.04321009305        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1931
 NPARAMETR:  9.9576E-01  1.0000E-02  1.0838E+00  1.6373E+00  6.3416E-01  9.4640E-01  6.5770E+00  6.0556E-01  7.4345E-01  8.5656E-01
             2.1216E+00
 PARAMETER:  9.5749E-02 -4.9463E+00  1.8043E-01  5.9303E-01 -3.5545E-01  4.4911E-02  1.9836E+00 -4.0160E-01 -1.9646E-01 -5.4830E-02
             8.5218E-01
 GRADIENT:  -1.3051E+00  0.0000E+00 -5.2692E-01 -7.9643E-01  2.9716E-01  3.4504E-01 -9.7106E-03  3.2638E-01  4.7882E-01  5.2046E-01
             5.9466E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1469.93160569739        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2112
 NPARAMETR:  9.9694E-01  1.0000E-02  9.7096E-01  1.6146E+00  5.9623E-01  9.4642E-01  1.1804E+01  2.5799E-01  7.3622E-01  8.3872E-01
             2.1343E+00
 PARAMETER:  9.6936E-02 -6.2874E+00  7.0531E-02  5.7910E-01 -4.1712E-01  4.4933E-02  2.5684E+00 -1.2548E+00 -2.0623E-01 -7.5882E-02
             8.5813E-01
 GRADIENT:   8.9639E-01  0.0000E+00 -3.6350E+00  1.3673E-01  6.5982E+00  5.2208E-01  4.1809E-01  1.2630E-01 -8.6187E-02 -6.4700E-01
             1.6165E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1469.94826020777        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2287
 NPARAMETR:  9.9631E-01  1.0000E-02  9.3495E-01  1.6070E+00  5.7965E-01  9.4633E-01  1.1484E+01  2.2020E-01  7.3893E-01  8.4291E-01
             2.1199E+00
 PARAMETER:  9.6300E-02 -6.1255E+00  3.2742E-02  5.7438E-01 -4.4533E-01  4.4834E-02  2.5409E+00 -1.4132E+00 -2.0255E-01 -7.0891E-02
             8.5136E-01
 GRADIENT:  -5.8894E-01  0.0000E+00 -2.9131E-01 -1.4728E+00  4.3490E-01  1.7283E-01  1.0386E-01  1.1818E-01  2.5028E-01  2.3001E-01
             2.7965E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1469.95948949466        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2465
 NPARAMETR:  9.9668E-01  1.0000E-02  9.2204E-01  1.6047E+00  5.7370E-01  9.4592E-01  1.1311E+01  1.2152E-01  7.3858E-01  8.4655E-01
             2.1189E+00
 PARAMETER:  9.6677E-02 -5.6630E+00  1.8839E-02  5.7291E-01 -4.5565E-01  4.4400E-02  2.5258E+00 -2.0077E+00 -2.0303E-01 -6.6588E-02
             8.5089E-01
 GRADIENT:   1.7962E-01  0.0000E+00  1.0245E+00 -1.1362E+00 -1.8756E+00 -1.4103E-02 -5.1014E-02  3.8413E-02 -9.8901E-02  2.7667E-01
            -3.3243E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1470.05340201561        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     2629
 NPARAMETR:  9.9605E-01  1.0000E-02  9.3286E-01  1.6054E+00  5.7886E-01  9.4543E-01  1.1766E+01  1.9338E-02  7.3827E-01  8.4913E-01
             2.1222E+00
 PARAMETER:  9.6047E-02 -4.6141E+00  3.0498E-02  5.7337E-01 -4.4669E-01  4.3883E-02  2.5652E+00 -3.8457E+00 -2.0344E-01 -6.3547E-02
             8.5247E-01
 GRADIENT:   8.4278E+01  0.0000E+00  1.1517E+00  2.1134E+02  1.3558E+01  7.6663E+00  2.0441E+01  2.1983E-03  4.5868E+00  2.2146E-01
             6.8430E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1470.05868076685        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:     2783
 NPARAMETR:  9.9570E-01  1.0000E-02  9.3248E-01  1.6041E+00  5.7870E-01  9.4513E-01  1.1815E+01  2.2025E-02  7.3786E-01  8.4885E-01
             2.1214E+00
 PARAMETER:  9.5654E-02 -4.6141E+00  3.0109E-02  5.7248E-01 -4.4696E-01  4.3552E-02  2.5697E+00 -3.6788E+00 -2.0403E-01 -6.3888E-02
             8.5209E-01
 GRADIENT:  -1.6592E+00  0.0000E+00  3.5788E-01 -7.6146E+01  6.6306E-01 -1.0869E-01  1.5087E+01  1.2218E-03 -1.6573E-01 -4.8039E-02
            -1.2706E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2783
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.0763E-04  5.6698E-03  1.8825E-05 -1.0164E-02 -1.8407E-02
 SE:             2.9317E-02  4.6846E-03  3.8569E-04  2.7271E-02  2.1895E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7802E-01  2.2616E-01  9.6107E-01  7.0936E-01  4.0052E-01

 ETASHRINKSD(%)  1.7858E+00  8.4306E+01  9.8708E+01  8.6382E+00  2.6648E+01
 ETASHRINKVR(%)  3.5397E+00  9.7537E+01  9.9983E+01  1.6530E+01  4.6195E+01
 EBVSHRINKSD(%)  1.8685E+00  8.8162E+01  9.8674E+01  8.5437E+00  2.6131E+01
 EBVSHRINKVR(%)  3.7020E+00  9.8599E+01  9.9982E+01  1.6358E+01  4.5433E+01
 RELATIVEINF(%)  9.5980E+01  9.9412E-01  8.7069E-04  3.6906E+01  2.7106E+00
 EPSSHRINKSD(%)  3.4386E+01
 EPSSHRINKVR(%)  5.6948E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1470.0586807668490     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -734.90785420311079     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1470.059       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.00E-02  9.32E-01  1.60E+00  5.79E-01  9.45E-01  1.18E+01  2.28E-02  7.38E-01  8.49E-01  2.12E+00
 


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
 #CPUT: Total CPU Time in Seconds,       79.199
Stop Time:
Sun Oct 24 02:09:15 CDT 2021
