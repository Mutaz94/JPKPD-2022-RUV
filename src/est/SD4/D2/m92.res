Sun Oct 24 04:46:02 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat92.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1524.90827953888        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1411E+02 -6.6286E+01 -3.9847E+01 -3.3952E+01  3.6633E+01 -5.9712E+00 -1.0045E+02 -2.2315E-01 -1.0671E+02 -9.4232E+00
             1.0821E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1566.70593130031        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      110
 NPARAMETR:  9.5937E-01  1.0394E+00  1.0737E+00  8.5056E-01  1.1560E+00  1.0058E+00  1.5739E+00  9.9999E-01  1.2827E+00  9.9018E-01
             9.3533E-01
 PARAMETER:  5.8526E-02  1.3868E-01  1.7113E-01 -6.1858E-02  2.4495E-01  1.0576E-01  5.5358E-01  9.9994E-02  3.4894E-01  9.0134E-02
             3.3148E-02
 GRADIENT:  -8.0451E+01 -1.4944E+02 -2.0636E+01 -9.7796E+01  4.4494E+01 -5.6210E+01 -6.6811E+01  6.1923E-01  7.3744E+00  1.4215E+00
            -6.3105E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1582.20608962633        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      288
 NPARAMETR:  1.0221E+00  1.3793E+00  1.0309E+00  9.6102E-01  1.0594E+00  1.1903E+00  1.8800E+00  1.1486E+00  1.1542E+00  9.2311E-01
             1.0069E+00
 PARAMETER:  1.2184E-01  4.2155E-01  1.3043E-01  6.0237E-02  1.5766E-01  2.7419E-01  7.3127E-01  2.3852E-01  2.4340E-01  1.9996E-02
             1.0686E-01
 GRADIENT:   3.9147E+01  7.8446E+01  9.9138E+00  6.7135E+01 -6.6459E+01  1.7901E+01  3.1351E-01  2.8837E+00  3.0011E+00  4.2234E+00
             2.2808E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1592.50123424594        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      466
 NPARAMETR:  9.9823E-01  1.2585E+00  1.2304E+00  9.0703E-01  1.1971E+00  1.1297E+00  1.8132E+00  1.3684E+00  1.1444E+00  1.0063E+00
             9.5172E-01
 PARAMETER:  9.8230E-02  3.2993E-01  3.0731E-01  2.4157E-03  2.7990E-01  2.2191E-01  6.9510E-01  4.1363E-01  2.3492E-01  1.0631E-01
             5.0511E-02
 GRADIENT:   4.1060E+00  2.4076E+00  1.1701E+00  2.8336E+00 -1.3617E+00 -9.0588E-01  1.0519E+00 -5.2547E-01  2.0229E-03  6.2032E-02
            -3.1454E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1592.56017589276        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      644
 NPARAMETR:  9.9337E-01  1.4398E+00  1.0371E+00  7.9228E-01  1.2168E+00  1.1319E+00  1.6435E+00  1.2659E+00  1.2442E+00  1.0038E+00
             9.5022E-01
 PARAMETER:  9.3352E-02  4.6451E-01  1.3640E-01 -1.3284E-01  2.9620E-01  2.2392E-01  5.9683E-01  3.3577E-01  3.1850E-01  1.0380E-01
             4.8938E-02
 GRADIENT:  -6.3591E+00  3.6104E+00 -7.4939E-01  5.5881E+00  2.2514E+00 -2.5229E-01  8.2616E-01  3.0571E-01  4.5280E-01 -1.2860E-01
            -5.8423E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1592.61163017286        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      827
 NPARAMETR:  9.9686E-01  1.4372E+00  1.0309E+00  7.8691E-01  1.2147E+00  1.1322E+00  1.6382E+00  1.2483E+00  1.2422E+00  1.0033E+00
             9.5165E-01
 PARAMETER:  9.6857E-02  4.6267E-01  1.3039E-01 -1.3964E-01  2.9450E-01  2.2414E-01  5.9361E-01  3.2179E-01  3.1689E-01  1.0331E-01
             5.0441E-02
 GRADIENT:  -3.2760E-01  2.4841E-02  5.2628E-01  1.0435E+00  2.8941E-01 -1.2441E-01  2.5017E-01  6.2852E-02 -8.9977E-02  7.9601E-02
             2.2337E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1592.63902914448        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1003
 NPARAMETR:  9.9647E-01  1.5215E+00  8.7342E-01  7.2906E-01  1.1879E+00  1.1302E+00  1.5658E+00  9.8559E-01  1.2805E+00  9.7874E-01
             9.5148E-01
 PARAMETER:  9.6464E-02  5.1973E-01 -3.5342E-02 -2.1600E-01  2.7218E-01  2.2243E-01  5.4841E-01  8.5490E-02  3.4722E-01  7.8508E-02
             5.0261E-02
 GRADIENT:  -2.2719E+00  2.0735E+00  3.2423E+00 -6.1104E-01 -8.0512E+00 -9.3068E-01 -1.6714E+00 -2.0315E-01 -1.2307E+00  3.6250E-01
             1.1247E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1592.70932101662        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1184             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9899E-01  1.5648E+00  7.7640E-01  6.9339E-01  1.1849E+00  1.1344E+00  1.5410E+00  8.1848E-01  1.3211E+00  9.6043E-01
             9.5081E-01
 PARAMETER:  9.8993E-02  5.4774E-01 -1.5309E-01 -2.6616E-01  2.6966E-01  2.2611E-01  5.3242E-01 -1.0030E-01  3.7846E-01  5.9630E-02
             4.9559E-02
 GRADIENT:   4.5622E+02  4.4470E+02 -2.4728E+00  1.1087E+02  2.3042E+01  1.2541E+02  9.5931E+01  6.3899E-01  2.4266E+01 -1.3233E-01
             1.1619E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1592.75140091735        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1360
 NPARAMETR:  9.9677E-01  1.5752E+00  7.7684E-01  6.8680E-01  1.1864E+00  1.1311E+00  1.5307E+00  7.4516E-01  1.3305E+00  9.7318E-01
             9.5063E-01
 PARAMETER:  9.6767E-02  5.5441E-01 -1.5253E-01 -2.7572E-01  2.7090E-01  2.2321E-01  5.2576E-01 -1.9416E-01  3.8553E-01  7.2817E-02
             4.9366E-02
 GRADIENT:  -2.5445E+00 -2.7967E+00  7.5378E-01  5.5345E-01 -6.6407E-01 -7.0478E-01 -5.1098E-01 -7.6451E-02 -2.6093E-01  2.0353E-01
            -3.0874E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1592.76078445322        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1519
 NPARAMETR:  9.9996E-01  1.5780E+00  7.6654E-01  6.8120E-01  1.1883E+00  1.1362E+00  1.5343E+00  7.5433E-01  1.3397E+00  9.7031E-01
             9.5063E-01
 PARAMETER:  9.9957E-02  5.5614E-01 -1.6587E-01 -2.8390E-01  2.7255E-01  2.2768E-01  5.2805E-01 -1.8192E-01  3.9241E-01  6.9864E-02
             4.9374E-02
 GRADIENT:   2.8578E+00 -6.2755E+00 -1.6128E+00  6.5650E-01  2.8613E+00  1.0829E+00  1.2379E+00  2.6086E-01  6.8405E-01  9.8558E-02
             2.8900E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1592.77392139230        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1698
 NPARAMETR:  9.9929E-01  1.5827E+00  7.6623E-01  6.8059E-01  1.1873E+00  1.1352E+00  1.5321E+00  7.3151E-01  1.3392E+00  9.7095E-01
             9.5054E-01
 PARAMETER:  9.9294E-02  5.5915E-01 -1.6627E-01 -2.8480E-01  2.7172E-01  2.2679E-01  5.2666E-01 -2.1264E-01  3.9209E-01  7.0518E-02
             4.9277E-02
 GRADIENT:   1.7270E+00 -3.7551E+00 -5.3504E-02  4.5953E-01  6.0783E-01  7.2628E-01  8.9449E-01  5.7927E-02  1.8229E-01  6.8303E-02
             6.0089E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1592.77788486178        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1876
 NPARAMETR:  9.9867E-01  1.5865E+00  7.6198E-01  6.7848E-01  1.1868E+00  1.1340E+00  1.5281E+00  7.2021E-01  1.3403E+00  9.7052E-01
             9.5055E-01
 PARAMETER:  9.8665E-02  5.6156E-01 -1.7183E-01 -2.8790E-01  2.7129E-01  2.2579E-01  5.2403E-01 -2.2821E-01  3.9289E-01  7.0074E-02
             4.9285E-02
 GRADIENT:   6.0813E-01 -3.3987E+00  1.7177E-01  5.0499E-01  1.5224E-01  3.2766E-01  5.6742E-01  3.9472E-02  2.9607E-02  6.8989E-02
             4.0180E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1592.78083714230        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2053
 NPARAMETR:  9.9933E-01  1.5874E+00  7.5793E-01  6.7693E-01  1.1869E+00  1.1352E+00  1.5287E+00  7.1311E-01  1.3433E+00  9.7023E-01
             9.5058E-01
 PARAMETER:  9.9333E-02  5.6212E-01 -1.7716E-01 -2.9019E-01  2.7133E-01  2.2680E-01  5.2444E-01 -2.3813E-01  3.9510E-01  6.9779E-02
             4.9315E-02
 GRADIENT:   1.7193E+00 -4.2633E+00 -4.3360E-01  6.6546E-01  1.0633E+00  7.2692E-01  9.6603E-01  1.0589E-01  3.1182E-01  1.2509E-01
             1.4277E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1592.78868412869        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2234             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9965E-01  1.5911E+00  7.5213E-01  6.7348E-01  1.1860E+00  1.1357E+00  1.5266E+00  6.8587E-01  1.3460E+00  9.6933E-01
             9.5051E-01
 PARAMETER:  9.9651E-02  5.6444E-01 -1.8485E-01 -2.9530E-01  2.7057E-01  2.2725E-01  5.2302E-01 -2.7707E-01  3.9717E-01  6.8850E-02
             4.9239E-02
 GRADIENT:   4.5786E+02  4.7514E+02  1.9667E+00  1.1282E+02  1.7049E+01  1.2652E+02  9.5816E+01  1.3867E-01  2.4468E+01  7.3398E-01
             8.6098E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1592.79219653521        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2410
 NPARAMETR:  9.9966E-01  1.5942E+00  7.4966E-01  6.7169E-01  1.1857E+00  1.1352E+00  1.5240E+00  6.7275E-01  1.3468E+00  9.6889E-01
             9.5044E-01
 PARAMETER:  9.9663E-02  5.6640E-01 -1.8814E-01 -2.9796E-01  2.7030E-01  2.2685E-01  5.2137E-01 -2.9638E-01  3.9770E-01  6.8395E-02
             4.9172E-02
 GRADIENT:   2.2386E+00 -3.9303E+00  8.6899E-01 -6.5786E-01 -6.2802E-01  7.4359E-01  7.7729E-01 -2.0165E-02 -1.3616E-01  4.4396E-03
            -4.3563E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1592.79595465763        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2589             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9941E-01  1.5971E+00  7.4618E-01  6.7011E-01  1.1854E+00  1.1352E+00  1.5222E+00  6.5501E-01  1.3485E+00  9.6840E-01
             9.5036E-01
 PARAMETER:  9.9407E-02  5.6818E-01 -1.9279E-01 -3.0031E-01  2.7010E-01  2.2682E-01  5.2014E-01 -3.2311E-01  3.9897E-01  6.7887E-02
             4.9087E-02
 GRADIENT:   4.5756E+02  4.8209E+02  2.8705E+00  1.1348E+02  1.6226E+01  1.2622E+02  9.5199E+01  4.9720E-02  2.4148E+01  5.7246E-01
             6.7209E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1592.79858812033        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2765
 NPARAMETR:  9.9892E-01  1.5990E+00  7.4360E-01  6.6973E-01  1.1855E+00  1.1343E+00  1.5202E+00  6.5869E-01  1.3496E+00  9.6854E-01
             9.5046E-01
 PARAMETER:  9.8918E-02  5.6935E-01 -1.9626E-01 -3.0088E-01  2.7020E-01  2.2606E-01  5.1887E-01 -3.1750E-01  3.9978E-01  6.8036E-02
             4.9192E-02
 GRADIENT:   8.9131E-01 -3.4813E+00  3.6506E-01  4.7414E-01  1.4277E-01  4.2437E-01  6.0821E-01  1.4752E-02 -9.7366E-03  1.5027E-02
            -2.8762E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1592.80327059133        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2943
 NPARAMETR:  9.9821E-01  1.6036E+00  7.3821E-01  6.6717E-01  1.1851E+00  1.1338E+00  1.5167E+00  6.4336E-01  1.3516E+00  9.6810E-01
             9.5046E-01
 PARAMETER:  9.8213E-02  5.7226E-01 -2.0353E-01 -3.0472E-01  2.6984E-01  2.2553E-01  5.1655E-01 -3.4105E-01  4.0129E-01  6.7584E-02
             4.9187E-02
 GRADIENT:  -3.7323E-01 -3.0882E+00  4.3674E-01  6.8865E-01 -8.1288E-02  2.1373E-01  4.6069E-01  1.4510E-02 -8.6789E-02  2.8075E-02
            -3.2449E-03

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1592.80498901351        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3118
 NPARAMETR:  9.9619E-01  1.6150E+00  7.2492E-01  6.5982E-01  1.1851E+00  1.1302E+00  1.5034E+00  6.0459E-01  1.3569E+00  9.6722E-01
             9.5043E-01
 PARAMETER:  9.6185E-02  5.7931E-01 -2.2169E-01 -3.1579E-01  2.6983E-01  2.2236E-01  5.0771E-01 -4.0320E-01  4.0522E-01  6.6674E-02
             4.9159E-02
 GRADIENT:  -4.0184E+00 -3.0637E+00  6.8360E-01  8.6279E-01 -2.6791E-01 -1.0711E+00 -6.3191E-01 -2.9761E-03 -4.5693E-01 -4.3761E-02
            -6.6437E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1592.80724289598        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3293
 NPARAMETR:  9.9289E-01  1.6397E+00  6.8662E-01  6.3963E-01  1.1835E+00  1.1246E+00  1.4762E+00  4.6978E-01  1.3749E+00  9.6378E-01
             9.5020E-01
 PARAMETER:  9.2863E-02  5.9454E-01 -2.7598E-01 -3.4687E-01  2.6848E-01  2.1747E-01  4.8946E-01 -6.5549E-01  4.1841E-01  6.3107E-02
             4.8913E-02
 GRADIENT:  -1.0197E+01 -5.8735E+00  5.9271E-01  1.8660E-01 -1.9220E-02 -3.1084E+00 -2.3213E+00  3.7861E-02 -7.8470E-01 -1.0861E-02
            -5.7588E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1592.80802983013        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3468
 NPARAMETR:  9.9134E-01  1.6546E+00  6.6202E-01  6.2639E-01  1.1820E+00  1.1222E+00  1.4607E+00  3.7401E-01  1.3869E+00  9.6038E-01
             9.4989E-01
 PARAMETER:  9.1297E-02  6.0358E-01 -3.1246E-01 -3.6777E-01  2.6723E-01  2.1527E-01  4.7895E-01 -8.8347E-01  4.2705E-01  5.9579E-02
             4.8587E-02
 GRADIENT:  -1.3193E+01 -8.4487E+00  1.9866E-01 -5.6156E-01  9.7005E-02 -4.0415E+00 -3.0545E+00  7.5260E-02 -8.3012E-01  4.1851E-02
            -7.2131E-03

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1592.80832619672        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3643
 NPARAMETR:  9.9064E-01  1.6633E+00  6.4786E-01  6.1813E-01  1.1813E+00  1.1212E+00  1.4520E+00  3.1553E-01  1.3945E+00  9.5793E-01
             9.4962E-01
 PARAMETER:  9.0596E-02  6.0879E-01 -3.3407E-01 -3.8106E-01  2.6665E-01  2.1438E-01  4.7295E-01 -1.0535E+00  4.3256E-01  5.7022E-02
             4.8308E-02
 GRADIENT:  -1.4560E+01 -1.0398E+01 -9.2251E-02 -1.3199E+00  3.7446E-02 -4.4208E+00 -3.3552E+00  8.5109E-02 -8.2162E-01  5.5421E-02
             1.8908E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1592.80845945360        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     3819
 NPARAMETR:  9.9033E-01  1.6686E+00  6.3918E-01  6.1263E-01  1.1811E+00  1.1208E+00  1.4468E+00  2.7744E-01  1.3999E+00  9.5628E-01
             9.4942E-01
 PARAMETER:  9.0285E-02  6.1199E-01 -3.4757E-01 -3.8999E-01  2.6647E-01  2.1406E-01  4.6933E-01 -1.1822E+00  4.3642E-01  5.5299E-02
             4.8097E-02
 GRADIENT:  -1.5177E+01 -1.1894E+01 -2.8284E-01 -1.9956E+00 -5.0260E-02 -4.5531E+00 -3.4659E+00  8.3748E-02 -7.9780E-01  5.3943E-02
             3.1138E-02

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1592.80906236903        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3994
 NPARAMETR:  9.9017E-01  1.6743E+00  6.2980E-01  6.0618E-01  1.1812E+00  1.1208E+00  1.4414E+00  2.3175E-01  1.4068E+00  9.5447E-01
             9.4917E-01
 PARAMETER:  9.0125E-02  6.1539E-01 -3.6235E-01 -4.0058E-01  2.6653E-01  2.1407E-01  4.6559E-01 -1.3621E+00  4.4131E-01  5.3398E-02
             4.7831E-02
 GRADIENT:  -1.5515E+01 -1.3916E+01 -4.9467E-01 -2.9890E+00 -1.7023E-01 -4.5488E+00 -3.4720E+00  7.3745E-02 -7.3941E-01  4.6152E-02
             4.2155E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1592.81212267721        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     4171
 NPARAMETR:  9.9038E-01  1.6796E+00  6.2054E-01  5.9900E-01  1.1818E+00  1.1216E+00  1.4369E+00  1.7400E-01  1.4156E+00  9.5275E-01
             9.4886E-01
 PARAMETER:  9.0337E-02  6.1858E-01 -3.7717E-01 -4.1249E-01  2.6703E-01  2.1478E-01  4.6249E-01 -1.6487E+00  4.4752E-01  5.1595E-02
             4.7505E-02
 GRADIENT:  -1.5163E+01 -1.6604E+01 -6.9789E-01 -4.4084E+00 -2.8713E-01 -4.2462E+00 -3.2423E+00  5.1345E-02 -6.1515E-01  3.3378E-02
             5.0152E-02

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1592.94237897985        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     4350
 NPARAMETR:  9.9378E-01  1.6663E+00  6.3475E-01  6.1218E-01  1.1800E+00  1.1274E+00  1.4611E+00  2.4045E-02  1.4082E+00  9.5839E-01
             9.4915E-01
 PARAMETER:  9.3763E-02  6.1059E-01 -3.5453E-01 -3.9073E-01  2.6549E-01  2.1992E-01  4.7917E-01 -3.6278E+00  4.4231E-01  5.7496E-02
             4.7806E-02
 GRADIENT:  -8.9309E+00 -1.2483E+01 -2.8604E-01 -2.1878E+00  8.8427E-01 -2.1061E+00 -1.1288E+00  7.5052E-04 -5.5823E-02  3.9386E-01
             1.4994E-01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1592.97883300225        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     4526
 NPARAMETR:  9.9683E-01  1.6934E+00  6.0671E-01  6.0113E-01  1.1698E+00  1.1306E+00  1.4491E+00  1.0000E-02  1.4057E+00  9.3836E-01
             9.4852E-01
 PARAMETER:  9.6825E-02  6.2676E-01 -3.9971E-01 -4.0894E-01  2.5682E-01  2.2276E-01  4.7097E-01 -9.3500E+00  4.4050E-01  3.6382E-02
             4.7145E-02
 GRADIENT:  -3.8909E+00 -6.8168E+00 -8.3948E-01 -8.9927E-01 -3.3246E+00 -8.9106E-01 -4.5080E-01  0.0000E+00 -4.2305E-01 -4.3089E-01
            -4.7320E-02

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1592.98292484736        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     4704
 NPARAMETR:  9.9874E-01  1.7233E+00  5.8947E-01  5.8295E-01  1.1769E+00  1.1335E+00  1.4325E+00  1.0000E-02  1.4312E+00  9.3667E-01
             9.4867E-01
 PARAMETER:  9.8740E-02  6.4423E-01 -4.2854E-01 -4.3966E-01  2.6288E-01  2.2534E-01  4.5942E-01 -1.4359E+01  4.5850E-01  3.4581E-02
             4.7301E-02
 GRADIENT:  -6.9882E-01 -6.3521E+00 -9.2731E-01 -1.3079E+00 -4.4343E+00  1.8227E-01  1.5897E-01  0.0000E+00 -3.1054E-01 -8.2416E-01
            -1.3806E-02

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1593.02573477162        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     4887
 NPARAMETR:  9.9945E-01  1.7348E+00  5.9207E-01  5.7805E-01  1.1868E+00  1.1342E+00  1.4279E+00  1.0000E-02  1.4432E+00  9.4997E-01
             9.4861E-01
 PARAMETER:  9.9453E-02  6.5088E-01 -4.2413E-01 -4.4810E-01  2.7130E-01  2.2593E-01  4.5618E-01 -1.5661E+01  4.6685E-01  4.8674E-02
             4.7245E-02
 GRADIENT:   5.5126E-01 -5.2751E+00 -8.3536E-01 -4.6657E-01 -3.2951E+00  4.1793E-01  3.3351E-01  0.0000E+00 -5.0593E-01 -2.6793E-01
            -5.9254E-02

0ITERATION NO.:  145    OBJECTIVE VALUE:  -1593.05043372071        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     5063
 NPARAMETR:  9.9603E-01  1.7434E+00  6.0624E-01  5.7287E-01  1.2053E+00  1.1292E+00  1.4182E+00  1.0000E-02  1.4713E+00  9.6869E-01
             9.4935E-01
 PARAMETER:  9.6024E-02  6.5586E-01 -4.0048E-01 -4.5710E-01  2.8671E-01  2.2153E-01  4.4936E-01 -1.5661E+01  4.8615E-01  6.8191E-02
             4.8021E-02
 GRADIENT:  -5.2335E+00 -6.0857E+00  3.8848E-01 -3.5649E-01 -8.7163E-01 -1.3854E+00 -8.2072E-01  0.0000E+00 -4.9963E-01 -1.7696E-01
            -9.1041E-02

0ITERATION NO.:  150    OBJECTIVE VALUE:  -1593.07527952573        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     5223
 NPARAMETR:  1.0006E+00  1.7447E+00  6.0594E-01  5.7205E-01  1.2067E+00  1.1361E+00  1.4247E+00  1.0000E-02  1.4782E+00  9.7021E-01
             9.4939E-01
 PARAMETER:  1.0058E-01  6.5659E-01 -4.0098E-01 -4.5852E-01  2.8787E-01  2.2757E-01  4.5399E-01 -1.5661E+01  4.9084E-01  6.9760E-02
             4.8068E-02
 GRADIENT:   2.6755E+00 -6.2501E+00 -2.0498E-01 -3.9313E-02 -2.6651E-01  1.0511E+00  7.5721E-01  0.0000E+00  1.2775E-01  2.1002E-02
             5.8513E-02

0ITERATION NO.:  151    OBJECTIVE VALUE:  -1593.07527952573        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     5251
 NPARAMETR:  9.9958E-01  1.7433E+00  6.0692E-01  5.7095E-01  1.2070E+00  1.1344E+00  1.4240E+00  1.0000E-02  1.4782E+00  9.6994E-01
             9.4914E-01
 PARAMETER:  1.0058E-01  6.5659E-01 -4.0098E-01 -4.5852E-01  2.8787E-01  2.2757E-01  4.5399E-01 -1.5661E+01  4.9084E-01  6.9760E-02
             4.8068E-02
 GRADIENT:   4.5795E-01  4.0796E-01 -1.9723E-01  5.0632E-01 -1.9258E-01  3.1168E-01  7.2243E-02  0.0000E+00  1.9004E-03  2.0141E-02
             5.8826E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     5251
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.6559E-04 -1.4719E-02 -3.4580E-04  2.1330E-02 -3.4660E-02
 SE:             2.9870E-02  2.7193E-02  1.1039E-04  2.0771E-02  2.0804E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8489E-01  5.8831E-01  1.7331E-03  3.0445E-01  9.5709E-02

 ETASHRINKSD(%)  1.0000E-10  8.9015E+00  9.9630E+01  3.0415E+01  3.0304E+01
 ETASHRINKVR(%)  1.0000E-10  1.7011E+01  9.9999E+01  5.1579E+01  5.1425E+01
 EBVSHRINKSD(%)  3.0091E-01  8.7283E+00  9.9712E+01  3.3080E+01  2.8102E+01
 EBVSHRINKVR(%)  6.0091E-01  1.6695E+01  9.9999E+01  5.5217E+01  4.8306E+01
 RELATIVEINF(%)  9.9328E+01  1.4383E+01  1.3881E-04  5.0712E+00  1.3851E+01
 EPSSHRINKSD(%)  4.4825E+01
 EPSSHRINKVR(%)  6.9557E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1593.0752795257301     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -857.92445296199196     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1593.075       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.74E+00  6.06E-01  5.72E-01  1.21E+00  1.14E+00  1.42E+00  1.00E-02  1.48E+00  9.70E-01  9.49E-01
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,      185.817
Stop Time:
Sun Oct 24 04:46:34 CDT 2021
