Wed Sep 29 15:34:18 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat6.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1684.82668127736        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3136E+02 -5.5759E+01 -8.7980E+01  4.4697E+01  1.4495E+02  5.2476E+01  4.4744E+00  1.2812E+01  6.7260E+00  6.2331E+00
             1.4737E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1696.63649921376        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  1.0107E+00  1.0794E+00  1.1687E+00  1.0001E+00  9.8976E-01  9.9621E-01  9.5440E-01  9.2664E-01  1.0148E+00  8.9688E-01
             9.4147E-01
 PARAMETER:  1.1068E-01  1.7639E-01  2.5589E-01  1.0010E-01  8.9708E-02  9.6203E-02  5.3327E-02  2.3813E-02  1.1471E-01 -8.8278E-03
             3.9687E-02
 GRADIENT:   2.3776E+01  9.9282E+00  1.5415E-01  1.3712E+01  1.7747E+01 -3.6534E-01  4.8100E+00  2.7320E+00 -6.5334E+00 -1.4002E+01
            -1.6619E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1698.54002081083        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  1.0180E+00  1.1977E+00  1.1457E+00  9.2639E-01  1.0241E+00  9.8914E-01  6.0125E-01  6.9190E-01  1.1302E+00  1.0010E+00
             9.8967E-01
 PARAMETER:  1.1780E-01  2.8038E-01  2.3601E-01  2.3537E-02  1.2384E-01  8.9080E-02 -4.0875E-01 -2.6831E-01  2.2241E-01  1.0096E-01
             8.9618E-02
 GRADIENT:   3.7888E+01  3.3969E+01  9.6156E+00  2.4491E+01 -4.6000E+00 -3.5245E+00 -2.9848E+00 -1.9901E+00 -1.7943E+01 -9.6382E+00
             4.4862E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.73459804878        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.9808E-01  1.2882E+00  1.1519E+00  8.4528E-01  1.0679E+00  9.9473E-01  3.9336E-01  8.8222E-01  1.3736E+00  1.0959E+00
             9.6294E-01
 PARAMETER:  9.8075E-02  3.5321E-01  2.4139E-01 -6.8091E-02  1.6569E-01  9.4713E-02 -8.3303E-01 -2.5317E-02  4.1745E-01  1.9162E-01
             6.2237E-02
 GRADIENT:  -6.3455E+00  2.9076E+00 -7.4909E-01  8.4213E+00 -3.1735E+00 -9.5422E-01  5.7126E-01  5.9371E-01  2.6631E+00  1.2417E+00
            -1.8909E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1700.90217919271        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  1.0020E+00  1.4027E+00  1.0867E+00  7.6174E-01  1.0963E+00  9.9761E-01  3.7022E-01  8.3343E-01  1.4823E+00  1.1065E+00
             9.6943E-01
 PARAMETER:  1.0202E-01  4.3839E-01  1.8317E-01 -1.7215E-01  1.9197E-01  9.7603E-02 -8.9365E-01 -8.2206E-02  4.9361E-01  2.0119E-01
             6.8951E-02
 GRADIENT:   1.4072E+00  4.4705E+00  1.0134E+00  2.6914E+00 -1.2642E+00  1.2925E-01 -8.6093E-01 -4.2211E-01 -1.8881E+00 -6.7379E-01
             2.4905E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1701.07376645591        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      892
 NPARAMETR:  1.0002E+00  1.5140E+00  9.9573E-01  6.9102E-01  1.1178E+00  9.9745E-01  4.7905E-01  7.9054E-01  1.5581E+00  1.1045E+00
             9.6554E-01
 PARAMETER:  1.0020E-01  5.1475E-01  9.5720E-02 -2.6959E-01  2.1138E-01  9.7451E-02 -6.3595E-01 -1.3504E-01  5.4344E-01  1.9935E-01
             6.4927E-02
 GRADIENT:  -3.4427E+00 -3.2364E-01  1.7999E+00  2.1170E+00 -4.5344E-01 -7.6117E-02 -2.9055E-01  6.5789E-02  3.7275E-01  3.3030E-01
            -1.2139E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1701.41498836515        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  1.0024E+00  1.5963E+00  7.7571E-01  6.4045E-01  1.0842E+00  9.9771E-01  6.5769E-01  4.3600E-01  1.4903E+00  1.0175E+00
             9.6887E-01
 PARAMETER:  1.0240E-01  5.6770E-01 -1.5398E-01 -3.4558E-01  1.8086E-01  9.7710E-02 -3.1902E-01 -7.3012E-01  4.9895E-01  1.1737E-01
             6.8380E-02
 GRADIENT:  -2.5212E-01 -8.5893E+00 -3.1280E+00  1.5332E+00  5.1942E+00 -1.6941E-01  1.8026E+00  6.3495E-01 -2.0282E+00  8.6968E-01
             1.0157E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1701.63330260187        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1246
 NPARAMETR:  1.0031E+00  1.7031E+00  6.6672E-01  5.6782E-01  1.0832E+00  9.9885E-01  6.2751E-01  1.7420E-01  1.6334E+00  9.8942E-01
             9.6851E-01
 PARAMETER:  1.0311E-01  6.3245E-01 -3.0538E-01 -4.6595E-01  1.7994E-01  9.8846E-02 -3.6600E-01 -1.6475E+00  5.9067E-01  8.9361E-02
             6.7999E-02
 GRADIENT:   8.0074E-01 -2.0516E+00 -1.7265E+00  2.4531E+00 -8.5200E+00  2.2258E-01 -8.9291E-01  1.3021E-01 -1.6015E-01 -1.5675E+00
             3.4225E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1701.76506759477        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1427             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0034E+00  1.7753E+00  6.3796E-01  5.1354E-01  1.1187E+00  9.9879E-01  6.1865E-01  5.2609E-02  1.7575E+00  1.0174E+00
             9.6840E-01
 PARAMETER:  1.0344E-01  6.7400E-01 -3.4949E-01 -5.6642E-01  2.1219E-01  9.8789E-02 -3.8022E-01 -2.8449E+00  6.6388E-01  1.1727E-01
             6.7894E-02
 GRADIENT:   4.9411E+02  8.4897E+02 -6.8947E-01  1.0552E+02  1.4193E+01  5.6973E+01  1.6691E+01  1.6070E-02  3.3112E+01  8.6513E-01
             1.3783E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1701.80622529252        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1606
 NPARAMETR:  1.0026E+00  1.7926E+00  6.5520E-01  5.0886E-01  1.1339E+00  9.9809E-01  6.1462E-01  2.7639E-02  1.7696E+00  1.0355E+00
             9.6806E-01
 PARAMETER:  1.0265E-01  6.8369E-01 -3.2281E-01 -5.7558E-01  2.2566E-01  9.8089E-02 -3.8676E-01 -3.4885E+00  6.7077E-01  1.3486E-01
             6.7537E-02
 GRADIENT:  -1.9732E-01 -5.3964E-01  4.6333E-01  1.6392E+00 -2.9036E+00 -4.9506E-02  1.5429E-01  2.4066E-03 -1.5037E+00 -4.1278E-01
            -1.6912E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1701.84147285366        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1793
 NPARAMETR:  1.0036E+00  1.7876E+00  6.5910E-01  5.0722E-01  1.1354E+00  9.9869E-01  6.0888E-01  1.0000E-02  1.7852E+00  1.0376E+00
             9.6781E-01
 PARAMETER:  1.0361E-01  6.8087E-01 -3.1688E-01 -5.7881E-01  2.2695E-01  9.8693E-02 -3.9613E-01 -4.8555E+00  6.7951E-01  1.3695E-01
             6.7283E-02
 GRADIENT:   2.1264E+00 -1.1309E+01 -7.3500E-02 -1.2620E+00 -8.7958E-01  2.0065E-01  4.6599E-02  0.0000E+00 -3.9389E-01 -4.4312E-01
            -1.3449E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1701.84873347549        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1987             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0036E+00  1.7852E+00  6.6252E-01  5.0902E-01  1.1367E+00  9.9869E-01  6.0738E-01  1.0000E-02  1.7896E+00  1.0422E+00
             9.6813E-01
 PARAMETER:  1.0357E-01  6.7952E-01 -3.1170E-01 -5.7527E-01  2.2811E-01  9.8687E-02 -3.9859E-01 -4.8555E+00  6.8202E-01  1.4129E-01
             6.7609E-02
 GRADIENT:   4.9591E+02  8.7099E+02  1.2110E+00  1.0604E+02  1.7165E+01  5.7017E+01  1.6809E+01  0.0000E+00  3.4268E+01  1.3497E+00
             9.4217E-01

0ITERATION NO.:   59    OBJECTIVE VALUE:  -1701.84988000708        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     2123
 NPARAMETR:  1.0036E+00  1.7849E+00  6.6540E-01  5.0951E-01  1.1363E+00  9.9869E-01  6.0711E-01  1.0000E-02  1.7885E+00  1.0418E+00
             9.6781E-01
 PARAMETER:  1.0357E-01  6.7939E-01 -3.0737E-01 -5.7431E-01  2.2782E-01  9.8688E-02 -3.9905E-01 -4.8555E+00  6.8136E-01  1.4099E-01
             6.7280E-02
 GRADIENT:   2.0743E+00 -1.0302E+01  2.7450E-01 -7.0504E-01 -1.1365E+00  1.9319E-01  9.8857E-02  0.0000E+00  9.3081E-02 -1.8323E-01
            -1.3063E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2123
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6963E-04 -4.5801E-02 -2.1743E-04  2.8579E-02 -3.4773E-02
 SE:             2.9854E-02  1.9679E-02  9.2882E-05  2.4209E-02  2.4252E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9547E-01  1.9945E-02  1.9238E-02  2.3780E-01  1.5163E-01

 ETASHRINKSD(%)  1.0000E-10  3.4072E+01  9.9689E+01  1.8897E+01  1.8752E+01
 ETASHRINKVR(%)  1.0000E-10  5.6535E+01  9.9999E+01  3.4223E+01  3.3988E+01
 EBVSHRINKSD(%)  4.1495E-01  3.2831E+01  9.9725E+01  1.9823E+01  1.6963E+01
 EBVSHRINKVR(%)  8.2818E-01  5.4883E+01  9.9999E+01  3.5717E+01  3.1048E+01
 RELATIVEINF(%)  9.9131E+01  4.1415E+00  1.8062E-04  7.2072E+00  2.1771E+01
 EPSSHRINKSD(%)  4.4145E+01
 EPSSHRINKVR(%)  6.8803E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1701.8498800070797     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -966.69905344334154     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.93
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1701.850       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.78E+00  6.65E-01  5.10E-01  1.14E+00  9.99E-01  6.07E-01  1.00E-02  1.79E+00  1.04E+00  9.68E-01
 


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
+        1.10E+03
 
 TH 2
+       -7.05E+00  4.64E+02
 
 TH 3
+        4.40E+00  1.38E+02  2.10E+02
 
 TH 4
+       -1.17E+01  4.02E+02 -1.13E+02  8.74E+02
 
 TH 5
+       -3.39E+00 -2.09E+02 -2.12E+02  1.27E+02  5.46E+02
 
 TH 6
+        2.76E-01 -9.90E-01  5.53E-01 -3.23E+00 -8.07E-01  1.97E+02
 
 TH 7
+        1.74E+00 -2.68E+01  1.63E+01 -2.61E+01 -2.99E+01  5.86E-02  1.21E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.21E+00 -2.37E+01 -1.67E+01  4.83E+01  2.69E+00 -3.06E-01  2.42E+01  0.00E+00  3.00E+01
 
 TH10
+       -2.51E-01 -1.70E+01 -3.32E+01 -2.43E+00 -5.17E+01  4.57E-02  2.26E+01  0.00E+00  2.28E+00  8.99E+01
 
 TH11
+       -8.09E+00 -2.23E+01 -2.06E+01 -1.51E-01  4.97E+00  2.68E+00  9.72E+00  0.00E+00  3.71E+00  1.72E+01  2.32E+02
 
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
 #CPUT: Total CPU Time in Seconds,       34.052
Stop Time:
Wed Sep 29 15:34:54 CDT 2021
