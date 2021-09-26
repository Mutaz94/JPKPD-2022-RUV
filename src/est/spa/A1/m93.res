Sat Sep 25 08:20:24 CDT 2021
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
$DATA ../../../../data/spa/A1/dat93.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m93.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1280.32702784695        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.8981E+01  1.4237E+01  2.2308E+01  2.5712E+00  1.2575E+02  2.4851E+01 -1.4158E+01 -7.3462E+00 -6.9709E+00 -8.3835E+01
            -6.4718E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1471.01049066427        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0082E+00  8.6013E-01  9.2486E-01  1.0922E+00  8.4690E-01  8.5191E-01  9.1017E-01  8.9314E-01  9.2349E-01  1.0610E+00
             1.8735E+00
 PARAMETER:  1.0818E-01 -5.0668E-02  2.1887E-02  1.8820E-01 -6.6173E-02 -6.0276E-02  5.8745E-03 -1.3012E-02  2.0401E-02  1.5924E-01
             7.2779E-01
 GRADIENT:   6.4540E+01  7.6215E+00 -2.2225E+00  1.8643E+00  1.6801E+01 -3.7160E+01 -1.8656E+00  7.3461E+00 -8.3954E+00  2.2098E+00
            -3.3549E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1480.88462346985        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0038E+00  6.4954E-01  6.6305E-01  1.2218E+00  6.2296E-01  9.1597E-01  4.5465E-01  2.2750E-01  1.0695E+00  8.7089E-01
             1.9095E+00
 PARAMETER:  1.0383E-01 -3.3149E-01 -3.1090E-01  3.0036E-01 -3.7328E-01  1.2227E-02 -6.8823E-01 -1.3806E+00  1.6717E-01 -3.8239E-02
             7.4687E-01
 GRADIENT:   3.5952E+01  1.9636E+01 -1.0039E+01  6.4339E+01  2.6710E+01 -8.7354E+00 -1.0649E+00  6.3797E-01  2.8822E+01 -1.4220E+00
            -1.2810E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1491.74961865184        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.8676E-01  5.3861E-01  3.5212E-01  1.1616E+00  3.8047E-01  9.5744E-01  1.0892E+00  6.6136E-02  7.9602E-01  6.3546E-01
             1.8557E+00
 PARAMETER:  8.6673E-02 -5.1876E-01 -9.4378E-01  2.4983E-01 -8.6636E-01  5.6508E-02  1.8543E-01 -2.6160E+00 -1.2813E-01 -3.5341E-01
             7.1827E-01
 GRADIENT:  -2.1726E+01  4.1911E+01  1.6734E+01  6.3066E+01 -4.2848E+01  1.9773E+00 -1.3442E+00  4.0081E-02 -2.9066E+01  5.4911E+00
             6.1366E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1494.27221381158        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      327
 NPARAMETR:  9.8904E-01  4.6922E-01  3.0172E-01  1.1284E+00  3.3697E-01  9.6197E-01  1.1185E+00  3.2006E-02  8.6097E-01  5.8188E-01
             1.8435E+00
 PARAMETER:  8.8975E-02 -6.5668E-01 -1.0982E+00  2.2084E-01 -9.8777E-01  6.1227E-02  2.1196E-01 -3.3418E+00 -4.9690E-02 -4.4149E-01
             7.1164E-01
 GRADIENT:  -2.4081E+01  4.7608E+00 -3.3631E+00 -5.3268E-01 -1.2919E+01  2.6508E+00  3.5114E-01  3.1652E-03 -9.2652E+00  4.2485E-01
             3.5824E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1496.98906678175        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      502
 NPARAMETR:  9.9070E-01  3.0967E-01  4.1593E-01  1.2809E+00  3.7608E-01  9.5176E-01  1.3903E+00  1.1187E-02  8.4922E-01  7.0890E-01
             1.8396E+00
 PARAMETER:  9.0656E-02 -1.0723E+00 -7.7725E-01  3.4759E-01 -8.7796E-01  5.0556E-02  4.2950E-01 -4.3930E+00 -6.3441E-02 -2.4404E-01
             7.0952E-01
 GRADIENT:  -1.2113E+01  4.7299E+00 -5.3504E+00  1.9355E+01 -2.6360E+00  2.6130E+00 -5.9440E-01  1.6496E-03  1.2769E+00  2.4026E+00
            -1.1650E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1498.74866703749        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      678
 NPARAMETR:  9.8460E-01  1.2736E-01  4.9106E-01  1.3910E+00  3.9639E-01  9.3950E-01  2.9561E+00  1.0000E-02  8.0059E-01  7.5052E-01
             1.8460E+00
 PARAMETER:  8.4477E-02 -1.9608E+00 -6.1119E-01  4.3001E-01 -8.2537E-01  3.7596E-02  1.1839E+00 -8.3998E+00 -1.2241E-01 -1.8698E-01
             7.1302E-01
 GRADIENT:  -3.4915E+00  1.5128E+00  4.1649E+00  7.3626E+00 -5.2554E+00  9.5956E-01  6.0257E-01  0.0000E+00  3.2625E-01 -1.2412E+00
            -1.5053E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1499.03742562781        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      853
 NPARAMETR:  9.8181E-01  4.8992E-02  5.0072E-01  1.4283E+00  3.9260E-01  9.3286E-01  4.8148E+00  1.0000E-02  7.8709E-01  7.6899E-01
             1.8512E+00
 PARAMETER:  8.1638E-02 -2.9161E+00 -5.9172E-01  4.5648E-01 -8.3497E-01  3.0495E-02  1.6717E+00 -1.3304E+01 -1.3942E-01 -1.6268E-01
             7.1584E-01
 GRADIENT:  -1.2053E+00  3.7662E-01  8.3015E-01  6.7445E+00 -3.1957E+00 -9.3078E-01  2.2444E-01  0.0000E+00 -4.6061E-01  4.4365E-01
             5.6664E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1499.11211105531        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  9.8030E-01  1.1717E-02  5.1164E-01  1.4470E+00  3.9453E-01  9.3371E-01  7.6226E+00  1.0000E-02  7.8208E-01  7.7788E-01
             1.8502E+00
 PARAMETER:  8.0100E-02 -4.3467E+00 -5.7013E-01  4.6953E-01 -8.3007E-01  3.1406E-02  2.1311E+00 -2.0775E+01 -1.4580E-01 -1.5118E-01
             7.1529E-01
 GRADIENT:  -3.1855E-02  1.0670E-02  5.8206E-01  1.3635E+00 -1.0781E+00 -1.7284E-01 -3.3324E-02  0.0000E+00  2.9749E-03  1.5075E-01
             2.2783E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1499.11855490827        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1203
 NPARAMETR:  9.8019E-01  1.0000E-02  5.1095E-01  1.4467E+00  3.9411E-01  9.3413E-01  9.6317E+00  1.0000E-02  7.8180E-01  7.7711E-01
             1.8489E+00
 PARAMETER:  7.9990E-02 -5.1506E+00 -5.7148E-01  4.6929E-01 -8.3112E-01  3.1859E-02  2.3651E+00 -2.5038E+01 -1.4616E-01 -1.5218E-01
             7.1458E-01
 GRADIENT:  -1.4609E-02  0.0000E+00 -9.2436E-02 -2.7860E-01  1.5832E-01  2.3882E-02 -7.0692E-03  0.0000E+00 -7.5922E-03 -2.0669E-02
            -1.0679E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1499.11866827563        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1381
 NPARAMETR:  9.8020E-01  1.0000E-02  5.1121E-01  1.4470E+00  3.9423E-01  9.3406E-01  9.8534E+00  1.0000E-02  7.8174E-01  7.7724E-01
             1.8493E+00
 PARAMETER:  8.0001E-02 -5.2268E+00 -5.7097E-01  4.6947E-01 -8.3081E-01  3.1785E-02  2.3878E+00 -2.5442E+01 -1.4624E-01 -1.5201E-01
             7.1480E-01
 GRADIENT:  -3.4060E-04  0.0000E+00  2.1733E-04 -4.2420E-03 -1.4551E-03 -1.2148E-03  3.5901E-04  0.0000E+00 -1.1118E-03  9.7011E-04
            -4.7032E-03

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1499.11866827563        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1403
 NPARAMETR:  9.8020E-01  1.0000E-02  5.1121E-01  1.4470E+00  3.9423E-01  9.3406E-01  9.8534E+00  1.0000E-02  7.8174E-01  7.7724E-01
             1.8493E+00
 PARAMETER:  8.0001E-02 -5.2268E+00 -5.7097E-01  4.6947E-01 -8.3081E-01  3.1785E-02  2.3878E+00 -2.5442E+01 -1.4624E-01 -1.5201E-01
             7.1480E-01
 GRADIENT:  -3.4060E-04  0.0000E+00  2.1733E-04 -4.2420E-03 -1.4551E-03 -1.2148E-03  3.5901E-04  0.0000E+00 -1.1118E-03  9.7011E-04
            -4.7032E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1403
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.5965E-05  5.5929E-04 -4.1952E-05 -7.1831E-03 -8.3790E-03
 SE:             2.9460E-02  1.6127E-03  2.4123E-04  2.7828E-02  2.4390E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9930E-01  7.2874E-01  8.6194E-01  7.9631E-01  7.3118E-01

 ETASHRINKSD(%)  1.3039E+00  9.4597E+01  9.9192E+01  6.7735E+00  1.8292E+01
 ETASHRINKVR(%)  2.5907E+00  9.9708E+01  9.9993E+01  1.3088E+01  3.3238E+01
 EBVSHRINKSD(%)  1.4529E+00  9.4847E+01  9.9192E+01  6.2764E+00  1.7309E+01
 EBVSHRINKVR(%)  2.8848E+00  9.9734E+01  9.9993E+01  1.2159E+01  3.1621E+01
 RELATIVEINF(%)  8.3506E+01  1.5298E-02  2.9576E-04  7.9413E+00  2.6103E+00
 EPSSHRINKSD(%)  3.7671E+01
 EPSSHRINKVR(%)  6.1151E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1499.1186682756340     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -763.96784171189586     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1499.119       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.00E-02  5.11E-01  1.45E+00  3.94E-01  9.34E-01  9.85E+00  1.00E-02  7.82E-01  7.77E-01  1.85E+00
 


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
+        1.29E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.70E+00  0.00E+00  2.81E+03
 
 TH 4
+       -2.90E+01  0.00E+00 -3.05E+02  7.64E+02
 
 TH 5
+        5.43E+01  0.00E+00 -4.44E+03 -1.58E+02  8.13E+03
 
 TH 6
+        2.06E+01  0.00E+00  8.63E+00 -7.28E+00 -3.01E+00  2.11E+02
 
 TH 7
+        1.25E-02  0.00E+00  1.39E-02 -9.39E-04 -1.43E-02 -4.88E-03  2.01E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.88E+00  0.00E+00  5.17E+01 -9.73E+00  3.47E+00  8.23E+00  2.56E-02  0.00E+00  2.53E+02
 
 TH10
+       -5.76E+00  0.00E+00 -4.46E+01  4.69E+00 -6.87E+01  1.64E+00 -5.15E-04  0.00E+00 -2.07E+00  1.40E+02
 
 TH11
+       -1.22E+01  0.00E+00 -2.28E+01 -9.59E+00  1.73E+01  1.72E+00  4.95E-03  0.00E+00  9.89E+00  2.99E+01  7.32E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.508
Stop Time:
Sat Sep 25 08:20:47 CDT 2021
