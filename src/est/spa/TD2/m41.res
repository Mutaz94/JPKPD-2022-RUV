Wed Sep 29 19:00:13 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat41.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1691.13853922884        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5199E+02 -6.6405E+01 -2.8068E+01 -2.5630E+01  9.1249E+01  4.8111E+01  6.9675E+00 -4.5997E-01  3.0784E+01 -1.2782E+01
            -3.0306E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1701.25613421641        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0323E+00  1.0640E+00  9.6374E-01  1.0408E+00  9.3495E-01  9.4853E-01  9.2856E-01  1.0236E+00  8.0780E-01  1.0434E+00
             1.0143E+00
 PARAMETER:  1.3175E-01  1.6199E-01  6.3062E-02  1.3996E-01  3.2736E-02  4.7161E-02  2.5879E-02  1.2335E-01 -1.1345E-01  1.4248E-01
             1.1419E-01
 GRADIENT:   4.4801E-01  1.2993E+01  2.4141E+00  9.0780E+00 -7.4311E+00 -1.0238E+01 -3.4035E+00  4.7616E-01 -1.1618E+01  4.7869E+00
             2.4201E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1701.54087097307        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0367E+00  9.9620E-01  8.2285E-01  1.0732E+00  8.3801E-01  9.5905E-01  1.0740E+00  8.7657E-01  7.6826E-01  9.1111E-01
             1.0015E+00
 PARAMETER:  1.3603E-01  9.6189E-02 -9.4987E-02  1.7061E-01 -7.6726E-02  5.8184E-02  1.7138E-01 -3.1741E-02 -1.6363E-01  6.9131E-03
             1.0151E-01
 GRADIENT:   8.4596E+00  9.7092E+00 -8.0366E+00  1.8643E+01  6.2705E+00 -6.2124E+00  1.0644E+00  2.6566E+00 -7.5543E+00  1.1121E+00
            -1.4925E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1702.08432321209        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      546
 NPARAMETR:  1.0318E+00  9.4241E-01  8.8272E-01  1.1027E+00  8.4562E-01  9.7407E-01  1.0559E+00  8.3261E-01  8.0352E-01  9.4190E-01
             1.0047E+00
 PARAMETER:  1.3129E-01  4.0690E-02 -2.4749E-02  1.9772E-01 -6.7683E-02  7.3731E-02  1.5438E-01 -8.3193E-02 -1.1876E-01  4.0140E-02
             1.0464E-01
 GRADIENT:  -8.5203E-01  2.0475E+00  1.1174E+00  1.4891E+00 -1.2506E+00  3.6400E-01  1.8888E-02 -3.4030E-01 -6.6553E-03 -2.2852E-01
             3.5402E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1702.08462778227        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  1.0316E+00  9.1491E-01  9.0628E-01  1.1218E+00  8.4547E-01  9.7332E-01  1.0783E+00  8.4933E-01  7.9440E-01  9.4660E-01
             1.0038E+00
 PARAMETER:  1.3110E-01  1.1067E-02  1.5968E-03  2.1490E-01 -6.7862E-02  7.2955E-02  1.7540E-01 -6.3307E-02 -1.3017E-01  4.5120E-02
             1.0376E-01
 GRADIENT:  -4.7421E-01  3.4079E+00  2.0080E+00  2.6462E+00 -2.8329E+00  2.0561E-01  5.1702E-02 -4.3865E-01 -1.6529E-01 -2.1966E-01
            -3.4749E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1702.08483895926        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      896
 NPARAMETR:  1.0315E+00  8.9519E-01  9.2337E-01  1.1353E+00  8.4557E-01  9.7263E-01  1.0949E+00  8.6263E-01  7.8805E-01  9.5010E-01
             1.0030E+00
 PARAMETER:  1.3101E-01 -1.0716E-02  2.0270E-02  2.2692E-01 -6.7745E-02  7.2251E-02  1.9062E-01 -4.7770E-02 -1.3820E-01  4.8809E-02
             1.0298E-01
 GRADIENT:  -7.0036E-02  4.1451E+00  2.5067E+00  3.3252E+00 -3.8361E+00  3.7789E-02  7.1751E-02 -4.5887E-01 -2.8134E-01 -1.7814E-01
            -3.7030E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1702.08493747695        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1072
 NPARAMETR:  1.0314E+00  8.8029E-01  9.3613E-01  1.1455E+00  8.4562E-01  9.7206E-01  1.1077E+00  8.7281E-01  7.8332E-01  9.5271E-01
             1.0023E+00
 PARAMETER:  1.3095E-01 -2.7500E-02  3.3999E-02  2.3584E-01 -6.7682E-02  7.1662E-02  2.0231E-01 -3.6035E-02 -1.4421E-01  5.1552E-02
             1.0231E-01
 GRADIENT:   2.9040E-01  4.6023E+00  2.8269E+00  3.7738E+00 -4.5392E+00 -1.1101E-01  8.5640E-02 -4.5423E-01 -3.7014E-01 -1.3303E-01
            -6.4759E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1702.08509892712        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1247
 NPARAMETR:  1.0314E+00  8.6385E-01  9.5003E-01  1.1566E+00  8.4565E-01  9.7139E-01  1.1223E+00  8.8402E-01  7.7820E-01  9.5555E-01
             1.0015E+00
 PARAMETER:  1.3091E-01 -4.6361E-02  4.8734E-02  2.4552E-01 -6.7648E-02  7.0969E-02  2.1541E-01 -2.3274E-02 -1.5078E-01  5.4533E-02
             1.0155E-01
 GRADIENT:   7.2553E-01  5.0000E+00  3.1193E+00  4.1869E+00 -5.2411E+00 -2.9089E-01  9.9680E-02 -4.3217E-01 -4.6636E-01 -7.2744E-02
            -9.6627E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1702.08533288264        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1422
 NPARAMETR:  1.0313E+00  8.4002E-01  9.6992E-01  1.1727E+00  8.4569E-01  9.7038E-01  1.1442E+00  9.0014E-01  7.7095E-01  9.5966E-01
             1.0004E+00
 PARAMETER:  1.3084E-01 -7.4328E-02  6.9457E-02  2.5928E-01 -6.7599E-02  6.9928E-02  2.3474E-01 -5.2045E-03 -1.6013E-01  5.8825E-02
             1.0042E-01
 GRADIENT:   1.3820E+00  5.4275E+00  3.4527E+00  4.6821E+00 -6.1381E+00 -5.6371E-01  1.1759E-01 -3.7861E-01 -6.0037E-01  2.5962E-02
            -1.4295E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1702.08569072346        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1597
 NPARAMETR:  1.0312E+00  8.0854E-01  9.9572E-01  1.1936E+00  8.4570E-01  9.6899E-01  1.1747E+00  9.2103E-01  7.6169E-01  9.6507E-01
             9.9889E-01
 PARAMETER:  1.3077E-01 -1.1253E-01  9.5706E-02  2.7701E-01 -6.7587E-02  6.8495E-02  2.6099E-01  1.7740E-02 -1.7222E-01  6.4448E-02
             9.8888E-02
 GRADIENT:   2.3050E+00  5.7701E+00  3.7598E+00  5.1523E+00 -7.1402E+00 -9.4476E-01  1.3800E-01 -2.7692E-01 -7.6869E-01  1.7291E-01
            -2.0528E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1702.09928525842        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1776             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0309E+00  7.5991E-01  1.0217E+00  1.2178E+00  8.5061E-01  9.7066E-01  1.2135E+00  9.4915E-01  7.5379E-01  9.6975E-01
             1.0034E+00
 PARAMETER:  1.3046E-01 -1.7455E-01  1.2150E-01  2.9703E-01 -6.1807E-02  7.0219E-02  2.9354E-01  4.7813E-02 -1.8264E-01  6.9279E-02
             1.0336E-01
 GRADIENT:   5.3578E+02  2.6602E+01 -2.6294E+00  2.7069E+02  1.7178E+01  3.4319E+01  6.1413E+00  1.0485E+00  9.2223E+00  5.5301E-01
             1.0675E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1702.14008597672        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1953
 NPARAMETR:  1.0307E+00  7.6213E-01  1.0282E+00  1.2195E+00  8.4824E-01  9.7080E-01  1.2130E+00  9.4113E-01  7.5380E-01  9.7237E-01
             1.0032E+00
 PARAMETER:  1.3023E-01 -1.7163E-01  1.2780E-01  2.9845E-01 -6.4593E-02  7.0368E-02  2.9309E-01  3.9324E-02 -1.8262E-01  7.1981E-02
             1.0324E-01
 GRADIENT:   2.5815E+00  1.6863E-01  7.6839E-01 -3.0969E+00 -1.3578E+00  9.1550E-02  5.8600E-02 -1.6876E-01  3.0403E-02 -6.8234E-02
            -2.9797E-02

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1702.14091919033        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     2010
 NPARAMETR:  1.0305E+00  7.6209E-01  1.0274E+00  1.2200E+00  8.4850E-01  9.7078E-01  1.2126E+00  9.4301E-01  7.5379E-01  9.7266E-01
             1.0033E+00
 PARAMETER:  1.3003E-01 -1.7169E-01  1.2699E-01  2.9886E-01 -6.4287E-02  7.0349E-02  2.9280E-01  4.1324E-02 -1.8265E-01  7.2276E-02
             1.0327E-01
 GRADIENT:   2.0683E+00  2.2388E-01 -2.5784E-01 -1.7411E+00 -2.2229E-01  7.7933E-02  3.9688E-02  2.3011E-02  3.6885E-02  1.6759E-02
             2.2771E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2010
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.2872E-04 -3.7812E-03 -3.0610E-02 -3.6196E-03 -2.8565E-02
 SE:             2.9843E-02  1.6916E-02  1.4968E-02  2.5009E-02  2.2052E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8586E-01  8.2312E-01  4.0845E-02  8.8492E-01  1.9519E-01

 ETASHRINKSD(%)  2.3190E-02  4.3330E+01  4.9856E+01  1.6215E+01  2.6125E+01
 ETASHRINKVR(%)  4.6374E-02  6.7885E+01  7.4856E+01  2.9801E+01  4.5424E+01
 EBVSHRINKSD(%)  4.5224E-01  4.3881E+01  5.3026E+01  1.6495E+01  2.3489E+01
 EBVSHRINKVR(%)  9.0244E-01  6.8506E+01  7.7934E+01  3.0269E+01  4.1460E+01
 RELATIVEINF(%)  9.7434E+01  8.4779E-01  2.4280E+00  2.1807E+00  6.7375E+00
 EPSSHRINKSD(%)  4.5128E+01
 EPSSHRINKVR(%)  6.9891E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1702.1409191903272     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -966.99009262658899     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.40
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
 





 #OBJV:********************************************    -1702.141       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  7.62E-01  1.03E+00  1.22E+00  8.48E-01  9.71E-01  1.21E+00  9.43E-01  7.54E-01  9.73E-01  1.00E+00
 


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
+       -1.51E+01  4.28E+02
 
 TH 3
+        1.08E+01  1.18E+02  2.77E+02
 
 TH 4
+       -9.91E+00  5.03E+02 -1.13E+02  9.11E+02
 
 TH 5
+        1.37E+00 -2.85E+02 -4.36E+02  9.49E+01  9.63E+02
 
 TH 6
+        2.92E-01 -3.10E+00  2.12E+00 -2.65E+00 -4.83E-01  2.08E+02
 
 TH 7
+        9.41E-01  1.39E+01  3.62E+00  1.17E-01 -7.38E+00  1.75E-01  1.66E+01
 
 TH 8
+       -2.60E-01 -1.39E+01 -4.49E+01  2.32E+00  7.38E+00 -1.66E-01  1.34E+00  2.93E+01
 
 TH 9
+        1.76E+00 -1.48E+01 -8.41E+00  5.67E+00  4.33E+00 -1.86E-01  3.01E+01  3.13E+00  1.77E+02
 
 TH10
+        2.48E-01  6.64E+00 -1.66E+01 -1.35E+01 -7.36E+01  2.19E-01  6.31E+00  1.91E+01  7.48E-01  7.54E+01
 
 TH11
+       -7.16E+00 -1.03E+01 -1.29E+01 -6.45E+00 -3.32E+00  1.87E+00  1.97E+00  5.07E+00  1.08E+01  1.32E+01  2.09E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       31.517
Stop Time:
Wed Sep 29 19:00:46 CDT 2021
