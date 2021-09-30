Thu Sep 30 00:11:00 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat34.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 RAW OUTPUT FILE (FILE): m34.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   142.135514913899        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8860E+02  9.9786E-01  9.2959E+01 -5.2431E+00  2.6397E+02 -1.9636E+01 -4.7625E+01 -1.0354E+02 -1.2426E+01 -2.0473E+02
            -3.9859E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1268.67494560164        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0407E+00  1.1018E+00  9.3228E-01  1.2245E+00  7.8121E-01  9.1014E-01  9.4486E-01  9.5525E-01  7.9892E-01  1.3472E+00
             7.3823E+00
 PARAMETER:  1.3989E-01  1.9690E-01  2.9876E-02  3.0253E-01 -1.4691E-01  5.8481E-03  4.3280E-02  5.4219E-02 -1.2449E-01  3.9801E-01
             2.0991E+00
 GRADIENT:  -5.2273E+01  3.0305E+01  3.9057E+00  3.6147E+01 -4.7608E+01 -1.6988E+01  1.2670E+01  4.7646E+00  2.3868E+01  3.3378E+01
             4.6982E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1401.76995833686        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.5893E-01  7.7480E-01  1.9492E-01  1.1420E+00  3.3094E-01  1.0328E+00  8.8630E-01  2.6659E-02  2.4888E-01  3.4454E-01
             5.2844E+00
 PARAMETER:  5.8059E-02 -1.5515E-01 -1.5352E+00  2.3275E-01 -1.0058E+00  1.3231E-01 -2.0702E-02 -3.5246E+00 -1.2908E+00 -9.6555E-01
             1.7648E+00
 GRADIENT:  -9.2052E+01  2.9196E+01 -6.7987E+01  1.2976E+02  2.1465E+01 -1.1750E+01  1.1626E+01  1.4617E-02  8.6462E-01  6.0745E+00
             2.5101E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1494.57479212716        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.8117E-01  1.1142E+00  6.4531E-01  1.0468E+00  7.5699E-01  9.8982E-01  8.7713E-01  1.0000E-02  9.3985E-01  6.3733E-01
             3.4112E+00
 PARAMETER:  8.0992E-02  2.0812E-01 -3.3803E-01  1.4569E-01 -1.7840E-01  8.9768E-02 -3.1099E-02 -7.2238E+00  3.7970E-02 -3.5046E-01
             1.3271E+00
 GRADIENT:   3.3180E+01  3.0618E+01 -3.6501E+00  4.7110E+01  7.4655E+00 -1.5383E+01  5.1116E+00  0.0000E+00  1.2886E+01  4.5248E+00
            -1.4671E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1501.19403375432        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.5658E-01  7.2367E-01  3.1283E-01  1.1321E+00  3.9365E-01  1.0743E+00  7.8057E-01  1.0000E-02  9.1471E-01  3.6251E-01
             3.4862E+00
 PARAMETER:  5.5606E-02 -2.2343E-01 -1.0621E+00  2.2406E-01 -8.3229E-01  1.7169E-01 -1.4773E-01 -5.8639E+00  1.0850E-02 -9.1471E-01
             1.3488E+00
 GRADIENT:  -2.5442E+01  3.0302E+01  1.0088E+01  3.7521E+01  5.1725E+00  9.6055E+00 -6.8443E+00  0.0000E+00 -9.2066E+00 -3.6757E+00
             4.9490E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1502.67545390067        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  9.8163E-01  6.7024E-01  3.2062E-01  1.1758E+00  3.8084E-01  1.0668E+00  9.1318E-01  1.0000E-02  8.6673E-01  3.4844E-01
             3.4137E+00
 PARAMETER:  8.1456E-02 -3.0011E-01 -1.0375E+00  2.6194E-01 -8.6539E-01  1.6465E-01  9.1763E-03 -5.8734E+00 -4.3033E-02 -9.5430E-01
             1.3278E+00
 GRADIENT:  -4.4979E+00  3.5033E+01  9.5996E+00  6.1048E+01 -2.8997E+01  3.9410E+00 -7.9874E+00  0.0000E+00 -1.7173E+01 -4.9697E+00
             1.4418E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1506.00254657875        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      639
 NPARAMETR:  9.7344E-01  6.3066E-01  2.7369E-01  1.1540E+00  3.3984E-01  1.0463E+00  1.0923E+00  1.0000E-02  1.0012E+00  2.4421E-01
             3.2212E+00
 PARAMETER:  7.3078E-02 -3.6099E-01 -1.1958E+00  2.4326E-01 -9.7928E-01  1.4531E-01  1.8830E-01 -5.8281E+00  1.0116E-01 -1.3097E+00
             1.2698E+00
 GRADIENT:  -1.2501E+01  3.5361E+01  1.3533E+01  5.5026E+01 -3.0247E+01 -4.5067E+00 -8.4711E+00  0.0000E+00 -7.2643E+00 -5.9346E+00
            -9.8496E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1510.94915936211        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  9.6082E-01  4.7381E-01  2.5575E-01  1.1429E+00  2.9297E-01  1.0848E+00  1.3034E+00  1.0000E-02  1.0258E+00  4.7834E-01
             3.3932E+00
 PARAMETER:  6.0029E-02 -6.4695E-01 -1.2635E+00  2.3355E-01 -1.1277E+00  1.8138E-01  3.6495E-01 -6.2520E+00  1.2548E-01 -6.3744E-01
             1.3218E+00
 GRADIENT:  -3.8676E+01  3.6934E+00 -2.7696E+00 -8.0102E+00  7.7615E+00  1.0048E+01  6.2815E+00  0.0000E+00 -8.8231E-01  2.1233E-01
             8.7278E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1518.87443711972        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      991
 NPARAMETR:  9.7494E-01  4.1047E-01  2.0561E-01  1.1071E+00  2.4732E-01  1.0639E+00  9.6830E-01  1.0000E-02  1.1849E+00  6.0665E-01
             2.9997E+00
 PARAMETER:  7.4616E-02 -7.9044E-01 -1.4818E+00  2.0178E-01 -1.2971E+00  1.6195E-01  6.7787E-02 -7.2516E+00  2.6962E-01 -3.9980E-01
             1.1985E+00
 GRADIENT:  -3.4676E+00  1.7159E+00 -1.0807E+00  1.1658E+01 -3.6660E+00  1.0281E+00  3.6472E+00  0.0000E+00  2.4382E+00 -1.4401E+00
             3.3390E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1520.74405023534        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1166
 NPARAMETR:  9.7564E-01  3.9154E-01  1.8153E-01  1.0631E+00  2.2961E-01  1.0596E+00  2.5439E-01  1.0000E-02  1.2385E+00  6.7800E-01
             2.9804E+00
 PARAMETER:  7.5338E-02 -8.3766E-01 -1.6064E+00  1.6119E-01 -1.3714E+00  1.5789E-01 -1.2689E+00 -7.4612E+00  3.1389E-01 -2.8861E-01
             1.1920E+00
 GRADIENT:  -2.4239E+00 -3.9144E+00 -6.3359E+00  1.0797E+00  1.0102E+01 -3.6018E-01  2.6100E-01  0.0000E+00 -2.7812E-01 -2.0838E-01
             3.4879E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1520.89191447241        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1341
 NPARAMETR:  9.7669E-01  3.9441E-01  1.8190E-01  1.0615E+00  2.2967E-01  1.0612E+00  5.2273E-02  1.0000E-02  1.2436E+00  6.8420E-01
             2.9687E+00
 PARAMETER:  7.6411E-02 -8.3036E-01 -1.6043E+00  1.5965E-01 -1.3711E+00  1.5937E-01 -2.8513E+00 -7.4195E+00  3.1799E-01 -2.7950E-01
             1.1881E+00
 GRADIENT:  -1.7706E-01  1.3622E-02 -1.8048E-01 -1.6737E-02  3.5652E-01  1.0522E-01  9.7030E-03  0.0000E+00  7.1595E-02 -1.0093E-01
            -2.4949E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1520.89736609145        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1522
 NPARAMETR:  9.7684E-01  3.9378E-01  1.8157E-01  1.0609E+00  2.2930E-01  1.0609E+00  1.0000E-02  1.0000E-02  1.2441E+00  6.8489E-01
             2.9692E+00
 PARAMETER:  7.6567E-02 -8.3195E-01 -1.6061E+00  1.5911E-01 -1.3727E+00  1.5909E-01 -5.2954E+00 -7.3765E+00  3.1844E-01 -2.7850E-01
             1.1883E+00
 GRADIENT:   6.7245E-02  1.0671E-01  6.9917E-02 -7.7161E-02 -2.8289E-01  2.7736E-02  0.0000E+00  0.0000E+00  5.5963E-02  1.5291E-02
             8.0232E-03

0ITERATION NO.:   56    OBJECTIVE VALUE:  -1520.89736609145        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1544
 NPARAMETR:  9.7684E-01  3.9378E-01  1.8157E-01  1.0609E+00  2.2930E-01  1.0609E+00  1.0000E-02  1.0000E-02  1.2441E+00  6.8489E-01
             2.9692E+00
 PARAMETER:  7.6567E-02 -8.3195E-01 -1.6061E+00  1.5911E-01 -1.3727E+00  1.5909E-01 -5.2954E+00 -7.3765E+00  3.1844E-01 -2.7850E-01
             1.1883E+00
 GRADIENT:   6.7245E-02  1.0671E-01  6.9917E-02 -7.7161E-02 -2.8289E-01  2.7736E-02  0.0000E+00  0.0000E+00  5.5963E-02  1.5291E-02
             8.0232E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1544
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4798E-03 -1.0010E-04  1.7988E-04 -1.1414E-02  1.6673E-03
 SE:             2.8982E-02  1.3343E-04  2.0668E-04  2.6212E-02  2.3910E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5928E-01  4.5310E-01  3.8411E-01  6.6324E-01  9.4440E-01

 ETASHRINKSD(%)  2.9069E+00  9.9553E+01  9.9308E+01  1.2185E+01  1.9900E+01
 ETASHRINKVR(%)  5.7292E+00  9.9998E+01  9.9995E+01  2.2886E+01  3.5840E+01
 EBVSHRINKSD(%)  2.5249E+00  9.9541E+01  9.9355E+01  9.9724E+00  2.0163E+01
 EBVSHRINKVR(%)  4.9861E+00  9.9998E+01  9.9996E+01  1.8950E+01  3.6260E+01
 RELATIVEINF(%)  9.4900E+01  3.1155E-04  3.0376E-04  3.7544E+01  2.3409E+00
 EPSSHRINKSD(%)  2.5664E+01
 EPSSHRINKVR(%)  4.4741E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1520.8973660914544     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -601.95883288678169     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.81
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1520.897       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.77E-01  3.94E-01  1.82E-01  1.06E+00  2.29E-01  1.06E+00  1.00E-02  1.00E-02  1.24E+00  6.85E-01  2.97E+00
 


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
+        9.86E+02
 
 TH 2
+       -3.19E+01  1.51E+03
 
 TH 3
+       -6.52E+01  2.78E+03  1.39E+04
 
 TH 4
+       -1.14E+01  1.15E+02 -5.68E+02  5.17E+02
 
 TH 5
+        1.69E+02 -5.50E+03 -1.74E+04 -3.79E+02  2.90E+04
 
 TH 6
+        5.50E+00 -1.59E+01  4.18E+01 -1.17E+01 -6.44E+00  1.56E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.57E+00 -3.19E+01  1.48E+02 -3.17E+00  5.64E+01 -1.05E+00  0.00E+00  0.00E+00  7.21E+01
 
 TH10
+       -4.44E+00 -2.83E+01 -1.26E+01  1.09E+01  7.62E+01  3.17E+00  0.00E+00  0.00E+00  3.26E+00  1.74E+02
 
 TH11
+       -1.64E+01 -1.03E+01 -4.36E+01 -9.19E+00  4.98E+01  3.10E+00  0.00E+00  0.00E+00  8.64E+00  1.82E+01  5.27E+01
 
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
 #CPUT: Total CPU Time in Seconds,       32.600
Stop Time:
Thu Sep 30 00:11:39 CDT 2021
