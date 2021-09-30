Wed Sep 29 09:41:36 CDT 2021
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
$DATA ../../../../data/int/D/dat75.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   29783.6680520758        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7474E+02  3.9587E+02 -6.5801E+01  7.3910E+01  3.1763E+02 -2.3564E+03 -1.1776E+03 -1.0123E+02 -2.3120E+03 -6.5245E+02
            -6.0248E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -972.282206904458        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.6166E+00  1.8559E+00  1.0597E+00  2.3483E+00  9.0243E-01  3.8374E+00  5.9197E+00  9.7657E-01  5.1257E+00  1.7003E+00
             1.1446E+01
 PARAMETER:  5.8035E-01  7.1837E-01  1.5796E-01  9.5369E-01 -2.6630E-03  1.4448E+00  1.8783E+00  7.6292E-02  1.7343E+00  6.3083E-01
             2.5377E+00
 GRADIENT:   4.7136E+01  2.1647E+01 -3.9280E+01  5.4887E+01 -5.9552E+01  1.2551E+02  1.2391E+02  4.5287E+00  1.3161E+02  3.3753E+01
             4.0708E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1033.04239560376        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      206
 NPARAMETR:  1.1350E+00  7.0225E+00  3.1841E+01  3.1078E+00  4.6222E+00  2.7275E+00  5.9499E+00  6.8932E-01  1.1008E+01  2.2036E+00
             1.1168E+01
 PARAMETER:  2.2660E-01  2.0491E+00  3.5607E+00  1.2339E+00  1.6309E+00  1.1034E+00  1.8834E+00 -2.7204E-01  2.4986E+00  8.9008E-01
             2.5131E+00
 GRADIENT:  -3.1964E+01  3.6360E+01 -3.0109E+00  3.0935E+01  4.7873E+01  4.7515E+01  4.5332E+01  7.4692E-02  2.1168E+01  4.4963E+01
             3.1540E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1196.08562212130        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      385
 NPARAMETR:  1.3925E+00  1.1846E+00  3.8042E+01  7.5978E-01  2.1171E+00  2.8385E+00  1.9617E+00  1.4535E+00  6.1691E+00  1.0233E+00
             9.4169E+00
 PARAMETER:  4.3111E-01  2.6943E-01  3.7387E+00 -1.7473E-01  8.5004E-01  1.1433E+00  7.7379E-01  4.7396E-01  1.9196E+00  1.2302E-01
             2.3425E+00
 GRADIENT:   1.4317E+01 -7.4874E+01 -9.3720E+00 -2.7915E+01 -1.0546E+01 -7.5845E+00  1.9428E+01  6.0730E-01  4.9905E+01  1.9123E+01
             2.5032E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1269.28160588829        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      561
 NPARAMETR:  1.1473E+00  1.6992E+00  1.1265E+02  7.8471E-01  2.3613E+00  2.2945E+00  1.2483E+00  3.7012E-01  6.0924E+00  3.9573E-01
             8.4571E+00
 PARAMETER:  2.3744E-01  6.3018E-01  4.8243E+00 -1.4244E-01  9.5920E-01  9.3050E-01  3.2181E-01 -8.9393E-01  1.9070E+00 -8.2703E-01
             2.2350E+00
 GRADIENT:  -1.0617E+01 -1.3275E+01 -1.2351E+00 -3.7429E+00 -6.7980E+00 -2.2628E+01  7.3909E+00  8.7051E-02  8.7581E+00  2.6559E+00
            -1.1066E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1272.40622608928        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      737
 NPARAMETR:  1.1884E+00  1.5225E+00  2.3392E+02  1.0192E+00  2.4550E+00  2.4555E+00  9.0428E-01  4.4774E-02  5.4708E+00  1.8741E-01
             8.5330E+00
 PARAMETER:  2.7259E-01  5.2034E-01  5.5550E+00  1.1902E-01  9.9813E-01  9.9835E-01 -6.1529E-04 -3.0061E+00  1.7994E+00 -1.5745E+00
             2.2439E+00
 GRADIENT:   4.1651E-01  6.2464E-01  5.1395E-01  4.7755E-01 -3.9865E-01  7.6815E-01  3.4427E-01  4.7411E-04 -4.1206E-01  5.7860E-01
            -2.6135E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1272.58668622422        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      919
 NPARAMETR:  1.1881E+00  1.5022E+00  1.9923E+02  1.0232E+00  2.4469E+00  2.4437E+00  8.3264E-01  3.7476E-02  5.4776E+00  1.4125E-01
             8.5598E+00
 PARAMETER:  2.7238E-01  5.0696E-01  5.3944E+00  1.2290E-01  9.9481E-01  9.9350E-01 -8.3156E-02 -3.1841E+00  1.8007E+00 -1.8572E+00
             2.2471E+00
 GRADIENT:   1.8187E-01  9.1200E-01  2.6688E-01  7.2399E-02 -6.8150E-01 -6.0460E-01 -7.0786E-01  4.2871E-04  1.7321E-01  3.2189E-01
             2.0753E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1272.61431560695        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1095
 NPARAMETR:  1.1857E+00  1.4794E+00  1.7106E+02  1.0500E+00  2.4422E+00  2.4483E+00  8.3026E-01  3.9152E-02  5.4151E+00  1.1358E-01
             8.5441E+00
 PARAMETER:  2.7033E-01  4.9163E-01  5.2420E+00  1.4880E-01  9.9288E-01  9.9538E-01 -8.6020E-02 -3.1403E+00  1.7892E+00 -2.0753E+00
             2.2452E+00
 GRADIENT:  -2.9168E-01  7.6066E-02 -9.5936E-02  3.8739E-01  2.8127E-01 -1.4305E-01 -3.0783E-01  5.8185E-04 -1.8417E-01  2.0715E-01
            -2.9036E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1272.63392911838        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1274
 NPARAMETR:  1.1868E+00  1.4655E+00  1.6091E+02  1.0573E+00  2.4366E+00  2.4546E+00  8.1260E-01  2.9706E-02  5.3899E+00  7.4549E-02
             8.5369E+00
 PARAMETER:  2.7130E-01  4.8220E-01  5.1809E+00  1.5571E-01  9.9059E-01  9.9798E-01 -1.0752E-01 -3.4164E+00  1.7845E+00 -2.4963E+00
             2.2444E+00
 GRADIENT:   1.4532E-01 -2.2272E-01 -2.5753E-01  1.8293E-01 -8.5919E-02  6.2456E-01 -4.0923E-01  3.6426E-04 -3.8503E-01  8.8380E-02
            -1.6658E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1272.84483148927        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1439
 NPARAMETR:  1.1864E+00  1.5245E+00  1.7739E+02  1.0086E+00  2.4438E+00  2.4532E+00  8.8291E-01  1.0000E-02  5.5798E+00  1.0000E-02
             8.5397E+00
 PARAMETER:  2.7092E-01  5.2165E-01  5.2784E+00  1.0859E-01  9.9357E-01  9.9741E-01 -2.4530E-02 -5.1224E+00  1.8192E+00 -4.6831E+00
             2.2447E+00
 GRADIENT:   1.9103E+01  1.2239E+01  1.1054E-01  2.5309E+00  4.0876E+00  2.9381E+01  2.2500E-01  0.0000E+00  1.0067E+02  0.0000E+00
             3.9324E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1272.86278393335        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1624
 NPARAMETR:  1.1885E+00  1.5373E+00  1.7579E+02  9.9261E-01  2.4431E+00  2.4553E+00  8.8520E-01  1.0000E-02  5.6009E+00  1.0000E-02
             8.5490E+00
 PARAMETER:  2.7267E-01  5.3003E-01  5.2693E+00  9.2585E-02  9.9326E-01  9.9823E-01 -2.1947E-02 -5.1224E+00  1.8229E+00 -4.6831E+00
             2.2458E+00
 GRADIENT:   4.6256E-01 -8.6692E-02  1.7109E-03  3.2269E-01 -1.5337E-01  9.3358E-01 -1.9885E-01  0.0000E+00  2.5504E+00  0.0000E+00
             3.3216E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1272.87886354012        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1804
 NPARAMETR:  1.1876E+00  1.5435E+00  1.7610E+02  9.8219E-01  2.4438E+00  2.4554E+00  9.0275E-01  1.0000E-02  5.6661E+00  1.0000E-02
             8.5464E+00
 PARAMETER:  2.7192E-01  5.3404E-01  5.2710E+00  8.2025E-02  9.9355E-01  9.9829E-01 -2.3114E-03 -5.1224E+00  1.8345E+00 -4.6831E+00
             2.2455E+00
 GRADIENT:   2.9298E-01 -2.5634E+00  1.5659E-02  2.2466E-01  4.1046E-01  9.1340E-01  3.7280E-01  0.0000E+00  5.0907E+00  0.0000E+00
             1.0558E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1272.89794417057        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1986             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1877E+00  1.5600E+00  1.7585E+02  9.6852E-01  2.4437E+00  2.4549E+00  9.1520E-01  1.0000E-02  5.7001E+00  1.0000E-02
             8.5471E+00
 PARAMETER:  2.7200E-01  5.4470E-01  5.2696E+00  6.8014E-02  9.9352E-01  9.9809E-01  1.1383E-02 -5.1224E+00  1.8405E+00 -4.6831E+00
             2.2456E+00
 GRADIENT:   1.9471E+01  1.2682E+01  1.0369E-01  1.9630E+00  4.2267E+00  2.9689E+01  4.1655E-01  0.0000E+00  1.0580E+02  0.0000E+00
             4.1034E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1272.90046814817        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     2159
 NPARAMETR:  1.1895E+00  1.5698E+00  1.7553E+02  9.6157E-01  2.4444E+00  2.4552E+00  9.1718E-01  1.0000E-02  5.6890E+00  1.0000E-02
             8.5534E+00
 PARAMETER:  2.7357E-01  5.5097E-01  5.2678E+00  6.0813E-02  9.9380E-01  9.9819E-01  1.3549E-02 -5.1224E+00  1.8385E+00 -4.6831E+00
             2.2463E+00
 GRADIENT:   7.0987E-01  1.8309E-01  1.1958E-03  1.8829E-01  4.1093E-02  1.0153E+00 -9.5302E-02  0.0000E+00  3.2509E+00  0.0000E+00
             1.0093E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1272.90832877719        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2339
 NPARAMETR:  1.1891E+00  1.5731E+00  1.7530E+02  9.5544E-01  2.4440E+00  2.4557E+00  9.1950E-01  1.0000E-02  5.7199E+00  1.0000E-02
             8.5508E+00
 PARAMETER:  2.7317E-01  5.5308E-01  5.2665E+00  5.4420E-02  9.9366E-01  9.9843E-01  1.6075E-02 -5.1224E+00  1.8439E+00 -4.6831E+00
             2.2460E+00
 GRADIENT:   6.2828E-01 -4.6708E-01  4.7177E-03  1.3219E-01  9.7286E-02  1.0739E+00 -2.4010E-02  0.0000E+00  4.1724E+00  0.0000E+00
             6.8453E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1272.91211986042        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2519
 NPARAMETR:  1.1876E+00  1.5762E+00  1.7513E+02  9.5081E-01  2.4438E+00  2.4549E+00  9.2449E-01  1.0000E-02  5.7405E+00  1.0000E-02
             8.5474E+00
 PARAMETER:  2.7192E-01  5.5502E-01  5.2655E+00  4.9564E-02  9.9354E-01  9.9807E-01  2.1486E-02 -5.1224E+00  1.8475E+00 -4.6831E+00
             2.2456E+00
 GRADIENT:   2.6701E-01 -9.5716E-01  7.1773E-03  9.0703E-02  1.4779E-01  9.3769E-01  8.7140E-02  0.0000E+00  4.7555E+00  0.0000E+00
             2.4702E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1272.91639437847        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2706             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1867E+00  1.5836E+00  1.7444E+02  9.4140E-01  2.4432E+00  2.4543E+00  9.3063E-01  1.0000E-02  5.7685E+00  1.0000E-02
             8.5449E+00
 PARAMETER:  2.7122E-01  5.5972E-01  5.2616E+00  3.9613E-02  9.9330E-01  9.9786E-01  2.8108E-02 -5.1224E+00  1.8524E+00 -4.6831E+00
             2.2453E+00
 GRADIENT:   1.9191E+01  1.3948E+01  9.2092E-02  1.6534E+00  4.1287E+00  2.9658E+01  2.7741E-01  0.0000E+00  1.0844E+02  0.0000E+00
             4.0160E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1272.91796710127        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2886
 NPARAMETR:  1.1877E+00  1.5877E+00  1.7441E+02  9.3896E-01  2.4435E+00  2.4552E+00  9.3286E-01  1.0000E-02  5.7698E+00  1.0000E-02
             8.5483E+00
 PARAMETER:  2.7201E-01  5.6231E-01  5.2614E+00  3.7018E-02  9.9343E-01  9.9822E-01  3.0500E-02 -5.1224E+00  1.8526E+00 -4.6831E+00
             2.2457E+00
 GRADIENT:   2.8382E-01 -5.9642E-01 -3.3063E-04 -3.3012E-03  7.1503E-02  1.0216E+00  2.7579E-02  0.0000E+00  4.7622E+00  0.0000E+00
             2.1195E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1272.91874790008        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3071
 NPARAMETR:  1.1874E+00  1.5902E+00  1.7429E+02  9.3490E-01  2.4436E+00  2.4549E+00  9.3677E-01  1.0000E-02  5.7810E+00  1.0000E-02
             8.5476E+00
 PARAMETER:  2.7173E-01  5.6386E-01  5.2607E+00  3.2685E-02  9.9347E-01  9.9807E-01  3.4682E-02 -5.1224E+00  1.8546E+00 -4.6831E+00
             2.2457E+00
 GRADIENT:   1.9855E-01 -9.2195E-01 -1.3392E-03 -9.8440E-02  1.7048E-01  9.7007E-01  1.1024E-01  0.0000E+00  4.9900E+00  0.0000E+00
             2.6190E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1272.91983291272        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3258
 NPARAMETR:  1.1873E+00  1.5930E+00  1.7427E+02  9.3396E-01  2.4434E+00  2.4549E+00  9.3751E-01  1.0000E-02  5.7860E+00  1.0000E-02
             8.5471E+00
 PARAMETER:  2.7171E-01  5.6560E-01  5.2606E+00  3.1677E-02  9.9337E-01  9.9807E-01  3.5476E-02 -5.1224E+00  1.8554E+00 -4.6831E+00
             2.2456E+00
 GRADIENT:   1.9857E-01 -5.1093E-01  2.1229E-04 -5.8435E-03  3.2598E-02  9.7430E-01  2.9232E-02  0.0000E+00  4.9145E+00  0.0000E+00
            -1.0134E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1272.92042681654        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     3452             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1873E+00  1.5961E+00  1.7419E+02  9.3255E-01  2.4432E+00  2.4549E+00  9.3775E-01  1.0000E-02  5.7921E+00  1.0000E-02
             8.5466E+00
 PARAMETER:  2.7168E-01  5.6755E-01  5.2601E+00  3.0163E-02  9.9329E-01  9.9807E-01  3.5729E-02 -5.1224E+00  1.8565E+00 -4.6831E+00
             2.2455E+00
 GRADIENT:   1.9335E+01  1.5347E+01  9.0225E-02  1.7335E+00  3.9171E+00  2.9756E+01  1.0542E-01  0.0000E+00  1.0883E+02  0.0000E+00
             3.9934E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1272.92082320343        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3637             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1874E+00  1.5976E+00  1.7428E+02  9.2994E-01  2.4435E+00  2.4546E+00  9.4123E-01  1.0000E-02  5.7984E+00  1.0000E-02
             8.5477E+00
 PARAMETER:  2.7180E-01  5.6850E-01  5.2606E+00  2.7360E-02  9.9345E-01  9.9798E-01  3.9429E-02 -5.1224E+00  1.8576E+00 -4.6831E+00
             2.2457E+00
 GRADIENT:   1.9370E+01  1.5072E+01  9.0257E-02  1.6465E+00  4.0549E+00  2.9725E+01  2.0300E-01  0.0000E+00  1.0923E+02  0.0000E+00
             4.0344E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1272.92104742608        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3822
 NPARAMETR:  1.1873E+00  1.5990E+00  1.7411E+02  9.2837E-01  2.4434E+00  2.4548E+00  9.4333E-01  1.0000E-02  5.8027E+00  1.0000E-02
             8.5473E+00
 PARAMETER:  2.7167E-01  5.6939E-01  5.2597E+00  2.5678E-02  9.9339E-01  9.9805E-01  4.1663E-02 -5.1224E+00  1.8583E+00 -4.6831E+00
             2.2456E+00
 GRADIENT:   1.7731E-01 -4.2219E-01 -7.7430E-04 -2.1039E-02  2.9410E-02  9.8280E-01  4.7770E-02  0.0000E+00  5.0295E+00  0.0000E+00
            -1.1581E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1272.92117426732        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     4016             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1873E+00  1.6011E+00  1.7412E+02  9.2801E-01  2.4432E+00  2.4548E+00  9.4209E-01  1.0000E-02  5.8048E+00  1.0000E-02
             8.5468E+00
 PARAMETER:  2.7169E-01  5.7070E-01  5.2597E+00  2.5293E-02  9.9331E-01  9.9803E-01  4.0350E-02 -5.1224E+00  1.8587E+00 -4.6831E+00
             2.2456E+00
 GRADIENT:   1.9334E+01  1.5658E+01  9.0206E-02  1.7180E+00  3.8975E+00  2.9750E+01  1.0020E-01  0.0000E+00  1.0926E+02  0.0000E+00
             3.9878E+01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1272.92139192396        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     4206             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1872E+00  1.6020E+00  1.7418E+02  9.2680E-01  2.4434E+00  2.4544E+00  9.4395E-01  1.0000E-02  5.8096E+00  1.0000E-02
             8.5467E+00
 PARAMETER:  2.7158E-01  5.7125E-01  5.2601E+00  2.3981E-02  9.9338E-01  9.9788E-01  4.2322E-02 -5.1224E+00  1.8595E+00 -4.6831E+00
             2.2455E+00
 GRADIENT:   1.9296E+01  1.5526E+01  9.1082E-02  1.6971E+00  3.9593E+00  2.9696E+01  1.4990E-01  0.0000E+00  1.0954E+02  0.0000E+00
             3.9931E+01

0ITERATION NO.:  124    OBJECTIVE VALUE:  -1272.92140961951        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     4344
 NPARAMETR:  1.1874E+00  1.6039E+00  1.7408E+02  9.2638E-01  2.4430E+00  2.4550E+00  9.4395E-01  1.0000E-02  5.8114E+00  1.0000E-02
             8.5471E+00
 PARAMETER:  2.7164E-01  5.7120E-01  5.2597E+00  2.3062E-02  9.9342E-01  9.9796E-01  4.3322E-02 -5.1224E+00  1.8596E+00 -4.6831E+00
             2.2456E+00
 GRADIENT:  -1.5300E-02 -1.6338E-01  1.9919E-04 -8.7571E-03  3.1401E-02 -1.7833E-02  9.0014E-03  0.0000E+00 -1.4320E-02  0.0000E+00
             1.3625E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4344
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1484E-02 -6.2797E-02  5.2140E-06  2.4853E-02 -4.4412E-06
 SE:             2.8905E-02  1.2016E-02  5.2816E-06  2.5479E-02  1.2549E-04
 N:                     100         100         100         100         100

 P VAL.:         6.9115E-01  1.7358E-07  3.2355E-01  3.2936E-01  9.7177E-01

 ETASHRINKSD(%)  3.1630E+00  5.9744E+01  9.9982E+01  1.4641E+01  9.9580E+01
 ETASHRINKVR(%)  6.2259E+00  8.3794E+01  1.0000E+02  2.7138E+01  9.9998E+01
 EBVSHRINKSD(%)  3.6692E+00  6.6858E+01  9.9953E+01  8.7956E+00  9.9497E+01
 EBVSHRINKVR(%)  7.2038E+00  8.9016E+01  1.0000E+02  1.6818E+01  9.9997E+01
 RELATIVEINF(%)  9.2670E+01  5.7523E+00  1.9637E-05  4.4139E+01  2.2290E-03
 EPSSHRINKSD(%)  6.7235E+00
 EPSSHRINKVR(%)  1.2995E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1272.9214096195083     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       381.16795014890249     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   155.55
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    16.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1272.921       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.19E+00  1.60E+00  1.74E+02  9.26E-01  2.44E+00  2.45E+00  9.45E-01  1.00E-02  5.81E+00  1.00E-02  8.55E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.48E+01
 
 TH 2
+       -1.48E+01  6.74E+01
 
 TH 3
+       -8.67E-04  2.75E-03  1.13E-07
 
 TH 4
+        3.92E+00  2.02E+01  7.96E-04  6.90E+00
 
 TH 5
+        3.26E+00 -1.77E+01 -7.19E-04 -5.36E+00  4.64E+00
 
 TH 6
+        7.25E+00 -6.81E+00 -2.94E-04 -1.44E+00  1.73E+00  1.11E+00
 
 TH 7
+        7.83E+00 -2.44E+01 -1.00E-03 -7.05E+00  6.38E+00  2.65E+00  8.92E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.51E+00 -7.59E+00 -3.15E-04 -2.08E+00  1.97E+00  8.97E-01  2.80E+00  0.00E+00  8.95E-01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.06E-01 -5.47E+00 -2.05E-04 -1.76E+00  1.40E+00  5.04E-01  1.96E+00  0.00E+00  5.82E-01  0.00E+00  8.36E-01
 
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
+        1.17E+02
 
 TH 2
+       -1.74E+00  7.42E+01
 
 TH 3
+       -4.21E-04 -2.19E-04  4.54E-05
 
 TH 4
+       -3.88E-01  2.37E+01  8.73E-04  3.21E+01
 
 TH 5
+       -9.29E-01 -1.02E+01 -1.43E-02 -4.68E+00  4.42E+01
 
 TH 6
+       -6.41E-01  9.26E-01 -1.67E-04  3.73E-02 -9.50E-02  2.69E+01
 
 TH 7
+        2.37E-02 -2.46E+01  1.25E-05 -8.72E-01  3.14E+00 -3.11E-01  1.63E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.33E-01 -8.57E+00  6.65E-04  2.87E+00  7.18E-01 -7.38E-02  2.31E+00  0.00E+00  3.64E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.42E+00 -7.73E+00  2.84E-05 -2.68E+00  5.06E-01  1.63E+00  3.41E+00  0.00E+00  5.25E-01  0.00E+00  1.40E+01
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.22E+02
 
 TH 2
+        4.23E+01  8.15E+01
 
 TH 3
+       -2.89E-03 -2.52E-04  1.10E-05
 
 TH 4
+        4.87E+01  2.48E+01  1.90E-03  3.25E+01
 
 TH 5
+       -1.02E+01 -1.31E+01 -6.32E-03 -3.54E+00  3.22E+01
 
 TH 6
+        1.45E+01 -1.46E+01 -1.98E-03  5.93E+00  5.08E+00  2.37E+01
 
 TH 7
+       -6.86E+00 -2.84E+01  2.73E-05 -3.44E+00  5.56E+00  5.91E+00  1.84E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.84E-01 -9.18E+00  1.66E-03  2.71E+00  1.07E+00  3.89E+00  2.66E+00  0.00E+00  3.69E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.41E+01 -7.71E+01 -3.28E-04 -2.23E+01  1.52E+01  1.99E+01  2.68E+01  0.00E+00  1.34E+01  0.00E+00  4.90E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,      171.862
Stop Time:
Wed Sep 29 09:44:30 CDT 2021
