Thu Sep 30 08:35:38 CDT 2021
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
$DATA ../../../../data/spa2/D/dat12.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m12.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1911.49781721057        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.5389E+02 -1.0042E+02 -6.1946E+01 -1.4716E+02  1.7172E+02 -2.9899E+02 -2.5117E+02 -1.7771E+01 -4.1967E+02 -9.6986E+01
            -8.0712E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2040.70606564081        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:       92
 NPARAMETR:  9.4715E-01  1.0353E+00  1.0273E+00  1.0369E+00  9.6529E-01  1.1041E+00  1.2032E+00  1.0218E+00  1.1553E+00  1.0560E+00
             1.0209E+00
 PARAMETER:  4.5706E-02  1.3466E-01  1.2691E-01  1.3621E-01  6.4671E-02  1.9906E-01  2.8501E-01  1.2156E-01  2.4433E-01  1.5453E-01
             1.2070E-01
 GRADIENT:   1.6958E+02 -4.0254E+01 -2.3639E+01  3.6815E+01  1.0450E+01 -8.0720E+01 -2.6212E+02 -7.4959E+00 -1.3225E+02 -4.5481E+00
            -6.1772E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2183.67578029019        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.2090E+00  1.1231E+00  1.0961E+00  8.0374E-01  1.1209E+00  1.3597E+00  1.7585E+00  1.0573E+00  1.1009E+00  1.0711E+00
             1.0418E+00
 PARAMETER:  2.8979E-01  2.1609E-01  1.9180E-01 -1.1848E-01  2.1417E-01  4.0729E-01  6.6445E-01  1.5573E-01  1.9614E-01  1.6872E-01
             1.4094E-01
 GRADIENT:   1.4540E+02 -1.5304E+02  1.1367E+01 -1.7464E+02 -1.1081E+02 -1.4004E+02 -1.9603E+02 -1.2844E+00 -8.6856E+01  8.0735E+00
            -6.1318E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2235.97689300638        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  1.1656E+00  1.5118E+00  9.6268E-01  6.8032E-01  1.3534E+00  1.4921E+00  1.8347E+00  8.5270E-01  1.2605E+00  9.6300E-01
             1.0793E+00
 PARAMETER:  2.5324E-01  5.1330E-01  6.1967E-02 -2.8519E-01  4.0261E-01  5.0016E-01  7.0690E-01 -5.9353E-02  3.3149E-01  6.2298E-02
             1.7635E-01
 GRADIENT:   8.1504E+01 -7.1539E+01 -1.3795E+01 -5.5026E+01  2.5968E+01 -7.3580E+01 -9.3431E+01 -4.4476E+00 -5.1377E+01 -3.7279E+01
            -4.0278E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2271.46361811851        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  1.0636E+00  1.1852E+00  1.4373E+00  9.6457E-01  1.2297E+00  1.5813E+00  2.4898E+00  1.3433E+00  1.2192E+00  9.7326E-01
             1.0751E+00
 PARAMETER:  1.6170E-01  2.6992E-01  4.6280E-01  6.3927E-02  3.0676E-01  5.5822E-01  1.0122E+00  3.9515E-01  2.9818E-01  7.2892E-02
             1.7241E-01
 GRADIENT:  -1.1821E+01 -1.4842E+01 -1.1725E+01 -2.4949E+01  6.2896E+00 -3.4833E+01 -2.4480E+01  1.8978E+00 -3.2108E+01 -6.6419E+00
            -1.3743E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2275.99028217498        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  1.0654E+00  1.1211E+00  1.5627E+00  1.0231E+00  1.2161E+00  1.6190E+00  2.6382E+00  1.3900E+00  1.2538E+00  9.6986E-01
             1.0760E+00
 PARAMETER:  1.6335E-01  2.1431E-01  5.4641E-01  1.2283E-01  2.9565E-01  5.8182E-01  1.0701E+00  4.2928E-01  3.2616E-01  6.9394E-02
             1.7328E-01
 GRADIENT:  -9.3179E+00 -7.7963E+00 -7.4250E+00 -1.5966E+01  2.3178E+00 -2.2932E+01 -1.6479E+01  1.3098E+00 -2.2511E+01 -4.6567E+00
            -9.3502E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2277.31093773288        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      950
 NPARAMETR:  1.0659E+00  1.1309E+00  1.5750E+00  1.0288E+00  1.2156E+00  1.6792E+00  2.6537E+00  1.3955E+00  1.2582E+00  9.7033E-01
             1.0804E+00
 PARAMETER:  1.6379E-01  2.2305E-01  5.5429E-01  1.2841E-01  2.9521E-01  6.1830E-01  1.0759E+00  4.3322E-01  3.2970E-01  6.9883E-02
             1.7729E-01
 GRADIENT:  -8.0554E+00 -3.8309E+00 -6.3757E+00 -1.2244E+01 -5.6103E-02 -5.3813E+00 -1.4791E+01  1.2848E+00 -2.1486E+01 -4.4607E+00
            -4.7907E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2279.36743844192        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1135
 NPARAMETR:  1.0672E+00  1.1290E+00  1.5817E+00  1.0278E+00  1.2097E+00  1.7562E+00  2.8029E+00  1.3909E+00  1.3990E+00  9.7107E-01
             1.0837E+00
 PARAMETER:  1.6503E-01  2.2137E-01  5.5849E-01  1.2744E-01  2.9038E-01  6.6312E-01  1.1307E+00  4.2995E-01  4.3577E-01  7.0639E-02
             1.8039E-01
 GRADIENT:  -6.0069E+00 -1.6899E+00 -7.3373E+00 -1.3788E+01 -1.7052E+00  1.4630E+01  4.4561E+00  1.8305E+00  3.6532E-01 -2.1079E+00
             7.1073E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2279.62249306278        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1314
 NPARAMETR:  1.0875E+00  1.1252E+00  1.5795E+00  1.0536E+00  1.2060E+00  1.7320E+00  2.8285E+00  1.3463E+00  1.3888E+00  9.8812E-01
             1.0838E+00
 PARAMETER:  1.8389E-01  2.1796E-01  5.5708E-01  1.5219E-01  2.8731E-01  6.4927E-01  1.1398E+00  3.9737E-01  4.2842E-01  8.8046E-02
             1.8043E-01
 GRADIENT:   8.5946E+00  2.1236E+00 -1.0523E+01  2.4908E-01  4.1185E+00  7.6941E+00  4.5574E+00  5.9768E-02  7.8213E-01  2.7490E-01
             5.5464E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2279.79441885851        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1493             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0819E+00  1.0728E+00  1.7016E+00  1.0724E+00  1.1850E+00  1.7716E+00  2.9061E+00  1.3311E+00  1.3742E+00  9.7464E-01
             1.0836E+00
 PARAMETER:  1.7869E-01  1.7029E-01  6.3157E-01  1.6988E-01  2.6973E-01  6.7189E-01  1.1668E+00  3.8601E-01  4.1790E-01  7.4310E-02
             1.8032E-01
 GRADIENT:   6.6774E+02  8.1823E+01  2.7452E+01  1.2966E+02  1.5230E+01  6.3931E+02  4.9535E+02 -2.2111E+00  5.3804E+01  6.0141E-01
             7.9883E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2280.76230108495        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1668
 NPARAMETR:  1.0808E+00  9.8519E-01  1.8555E+00  1.1437E+00  1.2164E+00  1.7097E+00  3.0212E+00  1.5227E+00  1.3416E+00  1.0009E+00
             1.0821E+00
 PARAMETER:  1.7774E-01  8.5084E-02  7.1818E-01  2.3425E-01  2.9587E-01  6.3631E-01  1.2057E+00  5.2050E-01  3.9384E-01  1.0091E-01
             1.7886E-01
 GRADIENT:   4.1908E+00 -4.6466E-01 -4.5068E+00  4.6269E+00  5.0562E+00  2.2532E+00  9.1653E-01  5.4675E-01  1.4135E+00  2.6197E-01
             1.1891E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2280.87102769268        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1855
 NPARAMETR:  1.0801E+00  9.8723E-01  1.8721E+00  1.1431E+00  1.2097E+00  1.7204E+00  3.0450E+00  1.5103E+00  1.3389E+00  1.0013E+00
             1.0828E+00
 PARAMETER:  1.7701E-01  8.7144E-02  7.2708E-01  2.3373E-01  2.9038E-01  6.4254E-01  1.2135E+00  5.1229E-01  3.9183E-01  1.0133E-01
             1.7956E-01
 GRADIENT:   3.6479E+00  1.5880E+00  8.7557E-02  4.1002E-01 -4.6365E+00  5.0744E+00  2.2921E+00 -5.0473E-01  1.0709E+00  7.7128E-01
             1.6360E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2280.99675614196        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     2016             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0812E+00  9.7160E-01  1.8795E+00  1.1439E+00  1.2137E+00  1.7472E+00  3.0964E+00  1.5103E+00  1.3294E+00  9.9607E-01
             1.0812E+00
 PARAMETER:  1.7805E-01  7.1185E-02  7.3102E-01  2.3448E-01  2.9366E-01  6.5803E-01  1.2302E+00  5.1233E-01  3.8472E-01  9.6061E-02
             1.7803E-01
 GRADIENT:   6.6825E+02  4.0789E+01  1.9286E+01  2.2400E+02  4.6876E+01  6.2337E+02  5.4885E+02  2.5969E+00  5.2440E+01  1.0875E+00
             1.7645E+00

0ITERATION NO.:   63    OBJECTIVE VALUE:  -2280.99872572976        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     2111
 NPARAMETR:  1.0817E+00  9.7257E-01  1.8795E+00  1.1462E+00  1.2134E+00  1.7437E+00  3.0989E+00  1.5103E+00  1.3290E+00  9.9692E-01
             1.0813E+00
 PARAMETER:  1.7857E-01  7.2185E-02  7.3102E-01  2.3643E-01  2.9342E-01  6.5602E-01  1.2311E+00  5.1233E-01  3.8441E-01  9.6911E-02
             1.7814E-01
 GRADIENT:   3.2507E-01 -2.5448E+00  5.0352E+01 -2.9458E+01  3.8586E-01 -1.2979E+04  6.9156E+03  1.6484E+04  3.1918E-01  1.3769E-01
             4.7785E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2111
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1788E-04  1.3852E-02 -4.9910E-02 -1.5568E-02 -2.4085E-02
 SE:             3.0038E-02  2.4771E-02  1.6490E-02  2.3416E-02  2.1211E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9687E-01  5.7602E-01  2.4726E-03  5.0615E-01  2.5617E-01

 ETASHRINKSD(%)  1.0000E-10  1.7015E+01  4.4757E+01  2.1553E+01  2.8941E+01
 ETASHRINKVR(%)  1.0000E-10  3.1135E+01  6.9482E+01  3.8461E+01  4.9507E+01
 EBVSHRINKSD(%)  1.2784E-01  1.5112E+01  4.9392E+01  2.3640E+01  2.6943E+01
 EBVSHRINKVR(%)  2.5552E-01  2.7941E+01  7.4388E+01  4.1691E+01  4.6627E+01
 RELATIVEINF(%)  9.9731E+01  3.4183E+01  1.2042E+01  2.4940E+01  2.6598E+01
 EPSSHRINKSD(%)  3.0442E+01
 EPSSHRINKVR(%)  5.1617E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2280.9987257297580     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1178.2724858841509     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    44.86
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.32
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2280.999       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.08E+00  9.73E-01  1.88E+00  1.15E+00  1.21E+00  1.74E+00  3.10E+00  1.51E+00  1.33E+00  9.97E-01  1.08E+00
 


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
+        3.11E+02
 
 TH 2
+        1.13E+07  2.25E+07
 
 TH 3
+        4.55E-01  1.43E+01  1.13E+05
 
 TH 4
+        4.07E+06 -8.07E+06 -5.73E+05  5.80E+06
 
 TH 5
+       -6.19E+06 -6.15E+06  1.87E+03  2.21E+06  4.10E+02
 
 TH 6
+       -9.63E+05 -1.91E+06  1.10E-01  6.87E+05 -2.64E-01  1.63E+05
 
 TH 7
+       -8.07E+01  5.74E+05 -3.21E+00  4.70E+01 -1.57E+05  4.88E+04  1.22E+01
 
 TH 8
+       -2.70E-01 -2.81E+06  2.02E+05 -1.02E+06 -7.66E+05  7.76E+00 -1.44E+05  1.74E+01
 
 TH 9
+       -2.15E+06 -1.54E+02 -6.92E+02 -1.25E+03  7.07E+00 -4.61E-02  4.94E+00  1.11E+00  1.63E+06
 
 TH10
+        6.90E-02 -2.20E+07 -1.56E+06  7.88E+06 -6.00E+06 -6.06E+01  5.60E+05  5.52E+06  4.17E+06  4.28E+07
 
 TH11
+       -5.72E+06  3.87E+02 -8.06E+05  1.69E+00  3.10E+06  9.66E+05 -8.16E+01  1.42E+06  2.16E+06  1.11E+07  1.15E+07
 
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
 #CPUT: Total CPU Time in Seconds,       55.250
Stop Time:
Thu Sep 30 08:36:35 CDT 2021
