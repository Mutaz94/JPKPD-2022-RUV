Wed Sep 29 20:03:51 CDT 2021
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
$DATA ../../../../data/spa/D/dat46.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12566.7436638266        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.2134E+02  2.4941E+02 -1.6048E+01  1.5307E+02  1.0701E+02 -1.5658E+03 -7.2465E+02 -4.8318E+01 -1.1156E+03 -4.1152E+02
            -2.4196E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -569.673451494742        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.3009E+00  1.3212E+00  9.9517E-01  1.8613E+00  1.1726E+00  1.7042E+00  1.2720E+00  9.8719E-01  1.3056E+00  1.0755E+00
             1.4733E+01
 PARAMETER:  3.6304E-01  3.7851E-01  9.5158E-02  7.2126E-01  2.5921E-01  6.3311E-01  3.4057E-01  8.7112E-02  3.6667E-01  1.7275E-01
             2.7901E+00
 GRADIENT:  -1.5440E+01  5.1830E+01 -5.8787E-01  9.9856E+01 -1.6516E+01  1.6299E+01 -2.2497E+00  3.9313E+00  3.2287E+00  3.0499E+00
             1.1753E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -590.925000518644        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2270E+00  9.0920E-01  1.6281E+00  1.8143E+00  2.0744E+00  1.4753E+00  3.5823E+00  4.9393E-01  8.6291E-01  4.8823E+00
             1.3409E+01
 PARAMETER:  3.0453E-01  4.8148E-03  5.8740E-01  6.9569E-01  8.2967E-01  4.8886E-01  1.3760E+00 -6.0537E-01 -4.7441E-02  1.6856E+00
             2.6959E+00
 GRADIENT:   3.0212E+00  3.0523E+01  5.6546E+00  5.4842E+01 -2.3969E+01 -1.0836E+01  1.2355E+01  2.1993E-01  3.3917E+00  2.2576E+01
             1.0382E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -626.022085618309        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0066E+00  2.6038E-01  1.3949E+00  1.7350E+00  2.5432E+00  1.4634E+00  2.8963E+00  4.6356E-01  5.9887E-01  2.2773E+00
             1.0788E+01
 PARAMETER:  1.0659E-01 -1.2456E+00  4.3285E-01  6.5099E-01  1.0334E+00  4.8073E-01  1.1634E+00 -6.6883E-01 -4.1270E-01  9.2301E-01
             2.4784E+00
 GRADIENT:  -5.1993E+01  8.5113E+00  7.4361E+00  7.0122E+01 -5.9253E+00  1.4081E+00  2.7033E-01  5.5568E-02 -2.3261E+00  4.0434E-01
             4.2574E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -655.663943555180        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  8.6166E-01  1.5900E-01  1.7632E-01  1.1343E+00  9.1294E+00  1.3117E+00  6.4332E-01  1.0000E-02  1.5096E-01  2.2850E+00
             9.4570E+00
 PARAMETER: -4.8898E-02 -1.7389E+00 -1.6355E+00  2.2603E-01  2.3115E+00  3.7135E-01 -3.4112E-01 -4.6884E+00 -1.7908E+00  9.2635E-01
             2.3468E+00
 GRADIENT:   5.9417E+01  6.4037E+01 -9.8369E+01  1.4437E+02 -1.8473E+01 -6.7379E+01  6.6273E+00  0.0000E+00  3.4481E-01  7.3387E+00
            -1.4246E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -697.804437156030        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  5.5881E-01  4.4178E-02  5.9390E-02  5.8740E-01  5.5549E+01  1.4154E+00  6.8212E-02  1.0000E-02  1.0000E-02  6.4056E+00
             9.9163E+00
 PARAMETER: -4.8195E-01 -3.0195E+00 -2.7236E+00 -4.3205E-01  4.1173E+00  4.4741E-01 -2.5851E+00 -1.0273E+01 -5.2394E+00  1.9572E+00
             2.3942E+00
 GRADIENT:  -2.0390E+01 -5.5918E+00 -1.8747E+01  1.1408E+02  1.2732E-01 -3.9353E+00  1.3913E-02  0.0000E+00  0.0000E+00  2.5235E-02
            -1.0391E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -703.094317413066        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  5.1115E-01  2.7511E-02  4.1338E-02  4.3253E-01  1.2746E+02  1.3752E+00  2.9290E-02  1.0000E-02  1.0000E-02  1.1313E+01
             1.0013E+01
 PARAMETER: -5.7109E-01 -3.4932E+00 -3.0860E+00 -7.3810E-01  4.9478E+00  4.1857E-01 -3.4305E+00 -1.1923E+01 -6.6500E+00  2.5260E+00
             2.4039E+00
 GRADIENT:  -5.2510E+00  1.8648E+00  9.4471E-01 -1.2542E+00 -1.8058E-02 -2.9468E+00  5.8192E-04  0.0000E+00  0.0000E+00  3.6678E-03
            -1.4538E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -703.531151084983        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      675
 NPARAMETR:  5.0073E-01  1.4334E-02  3.8674E-02  4.1390E-01  1.0322E+04  1.3790E+00  1.0000E-02  1.0000E-02  1.0000E-02  6.0422E-01
             1.0169E+01
 PARAMETER: -5.9168E-01 -4.1451E+00 -3.1526E+00 -7.8214E-01  9.3420E+00  4.2136E-01 -6.4337E+00 -1.2418E+01 -6.9488E+00 -4.0381E-01
             2.4193E+00
 GRADIENT:   3.9899E+01  9.1829E-03  3.1311E+01  2.6803E+01  4.6550E-05  2.5459E+00  0.0000E+00  0.0000E+00  0.0000E+00 -1.4077E-10
             1.4939E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -703.564568501475        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:      871             RESET HESSIAN, TYPE I
 NPARAMETR:  5.0078E-01  1.4874E-02  3.8787E-02  4.1314E-01  5.2657E+03  1.3845E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.3089E-01
             1.0200E+01
 PARAMETER: -5.9160E-01 -4.1082E+00 -3.1497E+00 -7.8396E-01  8.6690E+00  4.2537E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -8.5259E-02
             2.4224E+00
 GRADIENT:   3.8467E+01  2.4615E-02  3.7750E+01  1.7591E+01  4.5530E-05  4.3682E+00  0.0000E+00  0.0000E+00  0.0000E+00 -9.0949E-10
             1.8902E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -703.565456562818        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  5.0067E-01  1.4849E-02  3.8776E-02  4.1284E-01  3.1546E+03  1.3846E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.2923E-01
             1.0203E+01
 PARAMETER: -5.9182E-01 -4.1098E+00 -3.1499E+00 -7.8470E-01  8.1566E+00  4.2539E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -8.7259E-02
             2.4227E+00
 GRADIENT:   5.2541E-01 -6.2198E-04 -7.5523E-01  1.8878E-01  6.9213E-05  4.3698E-02  0.0000E+00  0.0000E+00  0.0000E+00 -2.8990E-09
            -9.8566E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -703.567255907061        NO. OF FUNC. EVALS.: 200
 CUMULATIVE NO. OF FUNC. EVALS.:     1257
 NPARAMETR:  5.0038E-01  1.4847E-02  3.8707E-02  4.1228E-01  1.3135E+03  1.3845E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.2509E-01
             1.0203E+01
 PARAMETER: -5.9239E-01 -4.1100E+00 -3.1517E+00 -7.8605E-01  7.2805E+00  4.2534E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -9.2259E-02
             2.4227E+00
 GRADIENT:   5.7660E-01 -7.9271E-04 -6.0057E-01 -9.3358E-02  1.6277E-04  6.3839E-02  0.0000E+00  0.0000E+00  0.0000E+00 -1.6371E-08
            -7.2898E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -703.568376387999        NO. OF FUNC. EVALS.: 207
 CUMULATIVE NO. OF FUNC. EVALS.:     1464             RESET HESSIAN, TYPE I
 NPARAMETR:  5.0019E-01  1.4869E-02  3.8652E-02  4.1205E-01  8.5603E+02  1.3843E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.2596E-01
             1.0202E+01
 PARAMETER: -5.9277E-01 -4.1085E+00 -3.1532E+00 -7.8661E-01  6.8523E+00  4.2520E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -9.1208E-02
             2.4226E+00
 GRADIENT:   3.8544E+01  2.5466E-02  3.8231E+01  1.7035E+01  2.6557E-04  4.3895E+00  0.0000E+00  0.0000E+00  0.0000E+00 -3.7630E-08
             1.9122E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -703.569115527769        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1662
 NPARAMETR:  5.0007E-01  1.4852E-02  3.8623E-02  4.1184E-01  6.5259E+02  1.3843E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.2596E-01
             1.0202E+01
 PARAMETER: -5.9301E-01 -4.1096E+00 -3.1539E+00 -7.8711E-01  6.5810E+00  4.2516E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -9.1208E-02
             2.4226E+00
 GRADIENT:   6.2919E-01  9.8783E-04 -1.3359E+00  9.0430E-01  3.6063E-04 -4.6731E-03  0.0000E+00  0.0000E+00  0.0000E+00 -6.6564E-08
            -2.7778E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -703.570207998391        NO. OF FUNC. EVALS.: 211
 CUMULATIVE NO. OF FUNC. EVALS.:     1873             RESET HESSIAN, TYPE I
 NPARAMETR:  4.9985E-01  1.4827E-02  3.8595E-02  4.1143E-01  4.4317E+02  1.3844E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.2431E-01
             1.0203E+01
 PARAMETER: -5.9346E-01 -4.1113E+00 -3.1546E+00 -7.8811E-01  6.1939E+00  4.2525E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -9.3208E-02
             2.4227E+00
 GRADIENT:   3.8433E+01  2.5236E-02  3.9058E+01  1.5966E+01  4.4708E-04  4.4714E+00  0.0000E+00  0.0000E+00  0.0000E+00 -1.3858E-07
             1.9403E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -703.570964859802        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2069
 NPARAMETR:  4.9973E-01  1.4811E-02  3.8567E-02  4.1122E-01  3.4677E+02  1.3843E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.2273E-01
             1.0203E+01
 PARAMETER: -5.9370E-01 -4.1124E+00 -3.1554E+00 -7.8862E-01  5.9487E+00  4.2523E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -9.5122E-02
             2.4227E+00
 GRADIENT:   4.3166E-01  8.9868E-04 -4.4615E-01 -3.1171E-01  5.7685E-04  8.4440E-02  0.0000E+00  0.0000E+00  0.0000E+00 -2.3181E-07
             4.1647E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -703.572168195873        NO. OF FUNC. EVALS.: 207
 CUMULATIVE NO. OF FUNC. EVALS.:     2276             RESET HESSIAN, TYPE I
 NPARAMETR:  4.9961E-01  1.4785E-02  3.8508E-02  4.1105E-01  2.4441E+02  1.3841E+00  1.0000E-02  1.0000E-02  1.0000E-02  7.9463E-01
             1.0201E+01
 PARAMETER: -5.9393E-01 -4.1141E+00 -3.1569E+00 -7.8904E-01  5.5988E+00  4.2503E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -1.2988E-01
             2.4225E+00
 GRADIENT:   3.8675E+01  2.5265E-02  3.8153E+01  1.7313E+01  9.6515E-04  4.3646E+00  0.0000E+00  0.0000E+00  0.0000E+00 -4.2130E-07
             1.9043E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -703.573025433833        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2473             RESET HESSIAN, TYPE I
 NPARAMETR:  4.9945E-01  1.4762E-02  3.8499E-02  4.1076E-01  1.8555E+02  1.3843E+00  1.0000E-02  1.0000E-02  1.0000E-02  7.9885E-01
             1.0203E+01
 PARAMETER: -5.9424E-01 -4.1157E+00 -3.1571E+00 -7.8974E-01  5.3233E+00  4.2518E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -1.2458E-01
             2.4227E+00
 GRADIENT:   3.8448E+01  2.5735E-02  3.9103E+01  1.6047E+01  1.0702E-03  4.4678E+00  0.0000E+00  0.0000E+00  0.0000E+00 -7.3340E-07
             1.9397E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -703.573885647672        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2669
 NPARAMETR:  4.9936E-01  1.4737E-02  3.8472E-02  4.1056E-01  1.5029E+02  1.3843E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.0029E-01
             1.0203E+01
 PARAMETER: -5.9442E-01 -4.1174E+00 -3.1578E+00 -7.9023E-01  5.1125E+00  4.2517E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -1.2278E-01
             2.4227E+00
 GRADIENT:   4.1750E-01  1.1649E-03 -5.0360E-01 -2.6772E-01  1.3196E-03  7.8377E-02  0.0000E+00  0.0000E+00  0.0000E+00 -1.1628E-06
             2.8010E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -703.575282487858        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:     2874             RESET HESSIAN, TYPE I
 NPARAMETR:  4.9927E-01  1.4693E-02  3.8415E-02  4.1039E-01  1.1126E+02  1.3840E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.0047E-01
             1.0201E+01
 PARAMETER: -5.9461E-01 -4.1204E+00 -3.1593E+00 -7.9065E-01  4.8119E+00  4.2500E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -1.2255E-01
             2.4225E+00
 GRADIENT:   3.8698E+01  2.5322E-02  3.8299E+01  1.7244E+01  2.1058E-03  4.3676E+00  0.0000E+00  0.0000E+00  0.0000E+00 -2.0574E-06
             1.9061E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -703.576282390077        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3071             RESET HESSIAN, TYPE I
 NPARAMETR:  4.9913E-01  1.4655E-02  3.8406E-02  4.1012E-01  8.7805E+01  1.3843E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.0398E-01
             1.0203E+01
 PARAMETER: -5.9488E-01 -4.1230E+00 -3.1595E+00 -7.9131E-01  4.5751E+00  4.2517E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -1.1817E-01
             2.4226E+00
 GRADIENT:   3.8466E+01  2.5830E-02  3.9239E+01  1.5999E+01  2.2464E-03  4.4688E+00  0.0000E+00  0.0000E+00  0.0000E+00 -3.3065E-06
             1.9410E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -703.577310399422        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     3270
 NPARAMETR:  4.9906E-01  1.4616E-02  3.8379E-02  4.0993E-01  7.3255E+01  1.3843E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.0686E-01
             1.0203E+01
 PARAMETER: -5.9503E-01 -4.1257E+00 -3.1602E+00 -7.9178E-01  4.3939E+00  4.2517E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -1.1461E-01
             2.4226E+00
 GRADIENT:   3.8749E-01  1.4022E-03 -4.6969E-01 -3.4916E-01  2.6218E-03  7.9395E-02  0.0000E+00  0.0000E+00  0.0000E+00 -4.9476E-06
             4.2756E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -703.578559625189        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     3468
 NPARAMETR:  4.9906E-01  1.4514E-02  3.8340E-02  4.0988E-01  5.9012E+01  1.3841E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.0878E-01
             1.0200E+01
 PARAMETER: -5.9503E-01 -4.1326E+00 -3.1613E+00 -7.9189E-01  4.1777E+00  4.2507E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -1.1222E-01
             2.4224E+00
 GRADIENT:   6.2180E-01 -4.2489E-05 -1.3914E+00  8.4746E-01  3.8566E-03 -1.2907E-02  0.0000E+00  0.0000E+00  0.0000E+00 -7.7085E-06
            -3.2920E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -703.579941211812        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3665             RESET HESSIAN, TYPE I
 NPARAMETR:  4.9914E-01  1.4392E-02  3.8290E-02  4.0976E-01  4.2615E+01  1.3842E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.1338E-01
             1.0198E+01
 PARAMETER: -5.9486E-01 -4.1411E+00 -3.1626E+00 -7.9217E-01  3.8522E+00  4.2514E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -1.0656E-01
             2.4222E+00
 GRADIENT:   3.9059E+01  2.1635E-02  3.7728E+01  1.8151E+01  6.3991E-03  4.3373E+00  0.0000E+00  0.0000E+00  0.0000E+00 -1.4514E-05
             1.8735E+01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -703.581047581558        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3840
 NPARAMETR:  4.9889E-01  1.4174E-02  3.8358E-02  4.0981E-01  3.0133E+01  1.3842E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.1598E-01
             1.0201E+01
 PARAMETER: -5.9537E-01 -4.1564E+00 -3.1608E+00 -7.9205E-01  3.5056E+00  4.2515E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -1.0336E-01
             2.4225E+00
 GRADIENT:  -4.2343E-01  3.1467E-03  2.1442E-01 -9.4136E-01  4.4305E-03 -2.3045E-02  0.0000E+00  0.0000E+00  0.0000E+00 -2.9832E-05
             1.5899E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -703.601285402804        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     4019
 NPARAMETR:  4.9965E-01  1.2677E-02  3.8502E-02  4.1091E-01  8.0210E+00  1.3855E+00  1.0000E-02  1.0000E-02  1.0000E-02  8.3802E-01
             1.0199E+01
 PARAMETER: -5.9384E-01 -4.2679E+00 -3.1570E+00 -7.8939E-01  2.1821E+00  4.2605E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -7.6714E-02
             2.4223E+00
 GRADIENT:  -4.2791E+00  4.3895E-02  7.5718E+00 -9.6215E+00 -4.8214E-02 -1.9384E-01  0.0000E+00  0.0000E+00  0.0000E+00  3.1810E-02
             1.9973E+00

0ITERATION NO.:  125    OBJECTIVE VALUE:  -703.754740827640        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     4197
 NPARAMETR:  4.9942E-01  1.0500E-02  3.8386E-02  4.1347E-01  6.5882E+00  1.3917E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.3475E-02
             1.0186E+01
 PARAMETER: -5.9431E-01 -4.4564E+00 -3.1601E+00 -7.8316E-01  1.9853E+00  4.3053E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -3.2970E+00
             2.4210E+00
 GRADIENT:  -6.5609E+00  1.2675E-02  6.1249E+00 -6.2301E+00  4.2301E-02  3.6075E-02  0.0000E+00  0.0000E+00  0.0000E+00  9.9791E-04
             1.5318E+00

0ITERATION NO.:  128    OBJECTIVE VALUE:  -703.814216682815        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     4289
 NPARAMETR:  5.0313E-01  1.1266E-02  3.8243E-02  4.1528E-01  6.2999E+00  1.3943E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.0187E+01
 PARAMETER: -5.8690E-01 -4.3859E+00 -3.1638E+00 -7.7879E-01  1.9405E+00  4.3240E-01 -6.1337E+00 -1.2418E+01 -6.9488E+00 -5.7412E+00
             2.4211E+00
 GRADIENT:  -4.2480E-01  1.5205E-02 -1.8516E-01 -3.1900E-01 -9.7275E-02  1.4294E-01  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -3.2004E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     4289
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3037E-03  2.5465E-06  1.1209E-04 -2.5538E-04  3.4085E-06
 SE:             2.8573E-02  3.1251E-06  2.6876E-04  3.6390E-04  3.5573E-05
 N:                     100         100         100         100         100

 P VAL.:         9.0795E-01  4.1516E-01  6.7665E-01  4.8282E-01  9.2367E-01

 ETASHRINKSD(%)  4.2758E+00  9.9990E+01  9.9100E+01  9.8781E+01  9.9881E+01
 ETASHRINKVR(%)  8.3688E+00  1.0000E+02  9.9992E+01  9.9985E+01  1.0000E+02
 EBVSHRINKSD(%)  4.3403E+00  9.9987E+01  9.9090E+01  9.8763E+01  9.9851E+01
 EBVSHRINKVR(%)  8.4922E+00  1.0000E+02  9.9992E+01  9.9985E+01  1.0000E+02
 RELATIVEINF(%)  5.7206E+00  1.3118E-07  5.9328E-05  1.1191E-04  1.7656E-05
 EPSSHRINKSD(%)  7.1012E+00
 EPSSHRINKVR(%)  1.3698E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -703.81421668281484     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       31.336609880923334     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    59.67
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -703.814       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.03E-01  1.13E-02  3.82E-02  4.15E-01  6.30E+00  1.39E+00  1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.02E+01
 


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
+        7.46E+01
 
 TH 2
+       -1.57E+02  3.32E+02
 
 TH 3
+       -4.14E+03  8.73E+03  2.30E+05
 
 TH 4
+        5.07E+02 -1.07E+03 -2.81E+04  3.45E+03
 
 TH 5
+        3.09E+00 -6.52E+00 -1.71E+02  2.10E+01  1.28E-01
 
 TH 6
+       -4.20E+00  8.88E+00  2.33E+02 -2.86E+01 -1.74E-01  2.37E-01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.97E+00  6.27E+00  1.65E+02 -2.02E+01 -1.23E-01  1.68E-01  0.00E+00  0.00E+00  0.00E+00  0.00E+00  1.18E-01
 
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
+        2.04E+03
 
 TH 2
+       -7.35E+02  2.84E+03
 
 TH 3
+       -6.74E+03  1.40E+04  3.68E+05
 
 TH 4
+       -3.31E+02 -1.57E+03 -4.51E+04  6.29E+03
 
 TH 5
+        5.33E+00 -2.53E+01 -2.75E+02  3.45E+01  5.13E-01
 
 TH 6
+        7.69E+00  2.80E+01  3.71E+02 -7.42E+01  2.16E-01  8.37E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.34E+01  1.73E+01  2.65E+02 -2.34E+01 -2.25E-01  1.73E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  4.30E+00
 
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
+        2.07E+03
 
 TH 2
+       -6.14E+02  9.44E+02
 
 TH 3
+       -8.05E+03  1.88E+04  5.92E+05
 
 TH 4
+       -2.94E+02 -2.11E+03 -7.01E+04  9.20E+03
 
 TH 5
+        5.79E+00 -1.75E+01 -4.30E+02  5.19E+01  3.82E-01
 
 TH 6
+        2.23E+02 -5.18E+01 -1.05E+03 -2.91E+01  1.27E+00  1.15E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.19E+01  3.66E+01  4.31E+02 -2.92E+01 -4.82E-01 -5.47E-02  0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.32E+01
 
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
 #CPUT: Total CPU Time in Seconds,       65.724
Stop Time:
Wed Sep 29 20:04:59 CDT 2021
