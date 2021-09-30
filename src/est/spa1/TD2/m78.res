Thu Sep 30 02:21:34 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2126.20048483338        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5670E+02 -7.4853E+01 -5.2644E+01 -2.9427E+01  3.7064E+01  3.3861E+01 -1.6203E+01  1.8248E+01  1.1090E+01  2.7653E+00
             2.7502E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2144.23810953918        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.7678E-01  1.2118E+00  1.3462E+00  9.7176E-01  1.2352E+00  9.8538E-01  1.1593E+00  8.4589E-01  9.1767E-01  1.0155E+00
             9.3610E-01
 PARAMETER:  7.6511E-02  2.9213E-01  3.9728E-01  7.1350E-02  3.1121E-01  8.5272E-02  2.4785E-01 -6.7364E-02  1.4079E-02  1.1537E-01
             3.3963E-02
 GRADIENT:   6.2433E+00  9.0607E+00 -4.5239E+00  2.0818E+01  2.7804E+01 -3.4310E+00 -5.3729E+00  2.1056E+00  6.7867E-01 -4.0718E+01
            -3.8908E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2147.88117367709        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.7274E-01  1.2108E+00  1.9363E+00  9.6726E-01  1.4081E+00  1.0078E+00  1.2135E+00  6.5950E-01  8.3271E-01  1.5341E+00
             9.9752E-01
 PARAMETER:  7.2363E-02  2.9126E-01  7.6080E-01  6.6716E-02  4.4222E-01  1.0775E-01  2.9355E-01 -3.1627E-01 -8.3067E-02  5.2793E-01
             9.7519E-02
 GRADIENT:  -2.4249E+00 -1.3803E+00  6.9680E+00 -1.2452E+01 -1.6505E+01  6.0252E+00 -2.5426E+00 -5.5505E-01 -6.3379E+00  1.6081E+01
             1.3384E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2150.15165776153        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.7102E-01  9.7247E-01  1.6995E+00  1.1267E+00  1.2538E+00  9.8799E-01  1.3896E+00  5.1399E-01  8.4349E-01  1.2843E+00
             9.7545E-01
 PARAMETER:  7.0596E-02  7.2080E-02  6.3032E-01  2.1929E-01  3.2618E-01  8.7917E-02  4.2898E-01 -5.6556E-01 -7.0211E-02  3.5024E-01
             7.5141E-02
 GRADIENT:  -2.8559E+00  6.5324E+00  5.3179E-01  7.1972E+00 -4.7566E-01 -1.0796E+00  1.3348E+00 -3.2930E-01  1.0369E-02 -2.1114E+00
            -2.0825E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2150.36747888304        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.7335E-01  7.9824E-01  2.1995E+00  1.2560E+00  1.3099E+00  9.8884E-01  1.3776E+00  7.2976E-01  8.6663E-01  1.3683E+00
             9.8161E-01
 PARAMETER:  7.2987E-02 -1.2535E-01  8.8821E-01  3.2791E-01  3.6993E-01  8.8776E-02  4.2034E-01 -2.1504E-01 -4.3149E-02  4.1357E-01
             8.1437E-02
 GRADIENT:   4.9668E+00  8.1120E+00  1.2144E+00  1.3092E+01 -1.2468E+00 -4.9661E-01 -1.9809E-02 -1.5228E-01 -3.8138E-01  1.1801E+00
            -7.7234E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2150.89959526627        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.6890E-01  5.9661E-01  2.7520E+00  1.3913E+00  1.3427E+00  9.8793E-01  1.3097E+00  1.3693E+00  8.6844E-01  1.4003E+00
             9.7418E-01
 PARAMETER:  6.8411E-02 -4.1649E-01  1.1123E+00  4.3021E-01  3.9471E-01  8.7853E-02  3.6978E-01  4.1433E-01 -4.1056E-02  4.3666E-01
             7.3836E-02
 GRADIENT:  -2.9334E+00  2.6886E+00 -2.1681E+00  4.8184E+00 -3.4005E-01 -7.9877E-01  7.1243E-02  1.1520E+00  1.7542E+00  2.0540E+00
            -4.7498E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2150.95338354191        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  9.6812E-01  5.0550E-01  3.1181E+00  1.4565E+00  1.3667E+00  9.8328E-01  1.2613E+00  1.5112E+00  8.5061E-01  1.4245E+00
             9.7271E-01
 PARAMETER:  6.7602E-02 -5.8220E-01  1.2372E+00  4.7603E-01  4.1239E-01  8.3141E-02  3.3214E-01  5.1289E-01 -6.1800E-02  4.5379E-01
             7.2326E-02
 GRADIENT:  -4.3850E+00  3.2311E+00 -3.6720E-01  5.4815E+00 -9.3857E-01 -2.7159E+00  1.5652E-03  4.6012E-01 -2.6598E-01  2.1902E+00
            -2.4643E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2151.01102368788        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  9.7039E-01  4.6328E-01  3.2417E+00  1.4736E+00  1.3739E+00  9.9005E-01  1.2267E+00  1.5516E+00  8.4658E-01  1.4164E+00
             9.7556E-01
 PARAMETER:  6.9939E-02 -6.6942E-01  1.2761E+00  4.8773E-01  4.1765E-01  9.0003E-02  3.0431E-01  5.3927E-01 -6.6552E-02  4.4812E-01
             7.5255E-02
 GRADIENT:   1.3245E+00 -1.7661E+00  1.8506E-01 -1.7093E+01  2.3399E+00  6.4300E-02  2.3491E-01  1.8756E-01  3.8853E-01  2.0450E-01
            -1.7721E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2151.09121209193        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1420
 NPARAMETR:  9.6866E-01  4.2090E-01  3.2219E+00  1.5128E+00  1.3550E+00  9.8967E-01  7.2233E-01  1.5151E+00  8.6190E-01  1.4033E+00
             9.7747E-01
 PARAMETER:  6.8163E-02 -7.6537E-01  1.2700E+00  5.1399E-01  4.0380E-01  8.9616E-02 -2.2527E-01  5.1549E-01 -4.8615E-02  4.3883E-01
             7.7216E-02
 GRADIENT:  -2.0807E+00  3.6631E+00  7.4241E-02  1.0696E+01 -2.3381E-01 -2.0726E-02  1.2746E-01 -1.7807E-01  1.6715E+00 -3.2850E-01
             8.8392E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2151.10985105478        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1595
 NPARAMETR:  9.6814E-01  3.7804E-01  3.2028E+00  1.5428E+00  1.3396E+00  9.8925E-01  4.6315E-01  1.4883E+00  8.5951E-01  1.3960E+00
             9.7759E-01
 PARAMETER:  6.7622E-02 -8.7276E-01  1.2640E+00  5.3361E-01  3.9241E-01  8.9195E-02 -6.6970E-01  4.9764E-01 -5.1396E-02  4.3361E-01
             7.7330E-02
 GRADIENT:  -2.4331E+00  3.9260E+00  1.9894E-01  1.5650E+01 -1.3662E+00 -4.6071E-02  8.8930E-02 -3.0180E-01  3.8908E+00 -3.5962E-01
             1.0639E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2151.11675120206        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1770
 NPARAMETR:  9.6760E-01  3.2163E-01  3.1739E+00  1.5818E+00  1.3203E+00  9.8868E-01  2.8910E-01  1.4552E+00  8.4268E-01  1.3867E+00
             9.7733E-01
 PARAMETER:  6.7062E-02 -1.0344E+00  1.2550E+00  5.5859E-01  3.7787E-01  8.8618E-02 -1.1410E+00  4.7515E-01 -7.1169E-02  4.2696E-01
             7.7066E-02
 GRADIENT:  -2.4587E+00  3.9776E+00  1.1043E-01  1.9274E+01 -2.0822E+00 -6.3262E-02  3.2815E-02 -3.3662E-01  4.1674E+00 -3.9846E-01
             1.0128E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2151.12228331385        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1945
 NPARAMETR:  9.6704E-01  2.5909E-01  3.1473E+00  1.6241E+00  1.2997E+00  9.8803E-01  1.6383E-01  1.4206E+00  8.2020E-01  1.3774E+00
             9.7687E-01
 PARAMETER:  6.6482E-02 -1.2506E+00  1.2465E+00  5.8493E-01  3.6217E-01  8.7956E-02 -1.7089E+00  4.5106E-01 -9.8212E-02  4.2022E-01
             7.6601E-02
 GRADIENT:  -2.3466E+00  3.6517E+00  4.5187E-02  2.0710E+01 -2.7344E+00 -7.7078E-02  7.9884E-03 -3.6700E-01  3.3587E+00 -3.9283E-01
             7.8450E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2151.14217974056        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2120
 NPARAMETR:  9.6662E-01  1.9646E-01  3.1310E+00  1.6648E+00  1.2814E+00  9.8742E-01  7.5716E-02  1.3940E+00  7.9724E-01  1.3706E+00
             9.7663E-01
 PARAMETER:  6.6052E-02 -1.5273E+00  1.2413E+00  6.0973E-01  3.4792E-01  8.7340E-02 -2.4808E+00  4.3215E-01 -1.2660E-01  4.1524E-01
             7.6355E-02
 GRADIENT:  -1.8597E+00  2.8817E+00 -1.6366E-01  1.8745E+01 -2.9697E+00 -6.2844E-02  1.2191E-03 -3.0446E-01  2.1464E+00 -2.1132E-01
             7.8706E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2151.16946622925        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2295
 NPARAMETR:  9.6638E-01  1.4288E-01  3.1449E+00  1.6991E+00  1.2700E+00  9.8694E-01  2.8853E-02  1.3847E+00  7.7804E-01  1.3661E+00
             9.7599E-01
 PARAMETER:  6.5803E-02 -1.8457E+00  1.2458E+00  6.3008E-01  3.3904E-01  8.6859E-02 -3.4455E+00  4.2551E-01 -1.5098E-01  4.1198E-01
             7.5697E-02
 GRADIENT:  -1.2358E+00  2.0371E+00 -2.4659E-01  1.4431E+01 -2.7308E+00 -5.5630E-02  1.2884E-04 -2.8590E-01  1.0014E+00 -2.6062E-01
             3.4246E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2151.18664199004        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2470
 NPARAMETR:  9.6618E-01  9.9712E-02  3.1989E+00  1.7274E+00  1.2672E+00  9.8668E-01  1.0000E-02  1.4018E+00  7.6341E-01  1.3669E+00
             9.7569E-01
 PARAMETER:  6.5594E-02 -2.2055E+00  1.2628E+00  6.4661E-01  3.3684E-01  8.6593E-02 -4.6107E+00  4.3774E-01 -1.6995E-01  4.1257E-01
             7.5390E-02
 GRADIENT:  -9.4656E-01  1.3593E+00 -2.9578E-01  1.0334E+01 -2.3794E+00 -5.1680E-02  0.0000E+00 -2.2764E-01  2.9145E-01 -1.6488E-01
             2.1544E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2151.28482448720        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2651             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6691E-01  6.2905E-02  3.2659E+00  1.7329E+00  1.2724E+00  9.8681E-01  1.0000E-02  1.4368E+00  7.5567E-01  1.3706E+00
             9.7544E-01
 PARAMETER:  6.6346E-02 -2.6661E+00  1.2835E+00  6.4982E-01  3.4094E-01  8.6727E-02 -5.5478E+00  4.6242E-01 -1.8015E-01  4.1522E-01
             7.5138E-02
 GRADIENT:   4.2792E+02  4.1735E+00  5.8916E+00  1.3361E+03  1.7647E+01  3.6559E+01  0.0000E+00  4.9129E-01  2.4346E+01  5.4148E+00
             1.7848E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2151.31149398632        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2839             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6696E-01  6.4277E-02  3.2703E+00  1.7383E+00  1.2709E+00  9.8674E-01  1.0000E-02  1.4362E+00  7.5104E-01  1.3723E+00
             9.7478E-01
 PARAMETER:  6.6405E-02 -2.6445E+00  1.2849E+00  6.5293E-01  3.3973E-01  8.6651E-02 -5.5478E+00  4.6198E-01 -1.8630E-01  4.1652E-01
             7.4458E-02
 GRADIENT:   4.2846E+02  4.7523E+00  5.9792E+00  1.3564E+03  1.5706E+01  3.6554E+01  0.0000E+00  4.2462E-01  2.2923E+01  5.8271E+00
             9.6737E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2151.31174917522        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3021
 NPARAMETR:  9.6696E-01  6.5272E-02  3.2701E+00  1.7380E+00  1.2711E+00  9.8674E-01  1.0000E-02  1.4363E+00  7.5198E-01  1.3724E+00
             9.7476E-01
 PARAMETER:  6.6406E-02 -2.6292E+00  1.2848E+00  6.5276E-01  3.3984E-01  8.6649E-02 -5.5478E+00  4.6204E-01 -1.8504E-01  4.1657E-01
             7.4431E-02
 GRADIENT:   1.5271E+00  1.4934E-01 -1.6505E-02 -2.0193E+01 -4.1350E-01  7.3442E-02  0.0000E+00 -1.4417E-02 -8.3950E-02  7.6313E-02
            -5.8403E-02

0ITERATION NO.:   86    OBJECTIVE VALUE:  -2151.31174917522        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     3049
 NPARAMETR:  9.6700E-01  6.3578E-02  3.2684E+00  1.7364E+00  1.2719E+00  9.8678E-01  1.0000E-02  1.4375E+00  7.5253E-01  1.3721E+00
             9.7486E-01
 PARAMETER:  6.6406E-02 -2.6292E+00  1.2848E+00  6.5276E-01  3.3984E-01  8.6649E-02 -5.5478E+00  4.6204E-01 -1.8504E-01  4.1657E-01
             7.4431E-02
 GRADIENT:  -5.8401E-02  2.2874E-02  4.4608E-02  2.2275E+00 -3.3721E-01 -1.1697E-02  0.0000E+00 -1.6198E-02 -1.7248E-01  3.1070E-02
            -5.7357E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3049
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1362E-04 -2.4856E-05 -1.9873E-02 -5.8323E-03 -3.9692E-02
 SE:             2.9820E-02  1.0787E-05  1.1067E-02  2.9463E-02  2.3101E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8091E-01  2.1211E-02  7.2536E-02  8.4308E-01  8.5759E-02

 ETASHRINKSD(%)  9.8392E-02  9.9964E+01  6.2925E+01  1.2959E+00  2.2609E+01
 ETASHRINKVR(%)  1.9669E-01  1.0000E+02  8.6255E+01  2.5750E+00  4.0107E+01
 EBVSHRINKSD(%)  3.6448E-01  9.9966E+01  6.6124E+01  1.4366E+00  1.8169E+01
 EBVSHRINKVR(%)  7.2763E-01  1.0000E+02  8.8524E+01  2.8526E+00  3.3037E+01
 RELATIVEINF(%)  9.8238E+01  5.9018E-07  3.1803E+00  4.9620E+00  1.8029E+01
 EPSSHRINKSD(%)  3.2013E+01
 EPSSHRINKVR(%)  5.3778E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2151.3117491752237     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1232.3732159705510     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    48.92
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2151.312       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.67E-01  6.53E-02  3.27E+00  1.74E+00  1.27E+00  9.87E-01  1.00E-02  1.44E+00  7.52E-01  1.37E+00  9.75E-01
 


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
+        1.21E+03
 
 TH 2
+       -1.73E+01  3.25E+02
 
 TH 3
+       -3.20E+00  1.25E+00  6.50E+00
 
 TH 4
+       -5.24E+00  4.18E+02 -1.04E+01  6.32E+02
 
 TH 5
+        2.05E+00 -6.73E+01 -2.92E+01 -1.59E+01  2.50E+02
 
 TH 6
+        1.54E+00 -3.15E+00 -4.07E-01 -2.30E+00 -1.11E+00  2.01E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        5.27E-01 -2.56E+00 -3.75E+00 -8.72E-01 -1.84E+00 -1.71E-01  0.00E+00  7.47E+00
 
 TH 9
+        2.92E+00 -9.48E+01  2.28E+00 -1.42E+00 -1.71E+00  3.14E-02  0.00E+00 -7.24E-01  3.35E+02
 
 TH10
+        1.45E+00  3.57E+00 -4.60E-01  8.61E-01 -3.56E+01  5.73E-01  0.00E+00  2.59E+00 -6.46E-01  5.16E+01
 
 TH11
+       -1.07E+01 -1.50E+01 -4.79E+00 -1.19E+01  4.01E+00  1.78E+00  0.00E+00  6.21E+00  7.32E+00  1.46E+01  4.48E+02
 
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
 #CPUT: Total CPU Time in Seconds,       56.636
Stop Time:
Thu Sep 30 02:22:32 CDT 2021
