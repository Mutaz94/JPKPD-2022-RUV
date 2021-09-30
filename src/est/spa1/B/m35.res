Wed Sep 29 21:03:42 CDT 2021
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
$DATA ../../../../data/spa1/B/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2150.90833414291        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4702E+02  1.2778E+01 -4.0686E+01  8.5235E+01  6.5707E+01  6.4105E+01  1.3061E+00  1.3953E+01  7.3921E+00 -1.0584E+01
            -1.6353E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2159.39185697907        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0594E+00  9.9567E-01  1.2238E+00  1.0419E+00  1.0639E+00  9.3639E-01  1.0105E+00  7.8260E-01  1.0011E+00  1.1585E+00
             9.9527E-01
 PARAMETER:  1.5774E-01  9.5665E-02  3.0197E-01  1.4101E-01  1.6197E-01  3.4275E-02  1.1049E-01 -1.4513E-01  1.0115E-01  2.4715E-01
             9.5257E-02
 GRADIENT:   6.8073E+00  2.2834E+01  1.4009E+00  2.8692E+01  2.0521E+00 -1.4906E+01  1.1401E+00 -3.9681E-02  2.8601E-01 -5.9712E+00
            -8.6466E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2160.15624183323        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0590E+00  8.8010E-01  1.2799E+00  1.1032E+00  1.0306E+00  9.5372E-01  9.6613E-01  7.0682E-01  9.7422E-01  1.1859E+00
             1.0057E+00
 PARAMETER:  1.5729E-01 -2.7721E-02  3.4677E-01  1.9825E-01  1.3018E-01  5.2611E-02  6.5546E-02 -2.4697E-01  7.3882E-02  2.7054E-01
             1.0572E-01
 GRADIENT:   8.3334E+00  1.3824E+01  1.0014E+01  1.2653E+01 -1.1482E+01 -6.6714E+00 -7.8505E-01 -3.3351E+00 -3.7351E-01 -1.2794E+00
            -6.5182E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2160.86568481159        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.0511E+00  7.1578E-01  1.5011E+00  1.2153E+00  1.0509E+00  9.7087E-01  1.0530E+00  9.5724E-01  9.1041E-01  1.2325E+00
             1.0045E+00
 PARAMETER:  1.4983E-01 -2.3438E-01  5.0623E-01  2.9502E-01  1.4961E-01  7.0440E-02  1.5168E-01  5.6301E-02  6.1440E-03  3.0901E-01
             1.0452E-01
 GRADIENT:  -4.5129E+00  1.1122E+01  4.0331E+00  1.5929E+01 -8.7004E+00  1.3495E+00  3.8904E-01  1.0939E-01 -9.5615E-01  2.9669E-01
            -2.5686E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2161.30362302469        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  1.0526E+00  5.1656E-01  1.6962E+00  1.3430E+00  1.0522E+00  9.6718E-01  8.3844E-01  1.0832E+00  8.7351E-01  1.2590E+00
             1.0049E+00
 PARAMETER:  1.5126E-01 -5.6057E-01  6.2837E-01  3.9488E-01  1.5086E-01  6.6631E-02 -7.6215E-02  1.7989E-01 -3.5237E-02  3.3028E-01
             1.0489E-01
 GRADIENT:   4.9282E+00  6.3083E+00  1.6029E+00  1.6032E+01 -2.7285E+00  9.4938E-01 -1.5075E-01 -2.7734E-01  3.1966E-01 -2.9852E-01
            -5.2998E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2161.40391037371        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  1.0487E+00  3.9308E-01  1.8022E+00  1.4215E+00  1.0437E+00  9.6257E-01  6.6148E-01  1.1723E+00  8.3964E-01  1.2673E+00
             1.0055E+00
 PARAMETER:  1.4752E-01 -8.3374E-01  6.8899E-01  4.5174E-01  1.4279E-01  6.1855E-02 -3.1328E-01  2.5898E-01 -7.4778E-02  3.3688E-01
             1.0545E-01
 GRADIENT:  -3.1341E-01  4.3812E+00  8.1356E-01  1.2593E+01 -2.9852E+00 -3.4149E-01 -6.3616E-02  2.5980E-01  6.9408E-01  6.9658E-01
             9.9638E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2161.50251264227        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1084             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0479E+00  3.0697E-01  1.8691E+00  1.4702E+00  1.0388E+00  9.6172E-01  7.6845E-01  1.2261E+00  8.1296E-01  1.2650E+00
             1.0056E+00
 PARAMETER:  1.4677E-01 -1.0810E+00  7.2545E-01  4.8537E-01  1.3806E-01  6.0972E-02 -1.6338E-01  3.0385E-01 -1.0707E-01  3.3506E-01
             1.0557E-01
 GRADIENT:   6.9048E+02  4.8013E+01  8.7512E+00  8.8804E+02  9.2633E+00  4.9994E+01  4.2389E-01  9.3198E-01  1.3301E+01  4.5724E+00
             1.4691E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2161.50779947337        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:     1206
 NPARAMETR:  1.0464E+00  3.0606E-01  1.8689E+00  1.4682E+00  1.0387E+00  9.6114E-01  6.7357E-01  1.2260E+00  8.1276E-01  1.2648E+00
             1.0056E+00
 PARAMETER:  1.4539E-01 -1.0840E+00  7.2538E-01  4.8406E-01  1.3801E-01  6.0363E-02 -2.9517E-01  3.0375E-01 -1.0731E-01  3.3493E-01
             1.0557E-01
 GRADIENT:  -2.5309E+00  5.4946E-02  2.2695E-01 -1.1280E+01  7.5486E-01 -5.4209E-01  5.0190E-02  2.8549E-01  2.4464E+00  2.7567E-01
             4.2753E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2161.51725073934        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1382
 NPARAMETR:  1.0467E+00  2.9578E-01  1.8690E+00  1.4779E+00  1.0361E+00  9.6280E-01  7.4103E-01  1.2218E+00  8.0263E-01  1.2613E+00
             1.0053E+00
 PARAMETER:  1.4566E-01 -1.1181E+00  7.2540E-01  4.9065E-01  1.3547E-01  6.2090E-02 -1.9971E-01  3.0034E-01 -1.1986E-01  3.3215E-01
             1.0528E-01
 GRADIENT:  -1.5504E+00  1.1640E+00 -1.1636E-01 -4.7517E+00  1.1759E+00  1.8477E-01  1.3055E-02 -6.4248E-02 -1.4116E-02 -2.0987E-02
            -3.7873E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2161.53721070870        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1565
 NPARAMETR:  1.0481E+00  2.8783E-01  1.8684E+00  1.4799E+00  1.0327E+00  9.6242E-01  6.0982E-01  1.2228E+00  8.0258E-01  1.2604E+00
             1.0053E+00
 PARAMETER:  1.4701E-01 -1.1454E+00  7.2506E-01  4.9200E-01  1.3220E-01  6.1692E-02 -3.9459E-01  3.0111E-01 -1.1993E-01  3.3142E-01
             1.0530E-01
 GRADIENT:   2.0542E+00  4.6654E-01  2.3064E-01 -1.0585E+01  4.1649E-01  9.1481E-02  1.1442E-02 -2.0306E-02  1.5153E-01  2.3221E-01
            -9.5873E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2161.53763408814        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1745
 NPARAMETR:  1.0471E+00  2.7844E-01  1.8689E+00  1.4870E+00  1.0302E+00  9.6219E-01  5.2235E-01  1.2241E+00  8.0096E-01  1.2583E+00
             1.0053E+00
 PARAMETER:  1.4604E-01 -1.1786E+00  7.2536E-01  4.9679E-01  1.2977E-01  6.1460E-02 -5.4941E-01  3.0224E-01 -1.2195E-01  3.2973E-01
             1.0533E-01
 GRADIENT:   1.9524E-02  8.0824E-01 -8.2355E-02 -7.6985E+00  4.8008E-01  4.0580E-02  7.5273E-03 -5.2966E-02  3.3337E-01  1.8561E-01
            -3.7292E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2161.53770972192        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1934
 NPARAMETR:  1.0467E+00  2.7201E-01  1.8704E+00  1.4915E+00  1.0286E+00  9.6209E-01  4.6976E-01  1.2258E+00  7.9950E-01  1.2570E+00
             1.0054E+00
 PARAMETER:  1.4567E-01 -1.2019E+00  7.2616E-01  4.9980E-01  1.2824E-01  6.1350E-02 -6.5553E-01  3.0359E-01 -1.2377E-01  3.2876E-01
             1.0535E-01
 GRADIENT:  -6.6517E-01  9.2413E-01 -1.7367E-01 -6.7326E+00  4.1133E-01  2.4319E-02  6.6192E-03 -7.4223E-02  4.1416E-01  1.7169E-01
            -4.5906E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2161.53914044632        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2130
 NPARAMETR:  1.0499E+00  2.6411E-01  1.8721E+00  1.4907E+00  1.0276E+00  9.6233E-01  3.7509E-01  1.2281E+00  7.9756E-01  1.2556E+00
             1.0054E+00
 PARAMETER:  1.4874E-01 -1.2314E+00  7.2707E-01  4.9924E-01  1.2725E-01  6.1607E-02 -8.8059E-01  3.0548E-01 -1.2620E-01  3.2763E-01
             1.0541E-01
 GRADIENT:   7.1732E+00 -6.3296E-01 -1.4164E-01 -1.9787E+01  1.5355E+00  1.8371E-01  9.6068E-03 -2.9610E-02  3.7913E-01  5.3619E-02
             7.7047E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2161.55190717862        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2313
 NPARAMETR:  1.0485E+00  2.6135E-01  1.8741E+00  1.4948E+00  1.0258E+00  9.6221E-01  2.0376E-01  1.2303E+00  7.9678E-01  1.2548E+00
             1.0055E+00
 PARAMETER:  1.4740E-01 -1.2419E+00  7.2811E-01  5.0202E-01  1.2547E-01  6.1473E-02 -1.4908E+00  3.0722E-01 -1.2718E-01  3.2694E-01
             1.0544E-01
 GRADIENT:   3.9369E+00  1.8454E-01  2.0250E-01 -1.4109E+01  6.7579E-02  1.4710E-01  3.2258E-03 -4.0394E-02  7.7928E-04  1.4994E-01
             8.3301E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2161.55332077647        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2497
 NPARAMETR:  1.0483E+00  2.5912E-01  1.8739E+00  1.4968E+00  1.0254E+00  9.6214E-01  9.4910E-02  1.2314E+00  7.9656E-01  1.2538E+00
             1.0055E+00
 PARAMETER:  1.4720E-01 -1.2504E+00  7.2802E-01  5.0336E-01  1.2508E-01  6.1401E-02 -2.2548E+00  3.0814E-01 -1.2745E-01  3.2619E-01
             1.0545E-01
 GRADIENT:   3.5266E+00  3.2017E-01 -6.2129E-02 -1.2686E+01  3.4193E-01  1.2820E-01  1.3992E-03 -1.9647E-03  1.1204E-01  6.8880E-02
             9.4382E-04

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2161.55343767104        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2677
 NPARAMETR:  1.0478E+00  2.5535E-01  1.8757E+00  1.5002E+00  1.0249E+00  9.6202E-01  7.2181E-02  1.2326E+00  7.9565E-01  1.2535E+00
             1.0055E+00
 PARAMETER:  1.4671E-01 -1.2651E+00  7.2900E-01  5.0560E-01  1.2462E-01  6.1279E-02 -2.5286E+00  3.0911E-01 -1.2860E-01  3.2595E-01
             1.0546E-01
 GRADIENT:   2.4624E+00  5.0824E-01 -2.0121E-01 -1.0782E+01  4.3299E-01  9.2561E-02  9.9416E-04 -3.8600E-02  3.1712E-01  5.8545E-02
             1.5041E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2161.55348329105        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2867
 NPARAMETR:  1.0475E+00  2.5089E-01  1.8784E+00  1.5036E+00  1.0243E+00  9.6194E-01  6.0218E-02  1.2346E+00  7.9422E-01  1.2532E+00
             1.0055E+00
 PARAMETER:  1.4644E-01 -1.2827E+00  7.3041E-01  5.0786E-01  1.2405E-01  6.1197E-02 -2.7098E+00  3.1075E-01 -1.3040E-01  3.2568E-01
             1.0547E-01
 GRADIENT:   1.9542E+00  6.0870E-01 -2.6829E-01 -9.8255E+00  4.4382E-01  7.7175E-02  7.7419E-04 -6.1061E-02  4.0909E-01  5.5972E-02
            -9.0908E-04

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2161.55645935387        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     3053
 NPARAMETR:  1.0477E+00  2.4784E-01  1.8798E+00  1.5042E+00  1.0237E+00  9.6198E-01  2.5548E-02  1.2360E+00  7.9278E-01  1.2528E+00
             1.0055E+00
 PARAMETER:  1.4656E-01 -1.2950E+00  7.3117E-01  5.0824E-01  1.2343E-01  6.1242E-02 -3.5672E+00  3.1188E-01 -1.3221E-01  3.2537E-01
             1.0546E-01
 GRADIENT:   2.3739E+00  2.9957E-01 -1.2954E-01 -1.3020E+01  4.2038E-01  1.1254E-01  2.2915E-04 -5.0069E-02  2.4389E-01  6.6322E-02
             1.6378E-02

0ITERATION NO.:   88    OBJECTIVE VALUE:  -2161.55716479234        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     3154
 NPARAMETR:  1.0483E+00  2.4605E-01  1.8752E+00  1.5054E+00  1.0231E+00  9.6199E-01  1.0000E-02  1.2402E+00  7.9162E-01  1.2484E+00
             1.0056E+00
 PARAMETER:  1.4718E-01 -1.2959E+00  7.3140E-01  5.0820E-01  1.2297E-01  6.1300E-02 -5.9129E+00  3.1263E-01 -1.3279E-01  3.2512E-01
             1.0544E-01
 GRADIENT:   4.9987E-03  3.6668E-02  1.0840E-01 -3.6032E-01  1.7048E-02  2.2378E-03  0.0000E+00 -2.3296E-02  4.1200E-02  7.5588E-02
            -1.3824E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3154
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.8154E-04 -1.0289E-04 -3.2209E-02 -5.0150E-03 -3.7655E-02
 SE:             2.9895E-02  4.6696E-05  1.5072E-02  2.9448E-02  2.2796E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8715E-01  2.7568E-02  3.2595E-02  8.6477E-01  9.8566E-02

 ETASHRINKSD(%)  1.0000E-10  9.9844E+01  4.9508E+01  1.3452E+00  2.3631E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  7.4506E+01  2.6723E+00  4.1678E+01
 EBVSHRINKSD(%)  3.3925E-01  9.9859E+01  5.4625E+01  1.6501E+00  1.8893E+01
 EBVSHRINKVR(%)  6.7734E-01  1.0000E+02  7.9411E+01  3.2729E+00  3.4216E+01
 RELATIVEINF(%)  9.8128E+01  1.4852E-05  6.5131E+00  7.7886E+00  1.8572E+01
 EPSSHRINKSD(%)  3.3635E+01
 EPSSHRINKVR(%)  5.5957E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2161.5571647923407     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1242.6186315876680     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    49.57
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2161.557       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  2.48E-01  1.88E+00  1.50E+00  1.02E+00  9.62E-01  1.00E-02  1.24E+00  7.92E-01  1.25E+00  1.01E+00
 


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
+        1.09E+03
 
 TH 2
+       -2.45E+01  3.83E+02
 
 TH 3
+       -3.20E+00  1.81E+01  4.64E+01
 
 TH 4
+       -7.29E+00  4.78E+02 -2.47E+01  7.57E+02
 
 TH 5
+        3.10E+00 -1.38E+02 -1.04E+02 -2.73E+01  4.51E+02
 
 TH 6
+        9.02E-02 -3.08E+00 -3.15E-01 -1.41E+00  6.27E-01  2.13E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.48E-01 -3.87E-01 -1.69E+01 -3.16E+00  1.56E-01 -1.04E-01  0.00E+00  2.30E+01
 
 TH 9
+        1.60E+00 -1.07E+02  3.26E+00  9.35E-01 -2.10E-02 -1.05E+00  0.00E+00 -1.82E-01  2.99E+02
 
 TH10
+        2.22E+00  3.27E+00 -6.00E+00  1.23E+00 -5.39E+01  2.97E-01  0.00E+00  7.51E+00 -6.90E-02  5.99E+01
 
 TH11
+       -8.28E+00 -1.21E+01 -4.86E+00 -1.23E+01  2.78E+00  2.09E+00  0.00E+00  5.13E+00  7.88E+00  9.96E+00  4.07E+02
 
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
 #CPUT: Total CPU Time in Seconds,       56.732
Stop Time:
Wed Sep 29 21:04:41 CDT 2021
