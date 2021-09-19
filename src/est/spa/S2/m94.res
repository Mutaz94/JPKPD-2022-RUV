Sat Sep 18 13:44:21 CDT 2021
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
$DATA ../../../../data/spa/S2/dat94.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1658.28461845980        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8601E+01 -6.0808E+01 -3.1256E+01 -7.6672E+01  4.6521E+01 -7.1652E+00 -7.6525E+00  8.5370E+00 -4.0932E+01  2.6010E+01
             1.1498E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1668.30643207743        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0017E+00  9.4642E-01  1.0111E+00  1.0817E+00  9.1361E-01  9.7757E-01  9.7313E-01  9.3200E-01  1.2029E+00  7.0034E-01
             1.0232E+00
 PARAMETER:  1.0169E-01  4.4927E-02  1.1108E-01  1.7850E-01  9.6456E-03  7.7316E-02  7.2767E-02  2.9581E-02  2.8476E-01 -2.5618E-01
             1.2295E-01
 GRADIENT:   5.1865E+01  7.9875E+00  1.5738E+01  8.2550E+00 -7.5800E+00 -1.7394E+01  1.2785E+00 -2.1362E+00  1.5030E+01 -4.1305E+00
             1.5467E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1669.99098692227        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0034E+00  9.6097E-01  7.5776E-01  1.0660E+00  7.9967E-01  1.0199E+00  1.1664E+00  7.6682E-01  1.1213E+00  5.6336E-01
             9.9027E-01
 PARAMETER:  1.0342E-01  6.0187E-02 -1.7739E-01  1.6392E-01 -1.2355E-01  1.1972E-01  2.5394E-01 -1.6550E-01  2.1449E-01 -4.7384E-01
             9.0225E-02
 GRADIENT:   5.3115E+01  4.6128E+00 -1.3335E+01  2.8350E+01  1.3428E+01  1.3688E-01  3.4856E+00  4.1320E+00  9.5851E+00 -5.1250E-01
             6.8369E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1670.39521977195        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.8570E-01  9.9448E-01  7.3680E-01  1.0275E+00  8.0531E-01  1.0149E+00  1.1275E+00  6.3731E-01  1.1088E+00  6.0355E-01
             9.7366E-01
 PARAMETER:  8.5595E-02  9.4467E-02 -2.0544E-01  1.2712E-01 -1.1653E-01  1.1478E-01  2.2002E-01 -3.5050E-01  2.0326E-01 -4.0492E-01
             7.3303E-02
 GRADIENT:   1.4898E+01 -1.4645E+00 -4.9287E+00  3.5513E+00  5.7247E+00 -2.0732E+00  3.6590E-01  1.7809E+00  2.6820E+00  3.8449E-01
             3.2187E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1670.39766205380        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.8400E-01  9.9066E-01  7.3006E-01  1.0283E+00  8.0005E-01  1.0161E+00  1.1358E+00  6.0858E-01  1.1018E+00  6.0295E-01
             9.7340E-01
 PARAMETER:  8.3870E-02  9.0612E-02 -2.1463E-01  1.2794E-01 -1.2308E-01  1.1600E-01  2.2732E-01 -3.9663E-01  1.9695E-01 -4.0592E-01
             7.3040E-02
 GRADIENT:   1.1159E+01 -1.2851E+00 -3.9916E+00  2.6080E+00  4.5773E+00 -1.6302E+00  3.5283E-01  1.4875E+00  2.0096E+00  4.2844E-01
             2.7313E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1670.39784346970        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.8365E-01  9.9028E-01  7.2715E-01  1.0281E+00  7.9831E-01  1.0164E+00  1.1377E+00  5.9888E-01  1.1002E+00  6.0218E-01
             9.7334E-01
 PARAMETER:  8.3518E-02  9.0231E-02 -2.1862E-01  1.2771E-01 -1.2526E-01  1.1631E-01  2.2902E-01 -4.1270E-01  1.9551E-01 -4.0720E-01
             7.2977E-02
 GRADIENT:   1.0370E+01 -1.2165E+00 -3.7411E+00  2.4172E+00  4.2825E+00 -1.5224E+00  3.3551E-01  1.3982E+00  1.8668E+00  4.1236E-01
             2.5665E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1670.39786764674        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  9.8353E-01  9.9023E-01  7.2591E-01  1.0279E+00  7.9761E-01  1.0166E+00  1.1384E+00  5.9492E-01  1.0997E+00  6.0182E-01
             9.7332E-01
 PARAMETER:  8.3396E-02  9.0181E-02 -2.2033E-01  1.2755E-01 -1.2613E-01  1.1642E-01  2.2963E-01 -4.1933E-01  1.9500E-01 -4.0779E-01
             7.2954E-02
 GRADIENT:   1.0092E+01 -1.1869E+00 -3.6453E+00  2.3516E+00  4.1716E+00 -1.4827E+00  3.2755E-01  1.3629E+00  1.8167E+00  4.0327E-01
             2.5017E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1670.98528506734        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      612
 NPARAMETR:  9.9966E-01  8.9817E-01  7.7176E-01  1.0952E+00  7.8508E-01  1.0356E+00  1.2274E+00  5.5849E-01  1.0554E+00  6.3141E-01
             9.7406E-01
 PARAMETER:  9.9664E-02 -7.3906E-03 -1.5909E-01  1.9095E-01 -1.4197E-01  1.3497E-01  3.0493E-01 -4.8252E-01  1.5394E-01 -3.5980E-01
             7.3719E-02
 GRADIENT:   2.3808E+00  3.0489E-01 -2.0803E+00  1.2738E+00  2.4020E+00  7.1140E-01  1.1225E-01  6.6125E-02  8.4212E-03  2.4992E-01
             1.4335E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1671.19139081092        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      787
 NPARAMETR:  9.9380E-01  7.1450E-01  9.0759E-01  1.2232E+00  7.8207E-01  1.0337E+00  1.3398E+00  7.3555E-01  1.0013E+00  6.5628E-01
             9.7264E-01
 PARAMETER:  9.3777E-02 -2.3617E-01  3.0425E-03  3.0150E-01 -1.4581E-01  1.3315E-01  3.9250E-01 -2.0713E-01  1.0134E-01 -3.2117E-01
             7.2256E-02
 GRADIENT:  -4.7215E+00  4.4794E+00  3.3540E+00  4.8175E+00 -5.2181E+00  8.6666E-01 -1.6349E-01 -3.2850E-01 -9.4919E-01 -3.3915E-01
            -8.0256E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1671.67097533082        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      965
 NPARAMETR:  9.9273E-01  4.3359E-01  1.0249E+00  1.3980E+00  7.5164E-01  1.0249E+00  1.6121E+00  8.5319E-01  9.2977E-01  6.8955E-01
             9.7550E-01
 PARAMETER:  9.2706E-02 -7.3566E-01  1.2461E-01  4.3502E-01 -1.8550E-01  1.2458E-01  5.7756E-01 -5.8768E-02  2.7177E-02 -2.7172E-01
             7.5198E-02
 GRADIENT:   2.5791E+00  3.1937E+00  4.0747E+00  3.0211E+00 -9.5381E+00 -1.0149E+00  9.5267E-01  9.8706E-01 -1.3512E-01  1.9099E+00
             1.2515E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1671.97112259378        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1140
 NPARAMETR:  9.8870E-01  2.8759E-01  1.0671E+00  1.4834E+00  7.3679E-01  1.0257E+00  1.6658E+00  8.9146E-01  8.9669E-01  6.8646E-01
             9.7357E-01
 PARAMETER:  8.8636E-02 -1.1462E+00  1.6496E-01  4.9431E-01 -2.0545E-01  1.2541E-01  6.1029E-01 -1.4893E-02 -9.0414E-03 -2.7621E-01
             7.3215E-02
 GRADIENT:   1.7804E-01  3.0714E-01  5.7807E-01  4.3974E-01  1.3876E+00  2.4934E-02  1.2093E-01 -6.0232E-01 -7.1203E-02 -5.8865E-01
            -2.3827E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1672.05721137664        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1316
 NPARAMETR:  9.8620E-01  1.8530E-01  1.0772E+00  1.5423E+00  7.1676E-01  1.0245E+00  1.5999E+00  9.2040E-01  8.7198E-01  6.8716E-01
             9.7332E-01
 PARAMETER:  8.6108E-02 -1.5858E+00  1.7435E-01  5.3325E-01 -2.3301E-01  1.2418E-01  5.6991E-01  1.7049E-02 -3.6992E-02 -2.7519E-01
             7.2955E-02
 GRADIENT:  -6.4432E-01  2.5000E-01  7.9993E-02  1.0356E+00 -1.0764E+00 -1.9185E-03 -4.6723E-02  2.5752E-01 -2.9460E-01  2.9472E-01
             5.1668E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1672.07397183400        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1491
 NPARAMETR:  9.8539E-01  1.2648E-01  1.0893E+00  1.5764E+00  7.0877E-01  1.0237E+00  1.4972E+00  9.4020E-01  8.5663E-01  6.8159E-01
             9.7337E-01
 PARAMETER:  8.5285E-02 -1.9677E+00  1.8551E-01  5.5512E-01 -2.4422E-01  1.2343E-01  5.0363E-01  3.8341E-02 -5.4747E-02 -2.8333E-01
             7.3011E-02
 GRADIENT:   4.2421E-01  6.9658E-02 -1.4375E-01  6.5232E-01  9.3765E-02 -8.9376E-03 -3.3507E-02  7.6632E-02 -2.4413E-01 -1.2879E-01
            -3.3651E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1672.07597771817        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1667
 NPARAMETR:  9.8447E-01  9.8711E-02  1.0906E+00  1.5919E+00  7.0321E-01  1.0234E+00  1.4195E+00  9.4094E-01  8.4950E-01  6.8261E-01
             9.7335E-01
 PARAMETER:  8.4353E-02 -2.2156E+00  1.8676E-01  5.6495E-01 -2.5209E-01  1.2312E-01  4.5034E-01  3.9124E-02 -6.3110E-02 -2.8183E-01
             7.2987E-02
 GRADIENT:  -1.8578E-01  4.1476E-02  2.2926E-01  4.1627E-01 -4.9133E-01  2.3959E-03 -1.8808E-02 -5.4285E-03 -2.0413E-02  7.4254E-02
             1.1149E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1672.07754565844        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1849
 NPARAMETR:  9.8446E-01  9.4509E-02  1.0914E+00  1.5942E+00  7.0250E-01  1.0233E+00  1.5400E+00  9.4259E-01  8.4809E-01  6.8159E-01
             9.7333E-01
 PARAMETER:  8.4334E-02 -2.2591E+00  1.8747E-01  5.6637E-01 -2.5311E-01  1.2306E-01  5.3177E-01  4.0874E-02 -6.4770E-02 -2.8332E-01
             7.2970E-02
 GRADIENT:  -1.0905E-02  2.3280E-02  4.2665E-01 -1.1352E-02 -6.3816E-01  1.9847E-03 -1.7627E-02 -7.9663E-03  3.3553E-02  6.2514E-02
             4.8377E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1672.08316786680        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2026
 NPARAMETR:  9.8443E-01  9.4265E-02  1.0904E+00  1.5944E+00  7.0197E-01  1.0233E+00  2.2813E+00  9.4281E-01  8.4579E-01  6.7836E-01
             9.7333E-01
 PARAMETER:  8.4311E-02 -2.2616E+00  1.8656E-01  5.6649E-01 -2.5387E-01  1.2304E-01  9.2473E-01  4.1105E-02 -6.7485E-02 -2.8807E-01
             7.2971E-02
 GRADIENT:  -1.7063E-02  1.0190E-03  3.2878E-01 -2.2256E-01 -5.4685E-02  2.3417E-03  1.7241E-04 -1.0734E-02  4.2391E-03 -3.7962E-03
            -4.5611E-06

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1672.08387234447        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2208             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8443E-01  9.3205E-02  1.0848E+00  1.5944E+00  6.9951E-01  1.0233E+00  2.3672E+00  9.3732E-01  8.4547E-01  6.7753E-01
             9.7327E-01
 PARAMETER:  8.4308E-02 -2.2730E+00  1.8144E-01  5.6648E-01 -2.5737E-01  1.2304E-01  9.6171E-01  3.5269E-02 -6.7863E-02 -2.8930E-01
             7.2911E-02
 GRADIENT:   4.2427E+01  1.0979E+00  7.1070E-01  1.0285E+02  2.5327E+00  5.5544E+00  1.3898E-01  1.5453E-02  1.4329E+00  1.1298E-01
             7.4072E-02

0ITERATION NO.:   83    OBJECTIVE VALUE:  -1672.08387369721        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     2303
 NPARAMETR:  9.8443E-01  9.3213E-02  1.0848E+00  1.5944E+00  6.9951E-01  1.0233E+00  2.3668E+00  9.3751E-01  8.4547E-01  6.7755E-01
             9.7328E-01
 PARAMETER:  8.4305E-02 -2.2729E+00  1.8143E-01  5.6650E-01 -2.5737E-01  1.2303E-01  9.6154E-01  3.5469E-02 -6.7868E-02 -2.8927E-01
             7.2915E-02
 GRADIENT:   1.1827E-02  2.1215E-04  7.6594E-03 -1.2201E-01 -1.6892E-02  1.8438E-03  1.8455E-04  3.7801E-03  1.0746E-03  3.3253E-03
             1.7862E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2303
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0119E-04 -2.7864E-03 -1.9482E-02 -3.2827E-03 -2.4806E-02
 SE:             2.9847E-02  3.5591E-03  1.8914E-02  2.9288E-02  1.9785E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9462E-01  4.3369E-01  3.0300E-01  9.1076E-01  2.0991E-01

 ETASHRINKSD(%)  8.8852E-03  8.8076E+01  3.6635E+01  1.8814E+00  3.3718E+01
 ETASHRINKVR(%)  1.7770E-02  9.8578E+01  5.9849E+01  3.7273E+00  5.6067E+01
 EBVSHRINKSD(%)  3.9410E-01  8.8319E+01  3.7976E+01  2.2535E+00  3.2747E+01
 EBVSHRINKVR(%)  7.8664E-01  9.8636E+01  6.1531E+01  4.4561E+00  5.4770E+01
 RELATIVEINF(%)  9.4285E+01  5.0488E-02  4.9121E+00  5.8928E+00  2.3736E+00
 EPSSHRINKSD(%)  4.4913E+01
 EPSSHRINKVR(%)  6.9655E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1672.0838736972098     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -936.93304713347163     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.07
 Elapsed covariance  time in seconds:     5.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1672.084       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  9.32E-02  1.08E+00  1.59E+00  7.00E-01  1.02E+00  2.37E+00  9.38E-01  8.45E-01  6.78E-01  9.73E-01
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.01E-02  3.03E-01  1.59E-01  1.70E-01  9.36E-02  7.27E-02  4.66E+00  1.45E-01  9.94E-02  1.05E-01  6.99E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.07E-04
 
 TH 2
+        9.14E-04  9.15E-02
 
 TH 3
+        6.49E-04 -1.48E-03  2.52E-02
 
 TH 4
+       -4.69E-04 -4.97E-02  4.09E-03  2.90E-02
 
 TH 5
+        4.01E-04  1.89E-02  1.01E-02 -8.84E-03  8.76E-03
 
 TH 6
+        2.30E-05  2.81E-03  5.43E-04 -1.04E-03  6.14E-04  5.29E-03
 
 TH 7
+       -1.08E-02 -9.38E-01  2.98E-01  5.32E-01 -7.37E-02 -3.43E-02  2.17E+01
 
 TH 8
+        2.53E-04 -8.87E-04  1.41E-02  2.29E-03  5.14E-03  1.19E-03  2.90E-01  2.12E-02
 
 TH 9
+        8.62E-05  2.22E-02 -2.62E-03 -1.20E-02  3.65E-03  8.40E-04 -2.93E-01 -2.46E-03  9.88E-03
 
 TH10
+        1.02E-04 -7.30E-04  5.68E-03  1.15E-03  3.26E-03 -5.19E-04  1.24E-02 -1.94E-03 -5.92E-04  1.11E-02
 
 TH11
+        6.32E-05 -3.83E-03  1.63E-03  1.94E-03 -1.53E-04  3.88E-04  1.13E-01  4.36E-04 -1.16E-03 -7.53E-04  4.89E-03
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.01E-02
 
 TH 2
+        1.00E-01  3.03E-01
 
 TH 3
+        1.36E-01 -3.08E-02  1.59E-01
 
 TH 4
+       -9.15E-02 -9.63E-01  1.51E-01  1.70E-01
 
 TH 5
+        1.42E-01  6.68E-01  6.82E-01 -5.54E-01  9.36E-02
 
 TH 6
+        1.05E-02  1.28E-01  4.70E-02 -8.41E-02  9.02E-02  7.27E-02
 
 TH 7
+       -7.71E-02 -6.66E-01  4.03E-01  6.70E-01 -1.69E-01 -1.01E-01  4.66E+00
 
 TH 8
+        5.78E-02 -2.02E-02  6.10E-01  9.25E-02  3.77E-01  1.12E-01  4.28E-01  1.45E-01
 
 TH 9
+        2.88E-02  7.37E-01 -1.66E-01 -7.08E-01  3.92E-01  1.16E-01 -6.33E-01 -1.70E-01  9.94E-02
 
 TH10
+        3.20E-02 -2.29E-02  3.39E-01  6.38E-02  3.30E-01 -6.77E-02  2.53E-02 -1.26E-01 -5.65E-02  1.05E-01
 
 TH11
+        3.00E-02 -1.81E-01  1.47E-01  1.63E-01 -2.34E-02  7.62E-02  3.48E-01  4.29E-02 -1.67E-01 -1.02E-01  6.99E-02
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.17E+03
 
 TH 2
+       -1.76E+01  4.00E+02
 
 TH 3
+       -1.06E+02  1.76E+02  4.69E+02
 
 TH 4
+        3.46E+01  3.93E+02 -5.01E+01  6.42E+02
 
 TH 5
+        1.42E+02 -6.39E+02 -9.50E+02 -1.27E+02  2.54E+03
 
 TH 6
+        1.24E+01 -4.61E+01 -2.27E+01 -5.40E+01  7.76E+01  2.08E+02
 
 TH 7
+        1.52E+00  2.04E+00  7.05E-02  6.31E-01 -2.35E+00  6.40E-01  1.55E-01
 
 TH 8
+       -4.12E-01 -1.47E+01 -6.14E+01  1.37E+00  2.14E+00 -2.08E+01 -1.46E+00  1.13E+02
 
 TH 9
+        3.07E+01 -6.89E+01  1.03E+01 -3.93E+01 -1.03E+01 -5.36E+00  1.03E+00  1.29E+01  2.48E+02
 
 TH10
+       -4.84E+00  7.46E+01  4.43E+01  2.03E+01 -2.86E+02 -6.12E+00  2.28E-01  5.38E+01  1.06E+01  1.66E+02
 
 TH11
+       -3.26E+01  3.22E+01 -1.19E+01  4.93E+01 -5.47E+01 -3.66E+01 -2.01E+00  4.51E+01 -6.64E+00  4.54E+01  2.64E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.04
 #CPUT: Total CPU Time in Seconds,       31.654
Stop Time:
Sat Sep 18 13:44:55 CDT 2021
