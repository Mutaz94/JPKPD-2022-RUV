Thu Sep 30 09:05:51 CDT 2021
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
$DATA ../../../../data/spa2/D/dat45.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17928.8417903208        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.3934E+02  3.7631E+02 -5.8573E+01  3.5330E+02  1.6665E+02 -1.7201E+03 -9.2269E+02 -5.6384E+01 -1.4950E+03 -3.7528E+02
            -3.5619E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -587.113722322629        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0058E+00  1.2663E+00  1.0310E+00  1.4477E+00  1.1083E+00  2.4598E+00  2.1108E+00  9.7692E-01  2.0431E+00  1.0760E+00
             1.3198E+01
 PARAMETER:  1.0576E-01  3.3613E-01  1.3052E-01  4.7000E-01  2.0283E-01  1.0001E+00  8.4705E-01  7.6653E-02  8.1448E-01  1.7323E-01
             2.6801E+00
 GRADIENT:  -7.6601E+01 -2.9664E+01 -2.8983E+01  3.7338E+01  2.8701E+01  5.9703E+01 -2.7423E+01  3.0649E+00  1.3377E+01  1.3334E+01
             3.3632E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -636.081518054447        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0073E+00  1.7224E+00  4.5753E+00  1.4068E+00  2.6426E+00  2.6410E+00  4.4851E+00  8.2859E-01  2.9290E+00  1.3772E+00
             1.2158E+01
 PARAMETER:  1.0729E-01  6.4371E-01  1.6207E+00  4.4132E-01  1.0717E+00  1.0712E+00  1.6008E+00 -8.8024E-02  1.1747E+00  4.2006E-01
             2.5980E+00
 GRADIENT:  -5.4168E+01  1.4734E+01 -1.1389E+01  1.5762E+01  1.3743E+01  7.9514E+01  6.0208E+01  9.1475E-02  4.2986E+01  9.5012E+00
             3.2841E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -714.023661002665        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0657E+00  7.1647E-01  1.8202E+01  1.5149E+00  2.3462E+00  2.2166E+00  4.7096E+00  9.0787E+00  1.8044E+00  6.2659E-01
             8.6532E+00
 PARAMETER:  1.6368E-01 -2.3342E-01  3.0015E+00  5.1532E-01  9.5278E-01  8.9596E-01  1.6496E+00  2.3059E+00  6.9024E-01 -3.6746E-01
             2.2579E+00
 GRADIENT:   1.1139E+01 -7.9630E+00 -3.3407E+00  2.4031E+01  2.1460E+01  3.8633E+01  1.3392E+01  7.3099E+00  1.6352E+01  3.3525E+00
             2.4948E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -731.365254270474        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0495E+00  1.6984E+00  4.4419E+00  7.0236E-01  1.7638E+00  2.1583E+00  2.9969E+00  3.6782E-01  1.3370E+00  4.4203E-01
             7.8653E+00
 PARAMETER:  1.4827E-01  6.2971E-01  1.5911E+00 -2.5331E-01  6.6748E-01  8.6931E-01  1.1976E+00 -9.0015E-01  3.9040E-01 -7.1638E-01
             2.1625E+00
 GRADIENT:  -1.8447E+00  7.9492E+00  1.3160E+00  3.2819E+00 -5.9569E+00  1.0383E+01  4.9097E+00  6.9635E-03  5.1452E+00  2.8845E+00
            -3.9923E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -741.464827230789        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      460
 NPARAMETR:  1.0814E+00  1.3758E+00  5.0205E+00  9.3638E-01  1.8298E+00  2.3207E+00  4.2302E+00  9.6197E-02  7.7269E-01  3.0872E-01
             8.4151E+00
 PARAMETER:  1.7825E-01  4.1905E-01  1.7135E+00  3.4261E-02  7.0423E-01  9.4185E-01  1.5422E+00 -2.2414E+00 -1.5788E-01 -1.0753E+00
             2.2300E+00
 GRADIENT:   4.9201E-02  2.9192E+00 -1.0564E-01  5.0425E-01  2.3714E+00  1.5256E+00  1.9619E+00  1.2666E-03 -1.3529E+00  1.2126E+00
            -3.0338E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -741.995944516064        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      635
 NPARAMETR:  1.0803E+00  1.1289E+00  5.1442E+00  1.0508E+00  1.7648E+00  2.3007E+00  4.5762E+00  1.0930E-01  9.4380E-01  1.6170E-01
             8.4235E+00
 PARAMETER:  1.7728E-01  2.2120E-01  1.7379E+00  1.4951E-01  6.6801E-01  9.3321E-01  1.6209E+00 -2.1136E+00  4.2156E-02 -1.7220E+00
             2.2310E+00
 GRADIENT:   8.5660E-02 -1.0573E+00  1.1944E-02 -5.7919E-01  1.4897E-02 -8.1765E-01 -1.8051E-01  2.1275E-03  1.2367E-01  3.3949E-01
             7.1587E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -742.261927937340        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      817             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0773E+00  1.2129E+00  4.2690E+00  1.0081E+00  1.7049E+00  2.3571E+00  4.3864E+00  2.2866E-02  8.9920E-01  1.8724E-02
             8.4144E+00
 PARAMETER:  1.7442E-01  2.9298E-01  1.5514E+00  1.0805E-01  6.3350E-01  9.5743E-01  1.5785E+00 -3.6781E+00 -6.2515E-03 -3.8779E+00
             2.2299E+00
 GRADIENT:   1.2986E+01  4.7880E+00  6.0400E-01  4.7309E+00 -3.0395E+00  6.3899E+01  4.8581E+01  1.3752E-04 -4.3642E-01  5.1110E-03
             2.6547E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -742.296927113316        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      951
 NPARAMETR:  1.0782E+00  1.2158E+00  4.1896E+00  1.0039E+00  1.7143E+00  2.3205E+00  4.3778E+00  1.9309E-02  9.0488E-01  1.5986E-02
             8.4112E+00
 PARAMETER:  1.7528E-01  2.9542E-01  1.5326E+00  1.0394E-01  6.3899E-01  9.4180E-01  1.5766E+00 -3.8472E+00  4.6493E-05 -4.0360E+00
             2.2296E+00
 GRADIENT:  -4.6831E-01  7.4944E-02  1.6078E-01  1.1286E+00 -2.0923E+00  2.3108E+00 -2.0408E+00  9.5193E-05 -1.7020E-01  3.4565E-03
            -6.0792E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -742.315031750216        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1092
 NPARAMETR:  1.0794E+00  1.2193E+00  4.1985E+00  9.9658E-01  1.7365E+00  2.3274E+00  4.3368E+00  1.1841E-02  9.0973E-01  1.0000E-02
             8.4063E+00
 PARAMETER:  1.7641E-01  2.9829E-01  1.5347E+00  9.6579E-02  6.5189E-01  9.4475E-01  1.5671E+00 -4.3361E+00  5.3888E-03 -5.4969E+00
             2.2290E+00
 GRADIENT:  -2.3559E-02 -1.2558E+00 -2.5509E-01 -2.9939E-01  9.1231E-01  3.2474E+00 -2.7972E+00  3.4851E-05  2.0529E-01  0.0000E+00
            -2.5942E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -742.329163797707        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1277
 NPARAMETR:  1.0798E+00  1.2305E+00  4.2915E+00  9.9368E-01  1.7389E+00  2.3377E+00  4.3494E+00  1.0000E-02  9.0797E-01  1.0000E-02
             8.4104E+00
 PARAMETER:  1.7674E-01  3.0740E-01  1.5566E+00  9.3662E-02  6.5325E-01  9.4917E-01  1.5700E+00 -4.5771E+00  3.4550E-03 -5.4927E+00
             2.2295E+00
 GRADIENT:   5.4486E-02 -5.7817E-01 -9.0607E-02 -8.0763E-01  1.6009E-01  4.6899E+00 -1.5121E+00  0.0000E+00  2.9487E-01  0.0000E+00
             7.1761E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -742.332740867204        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1465
 NPARAMETR:  1.0799E+00  1.2348E+00  4.3266E+00  9.9232E-01  1.7410E+00  2.3398E+00  4.3426E+00  1.0000E-02  9.0361E-01  1.0000E-02
             8.4096E+00
 PARAMETER:  1.7691E-01  3.1088E-01  1.5648E+00  9.2294E-02  6.5446E-01  9.5009E-01  1.5685E+00 -4.5766E+00 -1.3573E-03 -5.4927E+00
             2.2294E+00
 GRADIENT:   1.1431E-01 -4.2982E-01 -4.3740E-02 -4.9195E-01 -2.4016E-02  5.0055E+00 -1.6822E+00  0.0000E+00  1.7436E-01  0.0000E+00
             4.4566E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -742.336057954775        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1663             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0799E+00  1.2432E+00  4.3646E+00  9.9076E-01  1.7443E+00  2.3400E+00  4.3290E+00  1.0000E-02  8.9719E-01  1.0000E-02
             8.4072E+00
 PARAMETER:  1.7688E-01  3.1768E-01  1.5735E+00  9.0717E-02  6.5637E-01  9.5016E-01  1.5653E+00 -4.5766E+00 -8.4911E-03 -5.4927E+00
             2.2291E+00
 GRADIENT:   1.4035E+01  4.6358E+00  7.4831E-02  2.3710E+00  1.1328E+00  6.0430E+01  4.8254E+01  0.0000E+00  1.0874E-01  0.0000E+00
             2.6352E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -742.336642094622        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1855
 NPARAMETR:  1.0799E+00  1.2453E+00  4.3736E+00  9.8971E-01  1.7456E+00  2.3401E+00  4.3238E+00  1.0000E-02  8.9573E-01  1.0000E-02
             8.4069E+00
 PARAMETER:  1.7689E-01  3.1939E-01  1.5756E+00  8.9654E-02  6.5708E-01  9.5020E-01  1.5641E+00 -4.5766E+00 -1.0112E-02 -5.4927E+00
             2.2290E+00
 GRADIENT:   9.1331E-02 -3.0802E-02  3.0259E-03  4.7774E-01 -1.9280E-01  5.0493E+00 -2.1782E+00  0.0000E+00 -7.4092E-02  0.0000E+00
            -4.5670E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -742.337305763337        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2052             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0799E+00  1.2474E+00  4.3768E+00  9.8736E-01  1.7471E+00  2.3401E+00  4.3205E+00  1.0000E-02  8.9759E-01  1.0000E-02
             8.4076E+00
 PARAMETER:  1.7690E-01  3.2107E-01  1.5763E+00  8.7277E-02  6.5797E-01  9.5019E-01  1.5634E+00 -4.5766E+00 -8.0386E-03 -5.4927E+00
             2.2291E+00
 GRADIENT:   1.4047E+01  4.5541E+00  6.3200E-02  1.7522E+00  1.3580E+00  6.0419E+01  4.8325E+01  0.0000E+00  2.3586E-01  0.0000E+00
             2.6776E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -742.337634068432        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2233             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0799E+00  1.2492E+00  4.3781E+00  9.8658E-01  1.7478E+00  2.3401E+00  4.3173E+00  1.0000E-02  8.9694E-01  1.0000E-02
             8.4076E+00
 PARAMETER:  1.7689E-01  3.2251E-01  1.5766E+00  8.6484E-02  6.5835E-01  9.5019E-01  1.5626E+00 -4.5761E+00 -8.7684E-03 -5.4927E+00
             2.2291E+00
 GRADIENT:   1.4039E+01  4.6025E+00  6.1801E-02  1.7625E+00  1.3876E+00  6.0418E+01  4.8239E+01  0.0000E+00  2.3254E-01  0.0000E+00
             2.6778E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -742.337914133339        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2420
 NPARAMETR:  1.0799E+00  1.2506E+00  4.3793E+00  9.8599E-01  1.7482E+00  2.3401E+00  4.3145E+00  1.0000E-02  8.9599E-01  1.0000E-02
             8.4074E+00
 PARAMETER:  1.7690E-01  3.2364E-01  1.5769E+00  8.5891E-02  6.5860E-01  9.5021E-01  1.5620E+00 -4.5530E+00 -9.8274E-03 -5.4927E+00
             2.2291E+00
 GRADIENT:   9.6754E-02 -1.3917E-01 -1.0362E-02 -1.0816E-01  3.7127E-02  5.0252E+00 -1.9536E+00  0.0000E+00  5.6835E-02  0.0000E+00
            -1.8837E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -742.338170336340        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2618             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0799E+00  1.2530E+00  4.3842E+00  9.8563E-01  1.7483E+00  2.3402E+00  4.3098E+00  1.0000E-02  8.9300E-01  1.0000E-02
             8.4065E+00
 PARAMETER:  1.7690E-01  3.2551E-01  1.5780E+00  8.5526E-02  6.5863E-01  9.5024E-01  1.5609E+00 -4.5530E+00 -1.3174E-02 -5.4927E+00
             2.2290E+00
 GRADIENT:   1.4042E+01  4.8131E+00  8.5136E-02  2.2165E+00  1.2455E+00  6.0459E+01  4.7812E+01  0.0000E+00  1.0550E-01  0.0000E+00
             2.6381E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -742.338299181790        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2809
 NPARAMETR:  1.0799E+00  1.2539E+00  4.3820E+00  9.8510E-01  1.7485E+00  2.3402E+00  4.3080E+00  1.0000E-02  8.9267E-01  1.0000E-02
             8.4064E+00
 PARAMETER:  1.7690E-01  3.2628E-01  1.5775E+00  8.4991E-02  6.5876E-01  9.5025E-01  1.5605E+00 -4.5530E+00 -1.3533E-02 -5.4927E+00
             2.2290E+00
 GRADIENT:   9.4876E-02 -4.0056E-03  9.2843E-03  2.7260E-01 -8.4592E-02  5.0467E+00 -2.1673E+00  0.0000E+00 -4.6335E-02  0.0000E+00
            -3.6913E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -742.338463037376        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3006             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0799E+00  1.2547E+00  4.3773E+00  9.8413E-01  1.7489E+00  2.3402E+00  4.3068E+00  1.0000E-02  8.9366E-01  1.0000E-02
             8.4068E+00
 PARAMETER:  1.7690E-01  3.2686E-01  1.5764E+00  8.4005E-02  6.5896E-01  9.5024E-01  1.5602E+00 -4.5530E+00 -1.2432E-02 -5.4927E+00
             2.2290E+00
 GRADIENT:   1.4042E+01  4.7682E+00  7.2252E-02  1.8977E+00  1.3672E+00  6.0445E+01  4.7902E+01  0.0000E+00  1.8029E-01  0.0000E+00
             2.6617E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -742.338531309668        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3194
 NPARAMETR:  1.0799E+00  1.2554E+00  4.3753E+00  9.8382E-01  1.7489E+00  2.3402E+00  4.3054E+00  1.0000E-02  8.9326E-01  1.0000E-02
             8.4068E+00
 PARAMETER:  1.7690E-01  3.2743E-01  1.5760E+00  8.3687E-02  6.5898E-01  9.5024E-01  1.5599E+00 -4.5530E+00 -1.2876E-02 -5.4927E+00
             2.2290E+00
 GRADIENT:   9.6422E-02 -6.6883E-02 -9.6296E-04  1.3563E-03  1.3314E-02  5.0348E+00 -2.0533E+00  0.0000E+00  1.8544E-02  0.0000E+00
            -1.4931E-01

0ITERATION NO.:  103    OBJECTIVE VALUE:  -742.338563351648        NO. OF FUNC. EVALS.: 108
 CUMULATIVE NO. OF FUNC. EVALS.:     3302
 NPARAMETR:  1.0799E+00  1.2560E+00  4.3710E+00  9.8332E-01  1.7490E+00  2.3402E+00  4.3043E+00  1.0000E-02  8.9326E-01  1.0000E-02
             8.4070E+00
 PARAMETER:  1.7690E-01  3.2793E-01  1.5761E+00  8.3673E-02  6.5890E-01  9.5025E-01  1.5596E+00 -4.5530E+00 -1.3876E-02 -5.4927E+00
             2.2290E+00
 GRADIENT:   4.0888E-04  1.2865E-03  7.0977E-03  1.5137E-01 -3.8349E-02  4.2061E-03 -9.6412E-03  0.0000E+00 -1.1940E-02  0.0000E+00
            -7.9438E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3302
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.4638E-03  1.8686E-02  1.3718E-06 -5.2718E-02 -3.5512E-06
 SE:             2.8819E-02  2.5417E-02  9.0674E-06  1.0602E-02  6.8906E-05
 N:                     100         100         100         100         100

 P VAL.:         7.4262E-01  4.6223E-01  8.7974E-01  6.6297E-07  9.5890E-01

 ETASHRINKSD(%)  3.4511E+00  1.4850E+01  9.9970E+01  6.4480E+01  9.9769E+01
 ETASHRINKVR(%)  6.7831E+00  2.7495E+01  1.0000E+02  8.7384E+01  9.9999E+01
 EBVSHRINKSD(%)  3.1429E+00  1.0495E+01  9.9954E+01  7.0135E+01  9.9660E+01
 EBVSHRINKVR(%)  6.1871E+00  1.9889E+01  1.0000E+02  9.1081E+01  9.9999E+01
 RELATIVEINF(%)  9.3355E+01  3.2942E+01  2.9222E-06  3.7322E+00  1.5609E-04
 EPSSHRINKSD(%)  1.0191E+01
 EPSSHRINKVR(%)  1.9343E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -742.33856335164774     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       360.38767649395936     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    78.87
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    10.83
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -742.339       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.08E+00  1.26E+00  4.38E+00  9.84E-01  1.75E+00  2.34E+00  4.30E+00  1.00E-02  8.92E-01  1.00E-02  8.41E+00
 


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
+        1.98E+02
 
 TH 2
+       -5.10E+00  1.01E+01
 
 TH 3
+        1.15E+00 -2.33E-02  9.82E-03
 
 TH 4
+       -6.70E+01  4.06E+01 -4.36E-01  1.76E+02
 
 TH 5
+       -1.58E+01 -2.71E+00 -1.41E-01 -5.74E+00  2.95E+00
 
 TH 6
+       -3.75E+01 -4.14E-01 -4.63E-01  1.28E+01  7.12E+00  2.63E+01
 
 TH 7
+        8.51E+00 -3.84E+00  5.29E-02 -1.71E+01  3.66E-01 -1.56E+00  1.69E+00
 
 TH 8
+       -7.35E-22 -2.19E-24 -7.68E-24  2.44E-22  1.17E-22  4.06E-22 -3.02E-23 -2.43E-43
 
 TH 9
+        1.79E+01 -7.36E+00  1.29E-01 -3.37E+01  2.84E-01 -4.76E+00  3.33E+00  3.16E-21  6.69E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.31E+00 -2.30E+00 -1.61E-02 -8.40E+00  9.35E-01  9.19E-01  7.53E-01 -4.93E-22  1.40E+00  0.00E+00  9.24E-01
 
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
+        1.60E+02
 
 TH 2
+       -2.16E+00  3.08E+01
 
 TH 3
+        1.27E-01  3.00E-01  1.83E-01
 
 TH 4
+       -6.03E+00  4.09E+01 -3.23E-01  1.69E+02
 
 TH 5
+       -2.36E+00 -8.62E+00 -2.32E+00 -8.26E+00  3.97E+01
 
 TH 6
+       -1.01E-01 -6.71E-01  4.44E-02  2.20E+00 -1.17E+00  3.04E+01
 
 TH 7
+        2.29E-01  3.43E+00 -8.85E-02 -1.62E+01  1.33E+00 -8.74E-01  6.28E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  6.55E-01
 
 TH 9
+        7.02E-01 -2.86E+00 -1.59E-01 -2.81E+01  2.73E+00 -4.34E-01  2.99E+00  0.00E+00  1.56E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.57E+00 -2.84E+00 -1.18E-03 -1.13E+01  7.38E-01  1.13E+00  9.06E-01  0.00E+00  1.92E+00  0.00E+00  9.36E+00
 
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
+        1.63E+02
 
 TH 2
+        5.35E+01  3.07E+01
 
 TH 3
+        6.09E-01  2.98E-01  6.23E-02
 
 TH 4
+        5.53E+01  4.31E+01  1.35E-01  1.81E+02
 
 TH 5
+       -1.60E+01 -6.82E+00 -9.60E-01 -1.63E+01  1.80E+01
 
 TH 6
+        4.51E+01  1.31E+01 -8.45E-02 -1.15E+00 -9.90E-01  4.39E+01
 
 TH 7
+        6.04E+00  3.03E+00 -6.95E-02 -1.76E+01  3.75E+00  6.86E+00  6.53E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.86E+00 -3.97E+00  1.20E-01 -4.07E+01  1.87E+00 -1.14E+00  3.96E+00  0.00E+00  1.68E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.38E+01 -2.07E+01 -5.30E-01 -5.26E+01  1.17E+01  7.69E+00  5.90E+00  0.00E+00 -2.09E+00  0.00E+00  2.16E+02
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       89.791
Stop Time:
Thu Sep 30 09:07:22 CDT 2021
