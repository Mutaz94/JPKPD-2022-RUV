Wed Sep 29 21:00:37 CDT 2021
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
$DATA ../../../../data/spa1/B/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2042.22616920260        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7424E+02 -2.3231E+01 -4.7040E+01  3.6855E+01  6.7983E+01  2.3579E+01 -1.8184E+01  7.3160E+00 -1.4521E+01  1.4077E+01
            -9.2843E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2054.79992975510        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9072E-01  1.0533E+00  1.0816E+00  1.0256E+00  1.0151E+00  1.0705E+00  1.0789E+00  9.7363E-01  1.0784E+00  9.3045E-01
             1.1056E+00
 PARAMETER:  9.0673E-02  1.5193E-01  1.7847E-01  1.2526E-01  1.1503E-01  1.6815E-01  1.7598E-01  7.3278E-02  1.7547E-01  2.7915E-02
             2.0042E-01
 GRADIENT:   8.3523E+00  1.9032E+00 -1.7695E+01  2.5665E+01  2.0014E+01  1.8111E+00 -9.8033E+00  2.5031E+00 -2.0242E+00  3.4984E+00
            -6.0651E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2055.80508256854        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.8876E-01  9.9271E-01  1.1201E+00  1.0545E+00  9.9251E-01  1.0438E+00  1.3570E+00  9.3318E-01  1.0059E+00  8.9081E-01
             1.0966E+00
 PARAMETER:  8.8700E-02  9.2680E-02  2.1340E-01  1.5305E-01  9.2480E-02  1.4287E-01  4.0528E-01  3.0841E-02  1.0585E-01 -1.5622E-02
             1.9224E-01
 GRADIENT:   6.2959E+00  7.8400E+00 -7.0983E-01  5.9497E+00  9.8352E-01 -8.2916E+00  4.8674E+00  4.3968E-01 -9.6644E-01  3.6588E-01
            -1.3537E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2056.34440228671        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.8454E-01  8.6797E-01  1.2403E+00  1.1376E+00  9.8488E-01  1.0643E+00  1.3551E+00  9.8023E-01  9.9619E-01  9.0391E-01
             1.1154E+00
 PARAMETER:  8.4421E-02 -4.1594E-02  3.1534E-01  2.2888E-01  8.4763E-02  1.6234E-01  4.0387E-01  8.0034E-02  9.6187E-02 -1.0220E-03
             2.0919E-01
 GRADIENT:   8.4457E-02  6.4224E+00  2.4563E+00  6.1942E+00 -5.0323E+00  1.9720E-01  4.1591E-02 -2.4664E-01 -5.9517E-02 -9.5839E-03
            -2.5619E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2056.88229834296        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.8027E-01  5.0423E-01  1.5658E+00  1.3790E+00  9.7532E-01  1.0545E+00  1.5988E+00  1.1628E+00  9.1088E-01  9.4963E-01
             1.1120E+00
 PARAMETER:  8.0072E-02 -5.8472E-01  5.4841E-01  4.2137E-01  7.5010E-02  1.5303E-01  5.6928E-01  2.5083E-01  6.6604E-03  4.8317E-02
             2.0618E-01
 GRADIENT:   3.4322E-01  4.4540E+00  1.4318E+00  9.1710E+00 -1.6855E+00 -1.9605E+00 -1.0192E+00 -9.5489E-01 -1.7387E+00 -1.4976E-01
            -3.1955E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2056.99008439832        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  9.7864E-01  3.7263E-01  1.6951E+00  1.4731E+00  9.7046E-01  1.0527E+00  1.7856E+00  1.2725E+00  8.8119E-01  9.5868E-01
             1.1104E+00
 PARAMETER:  7.8408E-02 -8.8716E-01  6.2773E-01  4.8739E-01  7.0012E-02  1.5134E-01  6.7978E-01  3.4096E-01 -2.6483E-02  5.7803E-02
             2.0473E-01
 GRADIENT:   3.3726E-01  5.9253E+00  2.4836E+00  2.0318E+01 -6.4769E+00 -2.1950E+00 -7.4833E-01 -2.1983E-01 -1.3374E+00  8.2046E-01
            -4.0949E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2057.23407976595        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  9.7687E-01  2.5136E-01  1.7926E+00  1.5509E+00  9.6258E-01  1.0554E+00  2.1360E+00  1.3529E+00  8.5355E-01  9.5518E-01
             1.1133E+00
 PARAMETER:  7.6602E-02 -1.2809E+00  6.8364E-01  5.3886E-01  6.1860E-02  1.5390E-01  8.5891E-01  4.0223E-01 -5.8356E-02  5.4148E-02
             2.0734E-01
 GRADIENT:  -1.3639E-01  3.9088E+00  2.0482E+00  1.6840E+01 -6.4477E+00 -7.8141E-01 -4.4082E-01  2.3007E-01 -6.4869E-01  8.4894E-01
            -1.8011E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2057.43351991948        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1232
 NPARAMETR:  9.7535E-01  1.4199E-01  1.8821E+00  1.6210E+00  9.5690E-01  1.0571E+00  2.8598E+00  1.4095E+00  8.2807E-01  9.5021E-01
             1.1159E+00
 PARAMETER:  7.5038E-02 -1.8520E+00  7.3237E-01  5.8302E-01  5.5940E-02  1.5548E-01  1.1507E+00  4.4320E-01 -8.8653E-02  4.8931E-02
             2.0967E-01
 GRADIENT:  -3.8700E-01  2.1805E+00  1.9340E+00  1.2302E+01 -4.9168E+00  1.7953E-01 -1.8983E-01 -2.6184E-01 -4.6008E-01  2.5584E-01
            -1.6153E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2057.64915935338        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1413
 NPARAMETR:  9.7516E-01  9.3279E-02  1.9129E+00  1.6382E+00  9.5636E-01  1.0568E+00  3.5976E+00  1.4324E+00  8.2033E-01  9.5035E-01
             1.1161E+00
 PARAMETER:  7.4850E-02 -2.2722E+00  7.4861E-01  5.9361E-01  5.5380E-02  1.5529E-01  1.3803E+00  4.5937E-01 -9.8045E-02  4.9075E-02
             2.0981E-01
 GRADIENT:   7.3275E-01  4.0199E-01  1.1696E-01 -1.5928E+01  1.5187E+00  3.0758E-01 -2.7698E-02 -3.1085E-02  1.1567E+00 -2.2730E-02
             2.6083E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2057.67190220211        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1593
 NPARAMETR:  9.7460E-01  7.1332E-02  1.9122E+00  1.6542E+00  9.5097E-01  1.0558E+00  3.9217E+00  1.4316E+00  8.1448E-01  9.4854E-01
             1.1157E+00
 PARAMETER:  7.4270E-02 -2.5404E+00  7.4824E-01  6.0330E-01  4.9724E-02  1.5434E-01  1.4665E+00  4.5879E-01 -1.0520E-01  4.7173E-02
             2.0945E-01
 GRADIENT:   2.6031E-01  4.3572E-01 -5.7820E-01 -1.1894E+01  1.4225E+00  2.6598E-02 -8.9655E-02 -7.7198E-02  2.5804E-01  5.8443E-02
            -9.5023E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2057.67349688280        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1773
 NPARAMETR:  9.7431E-01  5.9052E-02  1.9155E+00  1.6636E+00  9.4875E-01  1.0556E+00  4.1288E+00  1.4336E+00  8.1256E-01  9.4729E-01
             1.1157E+00
 PARAMETER:  7.3979E-02 -2.7293E+00  7.5000E-01  6.0898E-01  4.7385E-02  1.5410E-01  1.5180E+00  4.6018E-01 -1.0757E-01  4.5855E-02
             2.0950E-01
 GRADIENT:   1.5826E-02  4.3483E-01 -8.3212E-01 -9.0030E+00  1.3107E+00 -3.1467E-02 -9.7736E-02 -1.3235E-01  3.5517E-01  2.5932E-02
            -1.1023E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2057.69579970087        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1938
 NPARAMETR:  9.7382E-01  4.8492E-02  1.9191E+00  1.6595E+00  9.4757E-01  1.0552E+00  5.4944E+00  1.4360E+00  8.1133E-01  9.4691E-01
             1.1158E+00
 PARAMETER:  7.3473E-02 -2.9263E+00  7.5184E-01  6.0652E-01  4.6147E-02  1.5372E-01  1.8037E+00  4.6185E-01 -1.0908E-01  4.5450E-02
             2.0957E-01
 GRADIENT:  -5.3625E-01  1.0812E-01 -4.7227E-01 -3.0067E+01  2.8085E+00 -1.0728E-01  1.6075E-01 -5.2309E-03  1.5383E+00  1.0947E-01
             2.0373E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2057.72426460091        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2118
 NPARAMETR:  9.7428E-01  3.6795E-02  1.9222E+00  1.6737E+00  9.4491E-01  1.0557E+00  5.3987E+00  1.4387E+00  8.0761E-01  9.4410E-01
             1.1158E+00
 PARAMETER:  7.3941E-02 -3.2024E+00  7.5349E-01  6.1502E-01  4.3336E-02  1.5420E-01  1.7862E+00  4.6374E-01 -1.1368E-01  4.2476E-02
             2.0954E-01
 GRADIENT:   6.7396E-01  1.9628E-01 -8.7158E-01 -1.7332E+01  1.8420E+00  1.0929E-01 -4.2547E-02 -1.4181E-01  4.1992E-01 -7.3922E-02
            -4.0233E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2057.74210020633        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     2281
 NPARAMETR:  9.7391E-01  2.7204E-02  1.9266E+00  1.6778E+00  9.4311E-01  1.0554E+00  7.1902E+00  1.4417E+00  8.0526E-01  9.4368E-01
             1.1158E+00
 PARAMETER:  7.3566E-02 -3.5044E+00  7.5574E-01  6.1747E-01  4.1425E-02  1.5396E-01  2.0727E+00  4.6582E-01 -1.1659E-01  4.2029E-02
             2.0955E-01
 GRADIENT:   2.2886E-01  1.8458E-01 -5.2524E-01 -2.1562E+01  1.3072E+00  5.1225E-02  6.3080E-02 -1.0254E-01  5.7504E-01  7.9892E-02
             2.6416E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2057.74723636813        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2464
 NPARAMETR:  9.7427E-01  2.2089E-02  1.9297E+00  1.6804E+00  9.4161E-01  1.0560E+00  7.0500E+00  1.4446E+00  8.0333E-01  9.4240E-01
             1.1158E+00
 PARAMETER:  7.3933E-02 -3.7127E+00  7.5735E-01  6.1905E-01  3.9833E-02  1.5445E-01  2.0530E+00  4.6783E-01 -1.1899E-01  4.0670E-02
             2.0953E-01
 GRADIENT:   1.1068E+00  1.1170E-01 -4.1039E-02 -2.2793E+01  2.2153E-01  2.7505E-01 -1.9548E-02 -1.1143E-01 -7.9651E-02  7.8247E-02
            -3.4300E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2057.75400356756        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     2629
 NPARAMETR:  9.7406E-01  1.7446E-02  1.9311E+00  1.6842E+00  9.4119E-01  1.0557E+00  8.8946E+00  1.4460E+00  8.0305E-01  9.4206E-01
             1.1158E+00
 PARAMETER:  7.3715E-02 -3.9486E+00  7.5808E-01  6.2130E-01  3.9385E-02  1.5416E-01  2.2854E+00  4.6882E-01 -1.1934E-01  4.0315E-02
             2.0956E-01
 GRADIENT:   7.9191E-01  1.3998E-01 -3.5661E-01 -2.1419E+01  6.1782E-01  1.7122E-01  3.7088E-02 -8.7480E-02  3.8515E-01  7.8330E-02
             1.3324E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2057.75632348963        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2810
 NPARAMETR:  9.7411E-01  1.2229E-02  1.9341E+00  1.6883E+00  9.4051E-01  1.0558E+00  8.9408E+00  1.4486E+00  8.0180E-01  9.4124E-01
             1.1158E+00
 PARAMETER:  7.3768E-02 -4.3040E+00  7.5965E-01  6.2372E-01  3.8662E-02  1.5426E-01  2.2906E+00  4.7061E-01 -1.2090E-01  3.9442E-02
             2.0956E-01
 GRADIENT:   1.0136E+00  8.2773E-02 -4.1743E-01 -2.0146E+01  4.8242E-01  2.1982E-01 -1.4622E-02 -1.0385E-01  9.8175E-02  1.3060E-02
            -1.9343E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2057.76141821716        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     2975
 NPARAMETR:  9.7398E-01  1.0009E-02  1.9352E+00  1.6894E+00  9.4037E-01  1.0556E+00  1.1883E+01  1.4495E+00  8.0161E-01  9.4121E-01
             1.1158E+00
 PARAMETER:  7.3635E-02 -4.5042E+00  7.6023E-01  6.2435E-01  3.8515E-02  1.5414E-01  2.5751E+00  4.7125E-01 -1.2113E-01  3.9412E-02
             2.0957E-01
 GRADIENT:   8.4595E-01  9.9540E-01 -4.7702E-01 -2.0991E+01  6.7004E-01  1.8478E-01  3.0563E-02 -8.1575E-02  3.8608E-01  3.6579E-02
             6.0173E-03

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2057.76259836048        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     3151
 NPARAMETR:  9.7423E-01  1.0000E-02  1.9387E+00  1.6885E+00  9.3957E-01  1.0560E+00  1.1553E+01  1.4530E+00  8.0095E-01  9.4051E-01
             1.1158E+00
 PARAMETER:  7.3899E-02 -4.5509E+00  7.6069E-01  6.2383E-01  3.8376E-02  1.5443E-01  2.5483E+00  4.7352E-01 -1.2197E-01  3.7664E-02
             2.0954E-01
 GRADIENT:   1.0491E-02  0.0000E+00 -3.1199E-01 -7.2014E-03  5.6800E-01 -5.4733E-03  2.6590E-04 -7.2830E-03 -4.8563E-03 -3.6713E-02
             8.2978E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3151
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2760E-04 -3.1189E-04 -2.9149E-02 -5.1445E-03 -3.8117E-02
 SE:             2.9827E-02  1.7270E-03  1.8163E-02  2.9419E-02  2.0173E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9391E-01  8.5668E-01  1.0852E-01  8.6118E-01  5.8818E-02

 ETASHRINKSD(%)  7.7442E-02  9.4214E+01  3.9152E+01  1.4430E+00  3.2419E+01
 ETASHRINKVR(%)  1.5482E-01  9.9665E+01  6.2975E+01  2.8651E+00  5.4329E+01
 EBVSHRINKSD(%)  3.9415E-01  9.4299E+01  4.2451E+01  1.7959E+00  2.9937E+01
 EBVSHRINKVR(%)  7.8675E-01  9.9675E+01  6.6881E+01  3.5595E+00  5.0912E+01
 RELATIVEINF(%)  9.5464E+01  8.0442E-03  6.9070E+00  2.7681E+00  8.2831E+00
 EPSSHRINKSD(%)  3.3344E+01
 EPSSHRINKVR(%)  5.5569E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2057.7625983604794     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1138.8240651558067     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    50.89
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:     7.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2057.763       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.74E-01  1.00E-02  1.94E+00  1.69E+00  9.40E-01  1.06E+00  1.16E+01  1.45E+00  8.01E-01  9.40E-01  1.12E+00
 


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
 
         3.09E-02  0.00E+00  3.65E-01  4.38E-02  9.21E-02  8.74E-02  1.67E+00  3.44E-01  5.31E-02  1.43E-01  5.51E-02
 


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
+        9.57E-04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.32E-04  0.00E+00  1.34E-01
 
 TH 4
+       -1.03E-04  0.00E+00  5.53E-03  1.92E-03
 
 TH 5
+       -2.74E-04  0.00E+00  3.00E-02  1.38E-03  8.49E-03
 
 TH 6
+       -3.79E-04  0.00E+00  2.62E-04  2.86E-04 -4.60E-04  7.64E-03
 
 TH 7
+       -1.08E-02  0.00E+00  4.50E-02 -7.60E-03  1.58E-02 -1.42E-02  2.78E+00
 
 TH 8
+       -2.25E-03  0.00E+00  8.38E-02  1.88E-03  1.69E-02  2.97E-03  6.42E-02  1.18E-01
 
 TH 9
+        2.79E-04  0.00E+00 -2.63E-03 -5.89E-04 -4.82E-04 -2.55E-04  1.65E-02 -2.97E-03  2.82E-03
 
 TH10
+       -1.59E-04  0.00E+00  1.66E-02  9.96E-04  6.29E-03  1.15E-03 -2.95E-02 -5.54E-03 -6.37E-04  2.03E-02
 
 TH11
+        2.82E-04  0.00E+00 -3.28E-04 -2.06E-04 -1.59E-04 -3.68E-04  1.36E-02 -1.44E-03  4.14E-04 -6.85E-04  3.03E-03
 
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
+        3.09E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -8.25E-02  0.00E+00  3.65E-01
 
 TH 4
+       -7.64E-02  0.00E+00  3.46E-01  4.38E-02
 
 TH 5
+       -9.60E-02  0.00E+00  8.90E-01  3.42E-01  9.21E-02
 
 TH 6
+       -1.40E-01  0.00E+00  8.19E-03  7.48E-02 -5.71E-02  8.74E-02
 
 TH 7
+       -2.10E-01  0.00E+00  7.39E-02 -1.04E-01  1.03E-01 -9.77E-02  1.67E+00
 
 TH 8
+       -2.12E-01  0.00E+00  6.67E-01  1.25E-01  5.32E-01  9.88E-02  1.12E-01  3.44E-01
 
 TH 9
+        1.70E-01  0.00E+00 -1.36E-01 -2.53E-01 -9.86E-02 -5.49E-02  1.87E-01 -1.62E-01  5.31E-02
 
 TH10
+       -3.60E-02  0.00E+00  3.18E-01  1.59E-01  4.78E-01  9.20E-02 -1.24E-01 -1.13E-01 -8.41E-02  1.43E-01
 
 TH11
+        1.65E-01  0.00E+00 -1.63E-02 -8.55E-02 -3.13E-02 -7.64E-02  1.48E-01 -7.61E-02  1.41E-01 -8.72E-02  5.51E-02
 
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
+        1.11E+03
 
 TH 2
+        5.48E-14  5.90E-29
 
 TH 3
+       -2.06E+01 -1.57E-14  3.52E+01
 
 TH 4
+       -2.05E+01  6.31E-15 -1.18E+01  6.23E+02
 
 TH 5
+        7.67E+01  6.93E-14 -1.58E+02 -6.72E+01  7.44E+02
 
 TH 6
+        5.67E+01  6.28E-14 -1.63E+01 -1.97E+01  6.25E+01  1.35E+02
 
 TH 7
+       -1.44E-02 -1.35E-17  3.26E-03  1.17E-02 -1.07E-02  6.58E-04  8.65E-06
 
 TH 8
+       -1.11E+00 -2.47E-15 -3.79E-01 -3.02E+00  4.81E-02  2.58E+00 -5.49E-04  4.93E-01
 
 TH 9
+       -9.97E+01 -8.79E-14  1.63E+01  1.01E+02 -4.95E+01  3.42E+00  5.80E-02 -3.96E+00  3.91E+02
 
 TH10
+       -7.44E+00 -9.04E-15  1.64E+01  1.61E+00 -8.00E+01  5.55E-01  3.80E-04  7.42E-01 -3.33E-01  9.83E+00
 
 TH11
+       -7.06E+01 -1.19E-13 -1.19E+01  3.56E+01 -1.91E+00  4.17E+00 -4.27E-03  1.19E+01 -3.09E+01  1.63E+01  3.50E+02
 
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
+        1.04E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        2.12E+01  0.00E+00  4.80E+01
 
 TH 4
+       -3.72E+01  0.00E+00 -8.71E+00  5.86E+02
 
 TH 5
+       -3.17E+01  0.00E+00 -1.52E+02 -4.52E+01  7.26E+02
 
 TH 6
+       -6.32E+01  0.00E+00  5.68E+00  2.83E+01 -8.77E+01  2.33E+02
 
 TH 7
+       -3.94E-03  0.00E+00  3.16E-04 -2.47E-02  1.21E-02 -3.59E-03  1.14E-05
 
 TH 8
+       -2.86E+01  0.00E+00 -1.27E+01 -2.07E+01 -3.35E+00  1.09E+01  2.39E-04  1.84E+01
 
 TH 9
+        8.32E+01  0.00E+00  4.86E+00 -8.39E+01  3.58E+01 -1.26E+01  4.09E-02 -6.33E+00  2.35E+02
 
 TH10
+       -1.47E+01  0.00E+00 -1.90E+00 -6.67E+00 -6.40E+01  2.85E+01 -3.48E-04  1.18E+01 -9.74E+00  4.99E+01
 
 TH11
+        6.42E+01  0.00E+00 -9.79E-01 -6.09E+01 -1.23E+01 -7.27E+00  1.15E-02  8.90E+00  4.15E+01  7.04E+00  3.16E+02
 
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
 #CPUT: Total CPU Time in Seconds,       58.362
Stop Time:
Wed Sep 29 21:01:37 CDT 2021
