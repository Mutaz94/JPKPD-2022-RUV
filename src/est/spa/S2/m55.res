Sat Sep 18 13:31:02 CDT 2021
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
$DATA ../../../../data/spa/S2/dat55.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m55.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1666.50155608319        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2375E+01 -9.3050E+01 -1.5451E+01 -9.8484E+01  6.7667E+01  5.3871E+01 -9.5573E+00 -3.3650E+00 -7.9857E-01 -1.6857E+01
             3.0134E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1674.86459922779        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0066E+00  1.0599E+00  9.2326E-01  1.0251E+00  9.2372E-01  7.8287E-01  1.0537E+00  1.0420E+00  9.5082E-01  1.0790E+00
             9.8350E-01
 PARAMETER:  1.0656E-01  1.5816E-01  2.0160E-02  1.2482E-01  2.0653E-02 -1.4479E-01  1.5234E-01  1.4112E-01  4.9574E-02  1.7602E-01
             8.3358E-02
 GRADIENT:   6.6736E+01 -2.0445E+00  2.8163E+00 -1.2209E+01 -1.8346E+01 -3.6605E+01 -1.1170E+00  3.6359E+00 -2.3933E+00  1.1227E+01
            -2.3061E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1676.05081005267        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0013E+00  1.0065E+00  7.0818E-01  1.0507E+00  7.8799E-01  8.0414E-01  1.1988E+00  8.0914E-01  8.7410E-01  8.5439E-01
             9.8601E-01
 PARAMETER:  1.0132E-01  1.0651E-01 -2.4506E-01  1.4948E-01 -1.3827E-01 -1.1798E-01  2.8128E-01 -1.1178E-01 -3.4558E-02 -5.7363E-02
             8.5907E-02
 GRADIENT:   3.7947E+01  3.9348E+00 -1.2486E+01  1.7159E+01  5.7676E+00 -2.5889E+01  2.4222E+00  5.2597E+00 -5.7661E+00  4.6446E+00
            -9.5013E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1677.61339773027        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.8960E-01  1.0356E+00  6.9050E-01  1.0239E+00  7.9203E-01  8.5676E-01  1.1449E+00  6.3930E-01  9.2591E-01  8.5830E-01
             9.8632E-01
 PARAMETER:  8.9547E-02  1.3494E-01 -2.7034E-01  1.2360E-01 -1.3316E-01 -5.4597E-02  2.3530E-01 -3.4738E-01  2.3026E-02 -5.2806E-02
             8.6229E-02
 GRADIENT:   2.1139E+00 -4.3597E-01 -2.9940E+00  3.0954E+00  1.8018E+00 -2.5516E-01 -3.2034E-02  1.6803E+00  1.4637E-01  1.2452E+00
            -4.3265E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1677.61997409172        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.8954E-01  1.0315E+00  6.4660E-01  1.0197E+00  7.6318E-01  8.5786E-01  1.1515E+00  5.4015E-01  9.2233E-01  8.2686E-01
             9.8660E-01
 PARAMETER:  8.9487E-02  1.3100E-01 -3.3602E-01  1.1953E-01 -1.7026E-01 -5.3315E-02  2.4105E-01 -5.1591E-01  1.9144E-02 -9.0123E-02
             8.6513E-02
 GRADIENT:   6.9314E-01 -4.3255E-01 -2.4927E+00  2.3127E+00  1.4499E+00  2.5044E-02  1.2412E-01  1.3212E+00  1.7807E-01  1.0414E+00
            -2.3593E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1677.62058343282        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.8960E-01  1.0306E+00  6.3666E-01  1.0186E+00  7.5658E-01  8.5798E-01  1.1524E+00  5.1499E-01  9.2159E-01  8.1964E-01
             9.8664E-01
 PARAMETER:  8.9549E-02  1.3017E-01 -3.5152E-01  1.1845E-01 -1.7894E-01 -5.3179E-02  2.4186E-01 -5.6361E-01  1.8341E-02 -9.8886E-02
             8.6550E-02
 GRADIENT:   6.0451E-01 -4.1367E-01 -2.3038E+00  2.1097E+00  1.3280E+00  3.1857E-02  1.2168E-01  1.2214E+00  1.6939E-01  9.6497E-01
            -2.1215E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1678.00089698176        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  1.0042E+00  1.0254E+00  6.6728E-01  1.0314E+00  7.7474E-01  8.6913E-01  1.1786E+00  4.5086E-01  9.2675E-01  8.5700E-01
             9.9268E-01
 PARAMETER:  1.0416E-01  1.2504E-01 -3.0454E-01  1.3090E-01 -1.5522E-01 -4.0257E-02  2.6432E-01 -6.9660E-01  2.3929E-02 -5.4315E-02
             9.2655E-02
 GRADIENT:   5.3189E+00  1.1051E+00 -3.1039E-01  2.2057E+00  6.2728E-01  3.0377E+00  1.1873E+00  4.3850E-03  7.9083E-01  4.6058E-01
             1.6677E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1678.11235839862        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      691
 NPARAMETR:  1.0047E+00  1.1486E+00  5.9267E-01  9.5222E-01  7.8423E-01  8.6406E-01  1.0655E+00  2.6665E-01  9.7694E-01  8.4705E-01
             9.9016E-01
 PARAMETER:  1.0470E-01  2.3856E-01 -4.2312E-01  5.1043E-02 -1.4305E-01 -4.6116E-02  1.6341E-01 -1.2218E+00  7.6669E-02 -6.5993E-02
             9.0113E-02
 GRADIENT:   5.3171E+00  5.0813E+00  3.5070E-01  5.9361E+00 -1.9346E+00  5.1995E-01 -1.3025E-02  8.1163E-02 -3.1352E-01 -7.5283E-02
             6.3282E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1678.24938994634        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  1.0020E+00  1.3119E+00  5.4311E-01  8.4709E-01  8.3834E-01  8.6277E-01  9.5503E-01  1.1580E-01  1.0679E+00  8.8005E-01
             9.8752E-01
 PARAMETER:  1.0199E-01  3.7147E-01 -5.1045E-01 -6.5954E-02 -7.6330E-02 -4.7604E-02  5.3985E-02 -2.0559E+00  1.6569E-01 -2.7777E-02
             8.7443E-02
 GRADIENT:  -2.5682E+00  1.9905E+00 -4.8500E-01  1.9889E+00 -1.5516E+00 -2.9176E-02  5.0722E-02  3.6665E-02  7.5729E-02  3.9139E-01
            -1.7325E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1678.32664222131        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1043
 NPARAMETR:  1.0023E+00  1.5009E+00  5.0803E-01  7.2939E-01  9.2865E-01  8.6257E-01  8.6112E-01  3.2576E-02  1.2004E+00  9.4848E-01
             9.8889E-01
 PARAMETER:  1.0228E-01  5.0608E-01 -5.7722E-01 -2.1555E-01  2.5979E-02 -4.7844E-02 -4.9517E-02 -3.3242E+00  2.8268E-01  4.7105E-02
             8.8833E-02
 GRADIENT:  -1.1781E+00  6.1292E-01  1.4626E-01  3.6673E-01 -8.0830E-02  1.7360E-02  9.9161E-02  2.3087E-03 -1.5980E-01 -1.0499E-02
            -9.1994E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1678.33474199264        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1218
 NPARAMETR:  1.0026E+00  1.5804E+00  4.8601E-01  6.7879E-01  9.6524E-01  8.6235E-01  8.2777E-01  1.5485E-02  1.2670E+00  9.7536E-01
             9.8965E-01
 PARAMETER:  1.0260E-01  5.5767E-01 -6.2153E-01 -2.8744E-01  6.4621E-02 -4.8091E-02 -8.9022E-02 -4.0679E+00  3.3662E-01  7.5048E-02
             8.9600E-02
 GRADIENT:  -1.1641E-02  1.3691E-01 -3.3724E-02  1.4127E-01  7.4401E-02 -2.6461E-02 -7.1675E-03  5.3813E-04 -1.5143E-02 -1.7684E-02
            -1.5838E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1678.33490510629        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1398
 NPARAMETR:  1.0026E+00  1.5809E+00  4.8588E-01  6.7833E-01  9.6551E-01  8.6240E-01  8.2756E-01  1.0000E-02  1.2677E+00  9.7553E-01
             9.8969E-01
 PARAMETER:  1.0260E-01  5.5797E-01 -6.2180E-01 -2.8811E-01  6.4904E-02 -4.8033E-02 -8.9272E-02 -4.5432E+00  3.3719E-01  7.5228E-02
             8.9640E-02
 GRADIENT:   2.0395E-02 -7.6521E-02 -7.7317E-03 -2.2738E-02  4.1548E-02 -2.3074E-03 -2.7364E-03  0.0000E+00  8.2280E-03 -1.7538E-02
             1.4337E-03

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1678.33490542824        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1455
 NPARAMETR:  1.0026E+00  1.5809E+00  4.8588E-01  6.7836E-01  9.6549E-01  8.6241E-01  8.2757E-01  1.0000E-02  1.2676E+00  9.7573E-01
             9.8969E-01
 PARAMETER:  1.0260E-01  5.5801E-01 -6.2179E-01 -2.8807E-01  6.4883E-02 -4.8028E-02 -8.9260E-02 -4.5476E+00  3.3715E-01  7.5431E-02
             8.9634E-02
 GRADIENT:  -3.4812E-03  4.8438E-02  2.6746E-03  2.0170E-02 -3.3605E-02 -6.3083E-04  5.5227E-03  0.0000E+00  3.7375E-03  1.2177E-02
             2.1117E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1455
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4442E-04 -2.8109E-02 -3.2966E-04  2.1504E-02 -3.1956E-02
 SE:             2.9800E-02  2.3224E-02  1.2305E-04  2.4031E-02  2.2551E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9613E-01  2.2614E-01  7.3811E-03  3.7086E-01  1.5646E-01

 ETASHRINKSD(%)  1.6689E-01  2.2198E+01  9.9588E+01  1.9494E+01  2.4452E+01
 ETASHRINKVR(%)  3.3350E-01  3.9468E+01  9.9998E+01  3.5188E+01  4.2925E+01
 EBVSHRINKSD(%)  5.5839E-01  2.2014E+01  9.9645E+01  2.0251E+01  2.3269E+01
 EBVSHRINKVR(%)  1.1137E+00  3.9181E+01  9.9999E+01  3.6400E+01  4.1124E+01
 RELATIVEINF(%)  9.8860E+01  3.6896E+00  1.2842E-04  4.0730E+00  8.6457E+00
 EPSSHRINKSD(%)  4.5265E+01
 EPSSHRINKVR(%)  7.0041E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1678.3349054282407     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -943.18407886450257     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1678.335       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.58E+00  4.86E-01  6.78E-01  9.65E-01  8.62E-01  8.28E-01  1.00E-02  1.27E+00  9.76E-01  9.90E-01
 


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
+        1.47E+03
 
 TH 2
+       -7.65E+00  4.14E+02
 
 TH 3
+        4.88E+00  2.27E+02  6.55E+02
 
 TH 4
+       -1.56E+01  3.04E+02 -3.92E+02  9.60E+02
 
 TH 5
+       -4.84E+00 -2.71E+02 -5.44E+02  3.24E+02  7.37E+02
 
 TH 6
+        1.62E+00 -1.34E+00  1.69E+00 -4.76E+00 -2.93E+00  2.64E+02
 
 TH 7
+        1.12E+00  8.54E+00 -3.86E+01 -2.74E+00 -3.66E+00  2.84E-03  1.17E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.06E+00 -2.21E+01 -3.76E+01  5.25E+01 -6.46E+00  2.97E-01  1.83E+01  0.00E+00  6.00E+01
 
 TH10
+       -4.94E+00 -1.63E+01 -5.13E+01 -1.04E+01 -6.41E+01 -1.03E+00  2.11E+01  0.00E+00  6.33E+00  7.97E+01
 
 TH11
+       -1.20E+01 -1.40E+01 -2.43E+01 -2.06E-01 -9.19E-01  2.38E+00  1.09E+01  0.00E+00  4.55E+00  1.60E+01  2.15E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.492
Stop Time:
Sat Sep 18 13:31:25 CDT 2021
