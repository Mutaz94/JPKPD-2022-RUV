Sat Sep 25 09:15:39 CDT 2021
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
$DATA ../../../../data/spa/A3/dat40.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   403.442528492339        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.4455E+01  3.4001E+01  1.1422E+02 -1.0268E+02  1.6023E+01  5.9884E+01 -5.7886E+01 -2.7928E+01 -9.9416E+01 -7.8926E+01
            -3.8800E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1228.85542335662        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0590E+00  9.9372E-01  9.5838E-01  1.1313E+00  1.1181E+00  7.1297E-01  9.9969E-01  9.6785E-01  1.0041E+00  9.4478E-01
             5.3146E+00
 PARAMETER:  1.5733E-01  9.3698E-02  5.7488E-02  2.2341E-01  2.1165E-01 -2.3831E-01  9.9695E-02  6.7319E-02  1.0405E-01  4.3197E-02
             1.7705E+00
 GRADIENT:   8.1495E+01 -3.7045E+01 -2.5778E+01 -3.4997E+01  1.0764E+01 -1.6030E+01  9.2661E+00  6.3766E+00  1.7625E+01  1.4396E+01
             1.8983E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1260.53045948482        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0129E+00  8.6975E-01  8.1230E-01  1.2245E+00  8.8348E-01  7.9316E-01  4.4577E-01  1.3430E-01  8.8725E-01  6.5690E-01
             4.6212E+00
 PARAMETER:  1.1279E-01 -3.9554E-02 -1.0788E-01  3.0257E-01 -2.3887E-02 -1.3173E-01 -7.0796E-01 -1.9077E+00 -1.9628E-02 -3.2022E-01
             1.6306E+00
 GRADIENT:  -4.6834E+01  1.3655E+01 -1.4071E+01  3.5759E+01 -5.0908E+00  1.3534E+00  1.1067E-01  2.5566E-01 -6.9050E-01  1.2265E+01
             1.0269E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1274.98404443393        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0095E+00  7.3189E-01  9.2978E-01  1.2601E+00  8.8936E-01  7.7467E-01  7.4047E-01  1.5595E-01  9.2074E-01  2.5659E-01
             3.9897E+00
 PARAMETER:  1.0950E-01 -2.1213E-01  2.7194E-02  3.3122E-01 -1.7252E-02 -1.5531E-01 -2.0047E-01 -1.7582E+00  1.7423E-02 -1.2603E+00
             1.4837E+00
 GRADIENT:   3.4636E+00  4.7570E-01 -2.8781E+00 -1.6790E+00  7.2896E-01 -2.5489E+00  5.5048E-01  3.1819E-01  3.7177E-01  1.6435E+00
            -8.2955E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1276.42886106730        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0042E+00  4.6893E-01  9.0837E-01  1.4086E+00  7.8090E-01  7.7837E-01  4.5225E-01  7.0724E-02  8.8575E-01  4.4439E-02
             3.9737E+00
 PARAMETER:  1.0420E-01 -6.5729E-01  3.8993E-03  4.4257E-01 -1.4731E-01 -1.5056E-01 -6.9353E-01 -2.5490E+00 -2.1321E-02 -3.0136E+00
             1.4797E+00
 GRADIENT:  -8.3289E-02  1.2660E+00 -9.5479E-03  4.9990E+00 -6.9945E-01 -6.9345E-01  7.4921E-03  8.4105E-02  1.9169E+00  5.1830E-02
            -4.7253E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1276.65714030140        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.9962E-01  3.0223E-01  8.5275E-01  1.4887E+00  6.9872E-01  7.7813E-01  2.0871E-01  1.0486E-02  8.3411E-01  1.0000E-02
             3.9542E+00
 PARAMETER:  9.9617E-02 -1.0966E+00 -5.9290E-02  4.9792E-01 -2.5850E-01 -1.5086E-01 -1.4668E+00 -4.4577E+00 -8.1390E-02 -5.2858E+00
             1.4748E+00
 GRADIENT:  -2.7397E+00  8.6406E-01  5.4510E-01  3.4890E+00 -1.5103E+00 -5.1096E-01 -4.8723E-03  2.0984E-03 -1.0808E+00  0.0000E+00
            -1.5080E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1276.75904860991        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  9.9517E-01  1.3314E-01  7.1210E-01  1.5473E+00  5.7276E-01  7.7791E-01  3.0427E-02  1.0000E-02  8.2222E-01  1.0000E-02
             3.9218E+00
 PARAMETER:  9.5157E-02 -1.9163E+00 -2.3954E-01  5.3653E-01 -4.5729E-01 -1.5115E-01 -3.3924E+00 -9.3186E+00 -9.5749E-02 -1.0306E+01
             1.4666E+00
 GRADIENT:  -1.8099E+00  6.1493E-01  9.7207E-01  6.7115E+00 -2.0938E+00 -1.4329E-01 -2.8317E-05  0.0000E+00 -1.2108E+00  0.0000E+00
            -1.7864E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1276.76066814131        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      519
 NPARAMETR:  9.9439E-01  1.1956E-01  6.8853E-01  1.5474E+00  5.5413E-01  7.7709E-01  2.3066E-02  1.0000E-02  8.2612E-01  1.0000E-02
             3.9173E+00
 PARAMETER:  9.4371E-02 -2.0239E+00 -2.7320E-01  5.3658E-01 -4.9035E-01 -1.5220E-01 -3.6694E+00 -1.0026E+01 -9.1014E-02 -1.1010E+01
             1.4654E+00
 GRADIENT:  -3.0461E+00  7.6018E-01  3.1065E+00  7.8116E+00 -4.9160E+00 -4.2883E-01 -1.6082E-05  0.0000E+00 -1.2619E+00  0.0000E+00
            -1.6442E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1276.83747070970        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.9523E-01  8.8177E-02  7.0816E-01  1.5707E+00  5.6088E-01  7.7779E-01  1.0378E-02  1.0000E-02  8.1880E-01  1.0000E-02
             3.9333E+00
 PARAMETER:  9.5222E-02 -2.3284E+00 -2.4508E-01  5.5150E-01 -4.7825E-01 -1.5130E-01 -4.4680E+00 -1.1944E+01 -9.9915E-02 -1.2961E+01
             1.4695E+00
 GRADIENT:   1.4559E-01  3.6548E-01  2.4169E+00  4.7576E+00 -3.7073E+00  9.3831E-02 -1.3273E-06  0.0000E+00  3.3833E-01  0.0000E+00
             7.9145E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1276.99421510141        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  9.9509E-01  3.8494E-02  5.1250E-01  1.5336E+00  4.2923E-01  7.8021E-01  1.0000E-02  1.0000E-02  8.7074E-01  1.0000E-02
             3.9051E+00
 PARAMETER:  9.5074E-02 -3.1572E+00 -5.6845E-01  5.2763E-01 -7.4575E-01 -1.4819E-01 -6.8138E+00 -1.8043E+01 -3.8416E-02 -1.8763E+01
             1.4623E+00
 GRADIENT:   2.8265E+00  6.0902E-01  2.6167E+00  2.4032E+01 -8.4746E+00 -6.9308E-01  0.0000E+00  0.0000E+00 -8.2673E-01  0.0000E+00
             3.5491E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1277.38940111223        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.9057E-01  1.0000E-02  4.7672E-01  1.5071E+00  4.0399E-01  7.8041E-01  1.0000E-02  1.0000E-02  8.8080E-01  1.0000E-02
             3.8702E+00
 PARAMETER:  9.0527E-02 -4.7433E+00 -6.4083E-01  5.1019E-01 -8.0636E-01 -1.4793E-01 -1.1412E+01 -2.9613E+01 -2.6921E-02 -2.9891E+01
             1.4533E+00
 GRADIENT:  -3.2866E-01  0.0000E+00  1.2210E-01 -1.5531E+00  4.9061E-01 -3.8760E-01  0.0000E+00  0.0000E+00  1.2447E-01  0.0000E+00
             5.0744E-01

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1277.39184640288        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1185
 NPARAMETR:  9.9073E-01  1.0000E-02  4.7072E-01  1.5048E+00  4.0014E-01  7.8193E-01  1.0000E-02  1.0000E-02  8.8298E-01  1.0000E-02
             3.8670E+00
 PARAMETER:  9.0691E-02 -4.7465E+00 -6.5349E-01  5.0866E-01 -8.1593E-01 -1.4599E-01 -1.1424E+01 -2.9659E+01 -2.4449E-02 -2.9925E+01
             1.4525E+00
 GRADIENT:   3.1291E-03  0.0000E+00  6.3791E-03 -5.5591E-03 -5.7063E-03 -4.1864E-04  0.0000E+00  0.0000E+00  5.9084E-03  0.0000E+00
            -3.9161E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1185
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.1070E-04 -3.3105E-06  1.2332E-04 -1.5541E-02  1.5324E-05
 SE:             2.7733E-02  1.5152E-06  1.8493E-04  2.5251E-02  2.6444E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7668E-01  2.8895E-02  5.0488E-01  5.3825E-01  9.5379E-01

 ETASHRINKSD(%)  7.0921E+00  9.9995E+01  9.9380E+01  1.5407E+01  9.9114E+01
 ETASHRINKVR(%)  1.3681E+01  1.0000E+02  9.9996E+01  2.8440E+01  9.9992E+01
 EBVSHRINKSD(%)  6.8228E+00  9.9996E+01  9.9272E+01  1.5059E+01  9.9041E+01
 EBVSHRINKVR(%)  1.3180E+01  1.0000E+02  9.9995E+01  2.7850E+01  9.9991E+01
 RELATIVEINF(%)  5.6475E+01  8.7906E-09  1.2452E-04  1.0891E+01  1.2633E-04
 EPSSHRINKSD(%)  2.0173E+01
 EPSSHRINKVR(%)  3.6277E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1277.3918464028836     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -542.24101983914545     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1277.392       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  1.00E-02  4.71E-01  1.50E+00  4.00E-01  7.82E-01  1.00E-02  1.00E-02  8.83E-01  1.00E-02  3.87E+00
 


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
+        1.64E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.97E+01  0.00E+00  2.47E+03
 
 TH 4
+       -1.09E+02  0.00E+00 -8.51E+01  4.54E+02
 
 TH 5
+        2.63E+02  0.00E+00 -3.73E+03 -2.92E+02  6.19E+03
 
 TH 6
+       -1.44E+01  0.00E+00  3.82E+01 -2.43E+01 -5.30E+00  2.60E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.75E+01  0.00E+00  4.78E+00 -1.59E+01  9.07E+01 -8.77E+00  0.00E+00  0.00E+00  1.74E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.73E+01  0.00E+00 -1.21E+01 -8.03E+00  1.99E+01  4.49E+00  0.00E+00  0.00E+00  6.23E+00  0.00E+00  3.06E+01
 
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
 #CPUT: Total CPU Time in Seconds,       18.894
Stop Time:
Sat Sep 25 09:15:59 CDT 2021
