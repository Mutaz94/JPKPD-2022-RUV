Sat Sep 25 07:18:35 CDT 2021
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
$DATA ../../../../data/spa/B/dat42.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1705.91940848043        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -6.4103E+01  1.4867E+01 -2.6224E+01  5.4165E+01  6.0649E+01 -5.5847E-01  4.2826E+00  2.4260E+00 -3.2427E+00  4.4938E+00
            -3.4834E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1711.33059036009        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0509E+00  8.8376E-01  9.4261E-01  1.0427E+00  8.7119E-01  1.0031E+00  8.8983E-01  9.7116E-01  1.0754E+00  8.0327E-01
             1.0460E+00
 PARAMETER:  1.4962E-01 -2.3571E-02  4.0900E-02  1.4185E-01 -3.7896E-02  1.0306E-01 -1.6724E-02  7.0733E-02  1.7265E-01 -1.1906E-01
             1.4495E-01
 GRADIENT:  -5.1701E+00  1.4241E+00  3.5617E+00  6.4021E+00 -4.7381E+00  9.5372E-01  2.1882E+00  1.0276E+00  1.7947E+01 -6.9534E+00
             1.4181E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1711.97201921883        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  1.0619E+00  7.6207E-01  9.1655E-01  1.1246E+00  8.1735E-01  1.0026E+00  7.9336E-01  8.6269E-01  9.9538E-01  8.2462E-01
             1.0252E+00
 PARAMETER:  1.6004E-01 -1.7172E-01  1.2864E-02  2.1740E-01 -1.0169E-01  1.0264E-01 -1.3148E-01 -4.7705E-02  9.5365E-02 -9.2829E-02
             1.2491E-01
 GRADIENT:   1.8621E+01  7.6844E+00 -7.2877E+00  2.6120E+01 -7.9944E-01  7.9310E-01 -1.1994E+00  9.8426E-01  6.8822E+00 -4.2425E+00
             6.5109E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1712.83904698320        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  1.0517E+00  7.2779E-01  1.0166E+00  1.1410E+00  8.5129E-01  9.9867E-01  9.0267E-01  9.1461E-01  9.4471E-01  8.7456E-01
             1.0044E+00
 PARAMETER:  1.5044E-01 -2.1775E-01  1.1646E-01  2.3190E-01 -6.1006E-02  9.8670E-02 -2.4010E-03  1.0744E-02  4.3120E-02 -3.4034E-02
             1.0436E-01
 GRADIENT:   9.3411E-01  4.4791E+00  3.0811E+00  4.4808E+00 -4.0440E+00 -1.0051E-01 -9.1810E-02 -1.1788E-01 -9.6945E-01 -8.0391E-01
            -2.0556E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1712.93334259693        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  1.0512E+00  6.6556E-01  9.5258E-01  1.1712E+00  8.0430E-01  9.9861E-01  1.1467E+00  7.7414E-01  8.9890E-01  8.5420E-01
             1.0084E+00
 PARAMETER:  1.4993E-01 -3.0713E-01  5.1424E-02  2.5807E-01 -1.1778E-01  9.8610E-02  2.3692E-01 -1.5601E-01 -6.5839E-03 -5.7595E-02
             1.0834E-01
 GRADIENT:  -6.4568E-02 -2.3115E-01 -4.6963E-02 -3.6704E-01  1.7139E-01 -6.5691E-03  5.2825E-02  2.5270E-02  1.2623E-01  3.4115E-02
             5.9954E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1712.95340976982        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      890
 NPARAMETR:  1.0529E+00  7.3090E-01  8.7036E-01  1.1278E+00  7.8716E-01  1.0003E+00  1.1986E+00  6.7946E-01  9.0157E-01  8.2817E-01
             1.0071E+00
 PARAMETER:  1.5151E-01 -2.1348E-01 -3.8843E-02  2.2027E-01 -1.3932E-01  1.0032E-01  2.8114E-01 -2.8646E-01 -3.6175E-03 -8.8542E-02
             1.0711E-01
 GRADIENT:   7.8970E-02  1.3888E+00 -5.2956E-01  3.5809E+00 -9.8094E-01  6.5106E-02  4.6571E-01  5.0662E-01 -7.7115E-02  2.2508E-01
            -1.1519E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1713.05905702087        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1067
 NPARAMETR:  1.0552E+00  8.8007E-01  7.4380E-01  1.0293E+00  7.7650E-01  1.0034E+00  1.1197E+00  4.5799E-01  9.4082E-01  8.1452E-01
             1.0119E+00
 PARAMETER:  1.5369E-01 -2.7749E-02 -1.9598E-01  1.2890E-01 -1.5296E-01  1.0336E-01  2.1304E-01 -6.8090E-01  3.8997E-02 -1.0516E-01
             1.1183E-01
 GRADIENT:  -4.0778E-01  7.0901E+00  4.0016E+00  4.5374E+00 -1.0768E+01  2.8740E-01  1.4625E+00  3.5334E-01  4.3384E-01  1.7668E+00
             2.3199E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1713.25228828192        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1245
 NPARAMETR:  1.0570E+00  9.9156E-01  6.6674E-01  9.5094E-01  7.8762E-01  1.0031E+00  1.0344E+00  2.7359E-01  9.8186E-01  8.0001E-01
             1.0027E+00
 PARAMETER:  1.5541E-01  9.1520E-02 -3.0535E-01  4.9696E-02 -1.3873E-01  1.0308E-01  1.3385E-01 -1.1961E+00  8.1696E-02 -1.2313E-01
             1.0271E-01
 GRADIENT:   1.1894E+00 -8.8059E-01 -1.9616E+00  3.6578E-01  2.3664E+00 -3.3542E-01  7.2422E-01  2.1109E-01 -1.9198E-01  3.5501E-01
            -1.1860E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1713.28783967995        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1421
 NPARAMETR:  1.0569E+00  1.0581E+00  6.3899E-01  9.0786E-01  8.0224E-01  1.0041E+00  9.7915E-01  1.8665E-01  1.0157E+00  7.9948E-01
             1.0049E+00
 PARAMETER:  1.5533E-01  1.5651E-01 -3.4786E-01  3.3349E-03 -1.2035E-01  1.0405E-01  7.8925E-02 -1.5785E+00  1.1561E-01 -1.2380E-01
             1.0485E-01
 GRADIENT:   4.1080E-01 -1.1875E+00 -9.7497E-01 -1.0832E+00  1.8217E+00 -8.3493E-02 -3.6849E-02  1.0333E-01 -1.8624E-01 -1.7967E-01
            -4.0495E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1713.34211552905        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1598
 NPARAMETR:  1.0566E+00  1.0268E+00  6.4117E-01  9.2686E-01  7.8895E-01  1.0042E+00  1.0048E+00  4.8142E-02  9.9949E-01  8.0054E-01
             1.0059E+00
 PARAMETER:  1.5505E-01  1.2646E-01 -3.4446E-01  2.4047E-02 -1.3705E-01  1.0421E-01  1.0479E-01 -2.9336E+00  9.9493E-02 -1.2247E-01
             1.0589E-01
 GRADIENT:  -1.7821E-01 -7.2608E-01 -6.4605E-01 -2.3429E-01  8.0142E-01  3.1535E-03  1.4573E-01  6.1338E-03  1.8597E-01  2.7430E-01
             1.7885E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1713.34542193192        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1773
 NPARAMETR:  1.0567E+00  1.0265E+00  6.4079E-01  9.2725E-01  7.8821E-01  1.0042E+00  1.0045E+00  1.0000E-02  9.9857E-01  7.9922E-01
             1.0057E+00
 PARAMETER:  1.5514E-01  1.2618E-01 -3.4506E-01  2.4472E-02 -1.3800E-01  1.0421E-01  1.0447E-01 -4.7136E+00  9.8565E-02 -1.2412E-01
             1.0569E-01
 GRADIENT:  -1.1004E-03 -6.7521E-03 -3.6385E-03 -2.7699E-03  4.8678E-03 -4.2036E-05  2.2658E-06  0.0000E+00  1.4182E-03  2.0702E-03
             1.2724E-03

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1713.34542193192        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1795
 NPARAMETR:  1.0567E+00  1.0265E+00  6.4079E-01  9.2725E-01  7.8821E-01  1.0042E+00  1.0045E+00  1.0000E-02  9.9857E-01  7.9922E-01
             1.0057E+00
 PARAMETER:  1.5514E-01  1.2618E-01 -3.4506E-01  2.4472E-02 -1.3800E-01  1.0421E-01  1.0447E-01 -4.7136E+00  9.8565E-02 -1.2412E-01
             1.0569E-01
 GRADIENT:  -1.1004E-03 -6.7521E-03 -3.6385E-03 -2.7699E-03  4.8678E-03 -4.2036E-05  2.2658E-06  0.0000E+00  1.4182E-03  2.0702E-03
             1.2724E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1795
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.7210E-06 -1.0416E-02 -4.0649E-04  4.1171E-03 -1.6001E-02
 SE:             2.9850E-02  2.0596E-02  1.8356E-04  2.5606E-02  2.2994E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9993E-01  6.1304E-01  2.6793E-02  8.7226E-01  4.8652E-01

 ETASHRINKSD(%)  1.0000E-10  3.1001E+01  9.9385E+01  1.4218E+01  2.2966E+01
 ETASHRINKVR(%)  1.0000E-10  5.2392E+01  9.9996E+01  2.6415E+01  4.0658E+01
 EBVSHRINKSD(%)  4.2799E-01  3.0709E+01  9.9443E+01  1.4379E+01  2.2610E+01
 EBVSHRINKVR(%)  8.5415E-01  5.1987E+01  9.9997E+01  2.6690E+01  4.0108E+01
 RELATIVEINF(%)  9.8962E+01  2.3915E+00  3.1041E-04  5.0722E+00  4.7444E+00
 EPSSHRINKSD(%)  4.4437E+01
 EPSSHRINKVR(%)  6.9127E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1713.3454219319153     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -978.19459536817715     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1713.345       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.03E+00  6.41E-01  9.27E-01  7.88E-01  1.00E+00  1.00E+00  1.00E-02  9.99E-01  7.99E-01  1.01E+00
 


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
+        9.78E+02
 
 TH 2
+       -6.70E+00  5.23E+02
 
 TH 3
+        1.44E+01  3.51E+02  8.92E+02
 
 TH 4
+       -1.35E+01  3.60E+02 -3.55E+02  9.54E+02
 
 TH 5
+       -2.14E+00 -5.72E+02 -1.06E+03  3.60E+02  1.66E+03
 
 TH 6
+        8.15E-03 -1.45E-01  2.79E+00 -2.43E+00 -5.68E-02  1.94E+02
 
 TH 7
+        9.60E-01  2.36E+01 -9.23E+00 -6.03E+00 -1.77E+01 -6.87E-01  4.75E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.88E-01 -2.66E+01 -3.22E+01  3.43E+01 -1.08E+00  4.74E-01  2.20E+01  0.00E+00  1.11E+02
 
 TH10
+       -1.05E+00 -9.76E+00 -8.28E+01 -1.99E+01 -5.14E+01 -1.75E-01  2.02E+01  0.00E+00  8.04E+00  1.21E+02
 
 TH11
+       -6.63E+00 -1.58E+01 -3.97E+01 -4.53E+00  7.07E+00  3.50E+00  8.22E+00  0.00E+00  9.24E+00  2.52E+01  2.06E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.600
Stop Time:
Sat Sep 25 07:19:04 CDT 2021
