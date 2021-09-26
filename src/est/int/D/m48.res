Sat Sep 25 05:52:40 CDT 2021
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
$DATA ../../../../data/int/D/dat48.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m48.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   28993.1621826354        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3273E+02  4.0348E+02  7.7242E+00  2.7803E+02  7.6000E+01 -2.0543E+03 -1.0175E+03 -1.1332E+02 -1.4953E+03 -5.8730E+02
            -6.0145E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -841.954365203272        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2922E+00  2.0995E+00  9.0755E-01  2.2777E+00  9.8679E-01  5.0214E+00  5.2157E+00  1.0117E+00  2.9875E+00  1.6210E+00
             1.2549E+01
 PARAMETER:  3.5637E-01  8.4171E-01  2.9963E-03  9.2319E-01  8.6705E-02  1.7137E+00  1.7517E+00  1.1161E-01  1.1944E+00  5.8306E-01
             2.6296E+00
 GRADIENT:  -7.1105E+00  2.8015E+01 -4.7964E+01  1.0959E+02 -1.0148E+01  1.4660E+02  6.4666E+01  4.1523E+00  3.8556E+01  3.3382E+01
             4.3731E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -978.520586659753        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.3710E+00  1.3613E+00  1.4986E+01  3.4289E+00  1.8501E+00  1.9590E+00  7.2827E+00  1.6379E+00  3.9751E+00  2.5458E+00
             1.2216E+01
 PARAMETER:  4.1553E-01  4.0845E-01  2.8071E+00  1.3322E+00  7.1525E-01  7.7242E-01  2.0855E+00  5.9341E-01  1.4800E+00  1.0345E+00
             2.6027E+00
 GRADIENT:   5.2691E+01  2.0718E+01 -7.6033E+00  7.2312E+01 -4.3614E+01  3.7870E+01  1.3882E+01  1.9010E+00  4.9216E+01  7.1666E+01
             3.7676E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1136.74728632993        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.3408E+00  1.2793E+00  9.2854E+00  9.2988E-01  2.6969E+00  1.9449E+00  4.4235E+00  1.1345E+00  1.7794E+00  7.4210E-01
             1.0129E+01
 PARAMETER:  3.9328E-01  3.4630E-01  2.3284E+00  2.7301E-02  1.0921E+00  7.6520E-01  1.5869E+00  2.2616E-01  6.7626E-01 -1.9827E-01
             2.4154E+00
 GRADIENT:   6.6917E+01 -3.7056E+01 -6.5225E+00 -4.9303E+01  5.4630E+01  9.3929E+00  7.8035E+00  2.4088E-02  2.1679E+01  7.2169E+00
             2.9021E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1189.91573196451        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.1246E+00  1.2734E+00  1.6509E+01  1.1148E+00  2.3652E+00  1.9621E+00  4.7911E+00  1.1102E-01  1.2028E+00  4.9482E-01
             8.3035E+00
 PARAMETER:  2.1747E-01  3.4166E-01  2.9039E+00  2.0863E-01  9.6085E-01  7.7404E-01  1.6668E+00 -2.0980E+00  2.8464E-01 -6.0357E-01
             2.2167E+00
 GRADIENT:   7.6646E+00  3.0164E+00  6.8156E-02 -2.8305E+00 -2.5085E-02  1.1827E+01  1.6746E+00 -3.4660E-04  8.4370E+00  1.4328E+00
            -4.7836E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1192.78225316207        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.1129E+00  1.3063E+00  1.2958E+01  1.0386E+00  2.3562E+00  1.8659E+00  4.7535E+00  1.0961E-01  7.2067E-01  3.8611E-01
             8.5514E+00
 PARAMETER:  2.0699E-01  3.6718E-01  2.6617E+00  1.3787E-01  9.5704E-01  7.2376E-01  1.6589E+00 -2.1108E+00 -2.2758E-01 -8.5163E-01
             2.2461E+00
 GRADIENT:  -1.2457E+00 -3.6124E-01  4.5347E-01  7.6316E-02 -3.1949E-01  3.2186E-01 -3.4169E-01 -3.4216E-04 -2.1287E-01  4.9772E-01
             1.7836E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1193.19094589281        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.1166E+00  1.0867E+00  1.6164E+01  1.1611E+00  2.3715E+00  1.8703E+00  5.3128E+00  8.2612E-02  8.7521E-01  3.5689E-01
             8.5733E+00
 PARAMETER:  2.1031E-01  1.8313E-01  2.8828E+00  2.4941E-01  9.6352E-01  7.2611E-01  1.7701E+00 -2.3936E+00 -3.3296E-02 -9.3034E-01
             2.2487E+00
 GRADIENT:  -9.4271E-01 -8.0050E-01  3.4516E-02 -2.0406E+00  4.1146E-01  1.9771E-01 -1.3855E-02 -1.9504E-04 -1.1580E-02  3.6041E-01
             1.4435E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1193.25877432137        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      726
 NPARAMETR:  1.1182E+00  1.0866E+00  1.4376E+01  1.1679E+00  2.3531E+00  1.8703E+00  5.3272E+00  5.0448E-02  8.9079E-01  2.4568E-01
             8.5760E+00
 PARAMETER:  2.1170E-01  1.8307E-01  2.7655E+00  2.5521E-01  9.5572E-01  7.2612E-01  1.7728E+00 -2.8868E+00 -1.5643E-02 -1.3037E+00
             2.2490E+00
 GRADIENT:  -2.1190E-01 -1.6427E-01 -1.3716E-01  1.0084E-01  7.9538E-02  2.6860E-02  3.3030E-02 -1.3608E-04 -1.9273E-02  5.8568E-02
             1.4895E-06

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1193.27412522703        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  1.1184E+00  1.0816E+00  1.6408E+01  1.1751E+00  2.3831E+00  1.8698E+00  5.3453E+00  3.1910E-02  9.0561E-01  1.7219E-01
             8.5802E+00
 PARAMETER:  2.1194E-01  1.7841E-01  2.8978E+00  2.6134E-01  9.6840E-01  7.2584E-01  1.7762E+00 -3.3448E+00  8.5243E-04 -1.6592E+00
             2.2495E+00
 GRADIENT:  -4.0558E-03  2.4035E-03  4.7957E-03  4.6649E-03 -1.0331E-02  7.4670E-03 -7.9243E-03 -5.0254E-05 -7.5773E-04  1.9380E-03
             2.2196E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1193.27415262956        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1080            RESET HESSIAN, TYPE II
 NPARAMETR:  1.1184E+00  1.0823E+00  1.6319E+01  1.1746E+00  2.3824E+00  1.8698E+00  5.3438E+00  3.0569E-02  9.0505E-01  1.6640E-01
             8.5802E+00
 PARAMETER:  2.1194E-01  1.7908E-01  2.8923E+00  2.6093E-01  9.6813E-01  7.2584E-01  1.7759E+00 -3.3878E+00  2.3945E-04 -1.6933E+00
             2.2495E+00
 GRADIENT:   1.5780E+00  2.3014E-01  8.8692E-03  8.2940E-01  3.7997E-01  1.2977E+00  8.3581E+00 -4.6587E-05  1.5299E-02  1.7566E-03
             4.0401E+00

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1193.27415262956        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:     1141
 NPARAMETR:  1.1184E+00  1.0823E+00  1.6319E+01  1.1746E+00  2.3824E+00  1.8698E+00  5.3438E+00  3.0569E-02  9.0505E-01  1.6640E-01
             8.5802E+00
 PARAMETER:  2.1194E-01  1.7908E-01  2.8923E+00  2.6093E-01  9.6813E-01  7.2584E-01  1.7759E+00 -3.3878E+00  2.3945E-04 -1.6933E+00
             2.2495E+00
 GRADIENT:   3.1545E-03  6.2137E-04 -5.5101E-05  1.4363E-02 -1.4671E-04 -2.0154E-03 -5.0011E-02  5.2878E-05 -1.3322E-02 -3.9419E-05
            -6.5027E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1141
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2956E-02  2.7397E-02 -1.3781E-05 -6.4550E-02  2.5940E-04
 SE:             2.8683E-02  2.5068E-02  3.4673E-05  1.2252E-02  2.5211E-03
 N:                     100         100         100         100         100

 P VAL.:         6.5149E-01  2.7444E-01  6.9102E-01  1.3794E-07  9.1805E-01

 ETASHRINKSD(%)  3.9069E+00  1.6018E+01  9.9884E+01  5.8953E+01  9.1554E+01
 ETASHRINKVR(%)  7.6611E+00  2.9471E+01  1.0000E+02  8.3151E+01  9.9287E+01
 EBVSHRINKSD(%)  5.1510E+00  1.0513E+01  9.9883E+01  6.4112E+01  9.1749E+01
 EBVSHRINKVR(%)  1.0037E+01  1.9921E+01  1.0000E+02  8.7120E+01  9.9319E+01
 RELATIVEINF(%)  8.9755E+01  3.8808E+01  2.5575E-05  6.3260E+00  1.2680E-01
 EPSSHRINKSD(%)  6.6131E+00
 EPSSHRINKVR(%)  1.2789E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1193.2741526295642     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       460.81520713884652     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.96
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1193.274       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.12E+00  1.08E+00  1.63E+01  1.17E+00  2.38E+00  1.87E+00  5.34E+00  3.06E-02  9.05E-01  1.66E-01  8.58E+00
 


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
+        2.16E+02
 
 TH 2
+        3.19E+00  3.00E+01
 
 TH 3
+        4.23E-02  3.64E-02  8.55E-03
 
 TH 4
+       -9.31E+00  3.40E+01 -6.37E-02  1.73E+02
 
 TH 5
+       -1.01E+00 -5.75E+00 -5.25E-01 -8.33E+00  4.33E+01
 
 TH 6
+       -5.24E+00  5.84E-03  2.11E-03  2.56E+00  9.02E-01  4.04E+01
 
 TH 7
+        7.58E-01  2.99E+00 -1.28E-02 -1.30E+01  9.54E-01 -4.06E-02  4.16E+00
 
 TH 8
+        9.39E+00 -1.07E+00 -1.06E-02  3.83E+00 -5.36E-01  1.14E+00 -3.05E-01 -8.06E-01
 
 TH 9
+        6.18E+00  2.38E+00 -2.96E-02 -4.49E+01  5.24E+00  1.96E+00  2.99E+00 -9.64E+00  1.69E+01
 
 TH10
+       -1.26E+00 -8.93E-01 -5.97E-04 -4.40E+00  1.87E+00  7.06E-01  1.84E-01  1.83E+00  5.02E+00  1.80E+00
 
 TH11
+       -7.09E+00 -2.91E+00 -1.18E-03 -1.12E+01  1.39E-01  2.49E+00  4.09E-01  9.24E-02  3.43E+00  5.40E-01  1.41E+01
 
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
 #CPUT: Total CPU Time in Seconds,       48.376
Stop Time:
Sat Sep 25 05:53:30 CDT 2021
