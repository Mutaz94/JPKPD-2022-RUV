Wed Sep 29 15:39:14 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat19.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1724.10025323237        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9772E+02 -1.1117E+01 -2.3752E+01  6.1962E+01  6.6519E+01  5.8041E+01  9.2185E+00  3.3052E+00  4.8944E+01 -1.0173E+01
             2.4015E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1729.72062133422        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0366E+00  1.0897E+00  9.5132E-01  9.5881E-01  9.9302E-01  9.4819E-01  9.4655E-01  9.8766E-01  6.8093E-01  1.0731E+00
             8.8348E-01
 PARAMETER:  1.3597E-01  1.8589E-01  5.0096E-02  5.7942E-02  9.2991E-02  4.6801E-02  4.5068E-02  8.7584E-02 -2.8429E-01  1.7051E-01
            -2.3892E-02
 GRADIENT:   2.1782E+01 -1.0962E+01 -2.0216E+00 -8.1895E+00  1.4310E+01 -1.2435E+01 -1.7891E+01  8.6675E-01 -2.2082E+01 -5.5205E-01
            -3.3328E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1730.84858400215        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.0307E+00  1.0267E+00  8.3922E-01  1.0025E+00  9.0605E-01  9.6798E-01  1.1416E+00  8.2883E-01  6.0710E-01  9.4954E-01
             9.1881E-01
 PARAMETER:  1.3023E-01  1.2637E-01 -7.5284E-02  1.0252E-01  1.3399E-03  6.7456E-02  2.3246E-01 -8.7740E-02 -3.9905E-01  4.8221E-02
             1.5327E-02
 GRADIENT:   4.0246E+00  1.3638E+01 -1.6081E+01  3.1589E+01  3.0867E+01 -3.9759E+00 -3.6170E+00  1.2556E+00 -2.3207E+01 -5.4666E+00
            -1.5239E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1734.85338460993        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0313E+00  9.3076E-01  7.2020E-01  1.0399E+00  7.8141E-01  9.8239E-01  1.2044E+00  3.7879E-01  7.1806E-01  9.0134E-01
             9.4576E-01
 PARAMETER:  1.3082E-01  2.8251E-02 -2.2823E-01  1.3917E-01 -1.4665E-01  8.2235E-02  2.8600E-01 -8.7078E-01 -2.3121E-01 -3.8699E-03
             4.4230E-02
 GRADIENT:   2.7661E+00  5.3847E+00 -1.4871E+00  3.5286E+00 -3.4711E+00  1.7190E+00  3.3480E+00  6.9917E-01  9.2713E-01  2.8052E+00
             2.7584E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1734.99482840393        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      720
 NPARAMETR:  1.0300E+00  8.6081E-01  7.1289E-01  1.0755E+00  7.5116E-01  9.7735E-01  1.2570E+00  2.6460E-01  7.0000E-01  8.8330E-01
             9.3995E-01
 PARAMETER:  1.2955E-01 -4.9883E-02 -2.3843E-01  1.7274E-01 -1.8614E-01  7.7094E-02  3.2876E-01 -1.2296E+00 -2.5667E-01 -2.4087E-02
             3.8070E-02
 GRADIENT:   2.3668E-01 -8.2176E-01 -6.3986E-01 -1.2137E-01  1.2815E-01 -9.3738E-02 -2.5298E-02  1.4300E-01  8.9999E-02  2.9005E-01
             2.3056E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1735.10037686981        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      899
 NPARAMETR:  1.0307E+00  9.4476E-01  6.6853E-01  1.0214E+00  7.6015E-01  9.7912E-01  1.1569E+00  8.9360E-02  7.2786E-01  8.8151E-01
             9.4046E-01
 PARAMETER:  1.3023E-01  4.3178E-02 -3.0267E-01  1.2119E-01 -1.7424E-01  7.8903E-02  2.4574E-01 -2.3151E+00 -2.1765E-01 -2.6116E-02
             3.8609E-02
 GRADIENT:  -6.2523E-02 -2.3106E+00 -2.1245E+00 -2.2125E+00  1.0068E+00  8.6632E-02  1.6719E-01  3.2209E-02  1.0171E-01  1.5377E+00
             5.5936E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1735.12227793227        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  1.0308E+00  9.6655E-01  6.6496E-01  1.0093E+00  7.6701E-01  9.7907E-01  1.1344E+00  3.2014E-02  7.3586E-01  8.7909E-01
             9.4021E-01
 PARAMETER:  1.3034E-01  6.5978E-02 -3.0802E-01  1.0922E-01 -1.6526E-01  7.8847E-02  2.2608E-01 -3.3416E+00 -2.0671E-01 -2.8862E-02
             3.8343E-02
 GRADIENT:   3.0447E-02 -9.3740E-01 -4.1613E-01 -1.1859E+00  4.5555E-01 -7.5898E-03  1.2741E-02  3.2800E-03  6.0903E-02  2.1876E-01
             9.6870E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1735.12309643030        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1251
 NPARAMETR:  1.0308E+00  9.6003E-01  6.6667E-01  1.0138E+00  7.6485E-01  9.7902E-01  1.1416E+00  1.4624E-02  7.3339E-01  8.7812E-01
             9.4011E-01
 PARAMETER:  1.3030E-01  5.9211E-02 -3.0545E-01  1.1370E-01 -1.6807E-01  7.8801E-02  2.3245E-01 -4.1251E+00 -2.1007E-01 -2.9973E-02
             3.8237E-02
 GRADIENT:   3.9911E-03  3.3512E-02 -7.0798E-03  5.2330E-02 -3.2899E-02  2.6490E-04 -1.8579E-03  6.5011E-04 -2.7541E-03  1.4777E-02
             3.3373E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1735.12515427770        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1438             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0330E+00  9.6223E-01  6.6632E-01  1.0124E+00  7.6557E-01  9.7939E-01  1.1396E+00  1.0000E-02  7.3433E-01  8.7831E-01
             9.4010E-01
 PARAMETER:  1.3244E-01  6.1494E-02 -3.0598E-01  1.1234E-01 -1.6714E-01  7.9173E-02  2.3072E-01 -5.7034E+00 -2.0879E-01 -2.9758E-02
             3.8235E-02
 GRADIENT:   6.6567E+02  3.7534E+01  1.1088E+01  1.0653E+02  2.1253E+01  4.9879E+01  1.0545E+01  0.0000E+00  1.3344E+01  9.2835E-01
             7.6912E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1735.12536178631        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1510
 NPARAMETR:  1.0316E+00  9.6214E-01  6.6630E-01  1.0123E+00  7.6555E-01  9.7910E-01  1.1394E+00  1.0000E-02  7.3420E-01  8.7830E-01
             9.4010E-01
 PARAMETER:  1.3112E-01  6.1404E-02 -3.0601E-01  1.1225E-01 -1.6717E-01  7.8882E-02  2.3054E-01 -5.7034E+00 -2.0897E-01 -2.9770E-02
             3.8231E-02
 GRADIENT:   1.9076E+00 -1.8295E-01  1.8345E-01 -5.8547E-01 -1.9930E-01  2.9406E-02  1.0983E-02  0.0000E+00  1.0174E-02  8.1773E-03
            -1.6561E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1510
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.3071E-04 -8.6248E-04 -4.6118E-04 -3.2868E-03 -9.6636E-03
 SE:             2.9855E-02  2.1857E-02  2.0070E-04  2.3654E-02  2.3943E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8582E-01  9.6852E-01  2.1566E-02  8.8949E-01  6.8650E-01

 ETASHRINKSD(%)  1.0000E-10  2.6775E+01  9.9328E+01  2.0758E+01  1.9788E+01
 ETASHRINKVR(%)  1.0000E-10  4.6381E+01  9.9995E+01  3.7207E+01  3.5660E+01
 EBVSHRINKSD(%)  3.9778E-01  2.6726E+01  9.9409E+01  2.1221E+01  1.8640E+01
 EBVSHRINKVR(%)  7.9398E-01  4.6310E+01  9.9997E+01  3.7939E+01  3.3805E+01
 RELATIVEINF(%)  9.8905E+01  3.1701E+00  3.6003E-04  3.9533E+00  5.5265E+00
 EPSSHRINKSD(%)  4.4446E+01
 EPSSHRINKVR(%)  6.9137E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1735.1253617863074     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -999.97453522256922     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.60
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1735.125       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  9.62E-01  6.66E-01  1.01E+00  7.66E-01  9.79E-01  1.14E+00  1.00E-02  7.34E-01  8.78E-01  9.40E-01
 


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
+        1.08E+03
 
 TH 2
+       -7.72E+00  5.15E+02
 
 TH 3
+        1.98E+01  2.51E+02  9.79E+02
 
 TH 4
+       -1.52E+01  4.71E+02 -4.85E+02  1.25E+03
 
 TH 5
+       -3.31E+00 -4.40E+02 -1.09E+03  4.79E+02  1.62E+03
 
 TH 6
+        1.36E+00 -1.95E+00  5.76E+00 -3.75E+00 -1.63E+00  2.05E+02
 
 TH 7
+        1.43E+00  3.51E+01 -1.75E+01 -1.05E+01  4.20E-01  3.58E-01  5.15E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.91E+00 -2.61E+01 -3.47E+01  2.72E+01  9.55E+00 -5.29E-01  2.28E+01  0.00E+00  1.62E+02
 
 TH10
+       -1.21E+00 -1.33E+01 -8.97E+01 -2.07E+01 -5.16E+01 -1.67E-01  1.41E+01  0.00E+00  1.40E+01  1.13E+02
 
 TH11
+       -7.41E+00 -1.48E+01 -3.92E+01 -2.20E+00  1.11E+01  2.04E+00  5.18E+00  0.00E+00  1.35E+01  2.12E+01  2.40E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.284
Stop Time:
Wed Sep 29 15:39:39 CDT 2021
