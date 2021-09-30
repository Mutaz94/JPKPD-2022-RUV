Wed Sep 29 17:10:27 CDT 2021
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
$DATA ../../../../data/spa/S2/dat16.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m16.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1711.27615989439        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5231E+02 -6.7851E+01 -6.0786E+01  2.0820E+01  1.1621E+02  7.7036E+01 -3.6625E+00  6.7396E+00  3.3787E+01 -1.1586E+01
             7.5521E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1724.62398007498        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0252E+00  1.1221E+00  1.0714E+00  9.7908E-01  9.8485E-01  8.3846E-01  1.0139E+00  9.7083E-01  8.5258E-01  1.0198E+00
             9.8439E-01
 PARAMETER:  1.2491E-01  2.1522E-01  1.6896E-01  7.8860E-02  8.4730E-02 -7.6183E-02  1.1378E-01  7.0396E-02 -5.9490E-02  1.1962E-01
             8.4272E-02
 GRADIENT:  -3.0771E+01  5.1693E+00  9.4785E+00 -3.1283E-01 -1.5088E+01 -2.5408E+01 -6.1823E+00  4.3634E-01 -2.7173E-01 -5.4767E+00
            -6.0403E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1726.40227959421        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0329E+00  1.0617E+00  9.8977E-01  1.0095E+00  9.5210E-01  8.7384E-01  1.1734E+00  5.2092E-01  7.8499E-01  1.0928E+00
             1.0076E+00
 PARAMETER:  1.3239E-01  1.5984E-01  8.9718E-02  1.0948E-01  5.0913E-02 -3.4864E-02  2.5992E-01 -5.5216E-01 -1.4208E-01  1.8876E-01
             1.0753E-01
 GRADIENT:  -7.2484E+00 -1.5929E-02  5.8542E-01 -2.6216E+00 -4.7458E-02 -7.2995E+00  2.0649E+00  1.9741E-01 -1.0669E+00  4.2157E+00
             4.3686E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1726.82016778542        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0378E+00  1.1776E+00  8.3502E-01  9.3188E-01  9.2434E-01  8.9773E-01  1.0724E+00  3.9602E-01  8.3711E-01  1.0086E+00
             9.9300E-01
 PARAMETER:  1.3712E-01  2.6344E-01 -8.0296E-02  2.9453E-02  2.1321E-02 -7.8836E-03  1.6994E-01 -8.2628E-01 -7.7804E-02  1.0854E-01
             9.2977E-02
 GRADIENT:   2.2564E+00  1.9737E+00  4.6577E-01  1.9490E+00 -2.0111E+00  2.4886E+00  5.5710E-01  2.7401E-01  6.3764E-01 -3.4831E-01
            -7.0956E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1726.94114667362        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  1.0376E+00  1.3294E+00  7.3989E-01  8.3336E-01  9.5313E-01  8.9210E-01  9.7073E-01  1.7277E-01  8.9748E-01  1.0178E+00
             9.9540E-01
 PARAMETER:  1.3690E-01  3.8471E-01 -2.0126E-01 -8.2293E-02  5.1993E-02 -1.4182E-02  7.0290E-02 -1.6558E+00 -8.1681E-03  1.1760E-01
             9.5388E-02
 GRADIENT:  -6.8037E-01  9.1635E-01 -5.6117E-01  1.6777E+00  3.6580E-01 -4.0135E-01 -1.9525E-01  5.8688E-02 -1.4178E-01  1.8707E-01
             1.9514E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1726.96048524188        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  1.0389E+00  1.3374E+00  7.3584E-01  8.2549E-01  9.5535E-01  8.9348E-01  9.6662E-01  6.7666E-02  9.0352E-01  1.0170E+00
             9.9496E-01
 PARAMETER:  1.3814E-01  3.9076E-01 -2.0674E-01 -9.1775E-02  5.4326E-02 -1.2628E-02  6.6046E-02 -2.5932E+00 -1.4544E-03  1.1686E-01
             9.4943E-02
 GRADIENT:   2.8952E+00 -1.7393E+00  1.2578E+00 -3.2269E+00 -1.2030E+00  2.1504E-01 -1.1785E-02  7.2591E-03 -2.2634E-01 -4.4884E-01
            -2.6379E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1726.96517345861        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  1.0381E+00  1.3383E+00  7.3430E-01  8.2673E-01  9.5520E-01  8.9294E-01  9.6625E-01  4.9701E-02  9.0393E-01  1.0194E+00
             9.9522E-01
 PARAMETER:  1.3735E-01  3.9143E-01 -2.0884E-01 -9.0277E-02  5.4167E-02 -1.3234E-02  6.5671E-02 -2.9017E+00 -1.0070E-03  1.1923E-01
             9.5208E-02
 GRADIENT:   5.4953E-01  1.6175E-01  6.1400E-02  2.1503E-01 -2.6662E-02 -4.0596E-02  7.0988E-03  4.5132E-03 -1.5157E-03  4.9564E-02
             7.9618E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1726.96565416788        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  1.0377E+00  1.3388E+00  7.3275E-01  8.2623E-01  9.5466E-01  8.9299E-01  9.6597E-01  2.6398E-02  9.0399E-01  1.0185E+00
             9.9517E-01
 PARAMETER:  1.3699E-01  3.9174E-01 -2.1095E-01 -9.0877E-02  5.3605E-02 -1.3175E-02  6.5373E-02 -3.5345E+00 -9.4084E-04  1.1836E-01
             9.5161E-02
 GRADIENT:  -4.7432E-01 -8.4863E-02 -2.5084E-02  9.4588E-02  7.1925E-02 -2.5265E-02 -1.1396E-02  1.3264E-03 -2.3433E-03  3.9533E-02
             1.7869E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1726.96605248945        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1362
 NPARAMETR:  1.0378E+00  1.3403E+00  7.3007E-01  8.2510E-01  9.5381E-01  8.9321E-01  9.6542E-01  1.0000E-02  9.0431E-01  1.0167E+00
             9.9507E-01
 PARAMETER:  1.3713E-01  3.9288E-01 -2.1462E-01 -9.2252E-02  5.2706E-02 -1.2929E-02  6.4804E-02 -4.6425E+00 -5.8354E-04  1.1658E-01
             9.5053E-02
 GRADIENT:  -1.4038E-01 -6.1225E-02 -1.1024E-02  2.5476E-02 -9.6727E-02  5.9309E-02  5.9654E-03  0.0000E+00 -9.0168E-04 -5.0862E-05
             2.2275E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1362
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.2889E-04 -1.2698E-02 -3.4690E-04  6.3257E-03 -2.1173E-02
 SE:             2.9819E-02  2.2971E-02  1.3527E-04  2.2067E-02  2.3796E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7515E-01  5.8043E-01  1.0334E-02  7.7438E-01  3.7357E-01

 ETASHRINKSD(%)  1.0286E-01  2.3043E+01  9.9547E+01  2.6072E+01  2.0281E+01
 ETASHRINKVR(%)  2.0561E-01  4.0776E+01  9.9998E+01  4.5346E+01  3.6449E+01
 EBVSHRINKSD(%)  5.1117E-01  2.2741E+01  9.9602E+01  2.7656E+01  1.8204E+01
 EBVSHRINKVR(%)  1.0197E+00  4.0310E+01  9.9998E+01  4.7664E+01  3.3093E+01
 RELATIVEINF(%)  9.8616E+01  2.3838E+00  1.3564E-04  1.9918E+00  8.3847E+00
 EPSSHRINKSD(%)  4.3784E+01
 EPSSHRINKVR(%)  6.8398E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1726.9660524894534     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -991.81522592571525     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1726.966       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.34E+00  7.30E-01  8.25E-01  9.54E-01  8.93E-01  9.65E-01  1.00E-02  9.04E-01  1.02E+00  9.95E-01
 


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
+        1.28E+03
 
 TH 2
+       -7.88E+00  4.10E+02
 
 TH 3
+        1.73E+01  1.25E+02  3.63E+02
 
 TH 4
+       -2.09E+01  4.21E+02 -2.88E+02  1.05E+03
 
 TH 5
+       -5.06E+00 -2.20E+02 -4.44E+02  3.22E+02  7.90E+02
 
 TH 6
+        2.89E+00 -1.41E+00  3.56E+00 -3.83E+00 -1.10E+00  2.45E+02
 
 TH 7
+        1.05E+00  2.14E+01 -4.03E+00 -1.61E+01 -3.53E+00 -1.90E-01  8.17E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.25E+00 -1.75E+01 -3.51E+01  3.66E+01  1.03E+01 -7.46E-01  2.83E+01  0.00E+00  7.73E+01
 
 TH10
+       -7.77E-01 -1.32E+01 -4.13E+01 -5.03E+00 -5.56E+01  1.83E-01  9.90E+00  0.00E+00  1.17E+01  8.49E+01
 
 TH11
+       -8.20E+00 -1.67E+01 -3.21E+01  3.80E+00  4.77E+00  2.70E+00  7.28E+00  0.00E+00  1.14E+01  1.83E+01  2.16E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.624
Stop Time:
Wed Sep 29 17:10:52 CDT 2021
