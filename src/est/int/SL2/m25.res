Sat Sep 25 01:04:14 CDT 2021
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
$DATA ../../../../data/int/SL2/dat25.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      999
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      899
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1903.33327314672        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.9769E+01 -4.5539E+01  1.0834E+02 -9.2195E+00  3.7270E+01  3.2215E+01 -6.7273E+01 -9.5771E+01 -8.2799E+01  1.1555E+01
            -3.7288E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3016.65681599634        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0060E+00  1.2842E+00  9.7133E-01  8.5220E-01  1.1709E+00  8.1650E-01  9.7587E-01  8.7007E-01  1.0539E+00  8.1430E-01
             1.9781E+00
 PARAMETER:  1.0598E-01  3.5012E-01  7.0911E-02 -5.9933E-02  2.5776E-01 -1.0273E-01  7.5570E-02 -3.9187E-02  1.5247E-01 -1.0542E-01
             7.8214E-01
 GRADIENT:   4.0751E+01 -1.7380E+01 -7.6185E+00 -4.8793E+01  2.2552E+00 -4.6018E+01  5.9584E-01  4.5476E+00 -2.2470E+01 -2.3884E+01
            -2.2606E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3028.61246127574        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0117E+00  1.4731E+00  8.5069E-01  7.8562E-01  1.2758E+00  9.3689E-01  7.1811E-01  1.6802E-01  1.2928E+00  1.0584E+00
             2.0471E+00
 PARAMETER:  1.1160E-01  4.8735E-01 -6.1711E-02 -1.4128E-01  3.4358E-01  3.4810E-02 -2.3113E-01 -1.6837E+00  3.5681E-01  1.5673E-01
             8.1640E-01
 GRADIENT:   4.6043E+01  1.6506E+01 -2.8627E+01  2.6142E+01  9.8450E+00  9.8358E+00 -1.1798E+01  2.3635E-02  2.3827E+00 -3.1334E+00
            -1.4453E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3035.02938112328        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      254
 NPARAMETR:  9.9507E-01  1.3824E+00  9.4896E-01  8.2792E-01  1.2330E+00  9.0775E-01  8.5312E-01  1.2184E-01  1.1745E+00  1.0053E+00
             2.1697E+00
 PARAMETER:  9.5058E-02  4.2381E-01  4.7615E-02 -8.8843E-02  3.0945E-01  3.2096E-03 -5.8860E-02 -2.0050E+00  2.6084E-01  1.0524E-01
             8.7458E-01
 GRADIENT:  -1.0553E+01 -8.0820E+00  5.9980E-01 -2.1596E+00 -8.3102E-01 -1.3770E+00 -1.2712E+00 -2.6783E-02 -1.8691E+00 -1.4730E-01
             5.0804E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3035.96279399654        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      430
 NPARAMETR:  9.9585E-01  1.6682E+00  8.9942E-01  6.5009E-01  1.4579E+00  9.0678E-01  7.8213E-01  1.7342E-01  1.3579E+00  1.1267E+00
             2.1725E+00
 PARAMETER:  9.5844E-02  6.1173E-01 -6.0069E-03 -3.3065E-01  4.7699E-01  2.1410E-03 -1.4573E-01 -1.6520E+00  4.0594E-01  2.1930E-01
             8.7589E-01
 GRADIENT:  -9.1973E+00 -1.4379E+00  2.4505E+00 -1.9742E+00 -8.7288E-01 -1.8360E+00 -9.2807E-01 -7.8231E-02 -2.0785E+00 -1.1186E+00
             8.6509E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3036.34519881355        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  9.9969E-01  1.9073E+00  7.6517E-01  4.9995E-01  1.6366E+00  9.1228E-01  7.4489E-01  2.0941E-01  1.6089E+00  1.2314E+00
             2.1577E+00
 PARAMETER:  9.9693E-02  7.4567E-01 -1.6766E-01 -5.9324E-01  5.9264E-01  8.1892E-03 -1.9452E-01 -1.4635E+00  5.7554E-01  3.0812E-01
             8.6905E-01
 GRADIENT:   9.9056E-01  6.1739E+00  1.2208E-02  3.7015E+00  8.1357E-01  2.5540E-01 -4.2607E-01 -2.9430E-03  2.9044E-01  3.5278E-01
            -1.4398E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3036.41377494652        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  9.9946E-01  2.0209E+00  6.8037E-01  4.2038E-01  1.7238E+00  9.1191E-01  7.3423E-01  2.1018E-01  1.7645E+00  1.2773E+00
             2.1563E+00
 PARAMETER:  9.9456E-02  8.0353E-01 -2.8512E-01 -7.6659E-01  6.4453E-01  7.7879E-03 -2.0893E-01 -1.4598E+00  6.6787E-01  3.4477E-01
             8.6837E-01
 GRADIENT:   4.1417E-01 -3.9949E-01 -1.2918E-01 -4.9251E-02 -1.9891E-01  1.0933E-01  2.2634E-01  4.7292E-02  1.1348E-01  9.7481E-02
             1.4116E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3036.43561798825        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      958
 NPARAMETR:  9.9955E-01  2.0348E+00  6.6651E-01  4.1135E-01  1.7336E+00  9.1205E-01  7.3301E-01  5.1552E-02  1.7891E+00  1.2816E+00
             2.1553E+00
 PARAMETER:  9.9549E-02  8.1040E-01 -3.0570E-01 -7.8830E-01  6.5021E-01  7.9355E-03 -2.1059E-01 -2.8652E+00  6.8170E-01  3.4807E-01
             8.6795E-01
 GRADIENT:   6.5434E-01  2.3593E-04 -2.6081E-01  2.5074E-01  1.8534E-02  1.4083E-01  3.2392E-01  3.2448E-03  2.3619E-01 -7.4431E-02
            -6.2185E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3036.43805242793        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1134
 NPARAMETR:  9.9932E-01  2.0319E+00  6.7011E-01  4.1311E-01  1.7317E+00  9.1170E-01  7.3235E-01  1.0000E-02  1.7828E+00  1.2808E+00
             2.1560E+00
 PARAMETER:  9.9323E-02  8.0897E-01 -3.0032E-01 -7.8405E-01  6.4912E-01  7.5594E-03 -2.1150E-01 -4.5488E+00  6.7820E-01  3.4752E-01
             8.6826E-01
 GRADIENT:   4.3890E-02 -6.8965E-03  2.9046E-02 -3.2332E-02 -6.9654E-02  8.6059E-03  8.6181E-03  0.0000E+00  1.0747E-02  5.4560E-03
            -1.8262E-02

0ITERATION NO.:   44    OBJECTIVE VALUE:  -3036.43805433708        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1261
 NPARAMETR:  9.9931E-01  2.0324E+00  6.6964E-01  4.1284E-01  1.7321E+00  9.1169E-01  7.3228E-01  1.0000E-02  1.7833E+00  1.2810E+00
             2.1560E+00
 PARAMETER:  9.9309E-02  8.0919E-01 -3.0102E-01 -7.8470E-01  6.4935E-01  7.5402E-03 -2.1159E-01 -4.5247E+00  6.7846E-01  3.4766E-01
             8.6826E-01
 GRADIENT:  -2.3928E-04  8.6695E-03  1.3284E-03  4.1328E-04 -2.8457E-03 -6.6982E-05 -4.4405E-04  0.0000E+00 -4.0076E-04 -1.8495E-04
             1.3663E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1261
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2606E-03 -2.8788E-02 -5.8246E-05  2.9115E-02 -2.2491E-02
 SE:             2.9511E-02  2.3782E-02  6.1711E-05  2.1124E-02  2.5525E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6593E-01  2.2610E-01  3.4525E-01  1.6811E-01  3.7824E-01

 ETASHRINKSD(%)  1.1354E+00  2.0326E+01  9.9793E+01  2.9231E+01  1.4489E+01
 ETASHRINKVR(%)  2.2580E+00  3.6521E+01  1.0000E+02  4.9918E+01  2.6879E+01
 EBVSHRINKSD(%)  1.3463E+00  1.9235E+01  9.9798E+01  3.2841E+01  1.2649E+01
 EBVSHRINKVR(%)  2.6744E+00  3.4771E+01  1.0000E+02  5.4896E+01  2.3698E+01
 RELATIVEINF(%)  9.7294E+01  8.7842E+00  2.0156E-04  5.8525E+00  2.3885E+01
 EPSSHRINKSD(%)  1.6811E+01
 EPSSHRINKVR(%)  3.0795E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          899
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1652.2514827020016     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3036.4380543370830     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1384.1865716350815     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3036.438       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  2.03E+00  6.70E-01  4.13E-01  1.73E+00  9.12E-01  7.32E-01  1.00E-02  1.78E+00  1.28E+00  2.16E+00
 


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
+        1.31E+03
 
 TH 2
+       -9.12E+00  3.39E+02
 
 TH 3
+        1.14E+00  4.97E+01  1.23E+02
 
 TH 4
+       -1.51E+01  3.49E+02 -1.22E+02  9.43E+02
 
 TH 5
+       -3.21E+00 -7.06E+01 -4.39E+01  1.38E+02  1.72E+02
 
 TH 6
+        4.92E+00 -2.06E+00  2.43E+00 -6.53E+00 -1.05E+00  2.29E+02
 
 TH 7
+        2.92E+00 -7.67E+00 -1.17E+01  3.63E+00 -6.23E+00 -1.56E+00  1.79E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.95E-01 -6.92E+00 -5.12E+00  4.20E+01  8.28E-01 -8.64E-01  1.60E+01  0.00E+00  2.00E+01
 
 TH10
+        9.48E-01 -1.06E+01  3.65E+00  1.39E+01 -1.32E+01 -3.64E-01  2.92E+00  0.00E+00  2.55E+00  7.03E+01
 
 TH11
+       -1.65E+01 -1.32E+01 -4.92E+00 -1.10E+01  7.57E-01  3.04E+00  5.62E+00  0.00E+00  2.91E+00  7.71E+00  2.52E+02
 
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
 #CPUT: Total CPU Time in Seconds,       40.840
Stop Time:
Sat Sep 25 01:04:56 CDT 2021
