Sat Sep 25 11:17:46 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat73.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m73.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1682.18562505738        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1066E+01 -1.6163E+01  3.7277E+01 -7.0691E+01 -7.1004E+01  5.5283E+01 -1.0946E+01 -1.4902E+01 -1.8000E+01  3.6220E+01
            -2.0681E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1699.53068562597        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0276E+00  1.1125E+00  9.8254E-01  9.0003E-01  1.1298E+00  8.4071E-01  1.1913E+00  1.2660E+00  1.0750E+00  6.1197E-01
             1.0441E+00
 PARAMETER:  1.2725E-01  2.0663E-01  8.2381E-02 -5.3222E-03  2.2203E-01 -7.3512E-02  2.7504E-01  3.3589E-01  1.7236E-01 -3.9107E-01
             1.4316E-01
 GRADIENT:   8.7335E+01 -7.0519E+01 -1.1765E+01 -6.0511E+01  7.4409E+01 -5.1484E+00  1.1623E+00 -1.0967E+01 -1.2014E+00  6.6924E-01
            -5.8877E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1705.65455055698        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0301E+00  1.2382E+00  9.3554E-01  8.6418E-01  1.0832E+00  8.8985E-01  1.0328E+00  1.5946E+00  1.1410E+00  4.4867E-01
             1.0506E+00
 PARAMETER:  1.2961E-01  3.1368E-01  3.3369E-02 -4.5972E-02  1.7989E-01 -1.6705E-02  1.3231E-01  5.6661E-01  2.3189E-01 -7.0146E-01
             1.4939E-01
 GRADIENT:   8.8177E+01 -3.2024E+00  1.1527E+01 -2.4504E+01 -1.8859E+01  1.6437E+01 -3.1033E+00 -2.1650E+00  2.4200E+00  2.1426E-01
            -1.5511E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1706.76046790664        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0186E+00  1.3089E+00  9.0034E-01  8.3666E-01  1.1081E+00  8.6437E-01  1.0159E+00  1.6810E+00  1.1312E+00  4.7831E-01
             1.0511E+00
 PARAMETER:  1.1840E-01  3.6920E-01 -4.9864E-03 -7.8336E-02  2.0269E-01 -4.5753E-02  1.1580E-01  6.1938E-01  2.2327E-01 -6.3750E-01
             1.4985E-01
 GRADIENT:   5.1610E+01  1.5174E+01  3.5904E+00 -3.7633E+00 -6.0394E+00  6.0878E+00 -1.1077E+00  3.0604E-01 -2.4635E+00  7.7410E-01
            -1.5421E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1706.76091504539        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0176E+00  1.2989E+00  9.0089E-01  8.4285E-01  1.1034E+00  8.6354E-01  1.0234E+00  1.6650E+00  1.1250E+00  4.7366E-01
             1.0513E+00
 PARAMETER:  1.1747E-01  3.6152E-01 -4.3710E-03 -7.0972E-02  1.9837E-01 -4.6720E-02  1.2313E-01  6.0982E-01  2.1775E-01 -6.4726E-01
             1.5003E-01
 GRADIENT:   4.8563E+01  1.4254E+01  3.4096E+00 -3.6221E+00 -5.7759E+00  5.7305E+00 -1.0554E+00  2.9004E-01 -2.3293E+00  7.5043E-01
            -1.4503E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1706.97851838144        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.0157E+00  1.3977E+00  7.9698E-01  7.8658E-01  1.1073E+00  8.5298E-01  9.8819E-01  1.5984E+00  1.2162E+00  4.6671E-01
             1.0575E+00
 PARAMETER:  1.1555E-01  4.3483E-01 -1.2693E-01 -1.4006E-01  2.0191E-01 -5.9020E-02  8.8124E-02  5.6898E-01  2.9569E-01 -6.6204E-01
             1.5591E-01
 GRADIENT:  -4.9407E+00  7.0289E+00  9.7785E-01  6.8043E+00 -4.0338E+00 -2.2362E+00  6.4797E-01 -4.6206E-01  1.7827E+00  5.0360E-01
             9.1220E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1707.14109701303        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      624
 NPARAMETR:  1.0186E+00  1.4740E+00  7.3194E-01  7.3275E-01  1.1193E+00  8.6016E-01  9.5132E-01  1.6206E+00  1.2554E+00  4.5502E-01
             1.0575E+00
 PARAMETER:  1.1846E-01  4.8800E-01 -2.1206E-01 -2.1095E-01  2.1274E-01 -5.0636E-02  5.0096E-02  5.8281E-01  3.2743E-01 -6.8741E-01
             1.5586E-01
 GRADIENT:   2.9357E+00  5.8219E+00  8.4940E-01  4.2692E+00 -7.2429E-01  9.2722E-01 -7.8404E-02 -5.6966E-01 -5.7796E-01 -2.7326E-01
            -2.8854E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1707.14324085713        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:      815            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0186E+00  1.4744E+00  7.3163E-01  7.3249E-01  1.1194E+00  8.5810E-01  9.5155E-01  1.6208E+00  1.2567E+00  4.5500E-01
             1.0582E+00
 PARAMETER:  1.1847E-01  4.8825E-01 -2.1248E-01 -2.1130E-01  2.1277E-01 -5.3032E-02  5.0338E-02  5.8295E-01  3.2848E-01 -6.8747E-01
             1.5654E-01
 GRADIENT:   4.8552E+01  3.7987E+01  1.0564E+00  1.1183E+01  2.1267E-01  2.8471E+00  5.5378E-01 -3.6839E-01  9.5184E-01 -1.7720E-01
             1.4240E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1707.14338598447        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  1.0187E+00  1.4744E+00  7.3164E-01  7.3251E-01  1.1194E+00  8.5810E-01  9.5142E-01  1.6209E+00  1.2571E+00  4.5496E-01
             1.0582E+00
 PARAMETER:  1.1847E-01  4.8825E-01 -2.1248E-01 -2.1130E-01  2.1277E-01 -5.3032E-02  5.0305E-02  5.8295E-01  3.2885E-01 -6.8747E-01
             1.5653E-01
 GRADIENT:  -4.8269E+05 -3.8336E+00 -2.6911E+05 -2.7064E+05  2.6877E+05 -4.4806E-04  3.9123E-02 -4.9097E+04  1.7388E+05  8.3182E+04
             2.2575E-02
 NUMSIGDIG:         3.3         8.1         3.3         3.3         3.3         4.9         2.3         3.3         3.3         3.3
                    3.5

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      915
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6580E-05 -1.6495E-02 -2.8606E-02  1.2272E-02 -3.3382E-02
 SE:             2.9823E-02  2.3764E-02  1.6546E-02  2.2973E-02  1.2734E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9929E-01  4.8760E-01  8.3833E-02  5.9321E-01  8.7534E-03

 ETASHRINKSD(%)  8.8953E-02  2.0389E+01  4.4569E+01  2.3037E+01  5.7340E+01
 ETASHRINKVR(%)  1.7783E-01  3.6621E+01  6.9274E+01  4.0767E+01  8.1801E+01
 EBVSHRINKSD(%)  6.4703E-01  2.0711E+01  4.3784E+01  2.3292E+01  5.8341E+01
 EBVSHRINKVR(%)  1.2899E+00  3.7132E+01  6.8397E+01  4.1159E+01  8.2645E+01
 RELATIVEINF(%)  9.8575E+01  2.8189E+00  1.7292E+00  2.8251E+00  3.7820E+00
 EPSSHRINKSD(%)  4.2795E+01
 EPSSHRINKVR(%)  6.7275E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1707.1433859844728     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -971.99255942073466     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1707.143       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.47E+00  7.32E-01  7.32E-01  1.12E+00  8.58E-01  9.52E-01  1.62E+00  1.26E+00  4.55E-01  1.06E+00
 


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
+        9.82E+08
 
 TH 2
+       -8.38E+00  5.52E+07
 
 TH 3
+       -2.83E+03  1.40E+02  5.91E+08
 
 TH 4
+        7.65E+08  2.77E+02  6.29E+04  5.97E+08
 
 TH 5
+        1.84E+03 -2.40E+02 -3.86E+08 -3.88E+08  5.04E+08
 
 TH 6
+       -2.67E+03 -1.33E+00 -2.06E+03 -2.08E+03  1.35E+03  2.61E+02
 
 TH 7
+       -1.25E+04  1.99E+01 -9.69E+03 -9.73E+03  6.31E+03  1.50E+00  9.60E+01
 
 TH 8
+        2.42E+05 -1.35E+01  1.88E+05  1.89E+05 -1.23E+05 -3.39E+02 -1.59E+03  1.60E+07
 
 TH 9
+        3.58E+04 -1.23E+01  2.78E+04  2.79E+04 -1.82E+04  7.76E+02  3.65E+03 -3.66E+07  8.36E+07
 
 TH10
+        1.00E+04  2.12E+00  7.78E+03  7.79E+03 -5.17E+03  1.03E+03  4.82E+03 -4.84E+07  1.11E+08  1.46E+08
 
 TH11
+        3.59E+04 -5.27E+00  2.79E+04  2.80E+04 -1.83E+04  4.47E+00  5.89E+00  4.60E+03 -1.05E+04 -1.38E+04  1.96E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.795
Stop Time:
Sat Sep 25 11:18:15 CDT 2021
