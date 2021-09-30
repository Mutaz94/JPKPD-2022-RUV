Wed Sep 29 16:13:48 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1620.64523252141        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2198E+02 -5.4011E+01 -5.4594E+01 -6.8643E+00  6.5732E+01  1.3103E+01 -2.7972E+00  1.5873E+01 -3.1665E+00  1.4123E+01
            -3.3091E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1633.29411162499        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.3639E-01  1.0652E+00  1.1447E+00  1.0186E+00  1.0317E+00  1.0524E+00  1.0079E+00  9.2195E-01  1.0398E+00  9.2330E-01
             1.0190E+00
 PARAMETER:  3.4277E-02  1.6320E-01  2.3519E-01  1.1846E-01  1.3122E-01  1.5112E-01  1.0788E-01  1.8736E-02  1.3904E-01  2.0195E-02
             1.1885E-01
 GRADIENT:  -6.6226E-01  1.7007E-01 -4.4351E+00 -4.5329E-01  4.5105E+00  7.8587E+00  1.2776E+00  6.0953E+00 -2.9130E+00 -6.6928E+00
            -5.5242E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1634.49085788617        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.3450E-01  9.7069E-01  1.2204E+00  1.0863E+00  1.0220E+00  1.0383E+00  8.5849E-01  6.7182E-01  1.0754E+00  1.0260E+00
             1.0386E+00
 PARAMETER:  3.2253E-02  7.0249E-02  2.9921E-01  1.8282E-01  1.2173E-01  1.3759E-01 -5.2582E-02 -2.9776E-01  1.7272E-01  1.2569E-01
             1.3783E-01
 GRADIENT:  -3.5211E+00  7.3078E+00  2.9554E+00  9.9303E+00 -8.0906E+00  3.2592E+00 -4.3310E-01 -2.0832E-02  2.8763E+00  1.1598E+00
             7.9001E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1634.93440782699        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  9.3675E-01  9.3804E-01  1.0829E+00  1.0966E+00  9.5832E-01  1.0287E+00  1.1164E+00  5.0636E-01  9.8121E-01  9.4077E-01
             1.0166E+00
 PARAMETER:  3.4664E-02  3.6040E-02  1.7968E-01  1.9223E-01  5.7428E-02  1.2827E-01  2.1008E-01 -5.8051E-01  8.1031E-02  3.8947E-02
             1.1648E-01
 GRADIENT:   6.8952E-01  6.9439E-01 -6.5607E-01  1.2815E+00  1.8276E+00 -6.8106E-01 -8.7967E-02  4.6839E-01 -1.5952E+00 -4.7379E-01
            -6.2958E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1635.03083223995        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.3699E-01  9.3490E-01  1.0054E+00  1.0916E+00  9.2318E-01  1.0309E+00  1.1530E+00  3.1674E-01  9.7336E-01  9.2030E-01
             1.0180E+00
 PARAMETER:  3.4917E-02  3.2679E-02  1.0535E-01  1.8766E-01  2.0074E-02  1.3047E-01  2.4234E-01 -1.0497E+00  7.2996E-02  1.6949E-02
             1.1785E-01
 GRADIENT:   2.1233E-01 -2.1121E+00 -2.0806E+00 -1.5707E+00  1.8946E+00  2.3903E-02  9.8260E-02  1.3780E-01 -4.6254E-01  1.2669E+00
             2.1554E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1635.12390709144        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.3786E-01  1.0931E+00  9.2985E-01  9.9093E-01  9.5520E-01  1.0324E+00  1.0371E+00  1.3329E-01  1.0453E+00  9.1548E-01
             1.0182E+00
 PARAMETER:  3.5841E-02  1.8904E-01  2.7265E-02  9.0893E-02  5.4164E-02  1.3190E-01  1.3645E-01 -1.9152E+00  1.4430E-01  1.1695E-02
             1.1808E-01
 GRADIENT:  -9.9472E-02  8.2503E-01 -7.9374E-02  7.3579E-01  1.8974E-01  4.0879E-02  3.2701E-01  2.0801E-02 -2.3791E-01 -1.7034E-01
            -1.6383E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1635.13136933370        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.3861E-01  1.1194E+00  9.2239E-01  9.7377E-01  9.6358E-01  1.0330E+00  1.0124E+00  6.8582E-02  1.0622E+00  9.2063E-01
             1.0185E+00
 PARAMETER:  3.6642E-02  2.1277E-01  1.9218E-02  7.3423E-02  6.2904E-02  1.3244E-01  1.1232E-01 -2.5797E+00  1.6035E-01  1.7298E-02
             1.1837E-01
 GRADIENT:   1.3375E+00  6.7947E-01  6.3631E-01  1.4268E-01 -6.3284E-01  2.0661E-01 -3.1648E-02  4.1002E-03 -5.4923E-02 -1.7307E-01
            -2.1116E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1635.13419531188        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  9.3864E-01  1.1203E+00  9.1843E-01  9.7279E-01  9.6273E-01  1.0330E+00  1.0126E+00  1.0000E-02  1.0622E+00  9.1957E-01
             1.0187E+00
 PARAMETER:  3.6676E-02  2.1364E-01  1.4908E-02  7.2412E-02  6.2017E-02  1.3245E-01  1.1250E-01 -4.5224E+00  1.6030E-01  1.6148E-02
             1.1854E-01
 GRADIENT:   1.3379E+00  2.2424E-02 -1.0397E-01  1.9820E-01  3.8566E-01  1.9583E-01 -5.7993E-02  1.4906E-04 -3.2963E-02 -8.1898E-02
            -4.3672E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1635.13423861766        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1290
 NPARAMETR:  9.3861E-01  1.1206E+00  9.1820E-01  9.7269E-01  9.6256E-01  1.0330E+00  1.0129E+00  1.0000E-02  1.0623E+00  9.1965E-01
             1.0188E+00
 PARAMETER:  3.6649E-02  2.1383E-01  1.4658E-02  7.2315E-02  6.1844E-02  1.3242E-01  1.1286E-01 -4.6630E+00  1.6046E-01  1.6243E-02
             1.1859E-01
 GRADIENT:   1.2777E+00  1.7494E-01 -2.3635E-02  2.1427E-01  1.1634E-01  1.8461E-01 -5.5751E-03  0.0000E+00  1.0903E-02 -2.1355E-02
            -3.8657E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1290
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2142E-05 -1.3643E-02 -3.3311E-04  3.9608E-03 -2.3306E-02
 SE:             2.9796E-02  1.9532E-02  1.5677E-04  2.4880E-02  2.3362E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9967E-01  4.8488E-01  3.3604E-02  8.7352E-01  3.1846E-01

 ETASHRINKSD(%)  1.8038E-01  3.4565E+01  9.9475E+01  1.6647E+01  2.1736E+01
 ETASHRINKVR(%)  3.6043E-01  5.7183E+01  9.9997E+01  3.0523E+01  3.8747E+01
 EBVSHRINKSD(%)  4.4758E-01  3.4049E+01  9.9516E+01  1.7021E+01  2.0388E+01
 EBVSHRINKVR(%)  8.9316E-01  5.6505E+01  9.9998E+01  3.1145E+01  3.6619E+01
 RELATIVEINF(%)  9.8819E+01  1.7929E+00  2.9786E-04  3.6641E+00  6.9470E+00
 EPSSHRINKSD(%)  4.2908E+01
 EPSSHRINKVR(%)  6.7405E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1635.1342386176582     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -899.98341205392001     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.04
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1635.134       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.39E-01  1.12E+00  9.18E-01  9.73E-01  9.63E-01  1.03E+00  1.01E+00  1.00E-02  1.06E+00  9.20E-01  1.02E+00
 


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
+        1.17E+03
 
 TH 2
+       -6.35E+00  3.95E+02
 
 TH 3
+        1.09E+01  1.70E+02  3.28E+02
 
 TH 4
+       -8.74E+00  3.37E+02 -1.44E+02  7.21E+02
 
 TH 5
+       -2.95E+00 -3.12E+02 -4.47E+02  1.84E+02  8.92E+02
 
 TH 6
+        1.07E+00 -1.46E+00  3.44E+00 -3.43E+00 -2.46E+00  1.83E+02
 
 TH 7
+        1.12E+00  1.78E+01  8.89E+00 -8.87E+00 -1.55E+01 -2.18E-01  4.23E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.07E+00 -2.26E+01 -2.04E+01  2.89E+01  4.94E+00 -5.10E-01  2.28E+01  0.00E+00  8.86E+01
 
 TH10
+       -7.60E-01 -7.31E+00 -3.89E+01 -1.27E+01 -5.66E+01  2.43E-01  1.15E+01  0.00E+00  5.83E+00  9.68E+01
 
 TH11
+       -9.13E+00 -1.69E+01 -3.65E+01 -2.84E+00  9.38E+00  2.81E+00  6.48E+00  0.00E+00  7.34E+00  2.56E+01  2.11E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.936
Stop Time:
Wed Sep 29 16:14:12 CDT 2021
