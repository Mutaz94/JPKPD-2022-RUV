Wed Sep 29 18:30:14 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat80.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1717.43041807304        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8678E+02  6.4301E+00  1.0729E+01  1.9149E+01 -3.1645E+01  7.9916E+01  1.5703E+01  7.1953E+00  4.0183E+01  2.0241E+01
            -3.1830E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1727.17285381988        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0268E+00  1.0560E+00  1.0312E+00  1.0173E+00  1.0580E+00  8.8203E-01  9.3293E-01  9.6106E-01  8.4098E-01  9.1087E-01
             1.1276E+00
 PARAMETER:  1.2640E-01  1.5449E-01  1.3077E-01  1.1717E-01  1.5636E-01 -2.5531E-02  3.0572E-02  6.0278E-02 -7.3186E-02  6.6470E-03
             2.2008E-01
 GRADIENT:  -5.8777E-01  1.5270E+01  2.0764E+00  1.9268E+01  2.8378E+00 -1.2060E+01  2.4377E+00  7.1012E-01  1.6659E-01  2.3021E+00
             1.1129E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1728.03249197500        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0323E+00  1.0358E+00  8.4891E-01  1.0107E+00  9.4320E-01  9.1046E-01  9.7463E-01  7.9856E-01  8.0814E-01  7.4144E-01
             1.1216E+00
 PARAMETER:  1.3182E-01  1.3518E-01 -6.3800E-02  1.1066E-01  4.1518E-02  6.1924E-03  7.4305E-02 -1.2494E-01 -1.1302E-01 -1.9916E-01
             2.1475E-01
 GRADIENT:   1.0436E+01  3.1837E+00 -2.1630E+00  9.4291E+00  6.0735E+00 -4.8915E-02 -4.2343E-01  1.3858E+00 -2.8787E+00 -3.8259E+00
             9.7263E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1729.04400698033        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0275E+00  1.1802E+00  6.7057E-01  9.0124E-01  9.1183E-01  9.1351E-01  8.9252E-01  4.1226E-01  8.7889E-01  7.7124E-01
             1.0751E+00
 PARAMETER:  1.2712E-01  2.6569E-01 -2.9963E-01 -3.9787E-03  7.6971E-03  9.5410E-03 -1.3703E-02 -7.8609E-01 -2.9096E-02 -1.5975E-01
             1.7242E-01
 GRADIENT:  -5.8669E+00 -3.2506E+00  1.1855E+00 -4.8268E+00 -4.0837E+00  1.1299E-01  2.1039E-01  5.3710E-01  7.7501E-01  1.7192E+00
            -4.4771E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1729.40669439685        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  1.0305E+00  1.3583E+00  5.3875E-01  7.8419E-01  9.2411E-01  9.1417E-01  8.0930E-01  1.5506E-01  9.5018E-01  7.2114E-01
             1.0881E+00
 PARAMETER:  1.3005E-01  4.0625E-01 -5.1851E-01 -1.4311E-01  2.1079E-02  1.0261E-02 -1.1158E-01 -1.7640E+00  4.8899E-02 -2.2692E-01
             1.8446E-01
 GRADIENT:  -5.4525E-01  4.8019E-01 -1.3600E+00  2.5474E+00  2.2540E-01 -1.2075E-01 -4.2367E-01  1.0202E-01  1.6122E-01  1.7502E-01
            -2.9733E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1729.44040085889        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      900             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0322E+00  1.3982E+00  5.2856E-01  7.5753E-01  9.4298E-01  9.1475E-01  7.9476E-01  4.7815E-02  9.7617E-01  7.3145E-01
             1.0902E+00
 PARAMETER:  1.3170E-01  4.3521E-01 -5.3761E-01 -1.7769E-01  4.1287E-02  1.0891E-02 -1.2971E-01 -2.9404E+00  7.5880E-02 -2.1273E-01
             1.8640E-01
 GRADIENT:   4.8906E+02  2.4594E+02  7.4838E+00  6.3326E+01  6.7183E+00  3.7052E+01  5.2320E+00  1.6890E-02  3.7783E+00  8.7976E-01
             1.6788E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1729.44289274850        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  1.0309E+00  1.3987E+00  5.2875E-01  7.5834E-01  9.4310E-01  9.1427E-01  7.9317E-01  3.4099E-02  9.7472E-01  7.3172E-01
             1.0899E+00
 PARAMETER:  1.3041E-01  4.3554E-01 -5.3724E-01 -1.7662E-01  4.1418E-02  1.0376E-02 -1.3172E-01 -3.2785E+00  7.4391E-02 -2.1236E-01
             1.8612E-01
 GRADIENT:   4.0701E-01  7.0738E-02 -1.2907E-01  2.2756E-01 -8.1049E-02 -8.5228E-02 -5.1992E-02  4.3250E-03 -3.8333E-02 -3.0484E-02
            -2.7733E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1729.44393090242        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  1.0304E+00  1.3987E+00  5.2998E-01  7.5832E-01  9.4425E-01  9.1451E-01  7.9335E-01  1.1735E-02  9.7549E-01  7.3372E-01
             1.0900E+00
 PARAMETER:  1.2997E-01  4.3554E-01 -5.3491E-01 -1.7665E-01  4.2639E-02  1.0630E-02 -1.3150E-01 -4.3452E+00  7.5180E-02 -2.0963E-01
             1.8618E-01
 GRADIENT:  -7.5064E-01 -2.6955E-01 -1.1908E-01 -3.6699E-02  6.7501E-03  1.7139E-02  3.3326E-02  5.2060E-04  2.0906E-02  3.3266E-02
             2.6994E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1729.44559896782        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1400
 NPARAMETR:  1.0321E+00  1.3984E+00  5.3146E-01  7.5823E-01  9.4561E-01  9.1475E-01  7.9302E-01  1.0000E-02  9.7638E-01  7.3552E-01
             1.0901E+00
 PARAMETER:  1.3159E-01  4.3531E-01 -5.3213E-01 -1.7677E-01  4.4071E-02  1.0895E-02 -1.3191E-01 -5.3976E+00  7.6096E-02 -2.0718E-01
             1.8623E-01
 GRADIENT:   3.6043E+00 -1.3423E+00 -1.2458E-01 -7.3054E-01  3.4298E-01  1.2704E-01  4.7648E-03  0.0000E+00  7.6531E-02  3.4555E-02
             2.1604E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1400
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.2201E-04 -1.5503E-02 -3.2842E-04  1.0359E-02 -2.3303E-02
 SE:             2.9808E-02  2.3263E-02  1.4903E-04  2.3760E-02  2.0828E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8870E-01  5.0514E-01  2.7544E-02  6.6284E-01  2.6323E-01

 ETASHRINKSD(%)  1.4020E-01  2.2065E+01  9.9501E+01  2.0400E+01  3.0222E+01
 ETASHRINKVR(%)  2.8021E-01  3.9261E+01  9.9998E+01  3.6638E+01  5.1311E+01
 EBVSHRINKSD(%)  5.7381E-01  2.1980E+01  9.9534E+01  2.0870E+01  3.0170E+01
 EBVSHRINKVR(%)  1.1443E+00  3.9129E+01  9.9998E+01  3.7384E+01  5.1237E+01
 RELATIVEINF(%)  9.8796E+01  2.6240E+00  1.3914E-04  2.9143E+00  4.5177E+00
 EPSSHRINKSD(%)  4.2767E+01
 EPSSHRINKVR(%)  6.7244E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1729.4455989678233     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -994.29477240408517     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.97
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1729.446       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.40E+00  5.31E-01  7.58E-01  9.46E-01  9.15E-01  7.93E-01  1.00E-02  9.76E-01  7.36E-01  1.09E+00
 


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
+        1.23E+03
 
 TH 2
+       -1.06E+01  5.66E+02
 
 TH 3
+        1.22E+01  3.32E+02  9.08E+02
 
 TH 4
+       -2.50E+01  3.99E+02 -5.50E+02  1.28E+03
 
 TH 5
+       -6.06E+00 -4.34E+02 -8.73E+02  4.60E+02  1.12E+03
 
 TH 6
+       -1.81E+00 -1.60E+00  3.02E+00 -4.77E+00 -1.09E+00  2.33E+02
 
 TH 7
+        6.20E-01  2.26E+01 -2.84E+01 -1.23E+01 -3.12E+00  3.74E-01  1.24E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.66E+00 -2.41E+01 -4.47E+01  5.82E+01 -7.03E+00 -4.74E-01  2.05E+01  0.00E+00  9.15E+01
 
 TH10
+       -4.70E-01 -1.29E+01 -5.52E+01 -2.08E+01 -6.91E+01 -6.98E-02  3.04E+01  0.00E+00  1.29E+01  9.33E+01
 
 TH11
+       -7.75E+00 -2.00E+01 -3.31E+01 -3.19E+00 -8.14E+00  2.69E+00  1.02E+01  0.00E+00  1.28E+01  2.51E+01  1.82E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.921
Stop Time:
Wed Sep 29 18:30:39 CDT 2021
