Sat Sep 25 12:42:49 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat14.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m14.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1679.98164888552        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.5529E+01 -9.3227E+01 -5.6954E+01 -5.6177E+01  6.9502E+01 -1.8224E+01  3.7148E+00  1.2924E+01  2.8969E+01  1.5929E+01
            -2.7192E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1689.70919065204        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.8891E-01  1.1346E+00  1.1448E+00  1.0094E+00  1.0464E+00  1.0766E+00  9.4945E-01  9.1590E-01  8.2197E-01  8.6360E-01
             1.1359E+00
 PARAMETER:  8.8844E-02  2.2629E-01  2.3527E-01  1.0939E-01  1.4538E-01  1.7384E-01  4.8132E-02  1.2146E-02 -9.6046E-02 -4.6642E-02
             2.2739E-01
 GRADIENT:  -5.9120E+00  4.6468E+01  3.6017E+00  5.8228E+01  1.2921E+01  1.3782E+01 -5.6598E+00 -1.0127E+00 -9.6274E+00 -9.4532E+00
             1.2996E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1691.68882526539        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0008E+00  1.0485E+00  9.0055E-01  1.0454E+00  9.1434E-01  1.0769E+00  1.1799E+00  4.9685E-01  7.4420E-01  8.1336E-01
             1.0971E+00
 PARAMETER:  1.0085E-01  1.4735E-01 -4.7459E-03  1.4439E-01  1.0444E-02  1.7405E-01  2.6542E-01 -5.9947E-01 -1.9545E-01 -1.0659E-01
             1.9265E-01
 GRADIENT:   1.7117E+01  3.2700E+01 -1.3667E+01  6.3337E+01  2.6617E+01  1.3904E+01  6.1860E+00  9.0357E-01 -8.5416E+00 -3.2524E+00
             6.3776E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1693.29816355747        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9298E-01  1.0488E+00  8.4722E-01  1.0135E+00  8.8822E-01  1.0480E+00  1.0903E+00  3.5234E-01  8.1807E-01  8.2665E-01
             1.0729E+00
 PARAMETER:  9.2959E-02  1.4769E-01 -6.5796E-02  1.1345E-01 -1.8539E-02  1.4685E-01  1.8642E-01 -9.4315E-01 -1.0080E-01 -9.0375E-02
             1.7037E-01
 GRADIENT:   1.8251E+00  7.5764E-02 -6.8720E+00  9.5014E+00  9.5918E+00  2.0669E+00  1.5978E+00  8.2380E-01  1.3948E+00  1.6114E+00
             1.4060E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1693.36304720545        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.9195E-01  1.0515E+00  8.0333E-01  1.0030E+00  8.6254E-01  1.0416E+00  1.0830E+00  2.4219E-01  8.1604E-01  7.9805E-01
             1.0689E+00
 PARAMETER:  9.1922E-02  1.5025E-01 -1.1900E-01  1.0295E-01 -4.7871E-02  1.4076E-01  1.7976E-01 -1.3180E+00 -1.0329E-01 -1.2558E-01
             1.6665E-01
 GRADIENT:  -1.1779E+00 -2.4976E+00 -1.2140E+00 -2.4137E+00 -2.0260E-01 -8.4305E-01  7.5706E-02  4.2047E-01  7.5633E-01  8.8341E-01
             2.4063E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1693.57941829815        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.9259E-01  1.0545E+00  8.0004E-01  1.0023E+00  8.6186E-01  1.0436E+00  1.0836E+00  5.0618E-02  8.1324E-01  7.9635E-01
             1.0698E+00
 PARAMETER:  9.2563E-02  1.5302E-01 -1.2310E-01  1.0229E-01 -4.8659E-02  1.4268E-01  1.8026E-01 -2.8835E+00 -1.0673E-01 -1.2772E-01
             1.6752E-01
 GRADIENT:   1.1380E-01  7.1876E-01  2.1666E+00 -8.5033E-01 -1.5674E+00  4.0198E-02 -2.2703E-01  1.4290E-02 -3.1506E-01 -8.6918E-01
            -1.9473E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1693.94751889415        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      491
 NPARAMETR:  1.0058E+00  1.0671E+00  7.8964E-01  9.9606E-01  8.6179E-01  1.0519E+00  1.0736E+00  1.0000E-02  8.2182E-01  7.9871E-01
             1.0693E+00
 PARAMETER:  1.0583E-01  1.6499E-01 -1.3618E-01  9.6056E-02 -4.8739E-02  1.5061E-01  1.7100E-01 -4.5393E+00 -9.6236E-02 -1.2476E-01
             1.6703E-01
 GRADIENT:  -9.9609E+00 -2.8003E+00 -7.6035E-02 -2.9811E+00 -1.5399E+00 -2.5595E+00 -4.0750E-01  0.0000E+00  4.6055E-01  2.4993E-01
            -1.6061E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1694.02479171656        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      666
 NPARAMETR:  1.0109E+00  1.1637E+00  7.7362E-01  9.3988E-01  8.9795E-01  1.0583E+00  1.0006E+00  1.0000E-02  8.6013E-01  8.1136E-01
             1.0722E+00
 PARAMETER:  1.1089E-01  2.5160E-01 -1.5667E-01  3.7998E-02 -7.6413E-03  1.5667E-01  1.0056E-01 -6.1463E+00 -5.0666E-02 -1.0904E-01
             1.6969E-01
 GRADIENT:  -4.9954E-01  9.3436E-01  3.0293E-01  7.5865E-01 -4.8570E-01 -2.5150E-01 -1.5349E-01  0.0000E+00 -1.1231E-01 -1.5373E-01
            -1.1396E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1694.02559692634        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      793
 NPARAMETR:  1.0112E+00  1.1574E+00  7.7546E-01  9.4328E-01  8.9635E-01  1.0590E+00  1.0059E+00  1.0000E-02  8.5777E-01  8.1171E-01
             1.0722E+00
 PARAMETER:  1.1111E-01  2.4614E-01 -1.5430E-01  4.1608E-02 -9.4286E-03  1.5728E-01  1.0585E-01 -6.1180E+00 -5.3414E-02 -1.0861E-01
             1.6968E-01
 GRADIENT:  -6.1433E-03  3.5361E-03  1.1415E-03  4.1022E-03 -1.8172E-03 -1.9607E-03 -1.2536E-03  0.0000E+00 -6.7420E-04 -5.0407E-04
             2.5434E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      793
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.9668E-05 -6.2179E-03 -3.4719E-04 -4.5000E-05 -1.7417E-02
 SE:             2.9823E-02  2.2007E-02  1.5876E-04  2.3224E-02  2.2517E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9760E-01  7.7752E-01  2.8751E-02  9.9845E-01  4.3922E-01

 ETASHRINKSD(%)  8.8628E-02  2.6275E+01  9.9468E+01  2.2195E+01  2.4566E+01
 ETASHRINKVR(%)  1.7718E-01  4.5646E+01  9.9997E+01  3.9464E+01  4.3097E+01
 EBVSHRINKSD(%)  4.3795E-01  2.6060E+01  9.9494E+01  2.2683E+01  2.3718E+01
 EBVSHRINKVR(%)  8.7399E-01  4.5328E+01  9.9997E+01  4.0221E+01  4.1810E+01
 RELATIVEINF(%)  9.8836E+01  2.3411E+00  2.1499E-04  2.8384E+00  5.0335E+00
 EPSSHRINKSD(%)  4.2443E+01
 EPSSHRINKVR(%)  6.6872E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1694.0255969263449     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -958.87477036260668     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1694.026       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.16E+00  7.75E-01  9.43E-01  8.96E-01  1.06E+00  1.01E+00  1.00E-02  8.58E-01  8.12E-01  1.07E+00
 


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
+        9.62E+02
 
 TH 2
+       -6.95E+00  4.64E+02
 
 TH 3
+        1.40E+01  2.00E+02  5.18E+02
 
 TH 4
+       -1.38E+01  4.08E+02 -3.01E+02  1.03E+03
 
 TH 5
+       -5.08E+00 -3.67E+02 -6.91E+02  3.60E+02  1.22E+03
 
 TH 6
+        1.14E+00 -1.30E+00  3.21E+00 -4.04E+00 -1.46E+00  1.76E+02
 
 TH 7
+        5.88E-01  2.77E+01 -5.81E+00 -1.19E+01 -7.33E+00 -2.20E-01  6.43E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.12E+00 -2.21E+01 -3.12E+01  3.46E+01  3.27E+00  4.03E-01  2.24E+01  0.00E+00  1.10E+02
 
 TH10
+       -7.89E-01 -8.60E+00 -5.08E+01 -2.08E+01 -6.08E+01  8.42E-01  2.23E+01  0.00E+00  1.05E+01  1.06E+02
 
 TH11
+       -7.24E+00 -1.78E+01 -3.20E+01 -1.27E+00  4.51E+00  1.49E+00  7.16E+00  0.00E+00  1.39E+01  2.53E+01  1.93E+02
 
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
 #CPUT: Total CPU Time in Seconds,       13.172
Stop Time:
Sat Sep 25 12:43:04 CDT 2021
