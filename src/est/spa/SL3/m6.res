Sat Sep 18 12:38:57 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat6.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1642.60157911045        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6344E+01 -9.9297E+01 -9.4605E+01 -1.8632E+01  1.6158E+02  5.7387E+00  7.3053E-01  1.0181E+01 -6.9670E-01 -3.1965E+00
            -5.5905E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1659.00552892961        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8939E-01  1.0631E+00  1.1505E+00  1.0085E+00  9.6830E-01  9.8057E-01  9.4874E-01  9.4094E-01  1.0183E+00  9.1582E-01
             1.1839E+00
 PARAMETER:  8.9334E-02  1.6117E-01  2.4017E-01  1.0843E-01  6.7791E-02  8.0375E-02  4.7384E-02  3.9123E-02  1.1810E-01  1.2059E-02
             2.6878E-01
 GRADIENT:  -3.0934E+00  2.5236E+00 -1.0093E+01  1.4581E+01  2.5851E+01 -2.5789E+00  5.7355E+00  3.1348E+00 -1.1894E+00 -4.3303E+00
             2.0273E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1660.66863080144        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.9949E-01  1.2129E+00  1.1600E+00  9.3218E-01  1.0395E+00  9.8063E-01  6.2856E-01  7.6194E-01  1.1652E+00  1.0673E+00
             1.1466E+00
 PARAMETER:  9.9494E-02  2.9301E-01  2.4841E-01  2.9775E-02  1.3869E-01  8.0444E-02 -3.6433E-01 -1.7188E-01  2.5289E-01  1.6516E-01
             2.3677E-01
 GRADIENT:   2.1227E+01  4.0269E+01 -4.8907E+00  4.4853E+01  2.2808E+01 -2.5087E+00  5.1677E-01 -1.6016E+00 -6.7634E+00  2.5239E+00
             8.1639E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1661.94001283197        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.8789E-01  1.1202E+00  1.2734E+00  9.5808E-01  1.0185E+00  9.8424E-01  4.4293E-01  1.0494E+00  1.2019E+00  1.0562E+00
             1.1221E+00
 PARAMETER:  8.7818E-02  2.1347E-01  3.4166E-01  5.7180E-02  1.1832E-01  8.4119E-02 -7.1434E-01  1.4820E-01  2.8388E-01  1.5472E-01
             2.1520E-01
 GRADIENT:  -3.8549E-01  3.2152E+00 -2.8056E+00  4.6108E+00  1.2248E+00 -8.7761E-01  1.6688E+00  1.6909E+00  1.7489E+00  2.2843E+00
             1.1303E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1662.00355204505        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.8668E-01  1.0731E+00  1.3332E+00  9.8171E-01  1.0163E+00  9.8625E-01  3.3587E-01  1.0509E+00  1.1935E+00  1.0623E+00
             1.1194E+00
 PARAMETER:  8.6593E-02  1.7052E-01  3.8762E-01  8.1542E-02  1.1613E-01  8.6150E-02 -9.9102E-01  1.4967E-01  2.7693E-01  1.6042E-01
             2.1277E-01
 GRADIENT:  -1.8080E+00 -3.3126E+00  4.9531E-01 -3.3795E+00 -1.6439E+00  2.6966E-02  1.1529E+00  5.7486E-01  2.3424E+00  1.5123E+00
            -1.0844E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1662.10873967974        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.8675E-01  1.0297E+00  1.3586E+00  1.0064E+00  1.0096E+00  9.8673E-01  1.5141E-01  1.0290E+00  1.1810E+00  1.0597E+00
             1.1192E+00
 PARAMETER:  8.6660E-02  1.2926E-01  4.0648E-01  1.0641E-01  1.0951E-01  8.6646E-02 -1.7878E+00  1.2862E-01  2.6638E-01  1.5796E-01
             2.1257E-01
 GRADIENT:  -1.0438E+00 -5.6841E+00  6.1766E-01 -4.9956E+00  4.8575E-02  3.1848E-01  2.7386E-01 -2.7332E-01  9.2769E-01 -5.3291E-02
            -5.4012E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1662.33460455745        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      469
 NPARAMETR:  9.8775E-01  1.0575E+00  1.3622E+00  9.9111E-01  1.0207E+00  9.8503E-01  2.6781E-02  1.0699E+00  1.2163E+00  1.0701E+00
             1.1213E+00
 PARAMETER:  8.7673E-02  1.5588E-01  4.0909E-01  9.1075E-02  1.2053E-01  8.4916E-02 -3.5201E+00  1.6759E-01  2.9585E-01  1.6772E-01
             2.1445E-01
 GRADIENT:  -3.4475E+01 -7.6115E+00 -8.3742E-01 -4.1188E+00  1.9752E-01 -4.7017E+00  5.0569E-03  3.3317E-01  9.6147E-01  6.4808E-01
             5.5774E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1662.92033424886        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  1.0035E+00  1.2041E+00  1.2892E+00  8.9224E-01  1.0535E+00  9.9550E-01  1.0000E-02  1.0921E+00  1.3576E+00  1.0878E+00
             1.1208E+00
 PARAMETER:  1.0346E-01  2.8573E-01  3.5400E-01 -1.4022E-02  1.5215E-01  9.5487E-02 -5.3271E+00  1.8814E-01  4.0570E-01  1.8413E-01
             2.1405E-01
 GRADIENT:  -4.7602E-01 -1.3666E+00 -9.1456E-01  1.6486E-01  1.3631E+00 -2.1951E-01  0.0000E+00  2.0833E-01  4.4382E-01  6.1470E-02
             3.7188E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1662.92326580447        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      773
 NPARAMETR:  1.0037E+00  1.2108E+00  1.2903E+00  8.8746E-01  1.0551E+00  9.9606E-01  1.0000E-02  1.0917E+00  1.3632E+00  1.0892E+00
             1.1201E+00
 PARAMETER:  1.0369E-01  2.9126E-01  3.5485E-01 -1.9397E-02  1.5364E-01  9.6052E-02 -5.3731E+00  1.8770E-01  4.0983E-01  1.8541E-01
             2.1346E-01
 GRADIENT:  -1.8046E-02  1.8726E-03 -7.8917E-03  1.0415E-02  4.3185E-03 -2.7735E-03  0.0000E+00  2.2592E-03  4.9852E-05  2.2941E-03
             2.6845E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      773
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1171E-04 -9.9885E-04 -1.9358E-02 -2.2707E-03 -2.4395E-02
 SE:             2.9817E-02  2.9154E-04  1.0991E-02  2.8995E-02  2.4037E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7828E-01  6.1241E-04  7.8194E-02  9.3758E-01  3.1015E-01

 ETASHRINKSD(%)  1.0936E-01  9.9023E+01  6.3179E+01  2.8631E+00  1.9472E+01
 ETASHRINKVR(%)  2.1860E-01  9.9990E+01  8.6442E+01  5.6442E+00  3.5153E+01
 EBVSHRINKSD(%)  5.3999E-01  9.9187E+01  6.5970E+01  3.1663E+00  1.7178E+01
 EBVSHRINKVR(%)  1.0771E+00  9.9993E+01  8.8420E+01  6.2324E+00  3.1405E+01
 RELATIVEINF(%)  9.8804E+01  8.6217E-04  4.1940E+00  1.4654E+01  2.2985E+01
 EPSSHRINKSD(%)  4.2644E+01
 EPSSHRINKVR(%)  6.7103E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1662.9232658044668     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -927.77243924072866     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1662.923       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.21E+00  1.29E+00  8.87E-01  1.06E+00  9.96E-01  1.00E-02  1.09E+00  1.36E+00  1.09E+00  1.12E+00
 


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
+        1.10E+03
 
 TH 2
+       -1.64E+01  5.70E+02
 
 TH 3
+        3.56E+00  6.80E+01  7.27E+01
 
 TH 4
+       -1.20E+01  5.01E+02 -4.84E+00  7.14E+02
 
 TH 5
+       -2.20E+00 -2.25E+02 -1.44E+02 -3.31E+01  5.85E+02
 
 TH 6
+       -1.02E-01 -3.01E+00  6.32E-01 -1.48E-01 -7.83E-01  1.94E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.49E+00 -1.35E+01 -1.64E+01 -1.66E+00  5.49E+00 -6.26E-01  0.00E+00  1.17E+01
 
 TH 9
+        3.07E+00 -9.89E+01  1.81E+00  7.03E+00  3.99E+00  6.65E-02  0.00E+00 -8.45E-01  9.28E+01
 
 TH10
+       -4.39E-02 -9.75E+00 -7.09E+00 -6.12E+00 -5.37E+01  1.34E+00  0.00E+00  7.74E+00  1.97E+00  8.25E+01
 
 TH11
+       -9.73E+00 -2.55E+01 -8.14E+00 -1.22E+01  5.62E+00  3.58E+00  0.00E+00  3.67E+00  8.62E+00  1.44E+01  1.78E+02
 
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
 #CPUT: Total CPU Time in Seconds,       12.867
Stop Time:
Sat Sep 18 12:39:11 CDT 2021
