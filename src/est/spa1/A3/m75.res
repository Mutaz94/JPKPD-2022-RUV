Thu Sep 30 00:32:13 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat75.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -237.645440665660        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1999E+02  2.2106E+01  1.8817E+02 -4.0875E+01  2.4105E+02  3.9409E+01 -8.6613E+01 -1.6342E+02 -8.6417E+01 -1.1862E+02
            -3.2555E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1460.29986715512        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1250E+00  9.5477E-01  9.3305E-01  1.2434E+00  9.3072E-01  8.8763E-01  1.0668E+00  8.3939E-01  1.3081E+00  7.4300E-01
             5.2761E+00
 PARAMETER:  2.1781E-01  5.3718E-02  3.0700E-02  3.1785E-01  2.8201E-02 -1.9200E-02  1.6464E-01 -7.5085E-02  3.6861E-01 -1.9706E-01
             1.7632E+00
 GRADIENT:   1.7027E+02 -1.3447E+01 -2.7349E+01  3.1849E+01  1.9566E+01 -2.2936E+01  1.2156E+01  7.5525E+00  4.4692E+01  1.4142E+01
             3.6398E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1501.96773052971        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0947E+00  7.0010E-01  3.1792E-01  1.2890E+00  3.5797E-01  1.0216E+00  1.3751E+00  1.2810E-01  1.1642E+00  1.4485E-01
             4.6325E+00
 PARAMETER:  1.9049E-01 -2.5653E-01 -1.0460E+00  3.5386E-01 -9.2731E-01  1.2136E-01  4.1851E-01 -1.9550E+00  2.5204E-01 -1.8320E+00
             1.6331E+00
 GRADIENT:   8.3331E+01  1.1333E+02  4.1737E+01  1.1827E+02 -1.2889E+02  1.3129E+01  9.4833E+00  3.2653E-01  6.2009E+00  1.0231E+00
             3.0775E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1602.17710827337        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      247
 NPARAMETR:  9.5789E-01  5.0825E-01  2.7106E-01  1.2240E+00  3.1455E-01  9.6849E-01  1.4411E+00  1.0000E-02  1.2176E+00  2.5740E-01
             2.8055E+00
 PARAMETER:  5.6983E-02 -5.7677E-01 -1.2054E+00  3.0209E-01 -1.0566E+00  6.7978E-02  4.6543E-01 -6.4440E+00  2.9688E-01 -1.2571E+00
             1.1316E+00
 GRADIENT:  -1.4416E+02  4.7777E+01 -8.7441E+00  6.8540E+01 -3.9288E+01 -2.1631E+01 -2.6049E+00  0.0000E+00 -2.8886E+00 -2.8941E+00
            -1.9055E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1616.49801610332        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      422
 NPARAMETR:  1.0105E+00  3.4162E-01  2.9204E-01  1.2537E+00  2.9724E-01  9.6672E-01  1.8282E+00  1.0000E-02  1.1685E+00  4.5110E-01
             2.7629E+00
 PARAMETER:  1.1048E-01 -9.7407E-01 -1.1309E+00  3.2613E-01 -1.1132E+00  6.6152E-02  7.0333E-01 -8.3479E+00  2.5570E-01 -6.9607E-01
             1.1163E+00
 GRADIENT:  -1.5118E+01  1.3231E+01 -2.1611E+00  3.1306E+01 -3.6475E+00 -1.0011E+01  1.1740E+00  0.0000E+00 -1.0507E+00  6.9685E-01
             6.1423E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1620.48954489653        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      597
 NPARAMETR:  1.0136E+00  2.0069E-01  2.4698E-01  1.1942E+00  2.4782E-01  9.8887E-01  2.1429E+00  1.0000E-02  1.1966E+00  5.3475E-01
             2.7221E+00
 PARAMETER:  1.1351E-01 -1.5060E+00 -1.2985E+00  2.7745E-01 -1.2951E+00  8.8811E-02  8.6215E-01 -1.0910E+01  2.7947E-01 -5.2595E-01
             1.1014E+00
 GRADIENT:   2.6812E+00  2.5398E+00  2.2896E+00 -1.2366E+01 -8.1502E+00  9.5919E-01 -3.2097E+00  0.0000E+00 -4.6556E+00 -9.5807E-01
            -4.9885E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1623.57492526943        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      773
 NPARAMETR:  1.0042E+00  1.0829E-01  2.7987E-01  1.2786E+00  2.6374E-01  9.7772E-01  4.0204E+00  1.0000E-02  1.1606E+00  5.0372E-01
             2.7402E+00
 PARAMETER:  1.0420E-01 -2.1229E+00 -1.1734E+00  3.4580E-01 -1.2328E+00  7.7471E-02  1.4914E+00 -1.3677E+01  2.4892E-01 -5.8573E-01
             1.1080E+00
 GRADIENT:   2.1198E+00  8.5981E+00 -1.2318E+01 -2.2411E+00  1.3186E+01 -1.4139E+00  1.0317E+01  0.0000E+00  4.9854E+00 -2.9238E+00
            -3.3621E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1624.06198017298        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      948
 NPARAMETR:  9.9902E-01  7.6879E-02  2.8377E-01  1.2912E+00  2.6637E-01  9.8183E-01  5.0258E+00  1.0000E-02  1.1084E+00  4.6639E-01
             2.7443E+00
 PARAMETER:  9.9015E-02 -2.4655E+00 -1.1596E+00  3.5557E-01 -1.2229E+00  8.1665E-02  1.7146E+00 -1.5321E+01  2.0294E-01 -6.6274E-01
             1.1095E+00
 GRADIENT:  -6.3658E+00 -5.0563E+00  9.8267E+00  9.6601E+00 -1.4168E+01  1.7087E+00 -8.4823E+00  0.0000E+00 -4.4673E+00  2.5010E+00
             2.6670E+00

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1624.19945713422        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1075
 NPARAMETR:  1.0007E+00  6.9505E-02  2.7873E-01  1.2828E+00  2.6319E-01  9.7797E-01  5.3268E+00  1.0000E-02  1.1281E+00  4.6201E-01
             2.7457E+00
 PARAMETER:  1.0068E-01 -2.5664E+00 -1.1775E+00  3.4907E-01 -1.2349E+00  7.7727E-02  1.7728E+00 -1.5877E+01  2.2052E-01 -6.7218E-01
             1.1100E+00
 GRADIENT:  -3.5191E-01 -2.7776E+00  4.1688E+00  1.7384E+00 -4.6337E+00  1.5784E-01 -4.2001E+00  0.0000E+00 -3.1875E-02  1.1846E+00
             2.1161E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1075
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0345E-03  3.1570E-02  4.1777E-05 -1.6256E-02  9.0959E-03
 SE:             2.9013E-02  1.2918E-02  2.4377E-04  2.7350E-02  1.7007E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4409E-01  1.4533E-02  8.6393E-01  5.5226E-01  5.9278E-01

 ETASHRINKSD(%)  2.8024E+00  5.6722E+01  9.9183E+01  8.3734E+00  4.3023E+01
 ETASHRINKVR(%)  5.5263E+00  8.1270E+01  9.9993E+01  1.6046E+01  6.7536E+01
 EBVSHRINKSD(%)  2.4566E+00  6.7713E+01  9.9101E+01  7.0168E+00  4.0745E+01
 EBVSHRINKVR(%)  4.8528E+00  8.9576E+01  9.9992E+01  1.3541E+01  6.4888E+01
 RELATIVEINF(%)  9.4100E+01  7.5161E+00  4.4809E-04  3.8404E+01  1.8956E+00
 EPSSHRINKSD(%)  2.5131E+01
 EPSSHRINKVR(%)  4.3946E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1624.1994571342229     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -705.26092392955024     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.74
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1624.199       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  6.95E-02  2.79E-01  1.28E+00  2.63E-01  9.78E-01  5.33E+00  1.00E-02  1.13E+00  4.62E-01  2.75E+00
 


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
+        1.24E+03
 
 TH 2
+       -8.29E+01  1.38E+04
 
 TH 3
+       -4.59E+01  7.72E+02  2.73E+04
 
 TH 4
+        3.10E+01  5.47E+01 -3.93E+02  9.02E+03
 
 TH 5
+        1.61E+02  2.74E+04 -1.60E+04 -3.83E+02  5.55E+04
 
 TH 6
+       -4.13E+01 -1.12E+01  1.35E+01 -2.30E+01 -5.92E-01  2.03E+02
 
 TH 7
+        3.72E-01 -7.06E+02 -1.75E+01 -1.37E+00  5.81E+02  1.16E-02  2.11E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.79E+01  3.52E+00  8.63E+01 -1.93E+01  1.31E+02 -3.39E+00  9.65E-01  0.00E+00  1.16E+02
 
 TH10
+        6.30E+01 -2.40E+02 -1.43E+02  1.24E+04  2.88E+02 -2.47E+01 -3.23E+00  0.00E+00 -5.24E+00  1.80E+04
 
 TH11
+       -1.32E+01  3.61E+03 -4.74E+01 -1.16E+01  2.79E+01  1.27E+00  6.71E+01  0.00E+00  3.84E+00  2.94E+01  2.51E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.028
Stop Time:
Thu Sep 30 00:32:40 CDT 2021
