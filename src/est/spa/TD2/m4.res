Sat Sep 18 14:25:45 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat4.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m4.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1637.47484137774        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.5233E+01 -3.1888E+01 -2.6071E+01  1.3695E+00  8.5146E+01 -2.4482E+01 -9.2757E-01 -1.2181E+00 -3.8464E+00 -1.9688E+01
             2.0440E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1642.89980083433        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  9.7827E-01  1.0032E+00  9.7573E-01  1.0095E+00  9.3629E-01  1.0653E+00  9.9531E-01  1.0098E+00  1.0187E+00  1.0423E+00
             9.3494E-01
 PARAMETER:  7.8032E-02  1.0316E-01  7.5431E-02  1.0948E-01  3.4172E-02  1.6324E-01  9.5298E-02  1.0975E-01  1.1852E-01  1.4144E-01
             3.2726E-02
 GRADIENT:   3.2467E+00  6.2979E-01 -5.5142E-02  4.2610E+00  6.1712E+00 -2.7252E+00  1.4450E+00  6.2202E-01 -8.3270E-01 -1.4771E+00
            -1.5540E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1643.26777400882        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      311
 NPARAMETR:  9.8093E-01  1.1166E+00  8.5231E-01  9.3566E-01  9.2335E-01  1.0703E+00  8.5688E-01  8.5676E-01  1.1055E+00  1.0300E+00
             9.4070E-01
 PARAMETER:  8.0750E-02  2.1027E-01 -5.9808E-02  3.3498E-02  2.0254E-02  1.6797E-01 -5.4452E-02 -5.4601E-02  2.0028E-01  1.2954E-01
             3.8873E-02
 GRADIENT:   5.8726E+00  6.2910E+00  2.6771E+00  1.0231E+01 -2.4361E+00 -1.5991E+00 -1.9214E+00 -5.5217E-01  1.2963E+00 -2.1857E+00
             5.0616E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1643.97592980830        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      489
 NPARAMETR:  9.7417E-01  1.2856E+00  6.2864E-01  8.0862E-01  8.8578E-01  1.0801E+00  9.0071E-01  5.2005E-01  1.1415E+00  9.5812E-01
             9.2992E-01
 PARAMETER:  7.3833E-02  3.5122E-01 -3.6420E-01 -1.1243E-01 -2.1292E-02  1.7708E-01 -4.5745E-03 -5.5384E-01  2.3233E-01  5.7215E-02
             2.7339E-02
 GRADIENT:  -1.0315E+01  5.3322E+00  4.5545E+00 -1.5712E+00 -1.1034E+01  9.9022E-01  1.2328E+00  3.1354E-01 -5.3130E-01  2.0183E-01
            -2.8656E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1644.45338828559        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  9.8061E-01  1.4721E+00  5.0676E-01  6.8343E-01  9.2420E-01  1.0779E+00  8.1576E-01  2.8797E-01  1.2740E+00  9.6003E-01
             9.3713E-01
 PARAMETER:  8.0424E-02  4.8672E-01 -5.7972E-01 -2.8063E-01  2.1179E-02  1.7497E-01 -1.0363E-01 -1.1449E+00  3.4215E-01  5.9214E-02
             3.5066E-02
 GRADIENT:   1.5806E+00  1.1083E-01 -1.3499E+00  1.7757E+00 -2.6114E-01 -1.4192E-01  2.6798E-01  2.1507E-01  2.3555E-01  9.4474E-01
             6.5184E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1644.49877489802        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      843
 NPARAMETR:  9.8049E-01  1.5366E+00  4.7659E-01  6.3902E-01  9.4610E-01  1.0779E+00  7.9179E-01  2.0619E-01  1.3331E+00  9.6758E-01
             9.3737E-01
 PARAMETER:  8.0294E-02  5.2954E-01 -6.4111E-01 -3.4782E-01  4.4594E-02  1.7501E-01 -1.3346E-01 -1.4790E+00  3.8752E-01  6.7047E-02
             3.5320E-02
 GRADIENT:   1.6322E+00 -2.1906E+00 -8.5875E-01 -6.8634E-01  6.7008E-01 -1.8260E-01  3.0846E-02 -3.1033E-02  2.1111E-01  6.5722E-01
             6.4998E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1644.52744946384        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1022
 NPARAMETR:  9.8012E-01  1.4896E+00  5.0103E-01  6.7031E-01  9.3289E-01  1.0781E+00  8.1058E-01  1.6129E-01  1.2933E+00  9.6599E-01
             9.3707E-01
 PARAMETER:  7.9919E-02  4.9853E-01 -5.9110E-01 -3.0002E-01  3.0529E-02  1.7517E-01 -1.1000E-01 -1.7245E+00  3.5723E-01  6.5396E-02
             3.5008E-02
 GRADIENT:   7.1257E-01 -1.6961E+00 -2.6121E-01 -4.6527E-01  2.0121E-01 -4.7539E-02  2.2160E-01  5.6445E-02  3.5429E-01  3.7617E-01
             4.0618E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1644.55813843745        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1198
 NPARAMETR:  9.7978E-01  1.5054E+00  4.8936E-01  6.6006E-01  9.3442E-01  1.0782E+00  8.0516E-01  1.5983E-02  1.3029E+00  9.6255E-01
             9.3630E-01
 PARAMETER:  7.9576E-02  5.0908E-01 -6.1466E-01 -3.1542E-01  3.2166E-02  1.7533E-01 -1.1671E-01 -4.0362E+00  3.6459E-01  6.1828E-02
             3.4178E-02
 GRADIENT:   6.1655E-02 -5.5539E-02  3.2426E-02 -6.8362E-02 -1.2341E-01 -1.0468E-02 -1.2571E-02  5.5838E-04  5.6252E-03  3.6264E-02
             2.4121E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1644.55831848998        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1325
 NPARAMETR:  9.7975E-01  1.5070E+00  4.8860E-01  6.5906E-01  9.3495E-01  1.0783E+00  8.0468E-01  1.0000E-02  1.3042E+00  9.6249E-01
             9.3625E-01
 PARAMETER:  7.9545E-02  5.1014E-01 -6.1620E-01 -3.1693E-01  3.2734E-02  1.7535E-01 -1.1731E-01 -4.9238E+00  3.6555E-01  6.1771E-02
             3.4130E-02
 GRADIENT:   2.5663E-03 -1.8608E-02 -2.3721E-03 -1.7460E-03  2.1725E-03 -3.6425E-03 -8.9033E-03  0.0000E+00 -2.1751E-03  3.7746E-03
             1.2803E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1325
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8425E-04 -2.9684E-02 -3.2459E-04  2.1243E-02 -3.1022E-02
 SE:             2.9884E-02  2.2574E-02  1.2349E-04  2.4486E-02  2.3042E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9508E-01  1.8851E-01  8.5762E-03  3.8565E-01  1.7819E-01

 ETASHRINKSD(%)  1.0000E-10  2.4376E+01  9.9586E+01  1.7967E+01  2.2806E+01
 ETASHRINKVR(%)  1.0000E-10  4.2810E+01  9.9998E+01  3.2707E+01  4.0411E+01
 EBVSHRINKSD(%)  3.3165E-01  2.4132E+01  9.9645E+01  1.8687E+01  2.1655E+01
 EBVSHRINKVR(%)  6.6220E-01  4.2440E+01  9.9999E+01  3.3883E+01  3.8620E+01
 RELATIVEINF(%)  9.9325E+01  3.4489E+00  1.2586E-04  4.3933E+00  9.2003E+00
 EPSSHRINKSD(%)  4.5803E+01
 EPSSHRINKVR(%)  7.0627E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1644.5583184899806     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -909.40749192624241     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.55
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1644.558       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.51E+00  4.89E-01  6.59E-01  9.35E-01  1.08E+00  8.05E-01  1.00E-02  1.30E+00  9.62E-01  9.36E-01
 


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
+        9.89E+02
 
 TH 2
+       -4.77E+00  4.57E+02
 
 TH 3
+        3.88E+00  2.45E+02  6.21E+02
 
 TH 4
+       -9.00E+00  3.32E+02 -3.79E+02  1.00E+03
 
 TH 5
+       -3.69E+00 -3.16E+02 -5.80E+02  3.22E+02  8.39E+02
 
 TH 6
+        8.75E-01 -1.02E+00  2.10E+00 -3.29E+00 -1.58E+00  1.70E+02
 
 TH 7
+        2.49E+00  6.00E+00 -2.36E+01 -6.85E+00 -1.59E+01  7.72E-01  1.05E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.23E+00 -2.25E+01 -3.64E+01  5.26E+01 -1.42E+00 -5.14E-01  2.34E+01  0.00E+00  5.65E+01
 
 TH10
+       -5.62E+00 -1.80E+01 -5.78E+01 -8.64E+00 -6.12E+01 -4.67E-01  2.51E+01  0.00E+00  5.97E+00  8.84E+01
 
 TH11
+       -5.80E+00 -1.69E+01 -3.22E+01 -1.94E+00 -2.51E+00  1.86E+00  1.12E+01  0.00E+00  5.32E+00  1.47E+01  2.33E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.631
Stop Time:
Sat Sep 18 14:26:09 CDT 2021
