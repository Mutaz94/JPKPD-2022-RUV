Wed Sep 29 19:11:30 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat66.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1643.07450824807        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8073E+02 -4.9468E+01 -5.0964E+00 -4.0842E+01  1.4400E+00  4.3529E+01 -5.1339E+00  8.8562E+00  2.4302E+01  2.5684E+00
            -6.6798E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1653.05335759029        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.5414E-01  1.1240E+00  1.0534E+00  1.0059E+00  1.0757E+00  9.3594E-01  1.0506E+00  9.4292E-01  8.4765E-01  1.0014E+00
             1.0455E+00
 PARAMETER:  5.3058E-02  2.1689E-01  1.5198E-01  1.0592E-01  1.7299E-01  3.3792E-02  1.4935E-01  4.1228E-02 -6.5284E-02  1.0144E-01
             1.4449E-01
 GRADIENT:  -2.6983E-01  1.1752E+01  7.4056E+00  1.5249E+01 -2.3440E+00 -1.2283E+01 -1.0719E+01  1.2789E+00 -1.3916E+00 -7.3515E+00
             4.6687E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1654.77979010679        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.5653E-01  8.8490E-01  1.0391E+00  1.1630E+00  9.7019E-01  9.9119E-01  1.4752E+00  6.0550E-01  6.8239E-01  1.0331E+00
             1.0694E+00
 PARAMETER:  5.5557E-02 -2.2284E-02  1.3834E-01  2.5099E-01  6.9739E-02  9.1155E-02  4.8881E-01 -4.0170E-01 -2.8215E-01  1.3255E-01
             1.6710E-01
 GRADIENT:   7.7728E+00  2.5137E+01  2.4198E+00  4.3263E+01 -1.2009E+00  1.0777E+01 -7.9004E-01  2.2846E-01 -5.9758E+00  3.0632E+00
             1.5295E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1657.11162997387        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.5391E-01  8.2676E-01  8.2764E-01  1.1602E+00  8.2808E-01  9.6341E-01  1.5302E+00  2.8180E-01  7.1479E-01  8.6908E-01
             1.0271E+00
 PARAMETER:  5.2819E-02 -9.0240E-02 -8.9183E-02  2.4863E-01 -8.8646E-02  6.2723E-02  5.2541E-01 -1.1666E+00 -2.3577E-01 -4.0319E-02
             1.2672E-01
 GRADIENT:  -1.0486E+00  5.0462E+00 -5.0930E+00  9.8541E+00  5.1533E+00 -6.4723E-01  7.2605E-01  4.2650E-01  1.4774E+00  3.1182E-01
             1.5658E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1657.66786322231        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  9.5249E-01  6.2142E-01  8.3515E-01  1.2720E+00  7.5499E-01  9.6341E-01  1.9119E+00  1.1868E-01  6.6287E-01  8.5647E-01
             1.0244E+00
 PARAMETER:  5.1322E-02 -3.7574E-01 -8.0147E-02  3.4057E-01 -1.8106E-01  6.2725E-02  7.4808E-01 -2.0313E+00 -3.1117E-01 -5.4932E-02
             1.2406E-01
 GRADIENT:   1.0801E-01  2.5377E+00  3.6199E+00  2.2803E+00 -5.7225E+00  1.0702E-01 -7.1301E-02  4.4691E-02 -2.7741E-01 -2.4769E-02
            -4.0594E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1657.69326255247        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  9.5272E-01  6.0556E-01  8.3828E-01  1.2771E+00  7.5589E-01  9.6318E-01  1.9462E+00  5.1056E-02  6.6133E-01  8.6045E-01
             1.0252E+00
 PARAMETER:  5.1567E-02 -4.0161E-01 -7.6407E-02  3.4457E-01 -1.7986E-01  6.2480E-02  7.6588E-01 -2.8748E+00 -3.1350E-01 -5.0305E-02
             1.2489E-01
 GRADIENT:   1.2266E+00 -1.3388E+00 -6.7586E-01 -3.7190E+00  2.0484E+00  1.0053E-01  9.9972E-04  8.1209E-03  2.0936E-01 -3.0406E-01
            -5.7544E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1657.69693915039        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.5219E-01  6.0688E-01  8.3801E-01  1.2783E+00  7.5517E-01  9.6292E-01  1.9438E+00  3.1028E-02  6.6090E-01  8.6292E-01
             1.0253E+00
 PARAMETER:  5.1013E-02 -3.9942E-01 -7.6731E-02  3.4556E-01 -1.8082E-01  6.2217E-02  7.6463E-01 -3.3729E+00 -3.1415E-01 -4.7438E-02
             1.2502E-01
 GRADIENT:  -1.5955E-01  9.8063E-02 -1.1137E-01  1.8062E-01  3.2411E-01 -1.2701E-02 -6.4508E-02  3.1877E-03 -9.7502E-03  5.6873E-02
             3.0283E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1657.69928052777        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1237             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5277E-01  6.0542E-01  8.3763E-01  1.2779E+00  7.5431E-01  9.6316E-01  1.9515E+00  1.0000E-02  6.6093E-01  8.6229E-01
             1.0253E+00
 PARAMETER:  5.1615E-02 -4.0183E-01 -7.7175E-02  3.4523E-01 -1.8195E-01  6.2461E-02  7.6859E-01 -4.9822E+00 -3.1411E-01 -4.8160E-02
             1.2494E-01
 GRADIENT:   3.4181E+02  4.1738E+01  4.3557E+00  3.1117E+02  1.6189E+01  3.5127E+01  3.1302E+01  0.0000E+00  1.1714E+01  7.2567E-01
             9.9364E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1657.69929907652        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1335
 NPARAMETR:  9.5336E-01  6.0656E-01  8.3661E-01  1.2762E+00  7.5435E-01  9.6331E-01  1.9530E+00  1.0000E-02  6.6027E-01  8.6136E-01
             1.0253E+00
 PARAMETER:  5.1346E-02 -4.0159E-01 -7.7393E-02  3.4558E-01 -1.8183E-01  6.2405E-02  7.6787E-01 -4.9276E+00 -3.1439E-01 -4.8280E-02
             1.2495E-01
 GRADIENT:  -4.9443E-01 -1.1863E-01  2.6745E-01  1.2723E+00  2.7686E-02 -1.9183E-02 -5.2126E-02  0.0000E+00  2.6331E-02  3.5937E-02
            -4.9346E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1335
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.6213E-04  2.3264E-02 -4.8442E-04 -2.6328E-02 -1.8797E-04
 SE:             2.9827E-02  2.1478E-02  2.1742E-04  2.2930E-02  2.3231E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8229E-01  2.7874E-01  2.5878E-02  2.5089E-01  9.9354E-01

 ETASHRINKSD(%)  7.7008E-02  2.8045E+01  9.9272E+01  2.3182E+01  2.2173E+01
 ETASHRINKVR(%)  1.5396E-01  4.8224E+01  9.9995E+01  4.0990E+01  3.9429E+01
 EBVSHRINKSD(%)  4.8565E-01  2.8994E+01  9.9311E+01  2.1943E+01  1.9840E+01
 EBVSHRINKVR(%)  9.6894E-01  4.9581E+01  9.9995E+01  3.9072E+01  3.5744E+01
 RELATIVEINF(%)  9.8424E+01  6.8673E+00  4.4083E-04  9.0472E+00  5.3121E+00
 EPSSHRINKSD(%)  4.2560E+01
 EPSSHRINKVR(%)  6.7007E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1657.6992990765236     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -922.54847251278545     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1657.699       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.53E-01  6.06E-01  8.37E-01  1.28E+00  7.54E-01  9.63E-01  1.95E+00  1.00E-02  6.61E-01  8.62E-01  1.03E+00
 


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
+       -1.35E+01  4.25E+02
 
 TH 3
+        2.22E+01  1.86E+02  8.19E+02
 
 TH 4
+       -8.30E+00  3.78E+02 -3.26E+02  9.77E+02
 
 TH 5
+       -6.67E+00 -3.60E+02 -1.07E+03  3.55E+02  1.69E+03
 
 TH 6
+        6.04E-01 -2.74E+00  4.32E+00 -2.99E+00 -2.92E+00  2.11E+02
 
 TH 7
+        1.65E+00  3.63E+01 -2.27E+00 -1.36E+01  7.03E-01  7.20E-02  1.94E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.60E+00 -2.13E+01 -4.15E+01 -1.91E+01  4.22E+01 -6.68E-01  1.55E+01  0.00E+00  1.82E+02
 
 TH10
+       -1.37E+00 -4.62E+00 -6.52E+01 -3.07E+01 -4.92E+01  2.06E-01  3.43E+00  0.00E+00  1.49E+01  1.08E+02
 
 TH11
+       -9.00E+00 -1.06E+01 -4.52E+01 -6.81E+00  1.83E+01  2.09E+00  2.02E+00  0.00E+00  1.63E+01  2.86E+01  2.10E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.129
Stop Time:
Wed Sep 29 19:11:55 CDT 2021
