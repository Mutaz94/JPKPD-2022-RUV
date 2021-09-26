Sat Sep 25 01:52:35 CDT 2021
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
$DATA ../../../../data/int/SL3/dat7.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      984
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      884
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -121.799719050927        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.4819E+01  6.4361E+01  1.9034E+02  9.8670E+01  1.7858E+02 -3.9982E+01 -1.8113E+02 -6.4095E+02 -1.8670E+02 -6.9704E+01
            -6.2344E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2679.31984736857        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9597E-01  1.2797E+00  9.4632E-01  8.3467E-01  1.1532E+00  1.1893E+00  1.0415E+00  1.0011E+00  7.7571E-01  1.0618E+00
             2.6892E+00
 PARAMETER:  9.5960E-02  3.4661E-01  4.4824E-02 -8.0717E-02  2.4252E-01  2.7338E-01  1.4062E-01  1.0111E-01 -1.5397E-01  1.5992E-01
             1.0892E+00
 GRADIENT:  -7.5001E-01  7.7297E+00 -2.2379E+01  3.0975E+01  4.5894E+01  3.0109E+01  1.1548E+01  1.0020E+00 -1.0468E+01 -2.1901E+01
             5.1613E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2683.20969937907        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0078E+00  1.4726E+00  1.1004E+00  7.2185E-01  1.3335E+00  1.0997E+00  9.5356E-01  8.1593E-01  7.3054E-01  1.3464E+00
             2.6962E+00
 PARAMETER:  1.0777E-01  4.8705E-01  1.9563E-01 -2.2594E-01  3.8781E-01  1.9501E-01  5.2444E-02 -1.0342E-01 -2.1397E-01  3.9745E-01
             1.0918E+00
 GRADIENT:   1.9112E+01  2.8808E+01  6.6855E-01  2.2572E+01  1.6663E+01  2.0459E+00  1.0310E+01 -1.1099E+00 -5.9177E+00 -3.7481E+00
             8.0030E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2686.77581937714        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9927E-01  1.4885E+00  9.7363E-01  6.9728E-01  1.2841E+00  1.0736E+00  7.4963E-01  8.8568E-01  1.0343E+00  1.3585E+00
             2.6626E+00
 PARAMETER:  9.9274E-02  4.9779E-01  7.3278E-02 -2.6058E-01  3.5004E-01  1.7098E-01 -1.8817E-01 -2.1398E-02  1.3369E-01  4.0640E-01
             1.0793E+00
 GRADIENT:   4.5584E+00  9.9119E+00  8.0630E-01  5.0724E+00 -9.1721E-01 -7.1146E+00 -9.4436E-01  2.8720E-02 -4.5368E+00  3.2367E+00
            -4.9192E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2687.12866997849        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      325
 NPARAMETR:  9.9734E-01  1.5037E+00  9.4129E-01  6.8466E-01  1.2846E+00  1.0894E+00  7.1786E-01  6.6866E-01  1.1085E+00  1.3495E+00
             2.6636E+00
 PARAMETER:  9.7335E-02  5.0792E-01  3.9499E-02 -2.7884E-01  3.5045E-01  1.8566E-01 -2.3148E-01 -3.0248E-01  2.0299E-01  3.9975E-01
             1.0797E+00
 GRADIENT:  -5.7642E+00 -3.0313E+00  6.8230E-01  1.4332E+00 -3.8340E+00 -3.3403E+00 -7.1010E-01 -1.8630E-02 -1.3728E+00 -4.0936E-01
            -6.0997E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2688.27722016583        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      503
 NPARAMETR:  9.9884E-01  1.8403E+00  7.4871E-01  4.6750E-01  1.5127E+00  1.0972E+00  5.8230E-01  1.1292E-01  1.5548E+00  1.5188E+00
             2.6608E+00
 PARAMETER:  9.8844E-02  7.0991E-01 -1.8940E-01 -6.6036E-01  5.1389E-01  1.9278E-01 -4.4078E-01 -2.0810E+00  5.4134E-01  5.1789E-01
             1.0786E+00
 GRADIENT:  -3.0745E+00  1.0279E+01 -2.6915E-01  4.4361E+00  2.1336E+00 -1.0465E+00 -6.9220E-01  2.9751E-03 -1.0567E+00 -1.2399E-01
            -4.8517E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2688.32341884915        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      678
 NPARAMETR:  1.0005E+00  1.8605E+00  7.3054E-01  4.4869E-01  1.5230E+00  1.1000E+00  5.7503E-01  8.1525E-02  1.6154E+00  1.5287E+00
             2.6651E+00
 PARAMETER:  1.0046E-01  7.2085E-01 -2.1397E-01 -7.0143E-01  5.2068E-01  1.9528E-01 -4.5333E-01 -2.4068E+00  5.7960E-01  5.2443E-01
             1.0802E+00
 GRADIENT:  -5.8893E-03 -4.7410E-03  1.6034E-04 -1.4912E-02  8.4333E-03 -1.4162E-03 -8.3342E-03  2.2268E-03 -3.3473E-03  1.8403E-03
             3.5703E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -2688.32341884915        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  1.0005E+00  1.8605E+00  7.3054E-01  4.4869E-01  1.5230E+00  1.1000E+00  5.7503E-01  8.1525E-02  1.6154E+00  1.5287E+00
             2.6651E+00
 PARAMETER:  1.0046E-01  7.2085E-01 -2.1397E-01 -7.0143E-01  5.2068E-01  1.9528E-01 -4.5333E-01 -2.4068E+00  5.7960E-01  5.2443E-01
             1.0802E+00
 GRADIENT:  -5.8893E-03 -4.7410E-03  1.6034E-04 -1.4912E-02  8.4333E-03 -1.4162E-03 -8.3342E-03  2.2268E-03 -3.3473E-03  1.8403E-03
             3.5703E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      700
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6534E-03 -3.4501E-02 -4.1617E-04  2.4100E-02 -1.9327E-02
 SE:             2.9480E-02  2.0016E-02  4.5244E-04  2.1854E-02  2.6483E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5527E-01  8.4760E-02  3.5766E-01  2.7014E-01  4.6551E-01

 ETASHRINKSD(%)  1.2380E+00  3.2945E+01  9.8484E+01  2.6785E+01  1.1280E+01
 ETASHRINKVR(%)  2.4606E+00  5.5037E+01  9.9977E+01  4.6396E+01  2.1287E+01
 EBVSHRINKSD(%)  1.3523E+00  3.2333E+01  9.8612E+01  2.8718E+01  1.0037E+01
 EBVSHRINKVR(%)  2.6862E+00  5.4212E+01  9.9981E+01  4.9189E+01  1.9066E+01
 RELATIVEINF(%)  9.7251E+01  4.7322E+00  1.0876E-02  5.6650E+00  2.7703E+01
 EPSSHRINKSD(%)  1.6199E+01
 EPSSHRINKVR(%)  2.9774E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          884
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1624.6833267058612     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2688.3234188491538     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1063.6400921432926     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2688.323       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.86E+00  7.31E-01  4.49E-01  1.52E+00  1.10E+00  5.75E-01  8.15E-02  1.62E+00  1.53E+00  2.67E+00
 


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
+        8.90E+02
 
 TH 2
+       -1.34E+01  4.63E+02
 
 TH 3
+        6.75E-01  4.26E+01  6.65E+01
 
 TH 4
+       -2.22E+01  5.30E+02 -6.68E+01  1.09E+03
 
 TH 5
+       -1.62E+00 -9.08E+01 -4.70E+01  7.78E+01  1.65E+02
 
 TH 6
+        3.62E+00 -4.70E+00  1.24E+00 -7.48E+00 -6.27E-01  1.55E+02
 
 TH 7
+        3.00E+00 -4.03E+01  5.58E+00 -1.80E+01 -8.02E+00  1.42E+00  1.19E+02
 
 TH 8
+        1.30E+00 -8.92E-02 -1.20E-01  2.64E-01  1.17E-01  1.19E+00 -3.84E-01 -3.46E-01
 
 TH 9
+        1.08E+00 -9.38E+00 -1.84E+00  4.62E+01  1.31E+00  4.30E-01  3.01E+01 -3.20E-02  2.45E+01
 
 TH10
+        4.29E-01 -8.21E+00 -4.56E+00  1.22E+01 -8.25E+00  1.77E-01  1.73E+01  1.97E-01  2.67E-01  5.16E+01
 
 TH11
+       -1.26E+01 -1.99E+01 -4.27E+00 -1.77E+01  1.64E+00  1.37E+00  6.11E+00  1.31E-01  3.10E+00  4.64E+00  1.64E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.888
Stop Time:
Sat Sep 25 01:53:03 CDT 2021
