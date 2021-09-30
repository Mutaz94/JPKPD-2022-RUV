Tue Sep 28 20:39:51 CDT 2021
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
$DATA ../../../../data/int/B/dat49.csv ignore=@
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
Current Date:       28 SEP 2021
Days until program expires : 201
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
 NO. OF DATA RECS IN DATA SET:     1000
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

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3205.14342768162        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0571E+02  2.0140E+01  1.0032E+02  1.9289E+01  8.0254E+01  2.1349E+01 -3.8058E+01 -2.5090E+02 -5.4743E+01  4.6192E+00
            -8.8513E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3493.58539044079        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:       97
 NPARAMETR:  7.1845E-01  9.8701E-01  9.3697E-01  9.8725E-01  9.4698E-01  9.9748E-01  1.0316E+00  1.1786E+00  1.0361E+00  9.9472E-01
             1.6749E+00
 PARAMETER: -2.3066E-01  8.6928E-02  3.4898E-02  8.7164E-02  4.5525E-02  9.7478E-02  1.3116E-01  2.6436E-01  1.3550E-01  9.4710E-02
             6.1577E-01
 GRADIENT:  -2.7577E+02 -2.6353E+01 -1.8715E+01 -1.7282E+01 -2.7730E+00 -1.3673E+02 -3.4135E+01  3.0706E+01  8.1507E+00  1.2064E+01
             7.1111E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3521.90150936542        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      281
 NPARAMETR:  7.3854E-01  8.5327E-01  8.6175E-01  1.0308E+00  9.4925E-01  9.9350E-01  1.8104E+00  4.7811E-01  8.4918E-01  7.1332E-01
             1.6544E+00
 PARAMETER: -2.0308E-01 -5.8675E-02 -4.8793E-02  1.3034E-01  4.7914E-02  9.3482E-02  6.9352E-01 -6.3791E-01 -6.3488E-02 -2.3782E-01
             6.0342E-01
 GRADIENT:  -5.7555E+02 -8.3136E+01 -7.2240E+01  3.6259E+01  1.5987E+02 -1.6208E+02  3.3637E+01  1.4396E+00 -2.6832E+01 -8.3050E+00
             6.6953E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3578.17865385428        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      458
 NPARAMETR:  8.0504E-01  6.5933E-01  6.2443E-01  1.1252E+00  6.2299E-01  9.5699E-01  1.6340E+00  2.2241E-01  9.8853E-01  1.6976E-01
             1.5004E+00
 PARAMETER: -1.1687E-01 -3.1652E-01 -3.7092E-01  2.1798E-01 -3.7323E-01  5.6039E-02  5.9103E-01 -1.4032E+00  8.8459E-02 -1.6734E+00
             5.0573E-01
 GRADIENT:  -4.0955E+02  1.4534E+00  2.1722E+02 -1.9197E+01 -3.9381E+02 -9.8976E+01 -6.2570E+01 -1.5978E+01  1.3361E+01 -2.1264E+01
             4.5880E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3785.46589765580        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      645             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5603E-01  9.0607E-01  8.4215E-01  1.0403E+00  8.6418E-01  1.0896E+00  1.2517E+00  4.6708E-01  1.0118E+00  8.6194E-01
             1.0730E+00
 PARAMETER:  5.5030E-02  1.3586E-03 -7.1803E-02  1.3955E-01 -4.5979E-02  1.8580E-01  3.2454E-01 -6.6126E-01  1.1168E-01 -4.8574E-02
             1.7045E-01
 GRADIENT:   3.8042E+02  4.8411E+01  2.6500E+00  1.1102E+02  5.5130E+01  1.3227E+02  1.8076E+01 -7.9336E+00  1.2968E+01  5.4981E-01
             6.5651E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3787.54858849824        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      780
 NPARAMETR:  9.5306E-01  9.1121E-01  8.4329E-01  1.0484E+00  8.6477E-01  1.0307E+00  1.3123E+00  4.6702E-01  1.0128E+00  8.6600E-01
             1.0433E+00
 PARAMETER:  5.1921E-02  7.0225E-03 -7.0439E-02  1.4726E-01 -4.5290E-02  1.3025E-01  3.7179E-01 -6.6138E-01  1.1275E-01 -4.3873E-02
             1.4239E-01
 GRADIENT:   1.9422E-01 -3.8943E+00 -3.5812E+00 -7.5535E-01  5.9348E+00  1.8339E+00 -1.6809E+00 -1.0577E+01  1.0123E-01 -3.6506E-01
             3.4472E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3788.22304275909        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  9.5301E-01  9.1398E-01  8.4470E-01  1.0487E+00  8.6300E-01  1.0255E+00  1.3242E+00  4.9724E-01  1.0112E+00  8.6554E-01
             1.0415E+00
 PARAMETER:  5.1871E-02  1.0059E-02 -6.8771E-02  1.4757E-01 -4.7341E-02  1.2515E-01  3.8081E-01 -5.9869E-01  1.1119E-01 -4.4398E-02
             1.4068E-01
 GRADIENT:   1.0838E-01  1.5389E-01 -3.2892E+00 -2.0951E-01 -1.5493E+00 -2.5449E-01  2.7326E-01 -9.8947E+00 -2.3561E-02  1.0316E+00
             3.0181E+00

0ITERATION NO.:   31    OBJECTIVE VALUE:  -3788.22304275909        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:      974
 NPARAMETR:  9.5308E-01  9.1397E-01  8.4495E-01  1.0487E+00  8.6307E-01  1.0257E+00  1.3242E+00  4.9713E-01  1.0113E+00  8.6468E-01
             1.0416E+00
 PARAMETER:  5.1871E-02  1.0059E-02 -6.8771E-02  1.4757E-01 -4.7341E-02  1.2515E-01  3.8081E-01 -5.9869E-01  1.1119E-01 -4.4398E-02
             1.4068E-01
 GRADIENT:  -9.6088E-01  1.3680E-01 -2.9185E+00  1.2514E-01 -1.0374E+00 -6.7377E-01  8.2466E-02  6.1543E+04 -6.5935E-02  9.8987E-01
            -2.6301E+05
 NUMSIGDIG:         2.0         2.8         1.4         3.3         1.9         1.6         3.0         2.3         2.7         0.9
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      974
 NO. OF SIG. DIGITS IN FINAL EST.:  0.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0662E-03 -1.9324E-02 -1.1228E-02  1.0986E-02 -2.3567E-02
 SE:             2.9942E-02  2.5290E-02  1.4212E-02  2.8066E-02  2.3829E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7159E-01  4.4480E-01  4.2948E-01  6.9547E-01  3.2266E-01

 ETASHRINKSD(%)  1.0000E-10  1.5276E+01  5.2389E+01  5.9745E+00  2.0171E+01
 ETASHRINKVR(%)  1.0000E-10  2.8218E+01  7.7332E+01  1.1592E+01  3.6273E+01
 EBVSHRINKSD(%)  2.7978E-01  1.4925E+01  5.8796E+01  6.9214E+00  2.1164E+01
 EBVSHRINKVR(%)  5.5877E-01  2.7623E+01  8.3022E+01  1.3364E+01  3.7849E+01
 RELATIVEINF(%)  9.9440E+01  3.6923E+01  7.6802E+00  6.2929E+01  1.9945E+01
 EPSSHRINKSD(%)  2.1099E+01
 EPSSHRINKVR(%)  3.7746E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3788.2230427590866     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2134.1336829906759     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.29
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3788.223       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.53E-01  9.14E-01  8.45E-01  1.05E+00  8.63E-01  1.03E+00  1.32E+00  4.97E-01  1.01E+00  8.66E-01  1.04E+00
 


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
+        1.02E+08
 
 TH 2
+       -7.57E-01  5.47E+02
 
 TH 3
+       -1.05E+00 -4.80E+01  1.04E+03
 
 TH 4
+       -6.26E+07  6.53E+07 -7.06E+07  3.85E+07
 
 TH 5
+       -9.02E-01 -3.65E+02 -6.74E+02  2.46E+02  1.22E+03
 
 TH 6
+        2.65E+00 -7.31E-01  1.28E+00 -7.57E-01 -1.31E+00  1.88E+02
 
 TH 7
+        5.54E-01  8.75E+00  1.50E-01  1.18E+07 -5.36E+00 -3.85E-01  6.08E+01
 
 TH 8
+        1.40E+03  6.45E+02  1.66E+04 -7.70E+03  3.11E+03  1.33E+03 -5.72E+02  1.04E+07
 
 TH 9
+        4.71E-01 -6.57E+00  2.61E+01  5.30E+07  2.33E+00 -5.30E-02  1.04E+01 -2.54E+04  1.45E+02
 
 TH10
+       -6.17E-01 -1.12E+01 -2.83E+01  6.89E+07 -2.17E+01  7.80E-01  2.75E+01 -1.96E+03 -6.04E+00  1.00E+02
 
 TH11
+        6.62E+07 -6.91E+07 -3.41E+04 -4.08E+07 -6.37E+03 -2.71E+03 -1.25E+07 -6.68E+04 -5.61E+07  4.04E+03  4.32E+07
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       38.161
Stop Time:
Tue Sep 28 20:40:30 CDT 2021
