Wed Sep 29 22:49:54 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat81.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1462.42430185077        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0904E+02  2.9856E+01  6.4438E+01  4.3886E+01  1.8675E+02  2.6052E+01 -3.1010E+01 -2.2230E+02 -2.7969E+01 -5.8852E+01
            -9.4603E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1759.71504863030        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0623E+00  1.0076E+00  8.9670E-01  1.0171E+00  8.4495E-01  1.2963E+00  1.0341E+00  1.3586E+00  9.1942E-01  7.9455E-01
             2.6486E+00
 PARAMETER:  1.6040E-01  1.0760E-01 -9.0318E-03  1.1691E-01 -6.8484E-02  3.5952E-01  1.3357E-01  4.0648E-01  1.5986E-02 -1.2998E-01
             1.0740E+00
 GRADIENT:   1.2900E+02 -8.1526E+00 -5.5206E-01 -1.2940E+01 -2.5978E+01  9.1715E+01  3.0682E+00  1.4455E+01  6.2608E+00  2.2529E+01
             1.5427E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1769.96362001853        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      202
 NPARAMETR:  1.0830E+00  1.2421E+00  5.9053E-01  8.9268E-01  7.7566E-01  1.1779E+00  1.0334E+00  1.4489E+00  1.0686E+00  4.7918E-01
             2.5552E+00
 PARAMETER:  1.7976E-01  3.1679E-01 -4.2674E-01 -1.3522E-02 -1.5404E-01  2.6371E-01  1.3288E-01  4.7078E-01  1.6631E-01 -6.3568E-01
             1.0381E+00
 GRADIENT:   5.2906E+01  3.3617E+01  1.7747E+00  2.8993E+01 -3.7517E+01  2.6751E+01  1.5992E+01  1.3345E+01  2.0214E+01  8.9140E+00
             1.2645E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1805.87758017232        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      381
 NPARAMETR:  1.0299E+00  1.2818E+00  1.4270E+00  8.7866E-01  1.1790E+00  1.0653E+00  9.8696E-01  2.8933E+00  1.0147E+00  8.2172E-01
             2.0466E+00
 PARAMETER:  1.2943E-01  3.4824E-01  4.5555E-01 -2.9354E-02  2.6466E-01  1.6321E-01  8.6872E-02  1.1624E+00  1.1461E-01 -9.6358E-02
             8.1619E-01
 GRADIENT:  -1.2813E+01 -2.4333E+01 -1.1340E+01  5.8971E+00  2.7674E+01 -8.0660E+00  5.2006E+00 -6.6034E+00  3.3618E+00 -1.3582E+00
             5.1732E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1811.56520958829        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      557
 NPARAMETR:  1.0399E+00  1.5467E+00  1.7763E+00  7.1889E-01  1.2746E+00  1.0871E+00  6.9122E-01  3.6277E+00  1.2962E+00  8.4108E-01
             2.0660E+00
 PARAMETER:  1.3913E-01  5.3613E-01  6.7453E-01 -2.3005E-01  3.4267E-01  1.8348E-01 -2.6930E-01  1.3886E+00  3.5947E-01 -7.3072E-02
             8.2564E-01
 GRADIENT:   1.8320E+00  4.7568E+00  1.8484E+00  9.8035E+00  3.8714E+00 -2.7866E-01 -2.6067E+00 -3.2950E+00 -1.0124E-01 -1.9743E+00
            -3.8481E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1813.05624783543        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  1.0398E+00  1.7299E+00  1.6561E+00  6.0097E-01  1.3034E+00  1.0896E+00  7.9612E-01  4.0503E+00  1.2549E+00  8.8888E-01
             2.0804E+00
 PARAMETER:  1.3901E-01  6.4805E-01  6.0448E-01 -4.0921E-01  3.6497E-01  1.8579E-01 -1.2800E-01  1.4988E+00  3.2705E-01 -1.7792E-02
             8.3254E-01
 GRADIENT:  -4.9720E-02  1.3880E+01  2.9025E+00  6.9928E+00 -9.0809E+00  2.7403E-01  1.6485E+00  4.5790E-01 -5.6263E-01 -9.8574E-01
            -1.1704E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1814.03882009027        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  1.0408E+00  2.0315E+00  7.2590E-01  3.8500E-01  1.3122E+00  1.0913E+00  7.1762E-01  3.5568E+00  1.7556E+00  9.0793E-01
             2.0803E+00
 PARAMETER:  1.3999E-01  8.0878E-01 -2.2034E-01 -8.5451E-01  3.7169E-01  1.8740E-01 -2.3182E-01  1.3689E+00  6.6282E-01  3.4135E-03
             8.3249E-01
 GRADIENT:   1.2574E-01  3.3691E+00  1.0361E+00  2.3732E+00  2.3604E+00  2.9971E-01 -1.2173E-01 -1.5560E+00  6.8727E-01  9.6432E-01
             6.4698E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1814.06195339741        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1087
 NPARAMETR:  1.0412E+00  2.0871E+00  5.9438E-01  3.4851E-01  1.3070E+00  1.0910E+00  7.1056E-01  3.4367E+00  1.8721E+00  8.9365E-01
             2.0822E+00
 PARAMETER:  1.4042E-01  8.3579E-01 -4.2024E-01 -9.5408E-01  3.6775E-01  1.8714E-01 -2.4170E-01  1.3345E+00  7.2705E-01 -1.2437E-02
             8.3342E-01
 GRADIENT:   3.6138E-01  8.2647E+00  6.0876E-01  3.6701E+00  1.2270E+00  6.0067E-02 -1.3238E-01 -1.4003E+00  7.6514E-01  7.7606E-01
             3.9839E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1814.11374220019        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1271            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0419E+00  2.0788E+00  5.9314E-01  3.4684E-01  1.3077E+00  1.0905E+00  7.0857E-01  3.4493E+00  1.8691E+00  8.7368E-01
             2.0898E+00
 PARAMETER:  1.4101E-01  8.3181E-01 -4.2232E-01 -9.5891E-01  3.6830E-01  1.8662E-01 -2.4451E-01  1.3382E+00  7.2543E-01 -3.5044E-02
             8.3708E-01
 GRADIENT:   1.5453E+02  3.1204E+02  1.3478E+00  2.3620E+01  8.7162E+00  2.9349E+01  4.0556E+00  2.2533E+01  5.8168E+00  1.2477E-01
             9.7379E+00

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1814.13595500611        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1369
 NPARAMETR:  1.0411E+00  2.0764E+00  5.9111E-01  3.4574E-01  1.3034E+00  1.0904E+00  7.0999E-01  3.4551E+00  1.8566E+00  8.7259E-01
             2.0915E+00
 PARAMETER:  1.4016E-01  8.3217E-01 -4.2655E-01 -9.6031E-01  3.6631E-01  1.8617E-01 -2.4398E-01  1.3423E+00  7.2111E-01 -3.5292E-02
             8.3634E-01
 GRADIENT:  -8.0844E-01  1.2192E+03 -2.3790E+03  1.0576E+03  1.6636E+00 -1.7777E-01 -3.4201E-01  7.3852E+02  2.4409E-01  5.9501E-02
            -1.2153E+03
 NUMSIGDIG:         2.5         2.3         2.3         2.3         1.9         2.3         1.8         2.3         2.1         1.6
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1369
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7096E-03 -3.3552E-02 -3.2057E-02  3.9429E-02 -4.0476E-02
 SE:             2.9643E-02  2.3177E-02  1.2554E-02  1.9675E-02  1.6567E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5401E-01  1.4771E-01  1.0662E-02  4.5064E-02  1.4556E-02

 ETASHRINKSD(%)  6.9179E-01  2.2355E+01  5.7944E+01  3.4087E+01  4.4500E+01
 ETASHRINKVR(%)  1.3788E+00  3.9713E+01  8.2313E+01  5.6555E+01  6.9197E+01
 EBVSHRINKSD(%)  1.1327E+00  2.0378E+01  6.1500E+01  4.1650E+01  4.3906E+01
 EBVSHRINKVR(%)  2.2525E+00  3.6603E+01  8.5177E+01  6.5953E+01  6.8534E+01
 RELATIVEINF(%)  9.7626E+01  7.8579E+00  4.3645E+00  4.0886E+00  1.3880E+01
 EPSSHRINKSD(%)  2.8865E+01
 EPSSHRINKVR(%)  4.9398E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1814.1359550061113     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -895.19742180143862     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.19
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1814.136       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.08E+00  5.91E-01  3.46E-01  1.31E+00  1.09E+00  7.09E-01  3.46E+00  1.86E+00  8.73E-01  2.09E+00
 


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
+        1.19E+06
 
 TH 2
+       -1.37E+01  2.14E+04
 
 TH 3
+        6.90E+05  2.36E+02  3.99E+05
 
 TH 4
+       -6.69E+00  2.78E+02  9.79E+02  2.30E+05
 
 TH 5
+        3.64E+05 -3.45E+01  2.10E+05  2.97E+02  1.11E+05
 
 TH 6
+        1.49E+00 -1.45E+00 -1.97E+01  1.03E+01 -7.39E-01  1.60E+02
 
 TH 7
+        4.10E+00 -1.25E+05  9.17E+01 -6.42E+01 -4.23E+00  2.25E-01  1.70E+02
 
 TH 8
+        1.01E+00  3.08E+03  6.76E+01 -2.80E+01 -1.64E+04  1.25E+00 -5.81E+00  1.13E+03
 
 TH 9
+        1.59E+02 -2.00E+01  8.63E+01 -3.51E+01 -3.96E+04  1.62E-01  1.59E+01 -3.17E+00  1.41E+04
 
 TH10
+       -2.44E-01  3.32E+00 -3.30E+01  3.01E+01 -2.73E+01 -3.25E-01 -5.72E+00  1.78E+00  6.36E-01  2.90E+01
 
 TH11
+       -1.23E+01 -2.00E+01 -2.06E+02 -4.37E+04 -6.06E+01 -9.77E-01  2.62E+01 -4.52E+03  1.60E+04  1.91E+01  8.44E+03
 
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
 #CPUT: Total CPU Time in Seconds,       34.705
Stop Time:
Wed Sep 29 22:50:30 CDT 2021
