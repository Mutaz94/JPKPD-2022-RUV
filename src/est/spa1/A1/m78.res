Wed Sep 29 22:48:32 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1858.25553621993        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6697E+02 -4.0899E+01 -2.4391E+01 -1.1896E+01  9.3700E+01  3.1761E+01 -3.7124E+01  5.6791E+00 -6.3027E+00 -2.2002E+01
            -4.2101E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1918.18841285812        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.5266E-01  1.0351E+00  1.1152E+00  1.0456E+00  1.0037E+00  9.7110E-01  1.2848E+00  8.3133E-01  1.0216E+00  9.4099E-01
             1.8719E+00
 PARAMETER:  5.1504E-02  1.3447E-01  2.0907E-01  1.4456E-01  1.0365E-01  7.0675E-02  3.5059E-01 -8.4734E-02  1.2133E-01  3.9179E-02
             7.2694E-01
 GRADIENT:   5.0023E+01  3.7623E+00 -1.2376E+01  1.5898E+01  4.1590E+00  3.4420E+00  8.9349E+00  8.2230E+00  1.4811E+01  4.7181E+00
             1.7871E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1927.28498375359        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      220
 NPARAMETR:  9.4697E-01  7.0179E-01  1.1535E+00  1.2871E+00  8.3705E-01  9.5773E-01  1.5797E+00  4.9099E-01  9.5402E-01  1.0416E+00
             1.7847E+00
 PARAMETER:  4.5517E-02 -2.5412E-01  2.4278E-01  3.5236E-01 -7.7868E-02  5.6806E-02  5.5723E-01 -6.1133E-01  5.2929E-02  1.4074E-01
             6.7927E-01
 GRADIENT:  -7.4160E+01  2.9098E+01  1.8671E+01  2.5470E+01 -6.4285E+01 -1.1900E+01  7.0313E-01  3.8908E+00  1.3159E+01  2.0269E+01
             1.5865E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1954.61654073868        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      397
 NPARAMETR:  9.7231E-01  5.0990E-01  6.9133E-01  1.3153E+00  5.8245E-01  9.7384E-01  2.2878E+00  1.9682E-01  8.1647E-01  6.9955E-01
             1.3609E+00
 PARAMETER:  7.1924E-02 -5.7354E-01 -2.6913E-01  3.7406E-01 -4.4051E-01  7.3492E-02  9.2759E-01 -1.5255E+00 -1.0277E-01 -2.5732E-01
             4.0816E-01
 GRADIENT:   3.4090E-02  1.9324E+01  1.4241E+01  2.2926E+01 -2.9149E+01 -4.3599E+00  3.7784E+00  9.5557E-01 -5.9789E+00  7.0582E-01
            -2.3989E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1957.46780872378        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      572
 NPARAMETR:  9.6497E-01  2.6519E-01  9.6075E-01  1.4878E+00  6.6234E-01  9.7277E-01  2.8820E+00  4.9534E-02  8.7010E-01  9.1793E-01
             1.3998E+00
 PARAMETER:  6.4344E-02 -1.2273E+00  5.9955E-02  4.9727E-01 -3.1198E-01  7.2397E-02  1.1585E+00 -2.9051E+00 -3.9146E-02  1.4367E-02
             4.3634E-01
 GRADIENT:  -3.1102E+00  5.5198E+00  7.6198E+00  4.2556E+00 -1.9057E+01 -2.0473E+00  2.7758E+00  4.4801E-02  8.0372E+00  6.3805E+00
            -3.6063E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1958.43279707974        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      747
 NPARAMETR:  9.6450E-01  1.4231E-01  1.1667E+00  1.5853E+00  7.3490E-01  9.8140E-01  3.5940E+00  1.0339E-02  8.1927E-01  9.8404E-01
             1.4163E+00
 PARAMETER:  6.3859E-02 -1.8498E+00  2.5422E-01  5.6077E-01 -2.0803E-01  8.1227E-02  1.3793E+00 -4.4718E+00 -9.9342E-02  8.3915E-02
             4.4807E-01
 GRADIENT:   1.2189E+00 -6.0458E-02 -1.8161E-01  1.0724E+01  1.5539E+00  1.8682E+00 -3.6449E+00  1.3586E-03 -2.1077E+00  1.6877E-01
             6.7465E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1959.00680180191        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      924
 NPARAMETR:  9.6157E-01  6.6782E-02  1.0627E+00  1.6008E+00  6.7448E-01  9.7242E-01  5.0162E+00  1.0000E-02  8.2913E-01  9.3006E-01
             1.4097E+00
 PARAMETER:  6.0808E-02 -2.6063E+00  1.6078E-01  5.7048E-01 -2.9381E-01  7.2032E-02  1.7127E+00 -6.6610E+00 -8.7376E-02  2.7498E-02
             4.4339E-01
 GRADIENT:  -1.0387E+00 -2.9076E+00  1.4082E+00 -2.5317E+00 -2.6337E-01 -9.3572E-01 -7.0030E+00  0.0000E+00  5.4680E+00  1.1971E+00
             3.8547E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1959.47526394363        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1101
 NPARAMETR:  9.6039E-01  4.2329E-02  1.0073E+00  1.6002E+00  6.4667E-01  9.7118E-01  6.1352E+00  1.0000E-02  8.2447E-01  8.9773E-01
             1.4075E+00
 PARAMETER:  5.9589E-02 -3.0623E+00  1.0725E-01  5.7016E-01 -3.3592E-01  7.0758E-02  1.9140E+00 -8.0335E+00 -9.3020E-02 -7.8838E-03
             4.4181E-01
 GRADIENT:  -1.8522E+00 -2.0819E+00 -1.3153E+00 -8.4442E+00  3.4420E+00 -1.2608E+00 -3.6969E+00  0.0000E+00  3.6896E+00 -3.6970E-01
             4.8374E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1959.48579339682        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1166
 NPARAMETR:  9.6093E-01  4.2791E-02  1.0072E+00  1.6041E+00  6.4562E-01  9.7334E-01  6.1181E+00  1.0000E-02  8.2333E-01  8.9899E-01
             1.4076E+00
 PARAMETER:  6.0065E-02 -3.0501E+00  1.0719E-01  5.7233E-01 -3.3739E-01  7.2758E-02  1.9119E+00 -7.9770E+00 -9.3401E-02 -6.4949E-03
             4.4191E-01
 GRADIENT:  -9.4400E-01  1.7731E+01 -1.2756E-01 -1.0223E+02  1.9023E+02 -4.3610E-01  2.4527E+01  0.0000E+00  2.9938E+00 -5.1974E-02
             4.1654E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1166
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5051E-03  1.8667E-02 -1.5209E-04 -1.1942E-02 -1.3372E-02
 SE:             2.9736E-02  1.0352E-02  1.7847E-04  2.8583E-02  2.4000E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5963E-01  7.1344E-02  3.9413E-01  6.7608E-01  5.7741E-01

 ETASHRINKSD(%)  3.8032E-01  6.5320E+01  9.9402E+01  4.2425E+00  1.9596E+01
 ETASHRINKVR(%)  7.5920E-01  8.7973E+01  9.9996E+01  8.3049E+00  3.5351E+01
 EBVSHRINKSD(%)  7.3058E-01  7.7947E+01  9.9334E+01  3.6342E+00  1.6506E+01
 EBVSHRINKVR(%)  1.4558E+00  9.5137E+01  9.9996E+01  7.1364E+00  3.0288E+01
 RELATIVEINF(%)  9.8351E+01  2.8556E+00  4.4025E-04  5.1407E+01  7.0835E+00
 EPSSHRINKSD(%)  3.0837E+01
 EPSSHRINKVR(%)  5.2165E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1959.4857933968162     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1040.5472601921435     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.47
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1959.486       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.61E-01  4.28E-02  1.01E+00  1.60E+00  6.46E-01  9.73E-01  6.12E+00  1.00E-02  8.24E-01  8.99E-01  1.41E+00
 


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
+        1.25E+03
 
 TH 2
+       -5.01E+01  1.67E+05
 
 TH 3
+       -1.83E+01  1.30E+02  5.67E+02
 
 TH 4
+       -2.40E+00  4.03E+02 -1.27E+02  2.45E+03
 
 TH 5
+       -1.09E+01 -7.52E+02 -6.89E+04  9.85E+01  3.58E+04
 
 TH 6
+        4.23E+00 -5.61E+00  1.90E+00 -3.92E+00  4.62E-01  2.03E+02
 
 TH 7
+        2.24E-01  9.13E+02 -3.97E+00 -1.27E+02  3.61E+00 -1.11E-01  3.56E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.45E+00  5.00E+02  1.16E+02 -8.94E+01 -8.89E+04  1.70E+00  7.52E+00  0.00E+00  4.51E+02
 
 TH10
+       -1.06E+01 -2.72E+02  8.78E+01 -1.89E+01 -8.17E+04  4.18E+00 -9.34E+00  0.00E+00  1.11E+02  2.56E+02
 
 TH11
+       -1.41E+01 -6.40E+01  6.71E-01 -1.51E+01 -1.18E+04  2.24E+00 -1.98E+00  0.00E+00  3.05E+01  5.42E+01  2.26E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.861
Stop Time:
Wed Sep 29 22:49:00 CDT 2021
