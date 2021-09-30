Wed Sep 29 14:46:33 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat1.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1666.05957175239        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3207E+02  2.5353E+01  3.4964E+01  2.1758E+01 -4.2823E+01  4.9700E+01  2.6839E+00 -1.0614E+00  1.2661E+01 -2.2872E+01
             1.0445E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1669.71995779271        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  9.9875E-01  1.0244E+00  9.5673E-01  9.9759E-01  1.0598E+00  9.8201E-01  1.0067E+00  9.8668E-01  9.2908E-01  1.2624E+00
             9.2763E-01
 PARAMETER:  9.8752E-02  1.2409E-01  5.5764E-02  9.7583E-02  1.5808E-01  8.1841E-02  1.0665E-01  8.6590E-02  2.6436E-02  3.3303E-01
             2.4873E-02
 GRADIENT:  -1.5597E+00 -2.4840E+01 -1.7236E+01 -2.5746E+01 -3.6870E-01 -4.3928E+00  4.5768E-01  6.0265E+00 -8.7973E+00  7.3956E+00
            -1.6280E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1671.48096642912        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.9823E-01  1.2140E+00  1.0342E+00  8.9559E-01  1.2196E+00  9.9336E-01  8.0223E-01  1.0336E+00  1.1080E+00  1.4893E+00
             9.5250E-01
 PARAMETER:  9.8225E-02  2.9391E-01  1.3366E-01 -1.0270E-02  2.9848E-01  9.3341E-02 -1.2036E-01  1.3308E-01  2.0256E-01  4.9831E-01
             5.1339E-02
 GRADIENT:  -2.9516E+00 -4.5831E+00  1.3913E+00 -2.5379E+00 -4.2477E+00  5.0952E-01  1.7178E+00 -3.5296E+00 -1.2419E-01  1.5150E+01
            -4.1561E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1672.75966165684        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.9980E-01  1.2876E+00  1.1006E+00  8.5884E-01  1.2751E+00  9.9251E-01  7.2852E-01  1.3720E+00  1.1742E+00  1.3927E+00
             9.6198E-01
 PARAMETER:  9.9798E-02  3.5279E-01  1.9589E-01 -5.2173E-02  3.4299E-01  9.2484E-02 -2.1674E-01  4.1626E-01  2.6060E-01  4.3125E-01
             6.1235E-02
 GRADIENT:   1.6880E-01  7.1347E+00  1.6881E+00  1.0374E+01 -1.6011E+00  6.5361E-02  6.5065E-02 -4.3912E-01  7.7780E-01 -5.2280E-01
            -3.3154E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1673.34638117096        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  1.0012E+00  1.5957E+00  8.7677E-01  6.5400E-01  1.3397E+00  9.9613E-01  6.8152E-01  1.3933E+00  1.4136E+00  1.4077E+00
             9.6183E-01
 PARAMETER:  1.0120E-01  5.6734E-01 -3.1515E-02 -3.2465E-01  3.9247E-01  9.6122E-02 -2.8343E-01  4.3168E-01  4.4614E-01  4.4196E-01
             6.1081E-02
 GRADIENT:   3.3832E-02  2.1694E+01  9.5527E+00  5.6333E+00 -1.2358E+01  8.1852E-01  1.7658E-01 -2.1537E+00 -1.6286E+00 -6.3973E-01
            -8.1994E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1673.79463495715        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  1.0016E+00  1.7871E+00  6.8172E-01  5.2703E-01  1.3832E+00  9.9525E-01  6.6272E-01  1.2842E+00  1.6379E+00  1.4409E+00
             9.5821E-01
 PARAMETER:  1.0160E-01  6.8061E-01 -2.8313E-01 -5.4050E-01  4.2441E-01  9.5243E-02 -3.1140E-01  3.5013E-01  5.9342E-01  4.6529E-01
             5.7309E-02
 GRADIENT:  -7.0181E-01  2.5383E+01  7.1741E+00  7.3200E+00 -1.1709E+01  6.6650E-02 -1.7586E-01 -1.3896E+00 -1.1362E+00  2.7065E+00
            -2.1006E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1674.03389633618        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  1.0021E+00  1.8821E+00  5.6783E-01  4.6377E-01  1.3964E+00  9.9580E-01  6.6192E-01  1.1248E+00  1.7684E+00  1.4304E+00
             9.5922E-01
 PARAMETER:  1.0210E-01  7.3241E-01 -4.6593E-01 -6.6836E-01  4.3393E-01  9.5790E-02 -3.1261E-01  2.1763E-01  6.7005E-01  4.5798E-01
             5.8367E-02
 GRADIENT:  -4.0151E-01  2.6100E+01  5.0644E+00  8.4804E+00 -1.0156E+01  3.1498E-02  2.9892E-02 -7.0426E-01 -8.8721E-01  2.6260E+00
            -1.8526E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1674.34926766781        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1243
 NPARAMETR:  1.0016E+00  1.8797E+00  5.4080E-01  4.4862E-01  1.3973E+00  9.9588E-01  6.6077E-01  1.1161E+00  1.7811E+00  1.4267E+00
             9.5977E-01
 PARAMETER:  1.0160E-01  7.3112E-01 -5.1471E-01 -7.0158E-01  4.3453E-01  9.5866E-02 -3.1435E-01  2.0986E-01  6.7721E-01  4.5534E-01
             5.8940E-02
 GRADIENT:  -1.4822E+00 -1.4007E+01  6.2591E-01 -7.6050E-01 -4.8864E+00  4.4954E-02  5.8824E-01  3.2241E-01 -7.0887E-01  3.5113E+00
            -7.7806E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1674.35033037041        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1413
 NPARAMETR:  1.0018E+00  1.8811E+00  5.4072E-01  4.4892E-01  1.3964E+00  9.9561E-01  6.5935E-01  1.0959E+00  1.7850E+00  1.4277E+00
             9.6081E-01
 PARAMETER:  1.0173E-01  7.3163E-01 -5.1408E-01 -7.0072E-01  4.3400E-01  9.5564E-02 -3.1643E-01  1.9354E-01  6.7962E-01  4.5590E-01
             5.9932E-02
 GRADIENT:  -3.3056E+04 -9.1941E+03  9.7081E-01  9.5843E+03  1.5492E+04 -1.3966E-01  1.7869E-01  1.4790E-01  9.8416E+03 -1.4746E+04
            -3.3944E-01
 NUMSIGDIG:         2.3         2.3         1.6         2.3         2.3         2.1         2.4         0.8         2.3         2.3
                    1.8

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1413
 NO. OF SIG. DIGITS IN FINAL EST.:  0.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7963E-03 -5.0274E-02 -2.5418E-02  3.7590E-02 -4.9908E-02
 SE:             2.9912E-02  2.1470E-02  8.2125E-03  2.3078E-02  2.2322E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5211E-01  1.9203E-02  1.9683E-03  1.0335E-01  2.5364E-02

 ETASHRINKSD(%)  1.0000E-10  2.8072E+01  7.2487E+01  2.2685E+01  2.5218E+01
 ETASHRINKVR(%)  1.0000E-10  4.8263E+01  9.2430E+01  4.0224E+01  4.4077E+01
 EBVSHRINKSD(%)  4.0460E-01  2.6037E+01  7.5339E+01  2.5133E+01  1.9912E+01
 EBVSHRINKVR(%)  8.0756E-01  4.5294E+01  9.3918E+01  4.3949E+01  3.5860E+01
 RELATIVEINF(%)  9.9140E+01  4.4785E+00  1.0081E+00  4.9580E+00  2.6501E+01
 EPSSHRINKSD(%)  4.5249E+01
 EPSSHRINKVR(%)  7.0023E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1674.3503303704072     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.19950380666899     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1674.350       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.88E+00  5.41E-01  4.49E-01  1.40E+00  9.96E-01  6.59E-01  1.10E+00  1.79E+00  1.43E+00  9.61E-01
 


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
+        1.62E+07
 
 TH 2
+        1.20E+06  8.91E+04
 
 TH 3
+       -3.26E+02  4.38E+05  2.16E+06
 
 TH 4
+       -5.24E+06 -3.89E+05 -1.92E+06  3.40E+06
 
 TH 5
+       -2.72E+06 -2.02E+05 -1.04E+02  8.81E+05  4.58E+05
 
 TH 6
+       -1.27E+03 -9.54E+01  2.29E+00  4.10E+02  2.14E+02  1.98E+02
 
 TH 7
+        1.76E+02 -2.86E+00  2.89E+06 -6.91E+01 -4.29E+01 -5.40E-01  3.86E+06
 
 TH 8
+        3.37E+04  2.49E+03 -2.13E+01 -1.09E+04 -5.67E+03 -9.44E-02  5.13E+00  3.64E+00
 
 TH 9
+        7.84E+01  1.05E+03 -2.08E+01  8.81E+05  5.15E+01  1.08E+02  3.72E+00 -2.80E+03  2.66E+01
 
 TH10
+       -2.54E+06 -1.88E+05 -1.95E+02 -1.64E+06  4.26E+05 -1.99E+02  3.44E+01  5.28E+03  2.25E+00  4.40E+01
 
 TH11
+        4.26E+03  2.99E+02  6.28E+06 -1.39E+03 -7.17E+02  2.61E+00  8.39E+06  1.99E+00 -3.53E+02 -2.69E+06  2.31E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.707
Stop Time:
Wed Sep 29 14:47:02 CDT 2021
