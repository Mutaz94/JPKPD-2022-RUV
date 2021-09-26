Sat Sep 25 12:19:43 CDT 2021
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
$DATA ../../../../data/spa/S2/dat49.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1587.35333180021        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4692E+02 -4.7491E+01 -1.5156E+01 -4.0188E+01  1.4438E+01 -1.9870E+01 -3.7326E+01  1.7989E+00 -2.9011E+01 -6.5602E+00
            -1.1996E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1596.10317972559        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:       91
 NPARAMETR:  9.7864E-01  1.1242E+00  1.0842E+00  9.1858E-01  1.0989E+00  1.0650E+00  1.4914E+00  9.8302E-01  1.2349E+00  1.0521E+00
             1.0607E+00
 PARAMETER:  7.8405E-02  2.1711E-01  1.8085E-01  1.5076E-02  1.9435E-01  1.6302E-01  4.9973E-01  8.2871E-02  3.1103E-01  1.5075E-01
             1.5897E-01
 GRADIENT:   8.8945E+01 -1.1711E+01  7.9776E+00 -2.5660E+01 -4.0689E+00  1.3178E+01  2.0850E+01 -2.5959E-01  2.6447E+01  2.4867E+00
             1.4499E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1598.56496808985        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      163
 NPARAMETR:  9.6668E-01  1.1144E+00  1.1370E+00  9.5963E-01  1.1332E+00  1.0812E+00  1.4176E+00  8.7373E-01  1.1081E+00  1.1901E+00
             1.0336E+00
 PARAMETER:  6.6114E-02  2.0831E-01  2.2840E-01  5.8797E-02  2.2502E-01  1.7811E-01  4.4899E-01 -3.4989E-02  2.0268E-01  2.7403E-01
             1.3307E-01
 GRADIENT:   6.6695E+01  4.3679E+00  1.0089E+01 -1.2686E+00  5.8810E+00  2.1755E+01  9.5237E+00 -4.7750E+00  1.0724E+01  8.8072E+00
             2.1458E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1600.87655759481        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  9.3515E-01  1.1585E+00  8.2177E-01  9.1261E-01  9.8913E-01  1.0388E+00  1.3641E+00  5.6952E-01  1.0539E+00  9.6755E-01
             1.0194E+00
 PARAMETER:  3.2952E-02  2.4714E-01 -9.6297E-02  8.5486E-03  8.9073E-02  1.3805E-01  4.1050E-01 -4.6295E-01  1.5248E-01  6.7010E-02
             1.1924E-01
 GRADIENT:  -2.6816E+00  1.9371E+00 -1.8729E+00  2.8730E+00  3.2129E+00  2.8069E+00  2.0773E+00  4.9748E-01  3.8756E+00  2.3941E+00
             2.5067E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1601.62102330289        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      394
 NPARAMETR:  9.5661E-01  1.2212E+00  7.2247E-01  8.7487E-01  9.5563E-01  1.0592E+00  1.3447E+00  4.5034E-01  1.0490E+00  8.9663E-01
             1.0207E+00
 PARAMETER:  5.5636E-02  2.9980E-01 -2.2508E-01 -3.3682E-02  5.4613E-02  1.5748E-01  3.9618E-01 -6.9775E-01  1.4786E-01 -9.1119E-03
             1.2045E-01
 GRADIENT:   1.8862E+00  5.8740E-01 -3.1695E+00  2.9403E+00 -8.0035E-01  1.0655E+00  1.8213E+00  9.8325E-01  6.3689E-01  5.6335E-01
             8.1594E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1601.93163958839        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      569
 NPARAMETR:  9.5563E-01  1.4210E+00  6.5217E-01  7.4829E-01  1.0324E+00  1.0554E+00  1.1813E+00  2.6873E-01  1.1727E+00  9.4535E-01
             1.0222E+00
 PARAMETER:  5.4616E-02  4.5133E-01 -3.2745E-01 -1.8996E-01  1.3193E-01  1.5393E-01  2.6664E-01 -1.2140E+00  2.5927E-01  4.3798E-02
             1.2192E-01
 GRADIENT:  -1.2193E+00  9.4604E-01 -6.8089E-01  5.8022E-01  6.1649E-01 -4.3357E-01 -5.0240E-02  2.5786E-01 -2.3721E-02 -8.9972E-02
             1.7669E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1602.04785279962        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      746
 NPARAMETR:  9.5673E-01  1.3902E+00  6.5052E-01  7.6378E-01  1.0147E+00  1.0572E+00  1.2029E+00  1.2381E-01  1.1534E+00  9.3410E-01
             1.0206E+00
 PARAMETER:  5.5764E-02  4.2945E-01 -3.2998E-01 -1.6947E-01  1.1457E-01  1.5565E-01  2.8473E-01 -1.9890E+00  2.4271E-01  3.1830E-02
             1.2040E-01
 GRADIENT:   1.0804E+00 -2.0556E+00 -7.0010E-01 -1.3119E+00  1.0804E+00  2.2892E-01  1.3569E-02  6.0243E-02  1.2474E-01  1.3985E-01
            -3.7517E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1602.05379688453        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:      919
 NPARAMETR:  9.5615E-01  1.3891E+00  6.5121E-01  7.6661E-01  1.0119E+00  1.0565E+00  1.2050E+00  1.2111E-01  1.1512E+00  9.3315E-01
             1.0206E+00
 PARAMETER:  5.5157E-02  4.2868E-01 -3.2893E-01 -1.6578E-01  1.1188E-01  1.5499E-01  2.8645E-01 -2.0110E+00  2.4081E-01  3.0813E-02
             1.2037E-01
 GRADIENT:  -1.9158E+01  1.8419E+01 -2.5209E+00  1.8937E+01 -2.4407E+00 -1.9171E+01 -7.8505E+00  5.4736E-02 -8.4896E+00  2.7525E-01
            -7.1920E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1602.05379688453        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  9.5615E-01  1.3891E+00  6.5121E-01  7.6661E-01  1.0119E+00  1.0565E+00  1.2050E+00  1.2111E-01  1.1512E+00  9.3315E-01
             1.0206E+00
 PARAMETER:  5.5157E-02  4.2868E-01 -3.2893E-01 -1.6578E-01  1.1188E-01  1.5499E-01  2.8645E-01 -2.0110E+00  2.4081E-01  3.0813E-02
             1.2037E-01
 GRADIENT:  -1.0808E-01  2.4794E+01  1.0816E+06 -1.6223E+04 -7.5545E-01  2.8540E+00 -1.1891E+06 -1.7825E+05  2.3551E+05 -2.0561E+03
             2.2215E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      944
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4409E-05 -1.2970E-02 -4.7523E-03  1.0256E-02 -2.5457E-02
 SE:             2.9842E-02  2.4766E-02  1.6810E-03  2.2695E-02  2.1485E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9774E-01  6.0048E-01  4.6971E-03  6.5132E-01  2.3607E-01

 ETASHRINKSD(%)  2.5559E-02  1.7031E+01  9.4368E+01  2.3969E+01  2.8021E+01
 ETASHRINKVR(%)  5.1111E-02  3.1161E+01  9.9683E+01  4.2192E+01  4.8190E+01
 EBVSHRINKSD(%)  4.0967E-01  1.6487E+01  9.5156E+01  2.5209E+01  2.6736E+01
 EBVSHRINKVR(%)  8.1765E-01  3.0256E+01  9.9765E+01  4.4063E+01  4.6324E+01
 RELATIVEINF(%)  9.9074E+01  6.3769E+00  3.1591E-02  4.6656E+00  8.6270E+00
 EPSSHRINKSD(%)  4.4213E+01
 EPSSHRINKVR(%)  6.8878E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1602.0537968845335     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -866.90297032079536     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.41
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1602.054       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.56E-01  1.39E+00  6.51E-01  7.67E-01  1.01E+00  1.06E+00  1.20E+00  1.21E-01  1.15E+00  9.33E-01  1.02E+00
 


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
+        7.88E+07
 
 TH 2
+       -1.48E+08  9.82E+02
 
 TH 3
+       -3.41E+07 -1.75E+05  1.98E+09
 
 TH 4
+        6.77E+08  6.74E+06  4.43E+07  6.87E+08
 
 TH 5
+        4.49E+08  1.08E+08  4.15E+06 -1.20E+08  5.09E+08
 
 TH 6
+        3.62E+05 -2.36E+07 -5.49E+07 -1.40E+03 -5.09E+08  3.47E+05
 
 TH 7
+        2.28E+08  1.29E+07  1.93E+06  1.86E+06  1.97E+08 -1.00E+06  2.11E+03
 
 TH 8
+        9.75E+06 -6.37E+03 -1.41E+09 -9.78E+06 -2.44E+06  8.06E+06  1.07E+06  1.55E+09
 
 TH 9
+       -3.00E+08 -3.79E+07  1.30E+09  3.87E+08  3.33E+08 -2.66E+07 -4.57E+07 -2.53E+07  1.14E+09
 
 TH10
+       -2.74E+09 -2.45E+06  3.72E+07 -3.59E+06  2.60E+08  3.32E+06  2.34E+08 -4.99E+05  7.62E+07  1.98E+08
 
 TH11
+        2.78E+03  1.64E+08  1.13E+07 -3.46E+08 -7.03E+08 -1.06E+08  2.67E+07  1.31E+07  1.33E+08  5.21E+08  6.32E+03
 
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
 #CPUT: Total CPU Time in Seconds,       21.245
Stop Time:
Sat Sep 25 12:20:06 CDT 2021
