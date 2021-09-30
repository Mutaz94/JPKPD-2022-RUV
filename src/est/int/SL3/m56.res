Wed Sep 29 04:28:34 CDT 2021
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
$DATA ../../../../data/int/SL3/dat56.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      970
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      870
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
 RAW OUTPUT FILE (FILE): m56.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -590.426908623787        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1365E+02 -3.4236E+01  1.9445E+02  1.6044E+02  1.4317E+02  4.0931E+01 -1.1697E+02 -2.6111E+02 -1.0484E+02 -2.1310E+01
            -5.8041E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2683.69676637820        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0271E+00  1.3942E+00  8.9642E-01  8.4857E-01  1.1434E+00  1.0755E+00  1.2256E+00  9.5584E-01  9.5366E-01  9.3954E-01
             2.6447E+00
 PARAMETER:  1.2674E-01  4.3232E-01 -9.3501E-03 -6.4208E-02  2.3397E-01  1.7282E-01  3.0345E-01  5.4836E-02  5.2550E-02  3.7631E-02
             1.0725E+00
 GRADIENT:   8.2638E+01  8.8494E+01 -9.2801E+00  3.6878E+01  4.4175E+00  2.6420E+01  2.9660E+01  4.9577E+00 -5.7437E+00 -1.8687E+01
             1.2724E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2691.36485745191        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      200
 NPARAMETR:  1.0455E+00  1.6031E+00  9.4847E-01  7.5359E-01  1.3020E+00  1.0516E+00  9.0340E-01  6.3397E-01  9.6891E-01  1.1336E+00
             2.6690E+00
 PARAMETER:  1.4445E-01  5.7197E-01  4.7094E-02 -1.8290E-01  3.6389E-01  1.5034E-01 -1.5908E-03 -3.5575E-01  6.8417E-02  2.2537E-01
             1.0817E+00
 GRADIENT:   3.4091E+01  6.2043E+01  4.1236E+00  6.4721E+01 -1.6197E+01  2.8256E+00 -1.7406E+01 -1.1531E-01 -1.4173E+01 -2.0411E+01
             1.6621E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2703.40314277977        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0233E+00  1.8680E+00  2.4490E+00  6.2062E-01  2.0001E+00  1.0359E+00  1.0274E+00  3.7162E+00  1.0636E+00  1.7582E+00
             2.5912E+00
 PARAMETER:  1.2306E-01  7.2485E-01  9.9566E-01 -3.7703E-01  7.9319E-01  1.3529E-01  1.2704E-01  1.4127E+00  1.6170E-01  6.6427E-01
             1.0521E+00
 GRADIENT:  -5.5165E+00  6.6233E+01 -1.7865E+01  6.7446E+01  6.2358E+01 -3.1713E+00  2.5537E+01 -5.6383E+00  5.5923E+00  1.8424E+01
             2.7758E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2717.43560198871        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      555
 NPARAMETR:  1.0219E+00  2.0816E+00  2.3941E+00  4.1884E-01  1.8584E+00  1.0570E+00  8.0752E-01  4.8083E+00  1.0486E+00  1.5413E+00
             2.5332E+00
 PARAMETER:  1.2168E-01  8.3311E-01  9.7303E-01 -7.7028E-01  7.1971E-01  1.5547E-01 -1.1379E-01  1.6703E+00  1.4749E-01  5.3264E-01
             1.0295E+00
 GRADIENT:  -4.9625E+00  2.4786E+01 -4.0623E-01  1.0091E+01 -1.2428E+01  5.4990E+00 -1.1186E+01 -2.7274E+00 -2.1508E+00 -1.6085E+01
            -2.4456E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2718.03140996952        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  1.0223E+00  2.0855E+00  2.4882E+00  4.1717E-01  1.8699E+00  1.0395E+00  8.2912E-01  4.9793E+00  1.1000E+00  1.5610E+00
             2.5316E+00
 PARAMETER:  1.2208E-01  8.3500E-01  1.0115E+00 -7.7426E-01  7.2589E-01  1.3871E-01 -8.7394E-02  1.7053E+00  1.9528E-01  5.4532E-01
             1.0289E+00
 GRADIENT:   8.2847E+01  2.9670E+02 -1.2975E-01  3.2436E+01  1.7189E+01  1.3571E+01 -4.3298E-01  1.6584E+00  4.2333E-01 -8.3703E+00
            -5.6848E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2718.16093642414        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  1.0226E+00  2.0807E+00  2.4847E+00  4.1709E-01  1.8705E+00  1.0423E+00  8.2777E-01  4.9956E+00  1.1423E+00  1.5613E+00
             2.5429E+00
 PARAMETER:  1.2231E-01  8.3270E-01  1.0102E+00 -7.7444E-01  7.2622E-01  1.4140E-01 -8.9019E-02  1.7085E+00  2.3308E-01  5.4550E-01
             1.0333E+00
 GRADIENT:  -3.9290E+00  1.5925E+01 -1.1886E+00  8.8050E+00 -1.0560E+01  2.8579E-01 -1.3448E+00  5.0086E-01  6.2633E-01 -1.3653E+01
            -1.3437E+01

0ITERATION NO.:   33    OBJECTIVE VALUE:  -2718.24905531735        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1001
 NPARAMETR:  1.0228E+00  2.0769E+00  2.4812E+00  4.1710E-01  1.8708E+00  1.0424E+00  8.3081E-01  5.0117E+00  1.1092E+00  1.5613E+00
             2.5520E+00
 PARAMETER:  1.2251E-01  8.3068E-01  1.0089E+00 -7.7459E-01  7.2650E-01  1.4159E-01 -8.5794E-02  1.7114E+00  2.0574E-01  5.4566E-01
             1.0371E+00
 GRADIENT:  -3.5346E+03 -2.4195E+03  8.5955E+02 -1.1108E+03  1.1826E+03  2.0150E-01 -1.2980E+00 -5.1976E+02  3.9519E-01  1.5733E+03
             8.2930E+02
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.3         2.1         1.0         2.3         0.6         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1001
 NO. OF SIG. DIGITS IN FINAL EST.:  0.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3331E-03 -1.7064E-02 -3.3422E-02  7.9417E-03 -2.1512E-02
 SE:             2.9487E-02  2.6606E-02  1.2966E-02  1.3226E-02  2.5788E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1000E-01  5.2127E-01  9.9489E-03  5.4821E-01  4.0417E-01

 ETASHRINKSD(%)  1.2144E+00  1.0868E+01  5.6562E+01  5.5690E+01  1.3607E+01
 ETASHRINKVR(%)  2.4140E+00  2.0555E+01  8.1131E+01  8.0366E+01  2.5363E+01
 EBVSHRINKSD(%)  1.3651E+00  1.1375E+01  6.4510E+01  6.1582E+01  1.3124E+01
 EBVSHRINKVR(%)  2.7116E+00  2.1455E+01  8.7404E+01  8.5241E+01  2.4525E+01
 RELATIVEINF(%)  9.7240E+01  5.2700E+00  4.0944E+00  8.9603E-01  4.7126E+01
 EPSSHRINKSD(%)  1.7005E+01
 EPSSHRINKVR(%)  3.1118E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          870
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1598.9530477761305     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2718.2490553173498     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1119.2960075412193     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.33
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.35
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2718.249       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.08E+00  2.48E+00  4.17E-01  1.87E+00  1.04E+00  8.30E-01  5.01E+00  1.11E+00  1.56E+00  2.55E+00
 


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
+        1.38E+06
 
 TH 2
+       -8.42E+00  7.60E+03
 
 TH 3
+       -3.26E+00 -7.19E+01  3.48E+03
 
 TH 4
+        1.03E+01  6.05E+02 -2.23E+02  2.09E+05
 
 TH 5
+       -8.77E+00 -4.82E+01  1.49E+01 -1.14E+02  1.18E+04
 
 TH 6
+        9.47E+01 -1.41E+00 -4.76E+00  3.39E+01 -1.01E+01  1.74E+02
 
 TH 7
+        3.11E+02  2.04E+01 -1.56E+01  7.01E+01 -2.89E+01 -2.43E+00  1.84E+02
 
 TH 8
+        1.34E+00  9.38E+01 -9.44E+00  4.69E+02 -7.51E+00  1.95E+00  4.53E+00  7.28E+02
 
 TH 9
+        1.28E+02  1.17E+01 -6.77E+00  4.81E+01 -1.12E+01  3.54E-01 -1.14E+06  3.29E+00  4.14E+05
 
 TH10
+       -9.57E+00 -3.18E+01  1.72E+01 -1.18E+02  2.59E+01 -1.32E+01 -4.61E+01 -6.63E+00 -1.75E+01  2.99E+04
 
 TH11
+       -1.62E+01  1.24E+01  8.64E-01 -1.74E+01  2.22E+00 -2.72E+00 -7.11E+00 -1.30E+03 -4.35E+00  9.61E+03  3.26E+03
 
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
 #CPUT: Total CPU Time in Seconds,       41.798
Stop Time:
Wed Sep 29 04:29:18 CDT 2021
