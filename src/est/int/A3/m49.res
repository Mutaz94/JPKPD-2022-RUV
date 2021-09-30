Wed Sep 29 00:06:13 CDT 2021
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
$DATA ../../../../data/int/A3/dat49.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -627.255815771601        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3780E+02  1.5998E+02  3.3485E+02 -2.3188E+01  2.4741E+02  2.3736E+01 -2.3765E+02 -3.9858E+02 -9.4009E+01 -1.5555E+02
            -5.4158E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2762.12578007137        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.8462E-01  9.4714E-01  8.0482E-01  1.0479E+00  8.1782E-01  9.5031E-01  1.5443E+00  8.6386E-01  1.0810E+00  1.0340E+00
             2.5715E+00
 PARAMETER:  8.4499E-02  4.5689E-02 -1.1714E-01  1.4683E-01 -1.0111E-01  4.9029E-02  5.3455E-01 -4.6341E-02  1.7787E-01  1.3347E-01
             1.0445E+00
 GRADIENT:   1.7228E+02  3.3643E+01  1.5233E+01 -3.6357E+01 -4.1118E+01 -5.0327E+00  1.3970E+01  1.7714E+01  1.2446E+01  1.3289E+01
             1.6904E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2775.44917231954        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      220
 NPARAMETR:  9.7590E-01  6.9157E-01  6.3796E-01  1.2407E+00  6.0548E-01  1.0535E+00  1.9299E+00  4.5692E-01  9.4542E-01  8.4183E-01
             2.5051E+00
 PARAMETER:  7.5603E-02 -2.6880E-01 -3.4948E-01  3.1565E-01 -4.0173E-01  1.5212E-01  7.5746E-01 -6.8324E-01  4.3874E-02 -7.2183E-02
             1.0183E+00
 GRADIENT:   7.5378E+01  3.4790E+01  3.7705E+01  6.9706E+01 -4.2050E+01  2.4582E+01  2.8331E+01  7.2600E+00 -2.9984E+01  4.3165E+00
             1.2824E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2799.28497703927        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      398
 NPARAMETR:  9.3153E-01  5.6374E-01  5.0317E-01  1.2852E+00  5.0138E-01  9.5725E-01  1.5943E+00  3.3064E-02  1.0652E+00  8.1521E-01
             2.2854E+00
 PARAMETER:  2.9073E-02 -4.7316E-01 -5.8682E-01  3.5095E-01 -5.9039E-01  5.6310E-02  5.6644E-01 -3.3093E+00  1.6313E-01 -1.0431E-01
             9.2656E-01
 GRADIENT:  -1.2334E+01  2.4218E+01  2.4076E+01  8.8188E+01  5.6168E+00 -8.0834E+00 -1.6327E+00  4.0274E-02 -1.7329E+01 -4.6285E+00
            -2.4154E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2846.04393637959        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      587             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3200E-01  2.8263E-01  2.1204E-01  1.1346E+00  2.5255E-01  9.9233E-01  1.5310E+00  1.0000E-02  1.2543E+00  8.7068E-01
             1.9486E+00
 PARAMETER:  2.9581E-02 -1.1636E+00 -1.4510E+00  2.2626E-01 -1.2762E+00  9.2304E-02  5.2592E-01 -1.5692E+01  3.2654E-01 -3.8483E-02
             7.6712E-01
 GRADIENT:   8.0271E+01  1.0500E+02  9.0041E+01  6.0918E+01  5.3848E+02  1.2369E+01  8.7708E+00  0.0000E+00 -1.4373E+01  2.3583E+00
            -1.5610E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2846.23439995971        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      744
 NPARAMETR:  9.3865E-01  2.7810E-01  2.1209E-01  1.1312E+00  2.5239E-01  9.7578E-01  1.5519E+00  1.0000E-02  1.2543E+00  8.7197E-01
             1.9506E+00
 PARAMETER:  3.6686E-02 -1.1798E+00 -1.4508E+00  2.2327E-01 -1.2768E+00  7.5485E-02  5.3948E-01 -1.5692E+01  3.2657E-01 -3.7004E-02
             7.6813E-01
 GRADIENT:   1.1517E+01 -4.4348E+00 -8.5216E+00 -3.2478E-01  6.8288E+01 -3.7370E+00  1.4092E+00  0.0000E+00 -3.7216E+01  2.0502E-01
            -1.6820E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2852.98488467705        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      923
 NPARAMETR:  9.3614E-01  2.7808E-01  2.1471E-01  1.1311E+00  2.4285E-01  9.8219E-01  1.5133E+00  1.0000E-02  1.2567E+00  8.5685E-01
             2.0761E+00
 PARAMETER:  3.4014E-02 -1.1798E+00 -1.4385E+00  2.2316E-01 -1.3153E+00  8.2033E-02  5.1426E-01 -1.5692E+01  3.2845E-01 -5.4494E-02
             8.3050E-01
 GRADIENT:   1.8748E+00  1.2872E+00  5.4974E+01 -2.0130E+00 -3.8291E+01 -2.6627E-01 -2.0668E-01  0.0000E+00 -3.3254E+01  1.6340E-01
            -3.2168E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2855.74185465850        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.3458E-01  2.7686E-01  2.0741E-01  1.1315E+00  2.4417E-01  9.8287E-01  1.5143E+00  1.0000E-02  1.3518E+00  8.5591E-01
             2.1015E+00
 PARAMETER:  3.2343E-02 -1.1843E+00 -1.4730E+00  2.2354E-01 -1.3099E+00  8.2723E-02  5.1493E-01 -1.5692E+01  4.0143E-01 -5.5587E-02
             8.4266E-01
 GRADIENT:  -1.4879E+00  5.2064E+00  4.4614E+00  2.8349E+00  3.0117E+01 -5.0761E-01  6.2125E-02  0.0000E+00 -9.6095E+00  2.1078E-01
             2.8962E+00

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2855.74185465850        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1088
 NPARAMETR:  9.3498E-01  2.7618E-01  2.0808E-01  1.1293E+00  2.4486E-01  9.8386E-01  1.5143E+00  1.0000E-02  1.3506E+00  8.5549E-01
             2.1054E+00
 PARAMETER:  3.2343E-02 -1.1843E+00 -1.4730E+00  2.2354E-01 -1.3099E+00  8.2723E-02  5.1493E-01 -1.5692E+01  4.0143E-01 -5.5587E-02
             8.4266E-01
 GRADIENT:  -1.1043E+00  4.0983E+00 -5.5601E+03  2.4752E+00 -6.2479E+03 -4.3517E-01  3.3616E-03  0.0000E+00  2.0488E+04  1.4496E-01
            -9.9258E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1088
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6282E-03  1.0760E-03  3.6690E-05 -6.8498E-03 -1.1512E-02
 SE:             2.9509E-02  2.5167E-02  2.5733E-04  2.9107E-02  2.6940E-02
 N:                     100         100         100         100         100

 P VAL.:         9.2903E-01  9.6590E-01  8.8662E-01  8.1395E-01  6.6915E-01

 ETASHRINKSD(%)  1.1415E+00  1.5687E+01  9.9138E+01  2.4879E+00  9.7490E+00
 ETASHRINKVR(%)  2.2701E+00  2.8914E+01  9.9993E+01  4.9138E+00  1.8547E+01
 EBVSHRINKSD(%)  1.3414E+00  1.3711E+01  9.9143E+01  3.7911E+00  1.0185E+01
 EBVSHRINKVR(%)  2.6647E+00  2.5543E+01  9.9993E+01  7.4384E+00  1.9333E+01
 RELATIVEINF(%)  9.7314E+01  1.7977E+01  5.9859E-04  6.0961E+01  6.9107E+00
 EPSSHRINKSD(%)  2.0518E+01
 EPSSHRINKVR(%)  3.6826E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2855.7418546585036     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1201.6524948900928     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2855.742       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.35E-01  2.77E-01  2.07E-01  1.13E+00  2.44E-01  9.83E-01  1.51E+00  1.00E-02  1.35E+00  8.56E-01  2.10E+00
 


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
+        1.28E+03
 
 TH 2
+        5.21E+00  4.77E+03
 
 TH 3
+       -1.72E+03 -2.43E+03  2.20E+06
 
 TH 4
+        4.01E-01 -4.14E+01  1.70E+03  4.41E+02
 
 TH 5
+       -1.68E+03 -2.22E+03  2.08E+06 -2.54E+06  2.03E+06
 
 TH 6
+        6.87E+00 -6.15E+00 -1.43E+03 -4.82E+00 -1.40E+03  1.96E+02
 
 TH 7
+        4.66E-02  2.36E+01 -1.06E+01 -1.39E+00 -4.10E+01 -6.48E-01  4.74E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.06E+06  1.16E+06 -1.23E+06  1.50E+06 -1.18E+06 -3.86E+06  4.93E+00  0.00E+00  6.98E+05
 
 TH10
+       -3.46E+00  1.39E+00 -7.19E+02  1.16E+01 -1.00E+03  1.88E+00  5.47E+00  0.00E+00 -4.43E+06  1.72E+02
 
 TH11
+       -3.22E+02 -2.04E+01  3.85E+05  3.54E+02  3.69E+05 -2.48E+02  6.46E+00  0.00E+00 -2.18E+05 -1.42E+02  6.80E+04
 
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
 #CPUT: Total CPU Time in Seconds,       43.523
Stop Time:
Wed Sep 29 00:06:58 CDT 2021
