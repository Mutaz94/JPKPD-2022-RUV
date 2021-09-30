Wed Sep 29 19:01:48 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat44.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m44.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1579.31084223941        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.6026E+02  6.3898E+01  1.7864E+01  1.0483E+02  2.1850E+01  2.2868E+01  3.6831E+00 -2.2073E+01  1.6831E+01 -7.6873E+00
            -3.3160E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1587.60605036030        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.3642E-01  9.6849E-01  8.5652E-01  9.7996E-01  8.9813E-01  1.0317E+00  9.5021E-01  1.2804E+00  8.2546E-01  8.9892E-01
             1.0092E+00
 PARAMETER:  3.4311E-02  6.7979E-02 -5.4872E-02  7.9757E-02 -7.4430E-03  1.3124E-01  4.8927E-02  3.4720E-01 -9.1814E-02 -6.5653E-03
             1.0912E-01
 GRADIENT:   4.0078E+02  2.4944E+01 -1.6709E+00  5.1410E+01  2.6798E+00  6.8944E+01 -1.1321E+01  4.3317E+00 -1.8484E+01  4.4719E+00
            -2.5879E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1588.83817464333        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      265
 NPARAMETR:  9.3853E-01  8.4712E-01  1.0171E+00  1.0769E+00  9.1635E-01  1.0540E+00  1.2898E+00  1.4756E+00  7.6578E-01  8.4934E-01
             1.1561E+00
 PARAMETER:  3.6555E-02 -6.5908E-02  1.1699E-01  1.7404E-01  1.2645E-02  1.5259E-01  3.5449E-01  4.8904E-01 -1.6686E-01 -6.3291E-02
             2.4503E-01
 GRADIENT:   1.0909E+01  9.9727E+00 -4.3791E+00  1.1366E+01 -7.5052E+00  7.2877E+00  1.1040E+00  8.7807E+00 -4.1963E+00 -6.8746E-01
             2.7464E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1591.50827445582        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  9.2936E-01  7.3505E-01  1.3351E+00  1.1663E+00  1.0095E+00  1.0282E+00  1.1756E+00  1.4924E+00  8.4320E-01  1.0468E+00
             1.0639E+00
 PARAMETER:  2.6745E-02 -2.0782E-01  3.8900E-01  2.5381E-01  1.0949E-01  1.2776E-01  2.6179E-01  5.0039E-01 -7.0551E-02  1.4574E-01
             1.6193E-01
 GRADIENT:  -1.7532E+00  1.5218E+01  5.0287E+00  1.9371E+01 -6.4896E+00 -1.7762E+00  2.2756E+00 -2.7922E+00  1.4176E+00  9.4175E-01
            -3.5646E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1594.23458235023        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      616
 NPARAMETR:  9.2754E-01  4.5250E-01  1.6544E+00  1.3366E+00  1.0123E+00  1.0365E+00  4.8408E-01  1.6901E+00  8.2043E-01  1.0772E+00
             1.0657E+00
 PARAMETER:  2.4784E-02 -6.9296E-01  6.0343E-01  3.9011E-01  1.1222E-01  1.3588E-01 -6.2550E-01  6.2481E-01 -9.7929E-02  1.7438E-01
             1.6364E-01
 GRADIENT:   1.9400E+00  3.1089E+00  8.0049E+00 -1.0914E+00 -1.0502E+01  2.8132E+00  5.3814E-01 -2.1748E+00  1.7774E+00 -6.0835E-01
            -3.2823E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1594.57033205549        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      791
 NPARAMETR:  9.2479E-01  3.4357E-01  1.6460E+00  1.4101E+00  9.9315E-01  1.0262E+00  2.4091E-01  1.6937E+00  7.8048E-01  1.0758E+00
             1.0700E+00
 PARAMETER:  2.1811E-02 -9.6838E-01  5.9834E-01  4.4367E-01  9.3123E-02  1.2589E-01 -1.3233E+00  6.2694E-01 -1.4784E-01  1.7304E-01
             1.6767E-01
 GRADIENT:  -1.3997E+00  2.0629E+00 -1.4085E+00  9.3478E+00  1.4532E+00 -6.9502E-01  7.5486E-02 -2.1466E-02  9.8506E-02  4.7681E-01
            -7.6212E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1594.62283962258        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      972
 NPARAMETR:  9.2537E-01  3.2928E-01  1.6597E+00  1.4147E+00  9.8985E-01  1.0280E+00  1.1693E-01  1.7021E+00  7.7594E-01  1.0701E+00
             1.0716E+00
 PARAMETER:  2.2434E-02 -1.0109E+00  6.0664E-01  4.4690E-01  8.9801E-02  1.2764E-01 -2.0461E+00  6.3184E-01 -1.5368E-01  1.6774E-01
             1.6918E-01
 GRADIENT:   4.0802E-01  3.3001E-01  1.5994E-01 -2.0546E+00 -1.3395E-01  9.4466E-02  1.8060E-02 -2.5822E-01  2.4728E-01  2.3245E-01
            -1.4127E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1594.63111295143        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1153
 NPARAMETR:  9.2542E-01  3.2756E-01  1.6600E+00  1.4148E+00  9.8927E-01  1.0281E+00  3.1114E-02  1.7067E+00  7.7490E-01  1.0674E+00
             1.0720E+00
 PARAMETER:  2.2489E-02 -1.0161E+00  6.0684E-01  4.4698E-01  8.9209E-02  1.2769E-01 -3.3701E+00  6.3458E-01 -1.5502E-01  1.6522E-01
             1.6952E-01
 GRADIENT:   5.8058E-01 -1.0257E-02 -4.7709E-03 -4.3993E+00  6.2699E-02  1.2355E-01  1.4988E-03  1.6188E-03 -1.4256E-02  7.0105E-02
            -6.8705E-03

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1594.63358778363        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1249
 NPARAMETR:  9.2606E-01  3.2857E-01  1.6597E+00  1.4132E+00  9.8913E-01  1.0290E+00  1.0000E-02  1.7076E+00  7.7507E-01  1.0662E+00
             1.0720E+00
 PARAMETER:  2.3187E-02 -1.0130E+00  6.0663E-01  4.4583E-01  8.9073E-02  1.2855E-01 -5.3218E+00  6.3508E-01 -1.5480E-01  1.6414E-01
             1.6955E-01
 GRADIENT:   2.0215E+00 -3.2855E-01  1.1381E-01 -6.7235E+00 -1.0488E-02  4.6683E-01  0.0000E+00  6.4488E-02 -1.0772E-01 -6.7921E-03
            -5.5104E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1249
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7981E-04 -1.3442E-04 -4.1306E-02 -6.4914E-03 -4.7791E-02
 SE:             2.9824E-02  6.1716E-05  1.9109E-02  2.9087E-02  1.9903E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9519E-01  2.9409E-02  3.0652E-02  8.2340E-01  1.6343E-02

 ETASHRINKSD(%)  8.7440E-02  9.9793E+01  3.5981E+01  2.5551E+00  3.3322E+01
 ETASHRINKVR(%)  1.7480E-01  1.0000E+02  5.9016E+01  5.0450E+00  5.5540E+01
 EBVSHRINKSD(%)  4.5421E-01  9.9798E+01  4.0448E+01  3.1099E+00  2.9201E+01
 EBVSHRINKVR(%)  9.0636E-01  1.0000E+02  6.4535E+01  6.1230E+00  4.9876E+01
 RELATIVEINF(%)  9.7866E+01  2.2730E-05  1.0849E+01  5.9560E+00  1.1297E+01
 EPSSHRINKSD(%)  4.5233E+01
 EPSSHRINKVR(%)  7.0005E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1594.6335877836254     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -859.48276121988727     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.92
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1594.634       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.26E-01  3.29E-01  1.66E+00  1.41E+00  9.89E-01  1.03E+00  1.00E-02  1.71E+00  7.75E-01  1.07E+00  1.07E+00
 


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
+        1.21E+03
 
 TH 2
+       -2.64E+01  4.53E+02
 
 TH 3
+        2.88E-01  2.90E+01  5.30E+01
 
 TH 4
+       -1.27E+01  5.62E+02 -1.35E+01  8.70E+02
 
 TH 5
+        1.21E+00 -1.77E+02 -1.19E+02 -5.83E+01  5.26E+02
 
 TH 6
+       -2.79E-01 -3.63E+00  1.80E-01 -3.36E+00 -6.31E-01  1.84E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -9.65E-02 -1.81E+00 -1.72E+01 -3.42E+00 -7.69E+00 -6.40E-02  0.00E+00  2.11E+01
 
 TH 9
+        4.10E+00 -1.05E+02  4.71E+00  1.32E-01  1.54E+00 -9.18E-01  0.00E+00  3.72E-02  2.93E+02
 
 TH10
+        1.16E+00  7.16E+00 -9.50E-01 -4.04E-01 -7.26E+01  5.91E-01  0.00E+00  9.23E+00  6.55E-01  5.52E+01
 
 TH11
+       -9.05E+00 -1.45E+01 -5.85E+00 -1.06E+01 -8.03E+00  2.57E+00  0.00E+00  3.80E+00  1.44E+01  1.07E+01  1.79E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.846
Stop Time:
Wed Sep 29 19:02:13 CDT 2021
