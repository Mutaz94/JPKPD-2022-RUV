Wed Sep 29 20:25:48 CDT 2021
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
$DATA ../../../../data/spa/D/dat93.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m93.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26272.9281570890        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.9693E+02  6.7664E+02  1.2392E+01  6.3287E+02  1.8758E+02 -3.0351E+03 -1.0947E+03 -6.6352E+01 -1.7108E+03 -8.8958E+02
            -4.8432E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -441.556956297461        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1316E+00  1.0319E+00  8.6247E-01  1.3401E+00  1.5008E+00  1.7923E+00  9.0552E-01  9.6020E-01  6.7068E-01  9.2110E-01
             1.5387E+01
 PARAMETER:  2.2363E-01  1.3140E-01 -4.7951E-02  3.9274E-01  5.0598E-01  6.8347E-01  7.5650E-04  5.9387E-02 -2.9947E-01  1.7817E-02
             2.8335E+00
 GRADIENT:  -3.0687E+01  1.0988E+01 -2.5135E+00  1.1949E+01 -3.4632E+00  4.4813E+01  2.4300E+00  3.3198E+00  7.3348E+00  8.0605E-01
             7.8228E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -457.467871532024        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.1852E+00  9.1053E-01  9.8762E-01  1.2970E+00  4.0108E+00  1.3268E+00  3.4438E-01  3.8539E-01  2.7169E-01  5.8716E+00
             1.6294E+01
 PARAMETER:  2.6995E-01  6.2670E-03  8.7544E-02  3.6006E-01  1.4890E+00  3.8278E-01 -9.6601E-01 -8.5349E-01 -1.2031E+00  1.8701E+00
             2.8908E+00
 GRADIENT:   1.0981E+01  8.7647E+00 -3.5913E-01  3.5698E+00 -8.6942E+00 -1.4361E+01  3.5493E-01  1.3970E-01  2.0229E+00 -6.0809E-01
             2.5360E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -468.031097539532        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.1452E+00  5.6661E-01  9.1013E-01  1.4561E+00  8.0767E+00  1.3411E+00  5.9670E-01  1.8099E-01  1.6788E-01  1.0752E+01
             1.5727E+01
 PARAMETER:  2.3561E-01 -4.6809E-01  5.8331E-03  4.7576E-01  2.1890E+00  3.9349E-01 -4.1635E-01 -1.6093E+00 -1.6845E+00  2.4751E+00
             2.8554E+00
 GRADIENT:   5.5499E+00  9.9031E+00  6.6408E+00  1.1284E+01  1.6124E+00 -1.6970E+01  4.6367E-01  3.7527E-02  1.2682E+00  5.4071E+00
             1.0397E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -473.510526667293        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0444E+00  1.8985E-01  7.2632E-01  1.5553E+00  5.8778E+00  1.4868E+00  1.5036E+00  8.0722E-02  6.8084E-02  8.8414E+00
             1.4752E+01
 PARAMETER:  1.4341E-01 -1.5615E+00 -2.1976E-01  5.4168E-01  1.8712E+00  4.9662E-01  5.0787E-01 -2.4167E+00 -2.5870E+00  2.2794E+00
             2.7914E+00
 GRADIENT:  -2.7709E+01  4.7487E+00  5.1124E+00  1.2317E+01  1.6347E+00  7.2850E+00  3.9934E-01  1.3816E-02  3.2797E-01 -4.6904E+00
            -5.4974E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -477.647661393977        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  1.0560E+00  2.2712E-02  6.3831E-01  1.6022E+00  5.3852E+00  1.4839E+00  5.4771E-01  3.7760E-02  1.0000E-02  8.6636E+00
             1.4791E+01
 PARAMETER:  1.5453E-01 -3.6849E+00 -3.4893E-01  5.7139E-01  1.7837E+00  4.9468E-01 -5.0201E-01 -3.1765E+00 -4.5848E+00  2.2591E+00
             2.7940E+00
 GRADIENT:  -4.9645E+00  5.5346E-01 -5.7274E+00  1.6377E+01 -2.0361E+00  2.4488E+00  7.5627E-04  3.5996E-03  0.0000E+00 -9.7623E-01
            -1.7379E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -478.380066158336        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  1.0616E+00  1.0000E-02  6.5279E-01  1.6070E+00  5.6179E+00  1.4669E+00  4.0603E-01  2.5721E-02  1.0000E-02  8.8425E+00
             1.4987E+01
 PARAMETER:  1.5978E-01 -1.0762E+01 -3.2650E-01  5.7438E-01  1.8260E+00  4.8316E-01 -8.0133E-01 -3.5604E+00 -5.1618E+00  2.2796E+00
             2.8072E+00
 GRADIENT:  -5.2675E+00  0.0000E+00  1.2705E+00 -9.1490E+00  2.5623E+00  1.2970E+00  8.7408E-05  1.7616E-03  0.0000E+00 -3.6068E+00
             1.8311E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -478.769964921227        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      679
 NPARAMETR:  1.0620E+00  1.0000E-02  6.4593E-01  1.6184E+00  5.6093E+00  1.4681E+00  4.0466E-01  1.1835E-02  1.0000E-02  9.0520E+00
             1.5132E+01
 PARAMETER:  1.6016E-01 -1.0686E+01 -3.3706E-01  5.8143E-01  1.8244E+00  4.8397E-01 -8.0471E-01 -4.3367E+00 -5.1618E+00  2.3030E+00
             2.8168E+00
 GRADIENT:  -7.4607E+01  0.0000E+00 -1.1549E+01 -8.6116E+01  3.5144E+01  3.6235E+01  1.0736E-04  1.9356E-04  0.0000E+00 -1.0975E+02
             9.3481E+01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -478.769964921227        NO. OF FUNC. EVALS.:  39
 CUMULATIVE NO. OF FUNC. EVALS.:      718
 NPARAMETR:  1.0617E+00  1.0000E-02  6.4544E-01  1.6169E+00  5.6301E+00  1.4697E+00  4.0142E-01  1.1333E-02  1.0000E-02  9.0241E+00
             1.5173E+01
 PARAMETER:  1.6016E-01 -1.0686E+01 -3.3706E-01  5.8143E-01  1.8244E+00  4.8397E-01 -8.0471E-01 -4.3367E+00 -5.1618E+00  2.3030E+00
             2.8168E+00
 GRADIENT:   1.9990E+02  0.0000E+00  9.8440E+01  4.7778E+01 -2.0381E+01 -7.5776E+01  7.7489E-04  2.3532E-04  0.0000E+00  1.3414E+01
            -1.3822E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      718
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8990E-02 -4.7057E-05 -1.1578E-04 -2.9409E-04 -3.1027E-02
 SE:             2.7863E-02  1.8555E-05  5.7265E-05  1.9043E-04  7.5825E-03
 N:                     100         100         100         100         100

 P VAL.:         4.9552E-01  1.1208E-02  4.3190E-02  1.2252E-01  4.2808E-05

 ETASHRINKSD(%)  6.6547E+00  9.9938E+01  9.9808E+01  9.9362E+01  7.4598E+01
 ETASHRINKVR(%)  1.2867E+01  1.0000E+02  1.0000E+02  9.9996E+01  9.3547E+01
 EBVSHRINKSD(%)  7.6487E+00  9.9914E+01  9.9719E+01  9.9147E+01  6.1506E+01
 EBVSHRINKVR(%)  1.4712E+01  1.0000E+02  9.9999E+01  9.9993E+01  8.5182E+01
 RELATIVEINF(%)  2.9347E+01  1.9410E-06  9.3133E-05  1.1970E-04  6.4307E+00
 EPSSHRINKSD(%)  2.2120E+00
 EPSSHRINKVR(%)  4.3751E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -478.76996492122669     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       256.38086164251149     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.13
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -478.770       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.00E-02  6.46E-01  1.62E+00  5.61E+00  1.47E+00  4.05E-01  1.18E-02  1.00E-02  9.05E+00  1.51E+01
 


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
+        3.48E+04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.13E+02  0.00E+00  1.94E+04
 
 TH 4
+        2.29E+02  0.00E+00 -3.08E+02  1.67E+03
 
 TH 5
+       -2.58E+01  0.00E+00  1.46E+01 -2.31E+01  1.29E+01
 
 TH 6
+       -5.25E+02  0.00E+00  1.04E+02 -2.08E+02  1.21E+01  2.16E+03
 
 TH 7
+       -6.33E-02  0.00E+00 -9.59E-02 -3.43E-02  6.23E-04 -4.48E-03 -1.01E-01
 
 TH 8
+       -1.37E+00  0.00E+00 -9.70E-01 -1.93E-01  5.76E-03 -2.14E-01  7.80E-01  9.51E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        7.42E+00  0.00E+00 -6.53E+00  7.12E+00 -1.19E+00 -3.61E+00  2.95E-04 -7.23E-04  0.00E+00  4.20E+00
 
 TH11
+       -1.45E+01  0.00E+00  8.37E+00 -1.76E+01  9.98E-01  3.80E+00  8.55E-04  6.79E-03  0.00E+00 -1.09E+00  3.56E+00
 
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
 #CPUT: Total CPU Time in Seconds,       19.091
Stop Time:
Wed Sep 29 20:26:09 CDT 2021
