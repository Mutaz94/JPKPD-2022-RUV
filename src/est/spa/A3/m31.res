Wed Sep 29 13:26:44 CDT 2021
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
$DATA ../../../../data/spa/A3/dat31.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -93.8382798347705        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4395E+02  1.3136E+02  2.1185E+01  1.5914E+02  3.6358E+02  2.2155E+01 -9.7111E+01 -2.4562E+01 -1.7624E+02 -2.4805E+02
            -2.4809E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1143.83858296630        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.4648E-01  7.7544E-01  8.8607E-01  1.1160E+00  7.6242E-01  9.2016E-01  1.0120E+00  9.3375E-01  1.4737E+00  9.2774E-01
             2.3700E+00
 PARAMETER:  4.4998E-02 -1.5432E-01 -2.0962E-02  2.0979E-01 -1.7126E-01  1.6793E-02  1.1192E-01  3.1450E-02  4.8776E-01  2.4995E-02
             9.6288E-01
 GRADIENT:   2.9998E+01  1.6950E+01 -5.7460E+00  6.0275E+01  1.3020E+02 -2.3782E+01 -2.4464E+00  4.0054E+00  1.7428E+01 -5.2190E+01
            -2.9274E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1165.63959918399        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.7181E-01  4.1257E-01  4.4171E-01  1.4533E+00  3.5873E-01  9.7504E-01  6.7422E-01  6.5691E-01  1.5141E+00  5.0883E-01
             2.2733E+00
 PARAMETER:  7.1400E-02 -7.8534E-01 -7.1710E-01  4.7387E-01 -9.2518E-01  7.4719E-02 -2.9420E-01 -3.2021E-01  5.1484E-01 -5.7565E-01
             9.2125E-01
 GRADIENT:   6.6442E+01  1.0722E+02  1.8197E+02  2.8393E+02 -1.3013E+02 -7.8676E+00 -2.7168E+00 -1.3851E+01  3.5389E+01 -2.3817E+01
            -2.8572E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1226.40453670653        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      318
 NPARAMETR:  9.9236E-01  4.0986E-01  5.3567E-01  1.4044E+00  4.5489E-01  9.9554E-01  1.9244E-01  3.1253E-01  1.4243E+00  5.6378E-01
             3.2330E+00
 PARAMETER:  9.2334E-02 -7.9194E-01 -5.2423E-01  4.3962E-01 -6.8770E-01  9.5530E-02 -1.5480E+00 -1.0631E+00  4.5367E-01 -4.7309E-01
             1.2734E+00
 GRADIENT:   2.3128E+01  9.4661E+00 -1.5856E+01  1.1276E+02  6.9517E+01  2.9369E+00  7.5460E-02  8.5813E-01  2.3926E+01 -5.3992E+00
            -4.7878E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1257.17659989244        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      493
 NPARAMETR:  9.9042E-01  2.8252E-01  2.6395E-01  1.1328E+00  2.6226E-01  1.0012E+00  4.8293E-01  4.9458E-02  1.3379E+00  5.1265E-01
             3.1969E+00
 PARAMETER:  9.0372E-02 -1.1640E+00 -1.2320E+00  2.2472E-01 -1.2384E+00  1.0119E-01 -6.2789E-01 -2.9066E+00  3.9110E-01 -5.6817E-01
             1.2622E+00
 GRADIENT:   3.8216E+00 -1.9041E+00 -5.5317E-01  7.8494E+00  3.0691E+00 -3.3639E-01  5.7971E-02 -9.0753E-03  2.3148E+00  3.4733E-01
            -4.5241E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1258.18046308547        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  9.9198E-01  3.7901E-01  2.1777E-01  1.0474E+00  2.4989E-01  1.0120E+00  6.0531E-01  2.8448E-02  1.4497E+00  5.1619E-01
             3.1620E+00
 PARAMETER:  9.1947E-02 -8.7018E-01 -1.4243E+00  1.4631E-01 -1.2867E+00  1.1192E-01 -4.0202E-01 -3.4597E+00  4.7137E-01 -5.6129E-01
             1.2512E+00
 GRADIENT:   7.1142E+00  7.9853E+00  1.7070E+01  8.1760E+00 -2.8311E+01 -3.4150E-01  6.3868E-01 -4.1624E-03  5.2579E+00  2.5464E+00
             6.5139E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1259.52052248560        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      846
 NPARAMETR:  9.8103E-01  4.4986E-01  1.7427E-01  9.6264E-01  2.4350E-01  1.0251E+00  5.2896E-01  1.0000E-02  1.5322E+00  4.7353E-01
             3.0885E+00
 PARAMETER:  8.0849E-02 -6.9881E-01 -1.6471E+00  6.1926E-02 -1.3126E+00  1.2476E-01 -5.3683E-01 -4.8011E+00  5.2672E-01 -6.4753E-01
             1.2277E+00
 GRADIENT:   2.9240E+00 -9.2394E-01 -1.2463E+00  3.1000E+00  1.9885E+00 -4.5829E-01  6.2426E-02  0.0000E+00 -1.8891E-01 -1.1131E+00
            -5.3565E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1259.58266628028        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1021
 NPARAMETR:  9.7929E-01  4.4700E-01  1.6917E-01  9.5052E-01  2.3943E-01  1.0269E+00  4.3060E-01  1.0000E-02  1.5536E+00  5.0505E-01
             3.0763E+00
 PARAMETER:  7.9070E-02 -7.0520E-01 -1.6769E+00  4.9258E-02 -1.3295E+00  1.2658E-01 -7.4257E-01 -4.9022E+00  5.4059E-01 -5.8309E-01
             1.2237E+00
 GRADIENT:   6.0127E-01  1.7167E-02  7.0231E-03 -2.0303E-01 -6.5809E-02  3.7213E-02  5.9104E-02  0.0000E+00  1.2144E-01  1.2810E-01
             5.3955E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1259.58377241746        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1113
 NPARAMETR:  9.7900E-01  4.4616E-01  1.6910E-01  9.5072E-01  2.3919E-01  1.0268E+00  4.1188E-01  1.0000E-02  1.5540E+00  5.0752E-01
             3.0740E+00
 PARAMETER:  7.8777E-02 -7.0708E-01 -1.6773E+00  4.9459E-02 -1.3305E+00  1.2646E-01 -7.8701E-01 -4.8935E+00  5.4085E-01 -5.7823E-01
             1.2230E+00
 GRADIENT:  -6.5701E-03  7.0574E-02  8.1760E-02  2.0286E-02 -1.8144E-01  6.5608E-04  5.6052E-03  0.0000E+00 -3.4299E-03  1.2096E-02
            -2.6952E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1113
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.8311E-04 -2.0441E-03  1.1593E-04 -9.8202E-03  8.3343E-03
 SE:             2.8780E-02  7.3345E-03  1.8267E-04  2.6643E-02  1.9174E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8106E-01  7.8047E-01  5.2566E-01  7.1244E-01  6.6381E-01

 ETASHRINKSD(%)  3.5818E+00  7.5428E+01  9.9388E+01  1.0742E+01  3.5763E+01
 ETASHRINKVR(%)  7.0354E+00  9.3962E+01  9.9996E+01  2.0331E+01  5.8736E+01
 EBVSHRINKSD(%)  3.2379E+00  7.5601E+01  9.9432E+01  8.0213E+00  3.6025E+01
 EBVSHRINKVR(%)  6.3710E+00  9.4047E+01  9.9997E+01  1.5399E+01  5.9073E+01
 RELATIVEINF(%)  8.8789E+01  2.9034E-01  2.6295E-04  4.1730E+01  1.0259E+00
 EPSSHRINKSD(%)  3.0768E+01
 EPSSHRINKVR(%)  5.2070E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1259.5837724174594     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -524.43294585372121     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1259.584       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  4.46E-01  1.69E-01  9.51E-01  2.39E-01  1.03E+00  4.12E-01  1.00E-02  1.55E+00  5.08E-01  3.07E+00
 


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
+        1.03E+03
 
 TH 2
+       -3.78E+01  1.95E+03
 
 TH 3
+       -3.28E+02  4.08E+03  1.24E+04
 
 TH 4
+       -2.12E+01  1.29E+02 -5.92E+02  4.33E+02
 
 TH 5
+        4.04E+02 -7.32E+03 -1.74E+04 -1.87E+02  3.10E+04
 
 TH 6
+        1.01E+00 -6.07E+00  6.08E+01 -1.06E+01  7.16E+00  1.64E+02
 
 TH 7
+       -4.67E-01 -1.51E+01 -3.31E+01 -1.73E+00  4.70E+01  1.61E-01  5.11E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.97E+00 -2.74E+01  1.18E+02 -4.68E+00  8.11E+01 -2.08E-01  2.40E+00  0.00E+00  5.46E+01
 
 TH10
+       -7.61E+00 -5.60E+01 -2.03E+02  2.08E+00  4.36E+02  6.14E+00  1.79E+01  0.00E+00  1.40E+00  1.33E+02
 
 TH11
+       -1.67E+01 -6.10E+00 -3.65E+01 -4.12E+00  1.96E+01  2.34E+00  4.34E+00  0.00E+00  4.35E+00  2.40E+01  3.24E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.987
Stop Time:
Wed Sep 29 13:27:05 CDT 2021
