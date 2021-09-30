Wed Sep 29 06:58:20 CDT 2021
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
$DATA ../../../../data/int/TD2/dat7.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2286.01033076424        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9442E+02  1.2570E+01  2.4723E+02  1.7369E+02  2.0616E+02  1.0181E+01 -9.1359E+01 -1.0001E+03 -2.9213E+02 -6.1008E+01
            -1.6540E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3285.63378606545        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  8.9476E-01  1.2061E+00  9.9693E-01  8.6697E-01  1.1082E+00  1.0167E+00  1.0257E+00  1.2649E+00  1.0571E+00  1.2547E+00
             1.7296E+00
 PARAMETER: -1.1198E-02  2.8740E-01  9.6921E-02 -4.2755E-02  2.0277E-01  1.1652E-01  1.2542E-01  3.3503E-01  1.5553E-01  3.2692E-01
             6.4792E-01
 GRADIENT:  -4.5654E+01  8.2723E+01 -1.4565E+01  1.0421E+02  7.3029E+01 -1.5740E+01  2.1749E+01 -1.5870E+00  4.3601E+00  1.8854E+00
             1.7534E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3291.14096844590        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  8.9621E-01  1.2571E+00  8.6759E-01  8.2692E-01  1.1014E+00  1.0318E+00  7.4922E-01  1.0356E+00  1.1179E+00  1.3556E+00
             1.7243E+00
 PARAMETER: -9.5763E-03  3.2880E-01 -4.2040E-02 -9.0050E-02  1.9658E-01  1.3133E-01 -1.8872E-01  1.3495E-01  2.1145E-01  4.0421E-01
             6.4484E-01
 GRADIENT:  -3.5216E+01  9.5401E+01 -2.4684E+01  8.3222E+01  6.2359E+01 -2.7348E+00 -4.8811E-01 -1.5850E+00  4.7459E+00  7.8958E+00
             1.6030E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3294.59012808660        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  8.9816E-01  1.3297E+00  9.8918E-01  7.7746E-01  1.1891E+00  1.0332E+00  7.4405E-01  1.2477E+00  1.4228E+00  1.2976E+00
             1.7027E+00
 PARAMETER: -7.4107E-03  3.8492E-01  8.9117E-02 -1.5172E-01  2.7317E-01  1.3269E-01 -1.9565E-01  3.2131E-01  4.5265E-01  3.6050E-01
             6.3223E-01
 GRADIENT:  -1.9563E+02 -6.5240E+01 -2.3656E+00  4.0208E+01  3.1088E+01 -4.5834E+01 -2.9372E+00 -5.2915E+00  3.7964E+01 -2.2106E+01
             1.4869E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3307.57175657410        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      537
 NPARAMETR:  9.8663E-01  1.3925E+00  9.5259E-01  7.7738E-01  1.1888E+00  1.1058E+00  6.8473E-01  1.2474E+00  1.4233E+00  1.4480E+00
             1.7017E+00
 PARAMETER:  8.6544E-02  4.3107E-01  5.1428E-02 -1.5183E-01  2.7298E-01  2.0059E-01 -2.7873E-01  3.2109E-01  4.5295E-01  4.7019E-01
             6.3165E-01
 GRADIENT:  -2.4315E+00 -2.4082E+00 -1.8200E+00  7.4819E+01  3.2191E+00  1.0955E+00  6.9622E-01 -5.5180E+00  3.0321E+01  8.2722E-01
             1.4225E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3309.13349806304        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      718
 NPARAMETR:  9.8749E-01  1.4268E+00  9.3672E-01  7.7727E-01  1.1898E+00  1.1132E+00  7.2801E-01  1.2488E+00  1.4198E+00  1.4616E+00
             1.6718E+00
 PARAMETER:  8.7408E-02  4.5545E-01  3.4627E-02 -1.5196E-01  2.7376E-01  2.0725E-01 -2.1744E-01  3.2216E-01  4.5052E-01  4.7955E-01
             6.1391E-01
 GRADIENT:  -5.3436E-01  3.9416E+01 -4.4980E+00  9.2796E+01 -6.2721E+00  3.4867E+00  5.2683E+00 -5.8963E+00  2.5433E+01  4.0865E+00
             1.0733E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3317.02147631865        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  9.8713E-01  1.4259E+00  9.3814E-01  6.9757E-01  1.2230E+00  1.1012E+00  6.9364E-01  1.2482E+00  1.2392E+00  1.4702E+00
             1.6655E+00
 PARAMETER:  8.7045E-02  4.5478E-01  3.6147E-02 -2.6015E-01  3.0134E-01  1.9644E-01 -2.6581E-01  3.2166E-01  3.1450E-01  4.8541E-01
             6.1010E-01
 GRADIENT:  -5.5788E-01 -4.7454E+01 -1.1815E+00  2.1274E-01  2.9305E+00 -4.0359E-01  2.9741E-01 -7.8185E+00 -1.2970E-01  1.0187E+00
             9.6112E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3317.20566663052        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  9.8759E-01  1.4279E+00  9.4124E-01  6.9768E-01  1.2196E+00  1.1025E+00  6.9344E-01  1.2689E+00  1.2401E+00  1.4663E+00
             1.6657E+00
 PARAMETER:  8.7510E-02  4.5621E-01  3.9444E-02 -2.5999E-01  2.9850E-01  1.9756E-01 -2.6610E-01  3.3818E-01  3.1521E-01  4.8277E-01
             6.1026E-01
 GRADIENT:   2.7728E-01 -4.3121E+01 -8.5136E-01  1.0609E+00 -8.2392E-01  3.7417E-02  2.8659E-01 -7.4841E+00  2.5719E-02  2.3591E-01
             9.7180E+01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -3317.24176513259        NO. OF FUNC. EVALS.:  67
 CUMULATIVE NO. OF FUNC. EVALS.:     1146
 NPARAMETR:  9.8755E-01  1.4269E+00  9.4129E-01  6.9766E-01  1.2196E+00  1.1029E+00  6.9154E-01  1.2723E+00  1.2401E+00  1.4648E+00
             1.6657E+00
 PARAMETER:  8.7476E-02  4.5646E-01  3.9495E-02 -2.6003E-01  2.9850E-01  1.9752E-01 -2.6618E-01  3.4156E-01  3.1519E-01  4.8277E-01
             6.1026E-01
 GRADIENT:   1.9852E+00  2.3987E+04  4.2453E+01 -1.9845E+01 -1.8039E+00 -1.9657E-01  2.2877E-01  3.2050E+04 -1.3366E+01  2.2727E+04
            -1.8770E+02
 NUMSIGDIG:         7.3         2.3         6.0         5.9         6.9         2.3         1.6         2.3         6.0         2.3
                    4.6

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1146
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0760E-03 -2.4687E-02 -1.2888E-02  2.1970E-02 -1.9808E-02
 SE:             2.9811E-02  1.9252E-02  1.3800E-02  2.5316E-02  2.6851E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7121E-01  1.9972E-01  3.5036E-01  3.8548E-01  4.6069E-01

 ETASHRINKSD(%)  1.2806E-01  3.5504E+01  5.3767E+01  1.5190E+01  1.0046E+01
 ETASHRINKVR(%)  2.5595E-01  5.8403E+01  7.8625E+01  2.8072E+01  1.9082E+01
 EBVSHRINKSD(%)  5.6601E-01  3.5766E+01  5.9699E+01  1.7029E+01  8.6579E+00
 EBVSHRINKVR(%)  1.1288E+00  5.8739E+01  8.3758E+01  3.1158E+01  1.6566E+01
 RELATIVEINF(%)  9.8862E+01  1.0230E+01  1.3633E+01  1.9688E+01  3.6962E+01
 EPSSHRINKSD(%)  2.2130E+01
 EPSSHRINKVR(%)  3.9363E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3317.2417651325886     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1663.1524053641779     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.41
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3317.242       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  1.43E+00  9.41E-01  6.98E-01  1.22E+00  1.10E+00  6.93E-01  1.27E+00  1.24E+00  1.47E+00  1.67E+00
 


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
+        5.63E+07
 
 TH 2
+        9.31E+01  1.29E+06
 
 TH 3
+        5.90E+07  1.82E+03  6.19E+07
 
 TH 4
+       -5.84E+00 -2.32E+06 -6.25E+01  1.67E+07
 
 TH 5
+       -1.17E+00  2.31E+06 -6.50E+01 -4.15E+06  2.07E+06
 
 TH 6
+       -8.47E-02  1.10E+02  8.95E-01 -6.16E-01 -4.77E-01  1.60E+02
 
 TH 7
+        3.37E-01  1.99E+02 -1.16E+01  2.16E+00 -8.15E+00  5.09E-02  7.57E+01
 
 TH 8
+        6.40E+06 -1.00E+01  2.62E+03  3.48E+06 -4.40E+00  1.67E+02  2.67E+02  1.44E+06
 
 TH 9
+       -1.42E+07 -4.29E+02 -1.49E+07  3.75E+01  1.10E+00  2.12E+00  1.90E+01 -1.62E+06  3.59E+06
 
 TH10
+        9.07E+01 -1.27E+01  1.64E+03 -2.14E+06  2.13E+06  1.03E+02  1.81E+02 -8.93E+05 -3.79E+02  5.47E+05
 
 TH11
+       -8.30E+00  4.17E+05 -1.23E+00 -1.48E+06 -7.35E+05  1.86E+00  7.43E+00 -6.14E+05  1.58E+01  3.84E+05  2.70E+05
 
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
 #CPUT: Total CPU Time in Seconds,       44.959
Stop Time:
Wed Sep 29 06:59:07 CDT 2021
