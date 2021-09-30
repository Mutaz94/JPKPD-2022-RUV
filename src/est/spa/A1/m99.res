Wed Sep 29 12:29:59 CDT 2021
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
$DATA ../../../../data/spa/A1/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1471.51789717135        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2282E+02 -5.7539E+00 -2.0523E+01  4.9440E+01  9.9747E+01  7.4256E+01 -1.9494E+01  1.5658E+00 -3.4015E+00 -2.2856E+01
            -3.5785E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1523.83728148917        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0087E+00  1.0376E+00  1.0454E+00  9.7447E-01  9.5303E-01  8.1262E-01  1.1021E+00  9.3334E-01  9.7070E-01  9.1994E-01
             2.4144E+00
 PARAMETER:  1.0867E-01  1.3692E-01  1.4442E-01  7.4143E-02  5.1893E-02 -1.0750E-01  1.9720E-01  3.1017E-02  7.0261E-02  1.6553E-02
             9.8144E-01
 GRADIENT:   6.5778E+01 -2.2016E+01 -5.0361E+00 -3.4811E+01 -1.9136E+01 -2.3791E+01  4.3213E+00  5.2775E+00  7.9323E+00  1.8423E+01
             1.4551E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1541.15193128017        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.9939E-01  8.3069E-01  6.2015E-01  1.1200E+00  6.7839E-01  8.5189E-01  1.4697E+00  2.5892E-01  8.0855E-01  7.2328E-01
             2.0359E+00
 PARAMETER:  9.9393E-02 -8.5493E-02 -3.7780E-01  2.1334E-01 -2.8803E-01 -6.0294E-02  4.8503E-01 -1.2512E+00 -1.1251E-01 -2.2396E-01
             8.1095E-01
 GRADIENT:   4.6146E+01  9.6340E+00 -5.0907E+01  1.1421E+02  5.9844E+01 -1.2786E+01  9.6426E+00  1.1674E+00  9.3667E-01  1.2716E+01
             1.0004E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1554.25277680247        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      272
 NPARAMETR:  9.9868E-01  9.0044E-01  7.2184E-01  1.0598E+00  7.5181E-01  8.9461E-01  1.2436E+00  1.9834E-01  8.8535E-01  7.7097E-01
             1.7745E+00
 PARAMETER:  9.8682E-02 -4.8686E-03 -2.2595E-01  1.5811E-01 -1.8527E-01 -1.1368E-02  3.1805E-01 -1.5178E+00 -2.1774E-02 -1.6011E-01
             6.7354E-01
 GRADIENT:  -3.4563E+01 -5.3608E+00 -7.3324E+00 -9.5698E-01  2.4579E+00 -3.6827E+00 -8.2739E+00  3.2907E-01 -2.1653E+00  2.5207E+00
             3.7998E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1557.59676490216        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  1.0063E+00  6.7356E-01  8.4790E-01  1.2146E+00  7.3636E-01  9.0593E-01  1.6684E+00  9.1001E-02  8.4442E-01  8.6219E-01
             1.6040E+00
 PARAMETER:  1.0626E-01 -2.9518E-01 -6.4987E-02  2.9442E-01 -2.0603E-01  1.2066E-03  6.1189E-01 -2.2969E+00 -6.9107E-02 -4.8284E-02
             5.7251E-01
 GRADIENT:  -1.9136E+00  7.5740E+00  1.7095E+00  1.2455E+01 -2.0687E+00  1.7943E+00  1.2574E+00  3.8413E-02 -2.0579E-02 -8.7422E-01
            -2.3377E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1558.07840127014        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  1.0031E+00  4.1684E-01  9.1087E-01  1.3685E+00  6.8844E-01  8.9308E-01  2.1009E+00  1.9796E-02  8.1029E-01  9.2507E-01
             1.6029E+00
 PARAMETER:  1.0310E-01 -7.7505E-01  6.6430E-03  4.1375E-01 -2.7333E-01 -1.3078E-02  8.4235E-01 -3.8223E+00 -1.1036E-01  2.2119E-02
             5.7182E-01
 GRADIENT:  -1.1927E-01  4.4260E+00  4.2247E+00  1.2741E+01 -8.3467E+00 -1.8382E+00 -6.2659E-01  1.7174E-03 -1.2467E+00 -3.8827E-01
            -6.1057E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1558.28961474785        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  9.9974E-01  2.6438E-01  9.1833E-01  1.4486E+00  6.5455E-01  8.9458E-01  2.6595E+00  1.0000E-02  8.0045E-01  9.6232E-01
             1.5941E+00
 PARAMETER:  9.9741E-02 -1.2304E+00  1.4806E-02  4.7059E-01 -3.2381E-01 -1.1406E-02  1.0781E+00 -5.6556E+00 -1.2258E-01  6.1596E-02
             5.6631E-01
 GRADIENT:  -8.8694E-01  1.7318E+00  4.1377E+00  3.7210E-01 -8.7550E+00 -1.1939E-01  7.9453E-01  0.0000E+00  1.4809E+00  9.2459E-01
             1.4431E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1558.32843720969        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      983
 NPARAMETR:  9.9912E-01  2.2480E-01  9.3893E-01  1.4707E+00  6.6049E-01  8.9403E-01  2.8305E+00  1.0000E-02  7.9205E-01  9.7477E-01
             1.5916E+00
 PARAMETER:  9.9119E-02 -1.3925E+00  3.6987E-02  4.8576E-01 -3.1478E-01 -1.2020E-02  1.1404E+00 -6.2831E+00 -1.3313E-01  7.4449E-02
             5.6475E-01
 GRADIENT:   2.8369E-01 -4.3511E-01 -2.7131E+00 -4.5047E+00  3.5179E+00  1.2386E-02 -1.7870E-01  0.0000E+00  2.5221E-01 -2.8078E-01
             4.6645E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1558.33348936294        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1040
 NPARAMETR:  9.9909E-01  2.2597E-01  9.4081E-01  1.4714E+00  6.5964E-01  8.9402E-01  2.8346E+00  1.0000E-02  7.9166E-01  9.7572E-01
             1.5912E+00
 PARAMETER:  9.9090E-02 -1.3874E+00  3.8987E-02  4.8620E-01 -3.1606E-01 -1.2030E-02  1.1419E+00 -6.2831E+00 -1.3363E-01  7.5417E-02
             5.6448E-01
 GRADIENT:   1.3594E-01  2.4912E-01  3.7799E-01 -2.5785E+00 -1.2902E+00  4.0115E-03  7.9478E-02  0.0000E+00  5.4587E-02 -1.5382E-02
             1.6241E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1040
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8626E-04  1.8089E-02 -1.6974E-04 -1.8523E-02 -1.3158E-02
 SE:             2.9536E-02  1.3732E-02  1.6312E-04  2.6356E-02  2.3147E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9497E-01  1.8775E-01  2.9807E-01  4.8218E-01  5.6972E-01

 ETASHRINKSD(%)  1.0495E+00  5.3997E+01  9.9454E+01  1.1704E+01  2.2453E+01
 ETASHRINKVR(%)  2.0881E+00  7.8837E+01  9.9997E+01  2.2038E+01  3.9865E+01
 EBVSHRINKSD(%)  1.2516E+00  6.2755E+01  9.9416E+01  9.8379E+00  1.8421E+01
 EBVSHRINKVR(%)  2.4876E+00  8.6128E+01  9.9997E+01  1.8708E+01  3.3449E+01
 RELATIVEINF(%)  9.6054E+01  1.8406E+00  2.5435E-04  1.3569E+01  4.6280E+00
 EPSSHRINKSD(%)  3.9170E+01
 EPSSHRINKVR(%)  6.2997E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1558.3334893629360     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -823.18266279919783     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1558.333       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  2.26E-01  9.41E-01  1.47E+00  6.60E-01  8.94E-01  2.83E+00  1.00E-02  7.92E-01  9.76E-01  1.59E+00
 


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
+        1.36E+03
 
 TH 2
+       -4.00E+01  4.59E+02
 
 TH 3
+        1.41E+01  1.39E+02  4.37E+02
 
 TH 4
+       -2.36E+01  3.52E+02 -9.34E+01  6.76E+02
 
 TH 5
+        2.35E+00 -4.16E+02 -8.37E+02  1.78E+00  1.79E+03
 
 TH 6
+        3.13E-01 -6.27E+00  3.54E+00 -6.06E+00 -2.00E+00  2.36E+02
 
 TH 7
+        5.35E-01  2.45E+01  1.69E-01 -2.69E+00 -5.10E-01  7.38E-02  4.28E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.27E+00 -2.95E+01  3.32E+00 -8.16E+00  5.66E-01 -1.05E+00  9.22E-01  0.00E+00  2.31E+02
 
 TH10
+       -1.87E+00  2.41E+01 -5.81E+00 -1.11E+01 -5.47E+01  8.57E-01  4.78E-01  0.00E+00 -1.30E+00  8.66E+01
 
 TH11
+       -1.31E+01 -6.69E+00 -2.74E+01 -1.06E+01  1.33E+01  3.30E+00  9.34E-01  0.00E+00  8.06E+00  2.27E+01  9.50E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.767
Stop Time:
Wed Sep 29 12:30:20 CDT 2021
