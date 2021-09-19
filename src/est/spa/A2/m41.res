Sat Sep 18 09:51:02 CDT 2021
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
$DATA ../../../../data/spa/A2/dat41.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1074.34403573879        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.6855E+00 -2.6873E+01 -3.3431E+00 -4.4191E+01  1.6683E+02  1.9411E+01 -2.7119E+01 -5.5681E+00 -5.4877E+01 -8.3389E+01
            -1.0460E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1429.16044207470        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0182E+00  9.1550E-01  9.7065E-01  1.0985E+00  8.4441E-01  8.9752E-01  9.8884E-01  9.3528E-01  1.0947E+00  9.6443E-01
             2.2130E+00
 PARAMETER:  1.1805E-01  1.1719E-02  7.0212E-02  1.9394E-01 -6.9117E-02 -8.1219E-03  8.8773E-02  3.3086E-02  1.9051E-01  6.3778E-02
             8.9436E-01
 GRADIENT:  -2.5915E+01  1.3576E+01  1.0825E+01  9.0162E-01  2.2683E+00 -1.7427E+01  6.1222E-01  4.5820E+00  8.7440E-01 -1.1359E+00
            -4.0692E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1439.14008429397        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0254E+00  6.8065E-01  3.7813E-01  1.1957E+00  4.3086E-01  9.5384E-01  7.1382E-01  3.1561E-01  1.0283E+00  3.8137E-01
             2.2840E+00
 PARAMETER:  1.2512E-01 -2.8470E-01 -8.7251E-01  2.7874E-01 -7.4197E-01  5.2739E-02 -2.3713E-01 -1.0533E+00  1.2794E-01 -8.6398E-01
             9.2593E-01
 GRADIENT:  -2.9477E+01  5.7907E+01  1.3992E+01  1.2532E+02 -1.4537E+01 -4.7438E-01 -1.8953E+01 -1.7919E+00 -8.3122E+00 -1.7136E+01
            -3.7998E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1456.50885580489        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0434E+00  5.4598E-01  2.3807E-01  1.1002E+00  3.0252E-01  9.6628E-01  1.0676E+00  1.0056E-01  1.0723E+00  4.1223E-01
             2.2374E+00
 PARAMETER:  1.4245E-01 -5.0517E-01 -1.3352E+00  1.9551E-01 -1.0956E+00  6.5702E-02  1.6544E-01 -2.1970E+00  1.6985E-01 -7.8618E-01
             9.0533E-01
 GRADIENT:   1.6068E+01  5.5556E+01  3.9088E+01  4.9228E+01 -8.2563E+01  8.4970E-01 -2.1188E+00 -6.6026E-02  1.0410E+01 -4.7652E+00
             1.9745E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1460.29638794119        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      346
 NPARAMETR:  1.0386E+00  4.4569E-01  2.2141E-01  1.0815E+00  2.7257E-01  9.6931E-01  1.0372E+00  7.1383E-02  1.0279E+00  5.5651E-01
             2.1148E+00
 PARAMETER:  1.3790E-01 -7.0813E-01 -1.4077E+00  1.7840E-01 -1.1999E+00  6.8825E-02  1.3651E-01 -2.5397E+00  1.2751E-01 -4.8608E-01
             8.4898E-01
 GRADIENT:  -4.2433E+00  1.4364E+01  6.8004E+00  7.8525E+00 -2.4701E+01  9.4082E-01  1.5903E+00  2.0391E-02 -6.3993E+00  4.8529E+00
             8.2531E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1460.72280168723        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  1.0404E+00  4.2395E-01  2.2558E-01  1.0880E+00  2.7254E-01  9.6614E-01  1.0821E+00  7.5038E-02  1.0449E+00  5.3752E-01
             2.0848E+00
 PARAMETER:  1.3961E-01 -7.5814E-01 -1.3891E+00  1.8434E-01 -1.2000E+00  6.5554E-02  1.7894E-01 -2.4898E+00  1.4393E-01 -5.2078E-01
             8.3467E-01
 GRADIENT:   2.5855E-02 -1.3234E-01 -2.0250E-01 -6.5413E-02  2.5204E-01  3.3633E-02  3.9724E-02  7.5871E-03  2.2075E-01  1.1909E-01
            -1.3545E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1460.72505355673        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      677
 NPARAMETR:  1.0405E+00  4.2490E-01  2.2667E-01  1.0890E+00  2.7342E-01  9.6591E-01  1.0852E+00  4.3323E-02  1.0428E+00  5.3498E-01
             2.0863E+00
 PARAMETER:  1.3967E-01 -7.5590E-01 -1.3843E+00  1.8528E-01 -1.1968E+00  6.5311E-02  1.8173E-01 -3.0391E+00  1.4192E-01 -5.2552E-01
             8.3539E-01
 GRADIENT:   1.0574E+01  1.5416E+00  3.7850E+00  3.2289E+00  1.3736E+01  5.7450E-01 -9.7744E-02  2.9052E-03  4.3246E-01 -1.7190E-01
             4.0982E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1460.72656667060        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      774
 NPARAMETR:  1.0404E+00  4.2467E-01  2.2659E-01  1.0890E+00  2.7332E-01  9.6587E-01  1.0825E+00  1.0000E-02  1.0431E+00  5.3724E-01
             2.0869E+00
 PARAMETER:  1.3962E-01 -7.5645E-01 -1.3846E+00  1.8529E-01 -1.1971E+00  6.5271E-02  1.7929E-01 -4.5635E+00  1.4223E-01 -5.2131E-01
             8.3567E-01
 GRADIENT:  -3.2942E-02 -3.2009E-02  1.2192E-01 -2.0495E-01 -6.7155E-02 -3.3242E-02 -8.8573E-02  0.0000E+00  7.5011E-03 -1.0056E-01
            -3.5417E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1460.72668977152        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  1.0404E+00  4.2463E-01  2.2659E-01  1.0892E+00  2.7332E-01  9.6598E-01  1.0850E+00  1.0000E-02  1.0430E+00  5.3745E-01
             2.0867E+00
 PARAMETER:  1.3964E-01 -7.5653E-01 -1.3846E+00  1.8541E-01 -1.1971E+00  6.5386E-02  1.8158E-01 -4.5407E+00  1.4214E-01 -5.2092E-01
             8.3559E-01
 GRADIENT:   1.3194E-02  2.4863E-02 -5.1036E-02 -1.0950E-02  3.5500E-02  9.0497E-03  7.2415E-03  0.0000E+00  4.9814E-03  1.9042E-02
            -2.7556E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      868
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.4388E-04  1.5433E-02 -1.4927E-04 -8.9416E-03  9.5522E-03
 SE:             2.9374E-02  1.5962E-02  2.1035E-04  2.7356E-02  2.0094E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9609E-01  3.3364E-01  4.7795E-01  7.4377E-01  6.3453E-01

 ETASHRINKSD(%)  1.5937E+00  4.6525E+01  9.9295E+01  8.3537E+00  3.2681E+01
 ETASHRINKVR(%)  3.1619E+00  7.1404E+01  9.9995E+01  1.6010E+01  5.4681E+01
 EBVSHRINKSD(%)  1.7045E+00  4.8046E+01  9.9290E+01  7.4597E+00  3.1966E+01
 EBVSHRINKVR(%)  3.3799E+00  7.3008E+01  9.9995E+01  1.4363E+01  5.3714E+01
 RELATIVEINF(%)  9.5287E+01  2.5228E+00  2.5145E-04  2.6241E+01  1.4079E+00
 EPSSHRINKSD(%)  3.8174E+01
 EPSSHRINKVR(%)  6.1775E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1460.7266897715240     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -725.57586320778580     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.43
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1460.727       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  4.25E-01  2.27E-01  1.09E+00  2.73E-01  9.66E-01  1.09E+00  1.00E-02  1.04E+00  5.37E-01  2.09E+00
 


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
+        1.06E+03
 
 TH 2
+       -1.59E+01  1.42E+03
 
 TH 3
+       -1.30E+02  3.00E+03  1.10E+04
 
 TH 4
+       -2.51E+01  1.83E+02 -1.01E+03  7.41E+02
 
 TH 5
+        2.04E+02 -5.36E+03 -1.48E+04  8.60E+01  2.42E+04
 
 TH 6
+        4.99E-01 -6.00E+00  3.06E+01 -8.83E+00  7.20E+00  2.00E+02
 
 TH 7
+        8.80E-01  3.96E+01 -1.36E+01 -1.60E+00 -5.97E+01 -3.66E-02  1.64E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.70E+00 -1.46E+01  1.18E+02 -7.47E+00  4.58E+01  1.36E+00  4.92E+00  0.00E+00  1.34E+02
 
 TH10
+       -3.99E+00  4.96E+01 -2.90E+02  1.38E+00  3.24E+02  2.55E+00  3.22E+01  0.00E+00 -1.85E+00  1.67E+02
 
 TH11
+       -1.30E+01 -3.02E+00 -8.43E+01 -6.85E+00  4.49E+01  1.66E+00  5.45E+00  0.00E+00  6.63E+00  2.04E+01  5.76E+01
 
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
 #CPUT: Total CPU Time in Seconds,       15.454
Stop Time:
Sat Sep 18 09:51:19 CDT 2021
