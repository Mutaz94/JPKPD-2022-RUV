Sat Sep 18 15:58:01 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/All/dat72.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK

ET1 = EXP(ETA(1)*THETA(4))
ET2 = EXP(ETA(2)*THETA(5))
ET3 = EXP(ETA(3)*THETA(6))


CL = 5.0 * THETA(1) * ET1
V = 85  * THETA(2) * ET2
KA = 0.7 * THETA(3) * ET3

SC = V
$ERROR
CVERR 	= 0.05
W  	= THETA(7)*F*CVERR
Y  	= F + W * ERR(1)
$THETA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvKA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvK
(0,1) ; RUV
$OMEGA
0.9 FIX ;     IIV CL
0.9 FIX  ;     IIV V
0.9 FIX ;      IIV KA
$SIGMA  1  FIX;        [P]
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
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
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
0INITIAL ESTIMATE OF OMEGA:
 0.9000E+00
 0.0000E+00   0.9000E+00
 0.0000E+00   0.0000E+00   0.9000E+00
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

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
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
 RAW OUTPUT FILE (FILE): m72.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   10714.1727323320        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:  -6.5893E+02  6.9113E+02 -5.6326E+02  3.6212E+03 -4.3862E+03  2.5036E+02 -2.3320E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -603.560423343686        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:       60
 NPARAMETR:  1.5615E+00  1.5634E+00  3.3490E+00  8.0704E-01  7.5716E-01  5.3832E-01  1.1053E+01
 PARAMETER:  5.4566E-01  5.4688E-01  1.3086E+00 -1.1439E-01 -1.7819E-01 -5.1930E-01  2.5027E+00
 GRADIENT:   2.8117E+01 -1.0142E+01 -1.9971E+00  5.3042E+00  1.2929E+01  2.5649E-01  5.8215E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -608.429725169838        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      114
 NPARAMETR:  1.3688E+00  1.4948E+00  8.6443E+00  7.4260E-01  6.2203E-01  2.7853E+00  1.0495E+01
 PARAMETER:  4.1393E-01  5.0197E-01  2.2569E+00 -1.9760E-01 -3.7476E-01  1.1243E+00  2.4509E+00
 GRADIENT:   3.7001E+00  4.2310E+00 -1.5637E-01  3.4737E+00 -4.7448E+00 -2.2527E-01  1.1929E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -608.822984232440        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  1.3546E+00  1.5141E+00  1.3066E+01  7.7808E-01  7.0910E-01  3.9332E+00  1.0160E+01
 PARAMETER:  4.0354E-01  5.1482E-01  2.6700E+00 -1.5093E-01 -2.4375E-01  1.4694E+00  2.4185E+00
 GRADIENT:  -1.1766E+00  4.6691E-01 -1.7378E-01  7.4564E-01  1.4155E+00  7.7171E-02  3.9413E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -608.825313665737        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      216
 NPARAMETR:  1.3549E+00  1.5060E+00  1.3686E+01  7.7403E-01  7.0229E-01  3.9680E+00  1.0116E+01
 PARAMETER:  4.0369E-01  5.0943E-01  2.7164E+00 -1.5614E-01 -2.5341E-01  1.4783E+00  2.4141E+00
 GRADIENT:  -1.0853E-01 -3.7687E-02 -1.6730E-01  4.1481E-02  1.5836E-01  6.3587E-02 -6.5790E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -608.831290913230        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      270
 NPARAMETR:  1.3564E+00  1.4979E+00  1.7572E+01  7.7070E-01  6.9734E-01  4.4087E+00  1.0089E+01
 PARAMETER:  4.0482E-01  5.0404E-01  2.9663E+00 -1.6046E-01 -2.6048E-01  1.5836E+00  2.4114E+00
 GRADIENT:   1.3938E+00 -6.2211E-01 -1.3121E-01 -5.4279E-01 -6.4273E-01  3.9052E-02 -3.9762E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -608.952632067332        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      324
 NPARAMETR:  1.3452E+00  1.4956E+00  6.0995E+02  7.6857E-01  6.9364E-01  2.4841E+01  1.0188E+01
 PARAMETER:  3.9656E-01  5.0253E-01  6.5134E+00 -1.6322E-01 -2.6580E-01  3.3125E+00  2.4212E+00
 GRADIENT:  -6.7619E-01  2.9883E-01 -2.6981E-03  2.8816E-01  3.8883E-01 -1.1936E-03  1.6331E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -608.980778631504        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.3622E+00  1.5080E+00  2.7793E+03  7.7408E-01  7.0138E-01  5.1228E+01  1.0181E+01
 PARAMETER:  4.0910E-01  5.1081E-01  8.0299E+00 -1.5608E-01 -2.5471E-01  4.0363E+00  2.4205E+00
 GRADIENT:   2.5170E-01 -1.1734E-01 -4.3867E-04 -2.9682E-02  1.5090E-02 -3.7732E-04 -1.3375E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -608.981520472526        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      558
 NPARAMETR:  1.3612E+00  1.5082E+00  2.8280E+05  7.7386E-01  7.0088E-01  4.7295E+02  1.0184E+01
 PARAMETER:  4.0834E-01  5.1094E-01  1.2652E+01 -1.5637E-01 -2.5541E-01  6.2590E+00  2.4209E+00
 GRADIENT:  -2.3748E-04  2.7850E-03 -4.3291E-06 -5.3052E-03 -4.6410E-03 -4.5460E-06 -1.1014E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -608.981527705717        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      674
 NPARAMETR:  1.3612E+00  1.5083E+00  2.6563E+07  7.7392E-01  7.0097E-01  4.1955E+03  1.0184E+01
 PARAMETER:  4.0835E-01  5.1096E-01  1.7195E+01 -1.5629E-01 -2.5529E-01  8.4418E+00  2.4208E+00
 GRADIENT:   7.7679E-04 -2.5339E-04 -5.0546E-08  1.7540E-04  1.4126E-04 -4.8118E-08  9.9653E-05

0ITERATION NO.:   47    OBJECTIVE VALUE:  -608.981527718831        NO. OF FUNC. EVALS.:  41
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  1.3612E+00  1.5083E+00  3.0558E+07  7.7392E-01  7.0097E-01  4.4876E+03  1.0184E+01
 PARAMETER:  4.0835E-01  5.1096E-01  1.7335E+01 -1.5629E-01 -2.5529E-01  8.5091E+00  2.4208E+00
 GRADIENT:  -3.6309E-03  2.3766E-03 -1.2179E-06 -9.3754E-04 -1.2988E-05  2.0224E-06  2.2634E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      715
 NO. OF SIG. DIGITS IN FINAL EST.:  4.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.8782E-03 -7.9410E-02 -1.8862E-06
 SE:             8.8787E-02  8.2703E-02  2.6930E-06
 N:                     100         100         100

 P VAL.:         9.5618E-01  3.3696E-01  4.8368E-01

 ETASHRINKSD(%)  5.9392E+00  1.2385E+01  9.9997E+01
 ETASHRINKVR(%)  1.1526E+01  2.3236E+01  1.0000E+02
 EBVSHRINKSD(%)  3.9816E+00  1.0061E+01  9.9998E+01
 EBVSHRINKVR(%)  7.8047E+00  1.9111E+01  1.0000E+02
 RELATIVEINF(%)  1.9388E+01  1.6093E+01  6.3727E-09
 EPSSHRINKSD(%)  1.4222E+01
 EPSSHRINKVR(%)  2.6422E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -608.98152771883088     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       126.16929884490730     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     6.47
 Elapsed covariance  time in seconds:     2.22
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -608.982       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.36E+00  1.51E+00  3.06E+07  7.74E-01  7.01E-01  4.49E+03  1.02E+01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-01
 
 ETA2
+        0.00E+00  9.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.49E-01
 
 ETA2
+        0.00E+00  9.49E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.49E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.69E-01  3.42E-01  9.19E+08  3.18E-01  5.25E-01  6.12E+04  3.09E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        2.86E-02
 
 TH 2
+        4.88E-02  1.17E-01
 
 TH 3
+        9.28E+07  2.82E+08  8.45E+17
 
 TH 4
+        3.31E-02  9.99E-02  2.91E+08  1.01E-01
 
 TH 5
+        5.32E-02  1.64E-01  4.70E+08  1.66E-01  2.75E-01
 
 TH 6
+        6.30E+03  1.90E+04  5.62E+13  1.94E+04  3.14E+04  3.74E+09
 
 TH 7
+       -1.57E-01 -7.32E-01 -2.49E+09 -8.62E-01 -1.45E+00 -1.67E+05  9.54E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        1.69E-01
 
 TH 2
+        8.45E-01  3.42E-01
 
 TH 3
+        5.97E-01  8.98E-01  9.19E+08
 
 TH 4
+        6.14E-01  9.20E-01  9.94E-01  3.18E-01
 
 TH 5
+        5.99E-01  9.17E-01  9.74E-01  9.92E-01  5.25E-01
 
 TH 6
+        6.09E-01  9.07E-01  9.99E-01  9.95E-01  9.77E-01  6.12E+04
 
 TH 7
+       -3.00E-01 -6.93E-01 -8.77E-01 -8.77E-01 -8.92E-01 -8.82E-01  3.09E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        7.42E+03
 
 TH 2
+       -2.32E+04  2.90E+05
 
 TH 3
+       -1.71E-04  8.81E-04  4.65E-12
 
 TH 4
+        5.27E+05  1.22E+05 -1.01E-02  5.52E+07
 
 TH 5
+       -1.61E+05 -1.42E+05  2.93E-03 -1.78E+07  5.79E+06
 
 TH 6
+        1.22E+00 -1.51E+01 -4.64E-08 -4.36E+00  6.78E+00  7.91E-04
 
 TH 7
+       -1.94E+03 -2.32E+04  2.81E-06 -3.98E+05  1.36E+05  1.20E+00  4.60E+03
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,        8.767
Stop Time:
Sat Sep 18 15:58:11 CDT 2021
