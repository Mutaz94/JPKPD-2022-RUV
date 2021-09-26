Sat Sep 25 15:04:14 CDT 2021
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
$DATA ../../../../data/spa/All/dat48.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m48.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   9117.69425227132        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   7.6672E+01  1.6487E+01 -1.1299E+02  4.9493E-01 -1.0690E+02 -4.7834E+01 -2.0181E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:   1817.83137361545        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:       60
 NPARAMETR:  1.0712E+00  6.3920E-01  3.1953E+00  6.6486E-01  9.8099E-01  1.3505E+00  1.9301E+00
 PARAMETER:  1.6877E-01 -3.4754E-01  1.2617E+00 -3.0818E-01  8.0809E-02  4.0049E-01  7.5758E-01
 GRADIENT:   7.7953E+01 -1.1101E+02  4.1007E+01 -1.7034E+02 -1.7632E+02  5.3564E+01 -5.3748E+03

0ITERATION NO.:   10    OBJECTIVE VALUE:  -504.723089180155        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      114
 NPARAMETR:  7.4997E-01  1.2955E+00  1.8538E+00  7.0632E-01  1.2506E+00  2.2570E+00  6.6607E+00
 PARAMETER: -1.8772E-01  3.5890E-01  7.1723E-01 -2.4769E-01  3.2360E-01  9.1403E-01  1.9962E+00
 GRADIENT:  -1.3991E+02 -1.0310E+01 -3.9350E+01 -1.4945E+00  6.1169E+01  7.3816E+01 -2.2592E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -630.157393129918        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      166
 NPARAMETR:  1.2295E+00  1.3767E+00  4.7566E+00  5.5909E-01  6.6743E-01  1.3449E+00  9.9022E+00
 PARAMETER:  3.0660E-01  4.1969E-01  1.6595E+00 -4.8145E-01 -3.0432E-01  3.9630E-01  2.3928E+00
 GRADIENT:   5.2036E+00 -2.7720E+00 -4.7418E-01 -1.4462E+00  3.4707E+01 -9.7131E-03  9.1041E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -648.654424010378        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      217
 NPARAMETR:  1.1546E+00  1.1672E+00  1.0389E+01  4.0938E-01  2.3543E-01  1.3542E+00  1.1347E+01
 PARAMETER:  2.4373E-01  2.5463E-01  2.4407E+00 -7.9312E-01 -1.3464E+00  4.0324E-01  2.5289E+00
 GRADIENT:   2.8525E+00  5.3063E-01 -4.4849E-01  2.3201E+00  1.0004E+00 -1.2000E-01  6.5869E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -648.700477251314        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  1.1409E+00  1.1484E+00  1.1154E+01  3.9983E-01  2.1650E-01  1.3751E+00  1.1311E+01
 PARAMETER:  2.3183E-01  2.3833E-01  2.5118E+00 -8.1671E-01 -1.4302E+00  4.1856E-01  2.5258E+00
 GRADIENT:   4.8275E-02 -2.6825E-01 -4.7153E-01  2.9881E-02 -2.7901E-01 -1.1442E-01  3.5571E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -648.757672966012        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      324
 NPARAMETR:  1.1328E+00  1.1385E+00  1.8270E+01  3.9612E-01  1.9883E-01  1.7330E+00  1.1296E+01
 PARAMETER:  2.2470E-01  2.2968E-01  3.0052E+00 -8.2603E-01 -1.5153E+00  6.4987E-01  2.5245E+00
 GRADIENT:   2.7296E-01  4.9460E+00 -2.7108E-01  7.2302E-02 -2.3683E+00 -7.1830E-02 -6.8542E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -649.194979775704        NO. OF FUNC. EVALS.:  53
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  1.1320E+00  1.1398E+00  8.0262E+02  3.9905E-01  2.1799E-01  1.1955E+01  1.1296E+01
 PARAMETER:  2.2399E-01  2.3089E-01  6.7879E+00 -8.1867E-01 -1.4233E+00  2.5812E+00  2.5245E+00
 GRADIENT:  -2.6386E-01 -6.7628E-01 -6.3296E-03  7.9062E-02  7.5679E-02 -1.6339E-03  2.4781E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -649.215619061532        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      458
 NPARAMETR:  1.1361E+00  1.1450E+00  2.6941E+04  4.0148E-01  2.2192E-01  7.0814E+01  1.1308E+01
 PARAMETER:  2.2759E-01  2.3539E-01  1.0301E+01 -8.1259E-01 -1.4054E+00  4.3601E+00  2.5255E+00
 GRADIENT:  -2.9128E-01 -5.3692E-01 -1.7809E-04 -1.6578E-01 -2.3892E-02 -5.0800E-05 -9.0167E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -649.217163342513        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      573
 NPARAMETR:  1.1376E+00  1.1469E+00  8.0213E+05  4.0208E-01  2.2274E-01  3.9322E+02  1.1322E+01
 PARAMETER:  2.2888E-01  2.3709E-01  1.3695E+01 -8.1111E-01 -1.4017E+00  6.0744E+00  2.5268E+00
 GRADIENT:   8.0090E-03  1.6743E-02 -5.8097E-06  3.1909E-03 -2.2113E-03 -1.7707E-06  1.1470E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -649.217170518223        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      689
 NPARAMETR:  1.1375E+00  1.1469E+00  7.5799E+07  4.0207E-01  2.2274E-01  3.9117E+03  1.1322E+01
 PARAMETER:  2.2885E-01  2.3706E-01  1.8244E+01 -8.1113E-01 -1.4017E+00  8.3717E+00  2.5267E+00
 GRADIENT:  -1.7235E-04 -2.1629E-04 -6.1369E-08 -1.2690E-04  4.8274E-05 -1.8944E-08 -5.8943E-04

0ITERATION NO.:   55    OBJECTIVE VALUE:  -649.217170644313        NO. OF FUNC. EVALS.: 113
 CUMULATIVE NO. OF FUNC. EVALS.:      802
 NPARAMETR:  1.1375E+00  1.1469E+00  1.1155E+08  4.0207E-01  2.2274E-01  4.7546E+03  1.1322E+01
 PARAMETER:  2.2885E-01  2.3706E-01  1.8630E+01 -8.1113E-01 -1.4017E+00  8.5669E+00  2.5267E+00
 GRADIENT:  -1.2381E-03  4.3440E-04  1.4370E-06 -2.3847E-03 -2.6000E-04 -5.5978E-07  1.7629E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      802
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7873E-02 -6.7576E-02  3.6664E-07
 SE:             8.4947E-02  5.2325E-02  1.2790E-06
 N:                     100         100         100

 P VAL.:         8.3336E-01  1.9654E-01  7.7437E-01

 ETASHRINKSD(%)  1.0007E+01  4.4567E+01  9.9999E+01
 ETASHRINKVR(%)  1.9012E+01  6.9272E+01  1.0000E+02
 EBVSHRINKSD(%)  9.1726E+00  4.5090E+01  9.9999E+01
 EBVSHRINKVR(%)  1.7504E+01  6.9849E+01  1.0000E+02
 RELATIVEINF(%)  1.0417E+01  3.2000E+00  1.0153E-09
 EPSSHRINKSD(%)  7.4452E+00
 EPSSHRINKVR(%)  1.4336E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -649.21717064431300     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       85.933655919425178     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     6.19
 Elapsed covariance  time in seconds:     1.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -649.217       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.14E+00  1.15E+00  1.12E+08  4.02E-01  2.23E-01  4.75E+03  1.13E+01
 


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
 
         8.74E-02  1.14E-01  9.03E+08  5.26E-02  1.28E-01  1.32E+04  1.57E+00
 


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
+        7.64E-03
 
 TH 2
+        8.08E-03  1.29E-02
 
 TH 3
+       -2.22E+07  1.66E+07  8.15E+17
 
 TH 4
+        1.55E-03  4.02E-03  2.75E+07  2.76E-03
 
 TH 5
+        3.77E-03  1.03E-02  9.19E+07  5.21E-03  1.63E-02
 
 TH 6
+        6.24E+01  7.44E+02  1.11E+13  4.83E+02  1.61E+03  1.74E+08
 
 TH 7
+        2.51E-02 -3.37E-02 -1.23E+09 -3.91E-02 -1.49E-01 -1.70E+04  2.45E+00
 
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
+        8.74E-02
 
 TH 2
+        8.13E-01  1.14E-01
 
 TH 3
+       -2.82E-01  1.61E-01  9.03E+08
 
 TH 4
+        3.37E-01  6.73E-01  5.79E-01  5.26E-02
 
 TH 5
+        3.38E-01  7.07E-01  7.98E-01  7.77E-01  1.28E-01
 
 TH 6
+        5.41E-02  4.95E-01  9.31E-01  6.95E-01  9.55E-01  1.32E+04
 
 TH 7
+        1.83E-01 -1.89E-01 -8.72E-01 -4.75E-01 -7.46E-01 -8.21E-01  1.57E+00
 
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
+        1.03E+07
 
 TH 2
+        5.99E+06  5.28E+06
 
 TH 3
+        1.61E-03  1.57E-03  4.76E-13
 
 TH 4
+       -7.99E+05 -2.38E+06 -8.00E-04  2.09E+06
 
 TH 5
+       -2.43E+07 -1.29E+07 -3.36E-03  5.66E+05  5.81E+07
 
 TH 6
+        8.06E+01 -3.72E+00 -5.22E-09  4.73E+01 -2.25E+02  2.04E-03
 
 TH 7
+       -1.44E+05 -4.30E+04 -8.18E-06 -3.19E+04  3.67E+05 -2.26E+00  2.92E+03
 
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
 #CPUT: Total CPU Time in Seconds,        7.973
Stop Time:
Sat Sep 25 15:04:26 CDT 2021
