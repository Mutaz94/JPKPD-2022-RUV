Sat Sep 18 15:49:01 CDT 2021
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
$DATA ../../../../data/spa/All/dat10.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m10.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1065.04967275994        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   4.7281E+01  1.8717E+01 -1.0468E+02  1.5920E+02  1.5568E+02 -3.0220E+01 -4.5831E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -819.301616906714        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:       60
 NPARAMETR:  9.8674E-01  1.0454E+00  4.0991E+00  1.5043E-01  8.9865E-02  4.5022E-01  1.0441E+01
 PARAMETER:  8.6655E-02  1.4441E-01  1.5108E+00 -1.7942E+00 -2.3094E+00 -6.9802E-01  2.4457E+00
 GRADIENT:  -2.1060E+02 -9.1386E+01 -4.6767E+00 -1.0006E+02  7.2216E+00 -8.8023E-02  2.4879E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -977.386533009802        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  1.0448E+00  1.0160E+00  5.2272E+01  3.3766E-01  1.2314E-02  5.7744E+00  6.7807E+00
 PARAMETER:  1.4383E-01  1.1585E-01  4.0565E+00 -9.8572E-01 -4.2970E+00  1.8534E+00  2.0141E+00
 GRADIENT:   9.0915E+01 -1.7387E+02  4.2064E-01  5.4480E+00 -1.3502E-01  1.0616E+00  1.6229E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -985.065956756675        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      166
 NPARAMETR:  1.0314E+00  1.0342E+00  3.1581E+02  3.0667E-01  1.1990E-01  2.7795E+00  6.6842E+00
 PARAMETER:  1.3092E-01  1.3359E-01  5.8551E+00 -1.0820E+00 -2.0211E+00  1.1223E+00  1.9997E+00
 GRADIENT:   4.4179E+01 -1.2840E+02 -4.6172E-02 -3.1695E+01 -4.5917E+00 -9.2818E-04  4.4188E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1007.55502042544        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      218
 NPARAMETR:  1.0356E+00  1.1400E+00  1.7293E+03  4.1457E-01  3.1928E-01  1.1300E+00  5.2927E+00
 PARAMETER:  1.3494E-01  2.3102E-01  7.5555E+00 -7.8052E-01 -1.0417E+00  2.2224E-01  1.7663E+00
 GRADIENT:  -7.0915E+00 -5.1811E+00 -3.0157E-03 -5.9248E-01  2.2889E+00 -4.4606E-06  6.1146E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1007.63485488483        NO. OF FUNC. EVALS.:  90
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  1.0449E+00  1.1471E+00  1.8657E+03  4.1815E-01  3.1726E-01  1.1456E+00  5.3084E+00
 PARAMETER:  1.4394E-01  2.3724E-01  7.6314E+00 -7.7192E-01 -1.0480E+00  2.3590E-01  1.7693E+00
 GRADIENT:   4.2647E-02  1.0685E-01 -2.0586E-03  2.8516E-03 -8.2163E-03 -3.8950E-06 -7.6141E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1007.63672974986        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      403
 NPARAMETR:  1.0445E+00  1.1466E+00  6.5539E+04  4.1784E-01  3.1707E-01  1.1463E+00  5.3099E+00
 PARAMETER:  1.4357E-01  2.3679E-01  1.1190E+01 -7.7266E-01 -1.0486E+00  2.3655E-01  1.7696E+00
 GRADIENT:  -2.6518E-01 -3.7198E-01 -5.9789E-05 -1.2366E-01  1.5871E-02 -2.8836E-09  5.5881E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1007.63692416460        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  1.0448E+00  1.1470E+00  6.6244E+05  4.1810E-01  3.1720E-01  1.1464E+00  5.3093E+00
 PARAMETER:  1.4387E-01  2.3713E-01  1.3504E+01 -7.7203E-01 -1.0482E+00  2.3659E-01  1.7695E+00
 GRADIENT:   1.6759E-02  9.6024E-03 -5.8240E-06  6.7951E-03  1.8836E-03 -4.8053E-10  2.6564E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1007.63693050257        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      637
 NPARAMETR:  1.0448E+00  1.1470E+00  1.9496E+07  4.1810E-01  3.1720E-01  1.1467E+00  5.3091E+00
 PARAMETER:  1.4386E-01  2.3712E-01  1.6886E+01 -7.7204E-01 -1.0482E+00  2.3693E-01  1.7694E+00
 GRADIENT:   2.8437E-03  2.3634E-03 -1.9799E-07  7.3309E-05  1.3877E-04  0.0000E+00  5.5315E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1007.63693085389        NO. OF FUNC. EVALS.: 112
 CUMULATIVE NO. OF FUNC. EVALS.:      749
 NPARAMETR:  1.0448E+00  1.1470E+00  1.3573E+08  4.1811E-01  3.1721E-01  1.1469E+00  5.3090E+00
 PARAMETER:  1.4385E-01  2.3712E-01  1.8826E+01 -7.7201E-01 -1.0482E+00  2.3703E-01  1.7694E+00
 GRADIENT:   2.7230E-03 -5.8280E-03 -2.8443E-08  9.0790E-03 -6.7554E-04  0.0000E+00  1.9191E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      749
 NO. OF SIG. DIGITS IN FINAL EST.:  4.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.0397E-04 -5.3747E-02  1.1402E-10
 SE:             8.9540E-02  7.9433E-02  2.9083E-10
 N:                     100         100         100

 P VAL.:         9.9284E-01  4.9864E-01  6.9501E-01

 ETASHRINKSD(%)  5.1406E+00  1.5849E+01  1.0000E+02
 ETASHRINKVR(%)  1.0017E+01  2.9185E+01  1.0000E+02
 EBVSHRINKSD(%)  4.4785E+00  1.5909E+01  1.0000E+02
 EBVSHRINKVR(%)  8.7565E+00  2.9287E+01  1.0000E+02
 RELATIVEINF(%)  9.0204E+01  6.9975E+01  1.0000E-10
 EPSSHRINKSD(%)  1.8984E+01
 EPSSHRINKVR(%)  3.4364E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1007.6369308538868     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -272.48610429014866     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     6.32
 Elapsed covariance  time in seconds:     1.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1007.637       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.04E+00  1.15E+00  1.36E+08  4.18E-01  3.17E-01  1.15E+00  5.31E+00
 


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
 
         5.21E-02  5.39E-02  7.92E+09  4.90E-02  5.75E-02  1.05E+00  5.25E-01
 


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
+        2.71E-03
 
 TH 2
+        1.79E-03  2.90E-03
 
 TH 3
+       -3.66E+08 -1.31E+08  6.28E+19
 
 TH 4
+       -3.76E-05  1.19E-03 -1.66E+06  2.40E-03
 
 TH 5
+       -4.07E-05  1.49E-03  8.35E+07  2.19E-03  3.31E-03
 
 TH 6
+        4.86E-02  1.71E-02 -8.35E+09 -2.50E-04 -1.13E-02  1.11E+00
 
 TH 7
+        1.23E-02 -3.54E-03 -2.61E+09 -1.23E-02 -2.17E-02  3.47E-01  2.76E-01
 
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
+        5.21E-02
 
 TH 2
+        6.39E-01  5.39E-02
 
 TH 3
+       -8.86E-01 -3.08E-01  7.92E+09
 
 TH 4
+       -1.47E-02  4.51E-01 -4.28E-03  4.90E-02
 
 TH 5
+       -1.36E-02  4.81E-01  1.83E-01  7.78E-01  5.75E-02
 
 TH 6
+        8.86E-01  3.02E-01 -1.00E+00 -4.84E-03 -1.86E-01  1.05E+00
 
 TH 7
+        4.51E-01 -1.25E-01 -6.27E-01 -4.78E-01 -7.17E-01  6.28E-01  5.25E-01
 
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
+        1.32E+08
 
 TH 2
+       -6.57E+07  3.33E+07
 
 TH 3
+       -4.22E-03  2.88E-03  1.04E-12
 
 TH 4
+        4.19E+07 -1.87E+07  1.10E-03  1.99E+07
 
 TH 5
+       -1.07E+07  4.07E+06 -1.08E-03 -7.22E+06  3.10E+06
 
 TH 6
+       -3.69E+07  2.42E+07  7.94E-03  6.58E+06 -7.66E+06  6.08E+07
 
 TH 7
+        8.25E+05 -3.75E+05  1.35E-05  3.68E+05 -1.29E+05  6.78E+04  6.90E+03
 
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
 #CPUT: Total CPU Time in Seconds,        8.287
Stop Time:
Sat Sep 18 15:49:11 CDT 2021
