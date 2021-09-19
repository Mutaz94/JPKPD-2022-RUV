Sat Sep 18 08:06:54 CDT 2021
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
$DATA ../../../../data/int/All/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26839.2017624821        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   7.7965E+01  9.6836E+00 -1.8484E+02 -1.9154E+02 -1.5940E+02 -3.0760E+01 -5.9764E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -998.488216932105        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:       63
 NPARAMETR:  1.5157E+00  2.5289E+00  3.6900E+00  6.4586E-01  1.1535E+00  4.9587E-01  1.4982E+01
 PARAMETER:  5.1588E-01  1.0278E+00  1.4056E+00 -3.3717E-01  2.4281E-01 -6.0145E-01  2.8068E+00
 GRADIENT:   7.5938E+01  3.2004E+01  3.5163E+01  2.0831E+01  1.0151E+01  2.1316E+01  5.3671E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1140.64153154614        NO. OF FUNC. EVALS.:  56
 CUMULATIVE NO. OF FUNC. EVALS.:      119
 NPARAMETR:  9.2383E-01  1.3173E+00  2.6072E+00  4.3863E-01  1.0443E+00  2.0328E-01  1.1851E+01
 PARAMETER:  2.0769E-02  3.7561E-01  1.0583E+00 -7.2411E-01  1.4337E-01 -1.4932E+00  2.5724E+00
 GRADIENT:  -6.4271E+01 -6.2800E+01 -5.5140E+00 -3.0796E+00 -5.4091E+01  4.1999E+00  4.6671E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1247.50472735146        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      170
 NPARAMETR:  1.2577E+00  1.5184E+00  2.3576E+00  8.6070E-01  1.0849E+00  1.1830E-01  8.1069E+00
 PARAMETER:  3.2927E-01  5.1766E-01  9.5765E-01 -5.0014E-02  1.8146E-01 -2.0345E+00  2.1927E+00
 GRADIENT:   1.2548E+01  5.8349E+00 -1.3717E+01 -1.6485E+01 -4.1477E+00 -6.1335E-02 -3.3307E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1248.23851737224        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      220
 NPARAMETR:  1.2127E+00  1.4752E+00  2.4112E+00  8.9786E-01  1.0994E+00  1.2529E-01  8.0895E+00
 PARAMETER:  2.9286E-01  4.8879E-01  9.8014E-01 -7.7424E-03  1.9480E-01 -1.9771E+00  2.1906E+00
 GRADIENT:   2.0873E+00 -1.0491E-01 -1.1742E-01 -1.0385E+00  4.3487E-01  1.7444E-01 -7.2896E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1248.23894390664        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      270
 NPARAMETR:  1.2084E+00  1.4755E+00  2.4114E+00  8.9892E-01  1.0990E+00  1.2459E-01  8.0893E+00
 PARAMETER:  2.8929E-01  4.8897E-01  9.8019E-01 -6.5626E-03  1.9439E-01 -1.9827E+00  2.1905E+00
 GRADIENT:   1.1406E+00 -6.4358E-02 -1.1475E-01 -6.2406E-01  2.6204E-01  1.6975E-01 -5.5252E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1248.23911548157        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      320
 NPARAMETR:  1.2053E+00  1.4757E+00  2.4119E+00  8.9969E-01  1.0987E+00  1.2269E-01  8.0897E+00
 PARAMETER:  2.8676E-01  4.8914E-01  9.8043E-01 -5.7096E-03  1.9409E-01 -1.9981E+00  2.1906E+00
 GRADIENT:   4.7502E-01 -3.0631E-02 -9.7650E-02 -3.2093E-01  1.3592E-01  1.5557E-01 -4.0683E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1248.23958539898        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.2014E+00  1.4761E+00  2.4135E+00  9.0068E-01  1.0982E+00  1.1753E-01  8.0910E+00
 PARAMETER:  2.8352E-01  4.8938E-01  9.8107E-01 -4.6091E-03  1.9367E-01 -2.0410E+00  2.1908E+00
 GRADIENT:  -3.7263E-01  1.3226E-02 -6.0774E-02  7.8998E-02 -3.0790E-02  1.1985E-01 -1.7943E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1248.24520662078        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  1.1974E+00  1.4767E+00  2.4200E+00  9.0164E-01  1.0975E+00  9.4252E-02  8.0977E+00
 PARAMETER:  2.8017E-01  4.8978E-01  9.8376E-01 -3.5365E-03  1.9306E-01 -2.2618E+00  2.1916E+00
 GRADIENT:  -1.2291E+00  6.3501E-02  6.9810E-02  5.4501E-01 -2.2615E-01  1.4125E-02  2.2808E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1248.30160717245        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  1.2128E+00  1.4942E+00  2.4231E+00  8.9984E-01  1.1010E+00  8.7901E-02  8.1180E+00
 PARAMETER:  2.9290E-01  5.0157E-01  9.8504E-01 -5.5348E-03  1.9619E-01 -2.3315E+00  2.1941E+00
 GRADIENT:  -2.1918E-01 -1.1900E+00 -5.1975E-01 -3.8473E-01 -1.4481E-01 -3.5934E-03 -5.0802E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1248.30662439294        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      622
 NPARAMETR:  1.2140E+00  1.5047E+00  2.4267E+00  9.0076E-01  1.1013E+00  8.7250E-02  8.1218E+00
 PARAMETER:  2.9392E-01  5.0859E-01  9.8652E-01 -4.5183E-03  1.9650E-01 -2.3390E+00  2.1946E+00
 GRADIENT:  -3.1605E-03 -4.8898E-03 -8.5183E-03  2.3452E-04 -2.0055E-03 -4.6335E-05 -1.6574E-02

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1248.30662439294        NO. OF FUNC. EVALS.:  14
 CUMULATIVE NO. OF FUNC. EVALS.:      636
 NPARAMETR:  1.2140E+00  1.5047E+00  2.4267E+00  9.0076E-01  1.1013E+00  8.7250E-02  8.1218E+00
 PARAMETER:  2.9392E-01  5.0859E-01  9.8652E-01 -4.5183E-03  1.9650E-01 -2.3390E+00  2.1946E+00
 GRADIENT:  -3.1605E-03 -4.8898E-03 -8.5183E-03  2.3452E-04 -2.0055E-03 -4.6335E-05 -1.6574E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      636
 NO. OF SIG. DIGITS IN FINAL EST.:  4.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0896E-02 -1.5618E-02 -9.8816E-03
 SE:             9.2926E-02  9.3706E-02  1.3446E-02
 N:                     100         100         100

 P VAL.:         8.2209E-01  8.6763E-01  4.6239E-01

 ETASHRINKSD(%)  1.5540E+00  7.2760E-01  8.5755E+01
 ETASHRINKVR(%)  3.0838E+00  1.4499E+00  9.7971E+01
 EBVSHRINKSD(%)  2.5461E+00  9.9024E-01  8.5645E+01
 EBVSHRINKVR(%)  5.0274E+00  1.9707E+00  9.7939E+01
 RELATIVEINF(%)  9.4924E+01  9.7150E+01  2.0413E+00
 EPSSHRINKSD(%)  8.6253E+00
 EPSSHRINKVR(%)  1.6507E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1248.3066243929436     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       405.78273537546715     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     9.36
 Elapsed covariance  time in seconds:     3.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1248.307       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.21E+00  1.50E+00  2.43E+00  9.01E-01  1.10E+00  8.72E-02  8.12E+00
 


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
 
         1.13E-01  1.71E-01  1.98E-01  7.14E-02  8.44E-02  4.58E-01  1.36E+00
 


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
+        1.28E-02
 
 TH 2
+        1.61E-02  2.92E-02
 
 TH 3
+        2.19E-03  7.20E-03  3.91E-02
 
 TH 4
+        4.15E-03  6.00E-03 -9.29E-05  5.10E-03
 
 TH 5
+        3.68E-03  9.64E-03  2.84E-03  3.00E-03  7.12E-03
 
 TH 6
+       -4.31E-03 -1.38E-02 -4.96E-02  3.36E-03 -2.93E-03  2.10E-01
 
 TH 7
+        4.57E-02  8.21E-02  1.49E-01  9.98E-03  2.31E-02 -5.05E-01  1.85E+00
 
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
+        1.13E-01
 
 TH 2
+        8.32E-01  1.71E-01
 
 TH 3
+        9.79E-02  2.13E-01  1.98E-01
 
 TH 4
+        5.13E-01  4.91E-01 -6.58E-03  7.14E-02
 
 TH 5
+        3.85E-01  6.68E-01  1.70E-01  4.98E-01  8.44E-02
 
 TH 6
+       -8.32E-02 -1.76E-01 -5.47E-01  1.03E-01 -7.58E-02  4.58E-01
 
 TH 7
+        2.96E-01  3.53E-01  5.53E-01  1.03E-01  2.01E-01 -8.10E-01  1.36E+00
 
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
+        3.66E+02
 
 TH 2
+       -2.34E+02  2.22E+02
 
 TH 3
+        1.13E+01 -8.00E+00  3.94E+01
 
 TH 4
+       -1.13E+02  3.80E+01  4.40E+00  3.36E+02
 
 TH 5
+        1.79E+02 -1.91E+02 -4.77E+00 -1.32E+02  3.63E+02
 
 TH 6
+       -1.08E+01  4.51E+00  4.44E+00 -1.34E+01 -3.83E+00  1.64E+01
 
 TH 7
+       -4.13E+00 -4.45E-02 -1.84E+00 -3.06E+00 -4.39E-01  4.28E+00  1.98E+00
 
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
 #CPUT: Total CPU Time in Seconds,       13.238
Stop Time:
Sat Sep 18 08:07:08 CDT 2021
