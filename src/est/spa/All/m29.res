Sat Sep 25 14:59:34 CDT 2021
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
$DATA ../../../../data/spa/All/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   8682.05395685858        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   5.5229E+02 -1.6821E+02 -1.3938E+01  3.4174E+02  2.4577E+02 -1.7678E+02 -1.9511E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -654.503312079309        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:       61
 NPARAMETR:  1.5670E+00  1.5189E+00  2.3496E+00  4.7498E-01  3.5924E-01  6.5466E-01  1.5118E+01
 PARAMETER:  5.4917E-01  5.1797E-01  9.5424E-01 -6.4449E-01 -9.2377E-01 -3.2364E-01  2.8159E+00
 GRADIENT:   8.4210E+01 -3.3590E+01 -8.6365E+00 -3.2067E+00  2.1660E+01  2.8920E+00  2.0357E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -710.627490889415        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      115
 NPARAMETR:  1.1245E+00  1.1874E+00  1.1057E+01  4.6185E-01  1.8286E-01  1.7519E+00  1.1100E+01
 PARAMETER:  2.1737E-01  2.7178E-01  2.5031E+00 -6.7251E-01 -1.5990E+00  6.6069E-01  2.5069E+00
 GRADIENT:  -1.9326E+01 -5.2703E+01 -9.0024E-01  2.7077E+01 -4.4216E+00 -2.0014E-01  8.2853E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -723.559996219843        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  1.1928E+00  1.3420E+00  3.4578E+01  5.2118E-01  4.0995E-01  4.5711E+00  9.3467E+00
 PARAMETER:  2.7633E-01  3.9419E-01  3.6432E+00 -5.5165E-01 -7.9171E-01  1.6197E+00  2.3350E+00
 GRADIENT:  -1.0495E-01 -5.2512E+00  4.5612E-03  2.7713E+00  9.3558E-01  2.1991E-01  8.8359E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -723.683078394403        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      217
 NPARAMETR:  1.1925E+00  1.3523E+00  3.5591E+01  5.2376E-01  4.2178E-01  4.4019E+00  9.2041E+00
 PARAMETER:  2.7605E-01  4.0182E-01  3.6721E+00 -5.4673E-01 -7.6328E-01  1.5820E+00  2.3196E+00
 GRADIENT:  -3.6371E-01 -7.1635E-01 -5.6862E-02  2.3450E-01  3.6379E-01 -1.8620E-02  4.9104E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -723.683297331757        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  1.1929E+00  1.3527E+00  3.6295E+01  5.2293E-01  4.2013E-01  4.4654E+00  9.2054E+00
 PARAMETER:  2.7637E-01  4.0208E-01  3.6917E+00 -5.4832E-01 -7.6720E-01  1.5964E+00  2.3198E+00
 GRADIENT:   1.0388E-02 -1.8912E-03 -8.0575E-02  1.7418E-02  1.7177E-03  3.5142E-02 -2.2412E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -723.684436514291        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      318
 NPARAMETR:  1.1932E+00  1.3527E+00  4.1020E+01  5.2152E-01  4.1743E-01  4.7453E+00  9.2092E+00
 PARAMETER:  2.7661E-01  4.0212E-01  3.8141E+00 -5.5100E-01 -7.7364E-01  1.6572E+00  2.3202E+00
 GRADIENT:   5.4460E-01  1.0228E+00 -9.1877E-02 -2.9895E-01 -5.1950E-01  7.6865E-02 -1.2234E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -723.695724236850        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.1924E+00  1.3509E+00  1.0816E+02  5.1853E-01  4.1199E-01  7.6689E+00  9.2216E+00
 PARAMETER:  2.7595E-01  4.0078E-01  4.7836E+00 -5.5676E-01 -7.8676E-01  2.1372E+00  2.3216E+00
 GRADIENT:   1.3188E+00  2.5892E+00 -3.0467E-02 -8.0580E-01 -1.3286E+00  2.3390E-02 -2.7000E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -723.736288710056        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      424
 NPARAMETR:  1.1900E+00  1.3486E+00  1.8946E+03  5.2167E-01  4.1851E-01  3.1933E+01  9.2186E+00
 PARAMETER:  2.7396E-01  3.9905E-01  7.6468E+00 -5.5073E-01 -7.7104E-01  3.5636E+00  2.3212E+00
 GRADIENT:  -9.3240E-02 -1.9557E-01 -8.2837E-04  5.0333E-02  9.4459E-02 -7.1792E-04  1.9093E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -723.736765764605        NO. OF FUNC. EVALS.:  53
 CUMULATIVE NO. OF FUNC. EVALS.:      477
 NPARAMETR:  1.1901E+00  1.3486E+00  2.0998E+04  5.2144E-01  4.1804E-01  1.0480E+02  9.2191E+00
 PARAMETER:  2.7400E-01  3.9908E-01  1.0052E+01 -5.5117E-01 -7.7217E-01  4.7520E+00  2.3213E+00
 GRADIENT:  -5.8238E-03 -2.9103E-03 -5.7518E-05  2.3709E-05  3.3126E-04 -9.7529E-05  7.8007E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -723.757345520296        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  1.1966E+00  1.3584E+00  5.2033E+05  5.2769E-01  4.2663E-01  5.1366E+02  9.2120E+00
 PARAMETER:  2.7948E-01  4.0633E-01  1.3262E+01 -5.3924E-01 -7.5185E-01  6.3416E+00  2.3205E+00
 GRADIENT:   1.9580E-02  4.2303E-03 -1.7062E-06  7.2090E-04 -3.6617E-03 -4.4126E-06  5.4448E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -723.757349704949        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      699
 NPARAMETR:  1.1966E+00  1.3584E+00  1.2334E+07  5.2770E-01  4.2665E-01  2.4556E+03  9.2118E+00
 PARAMETER:  2.7945E-01  4.0633E-01  1.6428E+01 -5.3923E-01 -7.5179E-01  7.9061E+00  2.3205E+00
 GRADIENT:   3.6426E-03  5.0097E-04 -7.3757E-08 -5.4051E-04 -8.5562E-04 -1.9005E-07 -2.5169E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -723.757350038054        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      814
 NPARAMETR:  1.1966E+00  1.3584E+00  6.9916E+07  5.2770E-01  4.2666E-01  5.7884E+03  9.2118E+00
 PARAMETER:  2.7945E-01  4.0633E-01  1.8163E+01 -5.3922E-01 -7.5178E-01  8.7636E+00  2.3205E+00
 GRADIENT:  -2.3304E-03  1.8111E-03 -2.4537E-06 -2.0142E-03 -4.5979E-04  3.2130E-06  1.1728E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      814
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3499E-03 -8.0370E-02 -1.7858E-06
 SE:             8.6803E-02  7.5119E-02  2.2397E-06
 N:                     100         100         100

 P VAL.:         9.6922E-01  2.8466E-01  4.2524E-01

 ETASHRINKSD(%)  8.0407E+00  2.0419E+01  9.9998E+01
 ETASHRINKVR(%)  1.5435E+01  3.6669E+01  1.0000E+02
 EBVSHRINKSD(%)  6.4886E+00  1.9398E+01  9.9998E+01
 EBVSHRINKVR(%)  1.2556E+01  3.5033E+01  1.0000E+02
 RELATIVEINF(%)  2.0545E+01  1.3718E+01  4.4116E-09
 EPSSHRINKSD(%)  1.3097E+01
 EPSSHRINKVR(%)  2.4479E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -723.75735003805437     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       11.393476525683809     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     7.75
 Elapsed covariance  time in seconds:     2.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -723.757       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.20E+00  1.36E+00  6.99E+07  5.28E-01  4.27E-01  5.79E+03  9.21E+00
 


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
 
         9.77E-02  1.89E-01  1.65E+09  1.50E-01  2.71E-01  1.50E+05  2.09E+00
 


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
+        9.55E-03
 
 TH 2
+        1.28E-02  3.57E-02
 
 TH 3
+       -3.03E+07 -2.51E+08  2.73E+18
 
 TH 4
+        4.97E-03  2.44E-02 -2.37E+08  2.26E-02
 
 TH 5
+        8.98E-03  4.43E-02 -4.40E+08  3.96E-02  7.34E-02
 
 TH 6
+        3.41E+03  2.32E+04 -2.47E+14  2.16E+04  4.03E+04  2.25E+10
 
 TH 7
+       -5.51E-03 -2.57E-01  3.27E+09 -2.74E-01 -5.11E-01 -2.95E+05  4.39E+00
 
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
+        9.77E-02
 
 TH 2
+        6.95E-01  1.89E-01
 
 TH 3
+       -1.88E-01 -8.03E-01  1.65E+09
 
 TH 4
+        3.39E-01  8.58E-01 -9.54E-01  1.50E-01
 
 TH 5
+        3.39E-01  8.65E-01 -9.84E-01  9.72E-01  2.71E-01
 
 TH 6
+        2.33E-01  8.18E-01 -9.98E-01  9.61E-01  9.93E-01  1.50E+05
 
 TH 7
+       -2.69E-02 -6.50E-01  9.46E-01 -8.69E-01 -9.01E-01 -9.39E-01  2.09E+00
 
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
+        8.85E+06
 
 TH 2
+       -6.35E+06  4.76E+06
 
 TH 3
+       -3.76E-03  3.02E-03  2.10E-12
 
 TH 4
+        6.04E+05 -1.99E+05  1.10E-04  3.11E+05
 
 TH 5
+        4.58E+04 -1.91E+06 -2.97E-03 -2.16E+06  1.73E+07
 
 TH 6
+       -3.51E+01  3.19E+01  2.54E-08  5.32E+00 -6.21E+01  3.60E-04
 
 TH 7
+        1.33E+05 -7.50E+04 -2.48E-05  3.22E+04 -1.85E+05  1.38E-01  3.98E+03
 
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
 #CPUT: Total CPU Time in Seconds,        8.881
Stop Time:
Sat Sep 25 14:59:50 CDT 2021
