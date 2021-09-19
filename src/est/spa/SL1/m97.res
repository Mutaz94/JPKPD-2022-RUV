Sat Sep 18 12:00:09 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat97.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1640.56252582434        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.9586E+02 -7.7879E+01 -4.0872E+01 -5.5036E+01  4.3417E+01 -2.5786E+01 -4.7158E+00  1.1320E+01  1.2597E+01  1.3307E+01
             2.8869E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1644.81326339932        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0013E+00  1.0900E+00  1.1017E+00  9.5930E-01  1.0796E+00  1.0853E+00  1.0559E+00  9.1508E-01  9.2247E-01  9.0942E-01
             8.9382E-01
 PARAMETER:  1.0130E-01  1.8616E-01  1.9685E-01  5.8450E-02  1.7660E-01  1.8189E-01  1.5440E-01  1.1259E-02  1.9302E-02  5.0470E-03
            -1.2256E-02
 GRADIENT:   1.8736E+02 -4.0987E+01 -1.7104E+00 -4.3590E+01  5.5221E+01  1.5779E+01 -5.2169E+00 -1.4925E+00 -6.8230E+00 -1.7988E+01
            -2.5852E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1648.26397475935        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  9.9264E-01  1.3225E+00  8.5827E-01  8.3667E-01  1.0624E+00  1.1209E+00  1.0226E+00  6.5399E-01  9.4192E-01  9.4199E-01
             9.0575E-01
 PARAMETER:  9.2608E-02  3.7949E-01 -5.2831E-02 -7.8321E-02  1.6052E-01  2.1410E-01  1.2231E-01 -3.2466E-01  4.0162E-02  4.0243E-02
             1.0065E-03
 GRADIENT:   1.0991E+02  3.1697E+00  4.4423E+00 -1.1819E+01  1.8288E+01  1.5706E+01  5.4884E+00 -9.3182E-02 -1.4923E+01 -5.5078E+00
            -1.8138E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1654.30131949788        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  9.2927E-01  1.3273E+00  6.8029E-01  8.2617E-01  9.5118E-01  1.0433E+00  9.4574E-01  3.4423E-01  9.8906E-01  8.4898E-01
             9.4133E-01
 PARAMETER:  2.6639E-02  3.8316E-01 -2.8523E-01 -9.0951E-02  4.9953E-02  1.4243E-01  4.4212E-02 -9.6645E-01  8.8996E-02 -6.3715E-02
             3.9539E-02
 GRADIENT:  -1.0137E+01  1.4950E+00 -5.1135E+00  1.3763E+00  8.9266E-01 -4.8644E+00 -4.4269E+00  8.2605E-01 -6.4730E+00  1.7796E+00
             2.2602E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1655.35510957570        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      624
 NPARAMETR:  9.3335E-01  1.1014E+00  7.6001E-01  9.6801E-01  8.7949E-01  1.0575E+00  1.1396E+00  1.7623E-01  9.0785E-01  8.2197E-01
             9.3450E-01
 PARAMETER:  3.1028E-02  1.9655E-01 -1.7443E-01  6.7490E-02 -2.8411E-02  1.5591E-01  2.3070E-01 -1.6360E+00  3.3244E-03 -9.6048E-02
             3.2254E-02
 GRADIENT:   2.8525E-01  3.7012E+00  5.9604E-01 -7.6147E-01 -3.2297E+00  4.0509E-01  8.9340E-01  1.7183E-01 -6.7811E-01  1.4071E+00
             1.5766E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1655.58083569669        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  9.3289E-01  9.3328E-01  7.7682E-01  1.0663E+00  8.1428E-01  1.0556E+00  1.3053E+00  8.0087E-02  8.4165E-01  7.8459E-01
             9.3371E-01
 PARAMETER:  3.0530E-02  3.0953E-02 -1.5254E-01  1.6416E-01 -1.0545E-01  1.5414E-01  3.6641E-01 -2.4246E+00 -7.2388E-02 -1.4260E-01
             3.1409E-02
 GRADIENT:   2.7725E-01 -4.4612E-01  8.3631E-01 -2.5318E+00 -1.3578E+00 -3.9091E-01  2.3227E-01  2.5446E-02 -6.5413E-02  2.1927E-01
            -2.0587E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1655.58603867338        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      978
 NPARAMETR:  9.3300E-01  9.2664E-01  7.7487E-01  1.0710E+00  8.1084E-01  1.0569E+00  1.3120E+00  6.3348E-02  8.3901E-01  7.8376E-01
             9.3374E-01
 PARAMETER:  3.0649E-02  2.3812E-02 -1.5506E-01  1.6858E-01 -1.0968E-01  1.5537E-01  3.7152E-01 -2.6591E+00 -7.5538E-02 -1.4365E-01
             3.1446E-02
 GRADIENT:   4.8185E-01 -4.0035E-01 -9.2381E-01  5.6424E-02  5.5715E-01  6.7630E-02  1.3207E-01  1.7520E-02  3.3889E-02  4.6962E-01
             1.2865E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1655.59522857338        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1154
 NPARAMETR:  9.3281E-01  9.3542E-01  7.7346E-01  1.0658E+00  8.1353E-01  1.0568E+00  1.3016E+00  1.0000E-02  8.4256E-01  7.8292E-01
             9.3378E-01
 PARAMETER:  3.0451E-02  3.3242E-02 -1.5689E-01  1.6371E-01 -1.0638E-01  1.5522E-01  3.6363E-01 -4.7229E+00 -7.1305E-02 -1.4472E-01
             3.1486E-02
 GRADIENT:   2.2821E-02 -1.6614E-03 -3.7189E-02 -3.8120E-03  4.2018E-03  6.1704E-03  1.7373E-02  0.0000E+00  3.4665E-03  2.0434E-02
             3.6461E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1655.59523276837        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1211
 NPARAMETR:  9.3280E-01  9.3555E-01  7.7361E-01  1.0657E+00  8.1367E-01  1.0568E+00  1.3014E+00  1.0000E-02  8.4264E-01  7.8300E-01
             9.3379E-01
 PARAMETER:  3.0438E-02  3.3381E-02 -1.5669E-01  1.6365E-01 -1.0620E-01  1.5520E-01  3.6341E-01 -4.7241E+00 -7.1218E-02 -1.4462E-01
             3.1496E-02
 GRADIENT:  -1.6624E-03 -3.3540E-04 -1.4253E-03 -3.3759E-03  3.8779E-04  6.9847E-04  2.0853E-03  0.0000E+00 -8.6382E-04  1.3050E-04
             2.2218E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1211
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.3327E-05  2.7580E-03 -4.7000E-04 -6.0857E-03 -9.2472E-03
 SE:             2.9860E-02  2.2149E-02  2.0054E-04  2.4211E-02  2.2598E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9751E-01  9.0090E-01  1.9095E-02  8.0154E-01  6.8240E-01

 ETASHRINKSD(%)  1.0000E-10  2.5799E+01  9.9328E+01  1.8889E+01  2.4292E+01
 ETASHRINKVR(%)  1.0000E-10  4.4943E+01  9.9995E+01  3.4209E+01  4.2684E+01
 EBVSHRINKSD(%)  3.4565E-01  2.5538E+01  9.9385E+01  1.9013E+01  2.3436E+01
 EBVSHRINKVR(%)  6.9010E-01  4.4554E+01  9.9996E+01  3.4411E+01  4.1380E+01
 RELATIVEINF(%)  9.9063E+01  4.0697E+00  4.1480E-04  5.6565E+00  4.9725E+00
 EPSSHRINKSD(%)  4.3910E+01
 EPSSHRINKVR(%)  6.8539E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1655.5952327683749     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -920.44440620463672     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1655.595       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.33E-01  9.36E-01  7.74E-01  1.07E+00  8.14E-01  1.06E+00  1.30E+00  1.00E-02  8.43E-01  7.83E-01  9.34E-01
 


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
+        1.14E+03
 
 TH 2
+       -5.00E+00  4.33E+02
 
 TH 3
+        1.52E+01  2.48E+02  7.72E+02
 
 TH 4
+       -7.63E+00  3.30E+02 -3.42E+02  9.14E+02
 
 TH 5
+       -3.81E+00 -4.31E+02 -9.71E+02  4.19E+02  1.60E+03
 
 TH 6
+        5.60E-01 -1.32E+00  3.09E+00 -2.84E+00 -2.64E+00  1.76E+02
 
 TH 7
+        1.03E+00  3.31E+01 -9.75E+00 -1.15E+01 -6.18E-02 -5.22E-02  4.23E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.97E+00 -2.75E+01 -3.70E+01  2.91E+01  4.72E+00  1.42E-01  1.49E+01  0.00E+00  1.35E+02
 
 TH10
+       -9.65E-01 -9.67E+00 -7.67E+01 -2.93E+01 -5.57E+01  2.58E-01  1.48E+01  0.00E+00  1.47E+01  1.15E+02
 
 TH11
+       -8.26E+00 -1.28E+01 -3.95E+01 -2.08E+00  1.35E+01  1.72E+00  4.69E+00  0.00E+00  1.01E+01  2.44E+01  2.49E+02
 
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
 #CPUT: Total CPU Time in Seconds,       20.274
Stop Time:
Sat Sep 18 12:00:31 CDT 2021
