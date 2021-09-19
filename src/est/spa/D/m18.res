Sat Sep 18 15:10:30 CDT 2021
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
$DATA ../../../../data/spa/D/dat18.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m18.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -412.066997270966        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.0740E+02 -1.3995E+02 -5.2539E+01 -2.5158E+02  4.1062E+02 -9.3778E+02 -3.4066E+02 -5.1601E+01 -6.4734E+02 -6.4019E+02
            -1.0929E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1269.74726170119        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      132
 NPARAMETR:  1.0253E+00  7.8431E-01  9.6787E-01  1.2906E+00  8.4934E-01  2.6573E+00  1.5942E+00  1.1266E+00  2.9669E+00  2.2263E+00
             1.0668E+00
 PARAMETER:  1.2500E-01 -1.4294E-01  6.7344E-02  3.5509E-01 -6.3298E-02  1.0773E+00  5.6640E-01  2.1917E-01  1.1875E+00  9.0032E-01
             1.6470E-01
 GRADIENT:  -4.1396E+01 -1.5322E+01 -6.2543E+01 -9.3328E+00  4.6864E+01  3.4968E+01  3.4654E+00 -5.6852E+00  4.1168E+01 -1.3873E+01
             1.5608E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1323.02812302722        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  8.5321E-01  1.1728E-01  1.5968E+01  2.1783E+00  1.7854E+00  3.2669E+00  3.8717E+00  5.1099E+00  1.9209E+00  1.7266E+00
             1.0138E+00
 PARAMETER: -5.8747E-02 -2.0432E+00  2.8706E+00  8.7855E-01  6.7966E-01  1.2838E+00  1.4537E+00  1.7312E+00  7.5281E-01  6.4613E-01
             1.1373E-01
 GRADIENT:  -6.1393E+01  4.0740E+00  6.0614E+00  5.5389E+01  5.7386E+01  9.4394E+01  2.6931E-01 -2.7259E+01 -5.4092E+01 -4.4418E+01
            -5.1646E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1354.49530600990        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      484
 NPARAMETR:  9.8535E-01  1.3173E-01  1.3113E+01  2.0732E+00  1.6533E+00  2.7852E+00  3.5450E+00  4.8988E+00  2.0792E+00  1.8378E+00
             1.0883E+00
 PARAMETER:  8.5241E-02 -1.9270E+00  2.6736E+00  8.2907E-01  6.0275E-01  1.1243E+00  1.3655E+00  1.6890E+00  8.3200E-01  7.0859E-01
             1.8459E-01
 GRADIENT:  -4.1540E+01  2.8084E+00  9.7917E+00  3.2088E+01  3.9921E+00  7.6586E+01  4.9922E-01 -3.1075E+01  2.3689E+00 -2.7168E+01
            -1.1379E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1357.04221925088        NO. OF FUNC. EVALS.: 214
 CUMULATIVE NO. OF FUNC. EVALS.:      698
 NPARAMETR:  9.9670E-01  1.3170E-01  1.3040E+01  2.0652E+00  1.6550E+00  2.7427E+00  3.5669E+00  4.9843E+00  2.0748E+00  1.8538E+00
             1.0886E+00
 PARAMETER:  9.6695E-02 -1.9272E+00  2.6680E+00  8.2521E-01  6.0382E-01  1.1090E+00  1.3717E+00  1.7063E+00  8.2988E-01  7.1723E-01
             1.8486E-01
 GRADIENT:  -3.9293E+01  2.8287E+00  9.0497E+00  3.1007E+01  3.0020E+00  7.4298E+01  4.6071E-01 -2.9810E+01  3.0813E+00 -2.5629E+01
            -1.0748E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1357.36519523253        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:      804
 NPARAMETR:  9.9665E-01  1.3183E-01  1.3054E+01  2.0654E+00  1.6545E+00  2.7418E+00  1.0000E-02  4.9791E+00  2.0755E+00  1.8532E+00
             1.0887E+00
 PARAMETER:  9.6642E-02 -1.9263E+00  2.6691E+00  8.2533E-01  6.0349E-01  1.1086E+00 -8.3931E+00  1.7052E+00  8.3022E-01  7.1689E-01
             1.8496E-01
 GRADIENT:  -3.0935E+00  3.8012E+00  1.4985E+01  1.6134E+02  5.2684E+00  2.5971E+02  0.0000E+00  1.7412E+01  4.7222E+01 -2.4318E+01
            -1.0676E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1357.37694339141        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      929
 NPARAMETR:  9.9665E-01  1.3183E-01  1.3053E+01  2.0653E+00  1.6545E+00  2.7414E+00  1.0000E-02  4.9789E+00  2.0755E+00  1.8532E+00
             1.0887E+00
 PARAMETER:  9.6642E-02 -1.9263E+00  2.6690E+00  8.2528E-01  6.0349E-01  1.1085E+00 -9.7324E+00  1.7052E+00  8.3021E-01  7.1690E-01
             1.8496E-01
 GRADIENT:  -3.9566E+01  2.4282E+00  8.9325E+00  3.1531E+01  2.7311E+00  7.3436E+01  0.0000E+00 -2.8653E+01 -3.0743E-01 -2.5900E+01
            -1.0825E+01

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1357.37694339141        NO. OF FUNC. EVALS.:  38
 CUMULATIVE NO. OF FUNC. EVALS.:      967
 NPARAMETR:  9.9655E-01  1.3208E-01  1.3088E+01  2.0670E+00  1.6535E+00  2.7444E+00  1.0000E-02  4.9705E+00  2.0772E+00  1.8518E+00
             1.0889E+00
 PARAMETER:  9.6642E-02 -1.9263E+00  2.6690E+00  8.2528E-01  6.0349E-01  1.1085E+00 -9.7324E+00  1.7052E+00  8.3021E-01  7.1690E-01
             1.8496E-01
 GRADIENT:   1.2405E+05 -6.4393E+03 -4.6387E+03 -1.5002E+04  1.0280E+04 -1.1120E+04  0.0000E+00  7.2462E+03 -7.4713E+03  1.7279E+04
            -3.3568E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      967
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.7816E-02 -6.1235E-05 -4.4036E-02 -3.5628E-02 -6.0140E-02
 SE:             2.3184E-02  1.2326E-05  1.6951E-02  2.9175E-02  2.5038E-02
 N:                     100         100         100         100         100

 P VAL.:         3.9165E-02  6.7764E-07  9.3829E-03  2.2202E-01  1.6308E-02

 ETASHRINKSD(%)  2.2330E+01  9.9959E+01  4.3210E+01  2.2591E+00  1.6119E+01
 ETASHRINKVR(%)  3.9674E+01  1.0000E+02  6.7750E+01  4.4672E+00  2.9640E+01
 EBVSHRINKSD(%)  2.1598E-01  9.9968E+01  6.3160E+01  1.3723E+00  1.9153E+01
 EBVSHRINKVR(%)  4.3150E-01  1.0000E+02  8.6428E+01  2.7258E+00  3.4637E+01
 RELATIVEINF(%)  9.9455E+01  3.3107E-06  1.1167E+01  3.2110E+01  5.1363E+01
 EPSSHRINKSD(%)  4.2458E+01
 EPSSHRINKVR(%)  6.6889E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1357.3769433914104     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -622.22611682767217     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.28
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0T MATRIX UNOBTAINABLE
 Elapsed covariance  time in seconds:     7.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1357.377       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.97E-01  1.32E-01  1.31E+01  2.07E+00  1.65E+00  2.74E+00  1.00E-02  4.98E+00  2.08E+00  1.85E+00  1.09E+00
 


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
+        3.12E+08
 
 TH 2
+        8.20E+03  4.81E+07
 
 TH 3
+       -8.93E+00 -2.35E+01  2.55E+03
 
 TH 4
+       -1.84E+02 -4.36E+02 -1.34E+01  1.07E+06
 
 TH 5
+        3.06E+02  8.13E+02  6.00E+01  4.63E+02  3.11E+06
 
 TH 6
+       -8.78E+01 -2.74E+02 -1.23E+01 -1.53E+02  4.26E+02  3.36E+05
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        3.68E+01  1.01E+02  1.05E+04  5.55E+01 -2.49E+02  5.09E+01  0.00E+00  4.30E+04
 
 TH 9
+       -1.80E+02 -5.29E+02 -3.28E+01 -2.79E+02  1.14E+03 -2.49E+02  0.00E+00  1.36E+02  1.04E+06
 
 TH10
+        2.33E+02  6.16E+02  3.47E+01  3.49E+02 -1.21E+03  3.21E+02  0.00E+00 -1.42E+02  6.96E+02  1.76E+06
 
 TH11
+       -1.23E+05  4.81E+04  3.51E+02  7.17E+03 -1.22E+04  4.02E+03  0.00E+00 -1.44E+03  7.09E+03 -9.19E+03  7.66E+07
 
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
 #CPUT: Total CPU Time in Seconds,       22.750
Stop Time:
Sat Sep 18 15:10:54 CDT 2021
