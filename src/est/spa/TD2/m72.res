Sat Sep 18 14:51:33 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat72.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1700.77998186782        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.9831E+01 -1.1983E+01  2.8753E+01 -4.7713E+01 -4.9688E+01  1.3903E+01  3.4528E+00  4.9079E-01  2.2481E+01 -2.0297E+00
            -3.1021E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1707.92813855886        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0447E+00  1.0326E+00  9.7524E-01  1.0136E+00  1.0604E+00  9.4401E-01  9.9069E-01  9.8076E-01  8.7667E-01  1.0240E+00
             1.0919E+00
 PARAMETER:  1.4369E-01  1.3212E-01  7.4930E-02  1.1352E-01  1.5868E-01  4.2383E-02  9.0644E-02  8.0571E-02 -3.1621E-02  1.2367E-01
             1.8789E-01
 GRADIENT:   4.1477E+01 -5.8481E+00 -1.1352E+01  9.9239E+00  2.0652E+01 -4.7884E+00 -3.5938E+00  2.5725E+00  2.0547E-01 -3.5289E+00
             6.3784E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1708.35520679619        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0407E+00  9.6159E-01  9.2659E-01  1.0592E+00  9.9236E-01  9.5904E-01  1.1271E+00  7.7812E-01  7.8687E-01  1.0350E+00
             1.0590E+00
 PARAMETER:  1.3987E-01  6.0829E-02  2.3757E-02  1.5752E-01  9.2328E-02  5.8174E-02  2.1967E-01 -1.5088E-01 -1.3970E-01  1.3441E-01
             1.5729E-01
 GRADIENT:   3.4408E+01  6.3165E+00 -6.6769E+00  2.4326E+01  9.4340E+00  1.5194E+00 -3.2915E+00  9.8794E-01 -7.7801E+00  2.6126E+00
            -5.7020E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1708.48522280143        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0343E+00  1.0048E+00  8.0934E-01  1.0188E+00  9.3889E-01  9.5944E-01  1.0996E+00  5.9408E-01  8.2084E-01  9.7594E-01
             1.0649E+00
 PARAMETER:  1.3376E-01  1.0481E-01 -1.1153E-01  1.1860E-01  3.6945E-02  5.8591E-02  1.9494E-01 -4.2075E-01 -9.7428E-02  7.5647E-02
             1.6284E-01
 GRADIENT:   1.3057E+01  5.8348E-01 -7.2052E+00  1.1199E+01  7.5114E+00  6.6138E-01 -1.6392E+00  9.6288E-01 -2.6715E+00  2.3281E+00
            -1.9466E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1708.48717876005        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      324
 NPARAMETR:  1.0334E+00  1.0059E+00  7.9861E-01  1.0162E+00  9.3129E-01  9.5941E-01  1.1021E+00  5.6828E-01  8.2329E-01  9.6616E-01
             1.0660E+00
 PARAMETER:  1.3286E-01  1.0591E-01 -1.2488E-01  1.1603E-01  2.8814E-02  5.8560E-02  1.9722E-01 -4.6514E-01 -9.4451E-02  6.5577E-02
             1.6388E-01
 GRADIENT:  -4.3059E+01 -3.0253E+00 -6.0533E+00  5.7470E-01  5.1912E+00 -4.1151E+00 -1.9070E+00  7.8132E-01 -2.4761E+00  1.8155E+00
            -1.6142E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1709.04470650909        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  1.0520E+00  9.0851E-01  7.9991E-01  1.0746E+00  8.8040E-01  9.6733E-01  1.2176E+00  4.5690E-01  7.9421E-01  9.3508E-01
             1.0701E+00
 PARAMETER:  1.5065E-01  4.0469E-03 -1.2325E-01  1.7191E-01 -2.7379E-02  6.6781E-02  2.9689E-01 -6.8328E-01 -1.3041E-01  3.2880E-02
             1.6777E-01
 GRADIENT:   5.9500E-01  1.5294E+00  5.0984E-01  1.8596E+00 -3.0440E-01  1.5643E-02 -2.6646E-01 -1.4608E-01 -5.6949E-01 -1.5960E-01
            -5.0364E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1709.07365437747        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      677
 NPARAMETR:  1.0504E+00  8.1631E-01  8.5990E-01  1.1326E+00  8.7598E-01  9.6600E-01  1.3227E+00  5.6189E-01  7.7075E-01  9.4450E-01
             1.0705E+00
 PARAMETER:  1.4913E-01 -1.0296E-01 -5.0935E-02  2.2455E-01 -3.2417E-02  6.5413E-02  3.7971E-01 -4.7644E-01 -1.6039E-01  4.2895E-02
             1.6812E-01
 GRADIENT:   1.9592E-01  2.4412E-01  2.6716E-01 -3.8024E-01 -1.0672E+00  2.7536E-02  1.6572E-01  8.6976E-02  2.7899E-01  1.8976E-01
             2.4575E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1709.07514165019        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      852
 NPARAMETR:  1.0500E+00  8.0504E-01  8.7203E-01  1.1404E+00  8.7927E-01  9.6571E-01  1.3347E+00  5.7570E-01  7.6708E-01  9.4883E-01
             1.0702E+00
 PARAMETER:  1.4881E-01 -1.1686E-01 -3.6929E-02  2.3138E-01 -2.8664E-02  6.5112E-02  3.8873E-01 -4.5216E-01 -1.6517E-01  4.7472E-02
             1.6782E-01
 GRADIENT:   6.6833E-03  2.7778E-02  3.1687E-02  2.4864E-02 -1.9095E-02  2.4493E-03  7.9300E-04 -3.9933E-03 -4.9119E-03 -1.2710E-02
            -4.0799E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1709.07514186151        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  1.0500E+00  8.0468E-01  8.7198E-01  1.1406E+00  8.7910E-01  9.6571E-01  1.3352E+00  5.7553E-01  7.6698E-01  9.4879E-01
             1.0702E+00
 PARAMETER:  1.4881E-01 -1.1732E-01 -3.6988E-02  2.3155E-01 -2.8860E-02  6.5105E-02  3.8907E-01 -4.5247E-01 -1.6529E-01  4.7437E-02
             1.6781E-01
 GRADIENT:   2.3101E-03  4.5527E-03  5.6912E-03  7.7105E-03  1.7460E-03  7.8330E-04 -4.2427E-04 -1.5966E-03 -1.5920E-03 -4.4713E-03
            -2.3282E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      909
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0076E-04  3.3153E-03 -2.3333E-02 -8.3801E-03 -2.0998E-02
 SE:             2.9829E-02  1.9744E-02  1.0596E-02  2.3931E-02  2.2261E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9730E-01  8.6666E-01  2.7664E-02  7.2621E-01  3.4553E-01

 ETASHRINKSD(%)  6.7927E-02  3.3854E+01  6.4501E+01  1.9827E+01  2.5423E+01
 ETASHRINKVR(%)  1.3581E-01  5.6247E+01  8.7398E+01  3.5722E+01  4.4383E+01
 EBVSHRINKSD(%)  5.1505E-01  3.4519E+01  6.7282E+01  1.9760E+01  2.3112E+01
 EBVSHRINKVR(%)  1.0274E+00  5.7122E+01  8.9295E+01  3.5616E+01  4.0882E+01
 RELATIVEINF(%)  9.8000E+01  2.2502E+00  1.2937E+00  3.9224E+00  6.3666E+00
 EPSSHRINKSD(%)  4.3579E+01
 EPSSHRINKVR(%)  6.8166E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1709.0751418615130     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -973.92431529777480     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.03
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1709.075       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  8.05E-01  8.72E-01  1.14E+00  8.79E-01  9.66E-01  1.34E+00  5.76E-01  7.67E-01  9.49E-01  1.07E+00
 


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
+        1.07E+03
 
 TH 2
+       -1.49E+01  4.32E+02
 
 TH 3
+        1.47E+01  1.73E+02  5.04E+02
 
 TH 4
+       -1.23E+01  4.42E+02 -2.17E+02  9.34E+02
 
 TH 5
+       -2.87E+00 -3.00E+02 -5.96E+02  2.00E+02  9.41E+02
 
 TH 6
+       -1.29E+00 -2.54E+00  1.12E+00 -1.43E+00  4.83E-01  2.03E+02
 
 TH 7
+        1.62E+00  2.73E+01  3.05E+00 -5.16E+00 -8.09E+00 -1.92E-01  2.69E+01
 
 TH 8
+       -1.06E+00 -1.34E+01 -5.23E+01  2.05E+00  1.17E+01 -1.48E-01  1.78E+00  1.91E+01
 
 TH 9
+        3.96E-01 -1.88E+01 -2.31E+01  5.13E+00  1.34E+01 -1.72E+00  2.40E+01  4.77E+00  1.50E+02
 
 TH10
+       -2.18E+00  7.02E-01 -3.12E+01 -1.80E+01 -6.21E+01 -2.42E-01  5.44E+00  1.74E+01  4.49E+00  7.99E+01
 
 TH11
+       -7.58E+00 -1.22E+01 -2.91E+01 -7.18E+00  1.00E+00  2.22E+00  3.58E+00  1.01E+01  9.45E+00  1.62E+01  1.88E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.017
Stop Time:
Sat Sep 18 14:51:50 CDT 2021
