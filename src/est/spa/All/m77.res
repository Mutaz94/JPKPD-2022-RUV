Sat Sep 25 15:11:53 CDT 2021
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
$DATA ../../../../data/spa/All/dat77.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12518.3977622009        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   1.9130E+02  7.8560E+00 -7.0461E+01 -3.3981E+02 -2.4664E+01 -6.6980E+01 -2.6647E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -516.095378868719        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:       60
 NPARAMETR:  1.6568E+00  1.8604E+00  2.2484E+00  6.3187E-01  8.1097E-01  5.7974E-01  1.5286E+01
 PARAMETER:  6.0491E-01  7.2081E-01  9.1020E-01 -3.5908E-01 -1.0952E-01 -4.4518E-01  2.8269E+00
 GRADIENT:   3.2164E+01  1.0338E+01 -9.0038E+00 -1.5707E+01  2.5339E+01  2.4428E+00  1.2404E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -527.795061558971        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:      119
 NPARAMETR:  1.4735E+00  1.8789E+00  3.5417E+00  6.5035E-01  8.4996E-01  7.8663E-01  1.4744E+01
 PARAMETER:  4.8766E-01  7.3066E-01  1.3646E+00 -3.3024E-01 -6.2570E-02 -1.4000E-01  2.7908E+00
 GRADIENT:  -3.0690E+00  1.0646E+01 -1.2423E+00 -1.5121E+00  4.1621E+00  2.7221E-01  1.3874E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -543.088455013954        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      176
 NPARAMETR:  1.3396E+00  1.4512E+00  6.6678E+00  6.4925E-01  6.7740E-01  2.1882E+00  1.1576E+01
 PARAMETER:  3.9237E-01  4.7238E-01  1.9973E+00 -3.3194E-01 -2.8949E-01  8.8307E-01  2.5489E+00
 GRADIENT:   1.4960E+00 -6.9432E+00 -1.1525E-02  2.7285E+00  1.2298E+01 -3.5514E-02 -2.5144E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -549.459480795953        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1859E+00  1.1880E+00  1.8722E+00  4.2156E-01  2.7466E-01  2.6461E-01  1.2311E+01
 PARAMETER:  2.7054E-01  2.7229E-01  7.2712E-01 -7.6379E-01 -1.1922E+00 -1.2295E+00  2.6105E+00
 GRADIENT:  -2.2387E+01  4.7871E+01 -1.3298E+01 -1.0896E+01 -3.4616E+00  1.0687E+00 -7.1120E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -557.636281281418        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      278
 NPARAMETR:  1.1684E+00  1.1471E+00  9.0111E+00  3.9211E-01  1.8563E-01  4.2797E-01  1.3609E+01
 PARAMETER:  2.5564E-01  2.3720E-01  2.2985E+00 -8.3622E-01 -1.5840E+00 -7.4870E-01  2.7108E+00
 GRADIENT:  -2.4312E+01  2.5608E+01 -8.9064E-01 -4.1317E+00 -1.4197E+00 -1.8567E-02  2.5032E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -558.768029506198        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      328
 NPARAMETR:  1.1718E+00  1.1130E+00  1.6128E+02  3.8847E-01  1.7844E-01  3.9039E+00  1.3594E+01
 PARAMETER:  2.5850E-01  2.0707E-01  5.1831E+00 -8.4555E-01 -1.6235E+00  1.4620E+00  2.7096E+00
 GRADIENT:   2.5692E+00 -2.5076E+00 -4.4461E-02 -1.0608E+00 -5.8791E-01 -4.9535E-03 -8.6504E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -558.867031417760        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  1.1742E+00  1.1222E+00  1.6901E+03  3.9366E-01  1.8957E-01  2.6043E+01  1.3573E+01
 PARAMETER:  2.6057E-01  2.1526E-01  7.5325E+00 -8.3228E-01 -1.5630E+00  3.3597E+00  2.7081E+00
 GRADIENT:   3.0806E-02  2.1741E-02 -2.0408E-03 -5.5738E-03 -1.1991E-02 -2.8307E-03 -6.3908E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -558.867900159420        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.1742E+00  1.1221E+00  2.1472E+03  3.9365E-01  1.8957E-01  3.1357E+01  1.3573E+01
 PARAMETER:  2.6056E-01  2.1524E-01  7.7719E+00 -8.3230E-01 -1.5630E+00  3.5455E+00  2.7081E+00
 GRADIENT:  -1.4755E+00 -6.1536E-01 -4.4259E-03 -7.6434E-01 -1.6037E-01  3.6261E-03 -2.2190E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -558.881933455032        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      560
 NPARAMETR:  1.1805E+00  1.1293E+00  2.0217E+03  3.9639E-01  1.9386E-01  3.0547E+01  1.3620E+01
 PARAMETER:  2.6596E-01  2.2159E-01  7.7117E+00 -8.2537E-01 -1.5406E+00  3.5193E+00  2.7115E+00
 GRADIENT:  -2.9944E-03  1.6768E-02 -1.9601E-02  7.5405E-03 -1.6034E-02  3.3884E-02  1.9289E-02

0ITERATION NO.:   46    OBJECTIVE VALUE:  -558.881933455032        NO. OF FUNC. EVALS.:  16
 CUMULATIVE NO. OF FUNC. EVALS.:      576
 NPARAMETR:  1.1805E+00  1.1293E+00  2.0217E+03  3.9639E-01  1.9386E-01  3.0547E+01  1.3620E+01
 PARAMETER:  2.6596E-01  2.2159E-01  7.7117E+00 -8.2537E-01 -1.5406E+00  3.5193E+00  2.7115E+00
 GRADIENT:   2.5213E+03  3.0249E+03 -4.3475E+01  5.4517E-04 -4.3473E+02  9.5672E+01 -2.4733E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      576
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3800E-02 -6.2547E-02 -3.9271E-04
 SE:             8.4613E-02  4.5287E-02  5.9576E-04
 N:                     100         100         100

 P VAL.:         7.7849E-01  1.6724E-01  5.0978E-01

 ETASHRINKSD(%)  1.0360E+01  5.2023E+01  9.9369E+01
 ETASHRINKVR(%)  1.9647E+01  7.6982E+01  9.9996E+01
 EBVSHRINKSD(%)  9.7450E+00  5.3096E+01  9.9623E+01
 EBVSHRINKVR(%)  1.8540E+01  7.8001E+01  9.9999E+01
 RELATIVEINF(%)  9.7954E+00  2.3400E+00  1.4184E-04
 EPSSHRINKSD(%)  5.1365E+00
 EPSSHRINKVR(%)  1.0009E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -558.88193345503214     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       176.26889310870604     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     4.78
 Elapsed covariance  time in seconds:     2.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -558.882       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.18E+00  1.13E+00  2.02E+03  3.96E-01  1.94E-01  3.05E+01  1.36E+01
 


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
 
         3.59E-01  7.67E-01  8.23E+01  1.37E+00  1.43E-03  5.64E-01  2.03E-01
 


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
+        1.29E-01
 
 TH 2
+        2.74E-01  5.88E-01
 
 TH 3
+       -2.95E+01 -6.31E+01  6.78E+03
 
 TH 4
+        4.92E-01  1.05E+00 -1.13E+02  1.89E+00
 
 TH 5
+       -5.11E-04 -1.09E-03  1.18E-01 -1.96E-03  2.04E-06
 
 TH 6
+        2.02E-01  4.32E-01 -4.64E+01  7.75E-01 -8.06E-04  3.18E-01
 
 TH 7
+       -7.20E-02 -1.54E-01  1.66E+01 -2.77E-01  2.88E-04 -1.14E-01  4.12E-02
 
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
+        3.59E-01
 
 TH 2
+        9.96E-01  7.67E-01
 
 TH 3
+       -9.97E-01 -1.00E+00  8.23E+01
 
 TH 4
+        9.98E-01  1.00E+00 -1.00E+00  1.37E+00
 
 TH 5
+       -9.97E-01 -1.00E+00  1.00E+00 -1.00E+00  1.43E-03
 
 TH 6
+        9.97E-01  1.00E+00 -1.00E+00  1.00E+00 -1.00E+00  5.64E-01
 
 TH 7
+       -9.88E-01 -9.93E-01  9.93E-01 -9.92E-01  9.93E-01 -9.93E-01  2.03E-01
 
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
+        7.11E+11
 
 TH 2
+        8.92E+11  1.12E+12
 
 TH 3
+       -5.69E+09 -7.14E+09  4.71E+07
 
 TH 4
+       -6.83E+11 -8.56E+11  5.46E+09  6.55E+11
 
 TH 5
+       -3.95E+10 -4.96E+10  4.21E+08  3.80E+10  2.00E+10
 
 TH 6
+       -8.31E+11 -1.04E+12  6.88E+09  7.98E+11  6.16E+10  1.01E+12
 
 TH 7
+        6.55E+06  8.22E+06 -3.29E+04 -6.29E+06 -2.63E+05 -4.81E+06  2.21E+03
 
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
 
 Elapsed finaloutput time in seconds:     0.10
 #CPUT: Total CPU Time in Seconds,        6.905
Stop Time:
Sat Sep 25 15:12:14 CDT 2021
