Sat Sep 25 10:51:05 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat95.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m95.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1605.90657311443        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0401E+01 -7.2484E+01 -3.6061E+01 -1.0298E+02  7.0578E+01  1.1083E+01 -1.7803E+01  1.6692E+01 -8.0694E+01  8.1813E-01
            -2.8085E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1623.92142521573        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0148E+00  9.7426E-01  1.0748E+00  1.0981E+00  9.7962E-01  9.7654E-01  1.0171E+00  8.4190E-01  1.3083E+00  9.4289E-01
             1.0345E+00
 PARAMETER:  1.1469E-01  7.3919E-02  1.7210E-01  1.9355E-01  7.9411E-02  7.6263E-02  1.1694E-01 -7.2091E-02  3.6871E-01  4.1191E-02
             1.3394E-01
 GRADIENT:   7.8257E+01  1.2243E+00 -7.6483E+00  2.5958E+01  3.0345E+01  9.1905E-01  2.7692E+00  9.6818E+00  9.7616E+00 -1.0652E+01
            -1.1426E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1629.19923062288        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0136E+00  9.1722E-01  1.0166E+00  1.1155E+00  9.3205E-01  9.7822E-01  9.0980E-01  3.3730E-01  1.2851E+00  1.1299E+00
             1.0497E+00
 PARAMETER:  1.1347E-01  1.3595E-02  1.1642E-01  2.0932E-01  2.9633E-02  7.7982E-02  5.4715E-03 -9.8678E-01  3.5082E-01  2.2215E-01
             1.4848E-01
 GRADIENT:   7.3584E+01 -3.2963E-01 -1.2755E+01  1.0610E+01  1.6275E+00  2.1854E+00  2.9117E+00  2.4423E+00  7.6999E+00  1.3876E+01
             3.9130E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1630.46478369371        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.8331E-01  8.5266E-01  9.2464E-01  1.1508E+00  8.5029E-01  9.7321E-01  1.2964E+00  2.4371E-01  1.1406E+00  9.2260E-01
             1.0427E+00
 PARAMETER:  8.3166E-02 -5.9391E-02  2.1646E-02  2.4048E-01 -6.2176E-02  7.2846E-02  3.5958E-01 -1.3118E+00  2.3155E-01  1.9444E-02
             1.4186E-01
 GRADIENT:  -2.0571E+00  5.8163E+00 -3.6139E+00  1.2672E+01  1.5954E+00  1.2517E-01  6.7534E-01  1.3948E+00 -4.8656E+00  3.9695E+00
            -1.3382E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1630.73978980794        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  9.8452E-01  8.7145E-01  8.7371E-01  1.1150E+00  8.3160E-01  9.7209E-01  1.2797E+00  1.0966E-01  1.1820E+00  8.4937E-01
             1.0488E+00
 PARAMETER:  8.4394E-02 -3.7594E-02 -3.5010E-02  2.0885E-01 -8.4400E-02  7.1693E-02  3.4659E-01 -2.1104E+00  2.6718E-01 -6.3255E-02
             1.4765E-01
 GRADIENT:  -7.4172E-01 -6.3533E+00 -2.9952E+00 -7.4869E+00  3.8401E+00 -5.9183E-01 -1.1343E+00  2.8527E-01  2.2094E+00 -8.4796E-01
             1.6059E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1631.72397641008        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      427
 NPARAMETR:  9.9428E-01  8.1252E-01  9.0528E-01  1.1689E+00  8.2435E-01  9.8023E-01  1.3752E+00  2.8470E-02  1.1393E+00  8.7945E-01
             1.0461E+00
 PARAMETER:  9.4260E-02 -1.0761E-01  4.8977E-04  2.5608E-01 -9.3163E-02  8.0036E-02  4.1862E-01 -3.4589E+00  2.3038E-01 -2.8463E-02
             1.4502E-01
 GRADIENT:  -1.1740E+01  1.5364E+00  1.0182E+00 -7.5273E+00 -1.0930E+00  4.7147E-01 -4.8474E-01  1.8262E-02 -1.4573E+00  1.8811E-01
             8.6580E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1632.86915067753        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  9.9476E-01  5.2862E-01  8.9670E-01  1.3378E+00  7.2329E-01  9.7084E-01  1.9229E+00  1.0000E-02  1.0222E+00  8.6148E-01
             1.0387E+00
 PARAMETER:  9.4750E-02 -5.3749E-01 -9.0388E-03  3.9102E-01 -2.2394E-01  7.0404E-02  7.5386E-01 -5.0213E+00  1.2192E-01 -4.9100E-02
             1.3794E-01
 GRADIENT:  -4.1009E+00  1.7725E+00  5.9284E+00 -1.2502E+00 -9.3383E+00 -1.4477E+00 -1.0120E+00  0.0000E+00  2.0875E+00 -1.8818E-01
             3.6644E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1633.05440863997        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  9.9456E-01  4.1454E-01  9.0601E-01  1.4037E+00  7.0137E-01  9.7183E-01  2.2860E+00  1.0000E-02  9.7523E-01  8.8441E-01
             1.0347E+00
 PARAMETER:  9.4547E-02 -7.8059E-01  1.2951E-03  4.3912E-01 -2.5472E-01  7.1422E-02  9.2682E-01 -5.5399E+00  7.4922E-02 -2.2838E-02
             1.3407E-01
 GRADIENT:  -1.2879E-01  1.2139E-01 -8.4399E-02  3.9862E-01  1.3745E-01  8.0815E-03  4.3005E-02  0.0000E+00 -1.3954E-01 -3.9595E-02
             5.4233E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1633.05459095123        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      873
 NPARAMETR:  9.9455E-01  4.1084E-01  9.0615E-01  1.4055E+00  7.0047E-01  9.7173E-01  2.2970E+00  1.0000E-02  9.7455E-01  8.8541E-01
             1.0344E+00
 PARAMETER:  9.4534E-02 -7.8954E-01  1.4459E-03  4.4039E-01 -2.5600E-01  7.1319E-02  9.3160E-01 -5.5484E+00  7.4225E-02 -2.1706E-02
             1.3381E-01
 GRADIENT:  -5.6720E-04  1.9993E-03  1.1040E-02 -9.0515E-03 -1.8572E-02  1.4528E-03  1.6360E-03  0.0000E+00 -3.6610E-03 -4.3642E-03
             2.1135E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      873
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8437E-04  1.9740E-02 -4.4064E-04 -1.3985E-02 -5.2401E-03
 SE:             2.9816E-02  1.6682E-02  2.0507E-04  2.6886E-02  2.3784E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9239E-01  2.3668E-01  3.1653E-02  6.0294E-01  8.2562E-01

 ETASHRINKSD(%)  1.1200E-01  4.4115E+01  9.9313E+01  9.9293E+00  2.0322E+01
 ETASHRINKVR(%)  2.2388E-01  6.8768E+01  9.9995E+01  1.8873E+01  3.6514E+01
 EBVSHRINKSD(%)  4.9510E-01  4.9632E+01  9.9264E+01  8.7556E+00  1.6076E+01
 EBVSHRINKVR(%)  9.8775E-01  7.4630E+01  9.9995E+01  1.6745E+01  2.9568E+01
 RELATIVEINF(%)  9.8254E+01  3.8317E+00  7.3007E-04  1.9377E+01  6.2573E+00
 EPSSHRINKSD(%)  4.3390E+01
 EPSSHRINKVR(%)  6.7953E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1633.0545909512282     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -897.90376438749001     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.43
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.27
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1633.055       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  4.11E-01  9.06E-01  1.41E+00  7.00E-01  9.72E-01  2.30E+00  1.00E-02  9.75E-01  8.85E-01  1.03E+00
 


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
+        1.18E+03
 
 TH 2
+       -2.09E+01  3.70E+02
 
 TH 3
+        1.63E+01  2.12E+02  7.72E+02
 
 TH 4
+       -4.14E+00  2.34E+02 -1.23E+02  4.97E+02
 
 TH 5
+       -3.02E+00 -4.66E+02 -1.12E+03  1.18E+02  2.02E+03
 
 TH 6
+       -3.06E+00 -4.72E+00  5.36E+00 -9.61E-01 -1.49E+00  2.08E+02
 
 TH 7
+        1.09E+00  2.73E+01  8.76E-01 -1.59E+00 -5.46E+00  1.77E-01  7.78E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.48E+00 -3.29E+01 -3.98E+00  1.05E+01 -1.10E+01 -1.65E+00  2.31E+00  0.00E+00  1.60E+02
 
 TH10
+       -5.42E-01  2.75E+01 -7.99E+01 -2.71E+01 -1.33E+01 -3.13E+00  6.28E+00  0.00E+00  2.41E+00  1.12E+02
 
 TH11
+       -1.18E+01 -8.33E+00 -3.92E+01 -7.61E+00  2.27E+01  1.91E+00  1.32E+00  0.00E+00  7.50E+00  2.53E+01  2.05E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.772
Stop Time:
Sat Sep 25 10:51:22 CDT 2021
