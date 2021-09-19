Sat Sep 18 14:03:43 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat43.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1663.36474061990        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.1421E+01 -2.6335E+01 -9.9453E+00 -1.0461E+01  3.7781E+01  8.7795E+00  7.9384E+00  6.5863E+00  2.7896E+01 -2.9084E+00
             2.9486E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1665.55379344958        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9084E-01  9.9247E-01  9.2157E-01  1.0021E+00  9.5079E-01  9.5861E-01  9.2680E-01  9.2956E-01  7.8685E-01  9.9086E-01
             9.9173E-01
 PARAMETER:  9.0800E-02  9.2437E-02  1.8318E-02  1.0210E-01  4.9541E-02  5.7728E-02  2.3982E-02  2.6955E-02 -1.3972E-01  9.0817E-02
             9.1691E-02
 GRADIENT:   5.0708E+01 -2.3272E+01 -9.8998E+00 -2.6786E+01  2.8686E+01 -7.3240E+00 -9.5607E+00  6.5508E+00 -1.8912E+01  2.6388E+00
            -3.2854E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1670.86629507141        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0106E+00  9.5232E-01  4.7297E-01  1.0028E+00  6.3077E-01  1.0035E+00  1.0045E+00  3.4600E-01  7.7412E-01  5.3077E-01
             1.0025E+00
 PARAMETER:  1.1058E-01  5.1145E-02 -6.4873E-01  1.0276E-01 -3.6082E-01  1.0346E-01  1.0451E-01 -9.6132E-01 -1.5603E-01 -5.3342E-01
             1.0251E-01
 GRADIENT:   8.5125E+01  2.7581E+01 -9.5450E+00  7.2224E+01  2.6370E+01  8.3353E+00 -1.6116E+01  1.5703E+00 -5.4483E+00 -7.1424E+00
             2.8406E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1675.29948133310        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      280
 NPARAMETR:  9.8796E-01  8.0187E-01  5.3452E-01  1.0772E+00  6.1011E-01  9.9080E-01  1.2654E+00  2.2231E-01  7.3209E-01  5.8155E-01
             9.9198E-01
 PARAMETER:  8.7891E-02 -1.2081E-01 -5.2639E-01  1.7435E-01 -3.9411E-01  9.0753E-02  3.3542E-01 -1.4037E+00 -2.1185E-01 -4.4206E-01
             9.1950E-02
 GRADIENT:  -5.0885E+00  1.8131E+01  1.2053E+01 -4.3146E-01 -1.4046E+01  2.8494E-01 -8.9025E-01  6.3847E-01 -5.6947E+00 -2.6950E+00
             1.0674E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1677.62592649985        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      455
 NPARAMETR:  9.8968E-01  6.0327E-01  5.3262E-01  1.1735E+00  5.4548E-01  9.8670E-01  1.5757E+00  6.1033E-02  7.0004E-01  5.8580E-01
             9.8669E-01
 PARAMETER:  8.9630E-02 -4.0539E-01 -5.2995E-01  2.5997E-01 -5.0610E-01  8.6611E-02  5.5468E-01 -2.6963E+00 -2.5661E-01 -4.3478E-01
             8.6596E-02
 GRADIENT:   6.4824E-01  3.9363E+00  7.5045E+00 -7.8583E+00 -9.3531E+00 -1.0198E+00  1.2264E-01  4.5204E-02  5.6561E-01 -2.1378E+00
            -1.9516E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1678.66990147995        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      630
 NPARAMETR:  9.8362E-01  4.0408E-01  6.3657E-01  1.3090E+00  5.5975E-01  9.8332E-01  2.0893E+00  1.0000E-02  6.5950E-01  7.3086E-01
             9.8777E-01
 PARAMETER:  8.3487E-02 -8.0615E-01 -3.5166E-01  3.6929E-01 -4.8026E-01  8.3180E-02  8.3685E-01 -4.9618E+00 -3.1628E-01 -2.1354E-01
             8.7690E-02
 GRADIENT:  -3.3338E+00  2.5011E+00 -3.8431E+00  1.0940E+01  2.6116E+00 -9.4887E-01  4.0508E-01  0.0000E+00  1.9719E-01 -9.7389E-01
            -5.0753E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1679.07635522567        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      806
 NPARAMETR:  9.7951E-01  2.4132E-01  6.9957E-01  1.4069E+00  5.6126E-01  9.8033E-01  2.8214E+00  1.0000E-02  6.3248E-01  8.1685E-01
             9.9438E-01
 PARAMETER:  7.9299E-02 -1.3217E+00 -2.5729E-01  4.4142E-01 -4.7757E-01  8.0139E-02  1.1372E+00 -8.5201E+00 -3.5810E-01 -1.0229E-01
             9.4367E-02
 GRADIENT:  -1.3969E+00  2.3410E+00  2.2800E+00  1.2589E+01 -5.3571E+00 -3.0892E-01  3.6282E-01  0.0000E+00 -7.6912E-01  8.8065E-02
            -4.4700E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1679.17278654874        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      982
 NPARAMETR:  9.7875E-01  1.9445E-01  7.0854E-01  1.4275E+00  5.5866E-01  9.8020E-01  3.1767E+00  1.0000E-02  6.2700E-01  8.2807E-01
             9.9644E-01
 PARAMETER:  7.8525E-02 -1.5376E+00 -2.4455E-01  4.5594E-01 -4.8221E-01  8.0000E-02  1.2558E+00 -1.0153E+01 -3.6680E-01 -8.8660E-02
             9.6432E-02
 GRADIENT:   2.5125E-01  7.6843E-02  8.6569E-01 -1.9403E-01 -1.1858E+00  9.3586E-02 -1.3613E-02  0.0000E+00 -2.5470E-02  7.7200E-02
             3.3033E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1679.17312642635        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  9.7861E-01  1.9330E-01  7.0873E-01  1.4281E+00  5.5882E-01  9.7994E-01  3.1874E+00  1.0000E-02  6.2688E-01  8.2823E-01
             9.9638E-01
 PARAMETER:  7.8376E-02 -1.5435E+00 -2.4427E-01  4.5637E-01 -4.8193E-01  7.9741E-02  1.2592E+00 -1.0200E+01 -3.6700E-01 -8.8469E-02
             9.6373E-02
 GRADIENT:  -6.3985E-03  5.8456E-04 -1.3494E-02 -8.0920E-02  7.1963E-03  8.6880E-04  9.7908E-03  0.0000E+00  1.3656E-02  1.0174E-02
             3.3241E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1074
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.6511E-04  3.7291E-02 -4.0853E-04 -2.4918E-02  5.9474E-03
 SE:             2.9867E-02  1.6604E-02  2.6406E-04  2.5991E-02  2.4901E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7689E-01  2.4707E-02  1.2184E-01  3.3769E-01  8.1123E-01

 ETASHRINKSD(%)  1.0000E-10  4.4375E+01  9.9115E+01  1.2928E+01  1.6577E+01
 ETASHRINKVR(%)  1.0000E-10  6.9059E+01  9.9992E+01  2.4184E+01  3.0407E+01
 EBVSHRINKSD(%)  4.3014E-01  5.3977E+01  9.9158E+01  1.0106E+01  1.2104E+01
 EBVSHRINKVR(%)  8.5843E-01  7.8819E+01  9.9993E+01  1.9191E+01  2.2743E+01
 RELATIVEINF(%)  9.8322E+01  5.2119E+00  5.1080E-04  2.2585E+01  5.4203E+00
 EPSSHRINKSD(%)  4.3935E+01
 EPSSHRINKVR(%)  6.8567E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1679.1731264263469     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -944.02229986260875     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1679.173       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.93E-01  7.09E-01  1.43E+00  5.59E-01  9.80E-01  3.19E+00  1.00E-02  6.27E-01  8.28E-01  9.96E-01
 


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
+        1.20E+03
 
 TH 2
+       -3.26E+01  8.87E+02
 
 TH 3
+        1.44E+01  2.68E+02  1.64E+03
 
 TH 4
+       -9.84E+00  4.39E+02 -3.31E+02  1.18E+03
 
 TH 5
+        3.08E-01 -6.33E+02 -2.40E+03  1.02E+02  4.10E+03
 
 TH 6
+       -1.33E+00 -5.50E-01  4.66E-01 -4.18E+00 -2.55E-02  2.10E+02
 
 TH 7
+        1.09E+00  5.57E+01 -3.85E+00 -1.20E+01  1.12E+01  4.57E-01  7.06E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.41E+00 -3.31E+01 -2.98E+01 -1.41E+01  4.10E+01 -4.64E+00  3.91E+00  0.00E+00  3.27E+02
 
 TH10
+       -5.88E+00 -1.73E+01 -8.06E+01 -6.71E+00 -6.78E+01 -3.75E-01 -5.49E+00  0.00E+00  1.06E+01  1.90E+02
 
 TH11
+       -5.13E+00 -1.14E+01 -3.83E+01 -1.25E+01  2.98E+01  6.26E+00 -6.54E-01  0.00E+00  2.00E+01  2.15E+01  2.17E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.590
Stop Time:
Sat Sep 18 14:04:03 CDT 2021
