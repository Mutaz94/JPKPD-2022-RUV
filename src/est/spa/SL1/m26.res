Sat Sep 25 10:25:21 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat26.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1710.45217098248        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6065E+01 -2.7542E+01 -1.6954E+01 -4.4856E+00 -5.6788E-01  2.5658E+01 -5.1427E+00  6.3268E+00  2.3371E+01 -5.3180E+00
             6.4567E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1719.62035470357        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:      118
 NPARAMETR:  9.9220E-01  1.0504E+00  1.0256E+00  9.6927E-01  1.0437E+00  9.5429E-01  1.0427E+00  9.8295E-01  9.3127E-01  1.0613E+00
             8.7904E-01
 PARAMETER:  9.2174E-02  1.4914E-01  1.2527E-01  6.8787E-02  1.4279E-01  5.3213E-02  1.4183E-01  8.2799E-02  2.8795E-02  1.5946E-01
            -2.8929E-02
 GRADIENT:  -1.3551E+01 -2.4421E+01 -8.4778E+00 -1.5522E+01 -5.7750E-01  2.6768E+00 -5.2932E+00  3.4036E+00  7.7523E+00 -4.5278E+00
             2.1551E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1719.92719807544        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.9093E-01  1.2175E+00  1.0723E+00  8.8414E-01  1.1558E+00  9.6099E-01  8.6802E-01  8.9618E-01  1.0421E+00  1.2045E+00
             8.7402E-01
 PARAMETER:  9.0886E-02  2.9681E-01  1.6980E-01 -2.3137E-02  2.4475E-01  6.0214E-02 -4.1537E-02 -9.6187E-03  1.4120E-01  2.8603E-01
            -3.4654E-02
 GRADIENT:  -1.6901E+01  1.1585E+00  2.5056E+00  9.1073E+00 -8.3816E-01  5.1459E+00 -9.0803E+00 -2.0974E+00  6.5447E+00 -2.7507E+00
             1.3903E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1721.63617979969        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      481
 NPARAMETR:  9.9486E-01  1.2367E+00  1.0124E+00  8.7271E-01  1.1346E+00  9.5245E-01  9.9048E-01  9.2855E-01  9.7927E-01  1.1731E+00
             8.3858E-01
 PARAMETER:  9.4851E-02  3.1247E-01  1.1229E-01 -3.6156E-02  2.2625E-01  5.1286E-02  9.0438E-02  2.5872E-02  7.9053E-02  2.5962E-01
            -7.6040E-02
 GRADIENT:  -8.0681E+00  8.1221E+00  2.6995E+00  1.0318E+01 -4.0255E+00  1.4691E+00 -1.4009E+00 -9.8557E-01  4.2355E+00 -1.2409E+00
             1.5604E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1721.70871952386        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:      676             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9674E-01  1.2400E+00  1.0053E+00  8.7011E-01  1.1353E+00  9.5045E-01  9.9767E-01  9.3345E-01  9.7530E-01  1.1738E+00
             8.3634E-01
 PARAMETER:  9.6738E-02  3.1511E-01  1.0531E-01 -3.9136E-02  2.2694E-01  4.9178E-02  9.7667E-02  3.1137E-02  7.4992E-02  2.6021E-01
            -7.8722E-02
 GRADIENT:   5.8107E+01  3.0744E+01  1.4863E+00  1.8149E+01 -2.6794E-01  8.1498E+00 -2.2308E-02 -6.7109E-01  4.6575E+00 -5.1250E-01
             8.8308E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1721.71793298913        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      855
 NPARAMETR:  9.9824E-01  1.2400E+00  1.0008E+00  8.7011E-01  1.1353E+00  9.4893E-01  1.0036E+00  9.3345E-01  9.7530E-01  1.1738E+00
             8.3425E-01
 PARAMETER:  9.8242E-02  3.1511E-01  1.0076E-01 -3.9136E-02  2.2694E-01  4.7576E-02  1.0354E-01  3.1137E-02  7.4992E-02  2.6021E-01
            -8.1221E-02
 GRADIENT:   5.0370E-02  6.4462E+00 -5.7937E-03  1.1171E+01 -3.4115E-01 -1.0995E-02 -7.8924E-03 -4.7011E-01  4.4827E+00 -6.4998E-01
             1.0863E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1721.71894261293        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  9.9817E-01  1.2400E+00  1.0008E+00  8.7011E-01  1.1353E+00  9.4894E-01  1.0036E+00  9.3555E-01  9.7530E-01  1.1738E+00
             8.3423E-01
 PARAMETER:  9.8173E-02  3.1511E-01  1.0077E-01 -3.9136E-02  2.2694E-01  4.7589E-02  1.0363E-01  3.3381E-02  7.4991E-02  2.6021E-01
            -8.1242E-02
 GRADIENT:  -1.2114E-01  6.3840E+00 -1.1513E-01  1.1199E+01 -3.2618E-01 -6.3447E-03  1.3609E-02 -4.2715E-01  4.5037E+00 -6.0750E-01
             3.9429E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1721.71894261293        NO. OF FUNC. EVALS.:  32
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  9.9818E-01  1.2400E+00  1.0009E+00  8.7012E-01  1.1353E+00  9.4894E-01  1.0036E+00  9.3554E-01  9.7531E-01  1.1737E+00
             8.3422E-01
 PARAMETER:  9.8173E-02  3.1511E-01  1.0077E-01 -3.9136E-02  2.2694E-01  4.7589E-02  1.0363E-01  3.3381E-02  7.4991E-02  2.6021E-01
            -8.1242E-02
 GRADIENT:  -1.7286E+06 -1.0971E+06 -1.0410E-01 -3.4572E+06  1.5234E+06 -1.0504E-02  1.1378E-02  3.4571E+06 -3.4572E+06  1.3286E+06
             2.9329E-02
 NUMSIGDIG:         3.3         3.3         2.2         3.3         3.3         3.6         2.8         3.3         3.3         3.3
                    2.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1060
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4038E-04 -1.9704E-02 -3.2566E-02  4.5769E-03 -3.5778E-02
 SE:             2.9879E-02  2.1151E-02  1.1758E-02  2.2523E-02  2.2903E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9625E-01  3.5155E-01  5.6100E-03  8.3897E-01  1.1825E-01

 ETASHRINKSD(%)  1.0000E-10  2.9143E+01  6.0610E+01  2.4544E+01  2.3273E+01
 ETASHRINKVR(%)  1.0000E-10  4.9793E+01  8.4485E+01  4.3064E+01  4.1130E+01
 EBVSHRINKSD(%)  3.3505E-01  2.8737E+01  6.6206E+01  2.4760E+01  1.9173E+01
 EBVSHRINKVR(%)  6.6897E-01  4.9216E+01  8.8580E+01  4.3389E+01  3.4671E+01
 RELATIVEINF(%)  9.8635E+01  1.4134E+00  1.1007E+00  1.6224E+00  1.3304E+01
 EPSSHRINKSD(%)  4.5818E+01
 EPSSHRINKVR(%)  7.0644E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1721.7189426129294     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -986.56811604919119     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.09
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1721.719       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.98E-01  1.24E+00  1.00E+00  8.70E-01  1.14E+00  9.49E-01  1.00E+00  9.36E-01  9.75E-01  1.17E+00  8.34E-01
 


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
+        8.67E+09
 
 TH 2
+       -2.90E+03  1.13E+09
 
 TH 3
+       -1.91E+05 -4.87E+04  1.61E+02
 
 TH 4
+       -1.30E+04  2.59E+04 -2.19E+05  1.14E+10
 
 TH 5
+        4.39E+03 -9.05E+03  7.39E+04 -3.99E+04  1.30E+09
 
 TH 6
+       -2.00E+04 -5.11E+03  7.55E-01 -2.30E+04  7.75E+03  2.22E+02
 
 TH 7
+       -6.73E+03 -1.71E+03  2.38E+00 -7.74E+03  2.60E+03  9.47E-01  6.25E+01
 
 TH 8
+       -9.25E+09  1.58E+05  2.04E+05  7.10E+05 -2.40E+05  2.13E+04  7.18E+03  9.87E+09
 
 TH 9
+       -1.16E+04  2.33E+04 -1.96E+05  1.05E+05 -3.57E+04 -2.05E+04 -6.86E+03  6.34E+05  9.09E+09
 
 TH10
+        3.71E+03  1.74E+03  6.24E+04 -3.26E+04  1.14E+04  6.54E+03  2.21E+03 -2.02E+05 -2.99E+04  9.27E+08
 
 TH11
+        3.31E+06  8.46E+05 -2.28E+01  3.80E+06 -1.28E+06 -1.89E+00  1.20E+01 -3.53E+06  3.39E+06 -1.08E+06  2.93E+02
 
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
 #CPUT: Total CPU Time in Seconds,       20.263
Stop Time:
Sat Sep 25 10:25:43 CDT 2021
