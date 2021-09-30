Wed Sep 29 13:34:05 CDT 2021
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
$DATA ../../../../data/spa/A3/dat46.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 NO. OF SIG. FIGURES REQUIRED:            2
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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   89.4676758156713        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7146E+02  2.6447E+01  1.3866E+02 -1.3267E+02  1.6023E+02  2.7498E+01 -6.3153E+01 -5.0110E+01 -1.4473E+02 -1.9582E+02
            -2.9346E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1176.31175874684        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9386E-01  1.0121E+00  8.4808E-01  1.2765E+00  8.9915E-01  8.5818E-01  1.0088E+00  1.0140E+00  1.0617E+00  1.1768E+00
             5.3260E+00
 PARAMETER:  9.3837E-02  1.1206E-01 -6.4779E-02  3.4411E-01 -6.3065E-03 -5.2936E-02  1.0877E-01  1.1388E-01  1.5987E-01  2.6283E-01
             1.7726E+00
 GRADIENT:  -2.9483E+01  5.7758E+00 -1.8874E+01  3.8347E+01 -1.4240E+01 -7.3921E+00  1.4713E+01  8.5999E+00  3.5559E+01  2.9209E+01
             1.8903E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1194.16655713033        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.8838E-01  7.1650E-01  3.5433E-01  1.3539E+00  3.8611E-01  9.2326E-01  4.3394E-01  2.8337E-01  9.1648E-01  6.0746E-01
             4.8618E+00
 PARAMETER:  8.8309E-02 -2.3338E-01 -9.3753E-01  4.0297E-01 -8.5163E-01  2.0156E-02 -7.3486E-01 -1.1610E+00  1.2787E-02 -3.9848E-01
             1.6814E+00
 GRADIENT:  -5.4847E+01  1.0958E+02  6.4888E+01  1.2841E+02 -1.4685E+02 -6.8294E+00  7.6096E-01  8.8399E-01  9.8812E+00  1.2493E+01
             1.2956E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1238.43947732427        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  9.5733E-01  6.0631E-01  1.8511E-01  1.1047E+00  2.9623E-01  9.1193E-01  1.5611E-01  1.0000E-02  1.2974E+00  2.6799E-01
             3.3590E+00
 PARAMETER:  5.6397E-02 -4.0036E-01 -1.5868E+00  1.9961E-01 -1.1166E+00  7.8087E-03 -1.7572E+00 -4.8179E+00  3.6036E-01 -1.2168E+00
             1.3116E+00
 GRADIENT:   2.1377E+01  3.8527E+01  4.0846E+01  2.3757E+01  9.0362E+00 -2.4333E+01 -9.7592E-01  0.0000E+00 -9.0454E-01 -9.8500E+00
            -4.9083E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1242.27866382298        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      390
 NPARAMETR:  9.6645E-01  6.1950E-01  1.7906E-01  1.0785E+00  2.9334E-01  1.0017E+00  1.4112E-01  1.0000E-02  1.2404E+00  3.1835E-01
             3.5849E+00
 PARAMETER:  6.5871E-02 -3.7884E-01 -1.6200E+00  1.7553E-01 -1.1264E+00  1.0169E-01 -1.8581E+00 -5.4554E+00  3.1543E-01 -1.0446E+00
             1.3767E+00
 GRADIENT:   9.2170E+00  3.5536E+01  2.9681E+01 -5.2114E+00 -4.1388E+01  7.8757E+00 -4.1111E-01  0.0000E+00 -8.0271E+00 -6.0260E+00
            -1.3060E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1249.99758064635        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      566
 NPARAMETR:  9.4689E-01  4.7746E-01  1.2698E-01  1.0393E+00  2.2223E-01  1.0004E+00  1.2659E-01  1.0000E-02  1.5330E+00  5.8277E-01
             3.1702E+00
 PARAMETER:  4.5424E-02 -6.3928E-01 -1.9637E+00  1.3857E-01 -1.4040E+00  1.0045E-01 -1.9668E+00 -8.8000E+00  5.2723E-01 -4.3997E-01
             1.2538E+00
 GRADIENT:  -3.5873E+00 -6.4846E+00 -4.4467E+00 -1.1413E+00  8.6768E+00  4.1288E-02  1.9410E-01  0.0000E+00  4.9291E-01  2.8867E-01
             4.2659E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1250.12801215910        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      748
 NPARAMETR:  9.4909E-01  4.8597E-01  1.2957E-01  1.0452E+00  2.2532E-01  9.9945E-01  5.5229E-02  1.0000E-02  1.5161E+00  5.8410E-01
             3.1563E+00
 PARAMETER:  4.7753E-02 -6.2160E-01 -1.9435E+00  1.4425E-01 -1.3902E+00  9.9452E-02 -2.7963E+00 -8.6141E+00  5.1617E-01 -4.3768E-01
             1.2494E+00
 GRADIENT:   1.8783E-01  8.3061E-02 -2.2194E-01  1.5734E+00 -6.7689E-01  7.5156E-02  3.3454E-02  0.0000E+00 -7.0532E-01  2.4460E-02
            -9.4210E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1250.14734367877        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      923
 NPARAMETR:  9.4843E-01  4.8696E-01  1.2901E-01  1.0421E+00  2.2520E-01  9.9937E-01  1.3162E-02  1.0000E-02  1.5235E+00  5.8565E-01
             3.1547E+00
 PARAMETER:  4.7056E-02 -6.1957E-01 -1.9478E+00  1.4122E-01 -1.3908E+00  9.9375E-02 -4.2304E+00 -8.6141E+00  5.2100E-01 -4.3503E-01
             1.2489E+00
 GRADIENT:  -3.8548E-01  3.1820E-01 -1.5183E-01  3.1915E-02 -6.1887E-01  7.3049E-04  1.9632E-03  0.0000E+00 -3.2488E-02  5.8648E-02
             1.6971E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1250.14787082213        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1015
 NPARAMETR:  9.4856E-01  4.8660E-01  1.2903E-01  1.0420E+00  2.2516E-01  9.9933E-01  1.0000E-02  1.0000E-02  1.5237E+00  5.8546E-01
             3.1531E+00
 PARAMETER:  4.7189E-02 -6.2032E-01 -1.9477E+00  1.4110E-01 -1.3909E+00  9.9328E-02 -4.5411E+00 -8.6141E+00  5.2115E-01 -4.3535E-01
             1.2484E+00
 GRADIENT:  -8.3065E-02 -6.5850E-02 -2.4227E-01 -1.2563E-01 -8.4174E-02 -1.3917E-02  1.2946E-04  0.0000E+00  1.2643E-02 -1.6550E-02
            -1.6484E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1015
 NO. OF SIG. DIGITS IN FINAL EST.:  2.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7958E-03 -9.7606E-05  2.5757E-04 -1.7555E-02  8.2876E-03
 SE:             2.8555E-02  1.9500E-04  1.8610E-04  2.5632E-02  2.1340E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4985E-01  6.1669E-01  1.6635E-01  4.9342E-01  6.9775E-01

 ETASHRINKSD(%)  4.3364E+00  9.9347E+01  9.9377E+01  1.4128E+01  2.8508E+01
 ETASHRINKVR(%)  8.4848E+00  9.9996E+01  9.9996E+01  2.6260E+01  4.8889E+01
 EBVSHRINKSD(%)  3.6366E+00  9.9310E+01  9.9477E+01  1.1095E+01  2.8689E+01
 EBVSHRINKVR(%)  7.1410E+00  9.9995E+01  9.9997E+01  2.0960E+01  4.9147E+01
 RELATIVEINF(%)  8.7075E+01  1.9847E-04  4.1554E-04  6.3832E+01  1.2584E+00
 EPSSHRINKSD(%)  3.0161E+01
 EPSSHRINKVR(%)  5.1226E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1250.1478708221327     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -514.99704425839457     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.21
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.92
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1250.148       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.49E-01  4.87E-01  1.29E-01  1.04E+00  2.25E-01  9.99E-01  1.00E-02  1.00E-02  1.52E+00  5.85E-01  3.15E+00
 


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
+        1.15E+03
 
 TH 2
+       -2.21E+01  2.26E+03
 
 TH 3
+       -7.09E+02  4.37E+03  1.86E+04
 
 TH 4
+       -7.31E+00  8.46E+01 -5.74E+02  3.52E+02
 
 TH 5
+        4.77E+02 -8.06E+03 -1.99E+04 -7.43E+01  3.27E+04
 
 TH 6
+        4.83E+00 -7.83E+00  1.02E+02 -1.42E+01  1.54E+01  1.70E+02
 
 TH 7
+       -8.49E-02 -5.17E-02 -3.01E-01  4.47E-02  1.71E-01 -9.44E-02  4.17E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.34E+01 -2.52E+01  2.25E+02 -9.00E+00  6.34E+01 -1.42E+00 -1.15E-02  0.00E+00  5.06E+01
 
 TH10
+       -8.82E+00 -7.43E+01 -1.45E+02  5.45E+00  4.13E+02  9.26E+00  5.48E-02  0.00E+00 -2.18E-01  1.35E+02
 
 TH11
+       -1.93E+01 -1.23E+01 -2.36E+01 -4.76E+00  2.69E+01  2.20E+00  1.24E-02  0.00E+00  5.60E+00  2.54E+01  3.06E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.188
Stop Time:
Wed Sep 29 13:34:27 CDT 2021
