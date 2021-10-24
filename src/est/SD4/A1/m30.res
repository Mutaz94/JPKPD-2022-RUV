Sun Oct 24 01:57:22 CDT 2021
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
$DATA ../../../../data/SD4/A1/dat30.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1357.12203843105        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3325E+02  3.8284E+01  5.4362E+01  5.0269E+01  4.3347E+00  5.1383E+01 -1.1047E+01 -2.4512E+01  2.0979E+01 -1.3271E+01
            -5.9312E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1477.35406808236        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1331E+00  1.1115E+00  8.9813E-01  9.6541E-01  1.0047E+00  1.2243E+00  1.0367E+00  1.0472E+00  7.3653E-01  9.1969E-01
             2.7344E+00
 PARAMETER:  2.2496E-01  2.0568E-01 -7.4433E-03  6.4803E-02  1.0464E-01  3.0240E-01  1.3601E-01  1.4615E-01 -2.0580E-01  1.6284E-02
             1.1059E+00
 GRADIENT:   2.7956E+02 -2.1845E+00 -2.7795E+00 -6.6917E+00 -1.4530E+01  6.6787E+01  2.8587E+00  5.0774E+00  4.1348E+00  1.1918E+01
             1.4593E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1496.78211374643        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0595E+00  8.6488E-01  3.0941E-01  1.0288E+00  4.8442E-01  1.0778E+00  1.2240E+00  6.2320E-01  6.3776E-01  3.4735E-01
             2.2039E+00
 PARAMETER:  1.5776E-01 -4.5165E-02 -1.0731E+00  1.2840E-01 -6.2480E-01  1.7495E-01  3.0211E-01 -3.7288E-01 -3.4979E-01 -9.5742E-01
             8.9023E-01
 GRADIENT:   1.8416E+02  5.9183E+01 -3.4771E+01  1.8291E+02  2.8266E+01  2.7753E+01  2.3094E+00  5.4688E+00 -1.5379E+01  2.1831E+00
             7.0251E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1514.41186483763        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      293
 NPARAMETR:  1.0203E+00  9.1016E-01  4.8967E-01  1.0037E+00  6.4818E-01  1.0434E+00  1.2048E+00  5.7360E-01  7.4618E-01  7.0723E-01
             1.9645E+00
 PARAMETER:  1.2014E-01  5.8668E-03 -6.1402E-01  1.0373E-01 -3.3359E-01  1.4245E-01  2.8630E-01 -4.5582E-01 -1.9278E-01 -2.4640E-01
             7.7524E-01
 GRADIENT:  -2.6777E+00 -1.6411E+01 -9.1201E+00 -2.3769E+01 -9.9196E+00  6.4347E+00 -1.9744E+00  5.0905E+00  2.2169E+00  1.2913E+01
             2.7858E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1518.22644626932        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      468
 NPARAMETR:  1.0209E+00  9.3367E-01  6.8527E-01  1.0356E+00  7.9833E-01  1.0203E+00  1.2740E+00  4.9908E-01  7.5804E-01  8.0657E-01
             1.8990E+00
 PARAMETER:  1.2068E-01  3.1368E-02 -2.7795E-01  1.3497E-01 -1.2523E-01  1.2005E-01  3.4218E-01 -5.9499E-01 -1.7702E-01 -1.1497E-01
             7.4132E-01
 GRADIENT:   3.3167E+00 -2.8077E+00  1.5614E+00 -7.7079E+00 -2.3785E+00  1.3301E-01 -6.4803E-01  4.1780E-01  8.5114E-01  9.2643E-01
             1.2262E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1518.34487776474        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  1.0205E+00  8.6760E-01  6.5413E-01  1.0720E+00  7.4860E-01  1.0202E+00  1.3601E+00  3.8190E-01  7.3341E-01  7.6959E-01
             1.8911E+00
 PARAMETER:  1.2029E-01 -4.2026E-02 -3.2446E-01  1.6956E-01 -1.8955E-01  1.2002E-01  4.0757E-01 -8.6259E-01 -2.1006E-01 -1.6190E-01
             7.3716E-01
 GRADIENT:   1.9982E+00  3.8357E-01 -5.3717E-01  1.4432E+00  1.3829E-01 -1.2198E-01 -7.9640E-02  1.6960E-01 -2.6874E-01  7.8995E-02
            -7.3339E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1518.37187531682        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      824
 NPARAMETR:  1.0201E+00  8.8530E-01  6.3200E-01  1.0583E+00  7.4191E-01  1.0204E+00  1.3357E+00  2.6060E-01  7.3947E-01  7.6897E-01
             1.8938E+00
 PARAMETER:  1.1988E-01 -2.1825E-02 -3.5887E-01  1.5670E-01 -1.9853E-01  1.2019E-01  3.8943E-01 -1.2448E+00 -2.0182E-01 -1.6270E-01
             7.3860E-01
 GRADIENT:   4.4063E-01  2.4921E-01 -2.4057E-01  1.6718E-01 -1.6573E-01 -1.9246E-01  1.2688E-01  2.9315E-02 -1.6831E-01  3.1610E-01
            -2.4210E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1518.37422879280        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      986
 NPARAMETR:  1.0199E+00  8.8521E-01  6.2814E-01  1.0576E+00  7.3944E-01  1.0209E+00  1.3331E+00  2.1750E-01  7.4107E-01  7.6648E-01
             1.8956E+00
 PARAMETER:  1.1972E-01 -2.1933E-02 -3.6499E-01  1.5605E-01 -2.0187E-01  1.2070E-01  3.8750E-01 -1.4256E+00 -1.9966E-01 -1.6595E-01
             7.3955E-01
 GRADIENT:   3.5928E-03  8.3629E-03 -1.4727E-02  1.7942E-02  1.7518E-03  3.4538E-04 -1.2935E-03  1.2858E-03 -8.2998E-03  1.0510E-02
            -2.7002E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      986
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.9008E-04  8.7243E-03 -6.1825E-03 -1.6581E-02 -9.6987E-03
 SE:             2.9550E-02  2.1670E-02  3.5915E-03  2.1534E-02  1.8934E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8407E-01  6.8724E-01  8.5174E-02  4.4132E-01  6.0848E-01

 ETASHRINKSD(%)  1.0045E+00  2.7403E+01  8.7968E+01  2.7857E+01  3.6569E+01
 ETASHRINKVR(%)  1.9988E+00  4.7296E+01  9.8552E+01  4.7954E+01  5.9765E+01
 EBVSHRINKSD(%)  1.3151E+00  2.7316E+01  8.8380E+01  2.7632E+01  3.5792E+01
 EBVSHRINKVR(%)  2.6129E+00  4.7170E+01  9.8650E+01  4.7629E+01  5.8774E+01
 RELATIVEINF(%)  9.6837E+01  3.9235E+00  9.2818E-02  4.6300E+00  2.3983E+00
 EPSSHRINKSD(%)  3.6236E+01
 EPSSHRINKVR(%)  5.9342E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1518.3742287927989     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -783.22340222906075     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     4.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1518.374       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  8.85E-01  6.28E-01  1.06E+00  7.39E-01  1.02E+00  1.33E+00  2.18E-01  7.41E-01  7.66E-01  1.90E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       27.493
Stop Time:
Sun Oct 24 01:57:29 CDT 2021
