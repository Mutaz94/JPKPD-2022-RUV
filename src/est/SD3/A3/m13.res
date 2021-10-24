Sat Oct 23 22:29:26 CDT 2021
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
$DATA ../../../../data/SD3/A3/dat13.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m13.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -269.979976479127        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5676E+02  5.8909E+01  1.3387E+02  3.6190E+01  2.4466E+02  3.3241E+01 -6.6530E+01 -1.1629E+02 -5.4151E+01 -1.1027E+02
            -3.2630E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1549.42388320439        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0480E+00  1.0054E+00  9.9342E-01  1.1517E+00  9.5470E-01  9.4175E-01  9.8525E-01  9.4537E-01  1.0703E+00  8.2399E-01
             3.7864E+00
 PARAMETER:  1.4688E-01  1.0536E-01  9.3403E-02  2.4123E-01  5.3639E-02  3.9989E-02  8.5141E-02  4.3824E-02  1.6794E-01 -9.3596E-02
             1.4314E+00
 GRADIENT:   1.1063E+02  2.7346E+01 -2.2020E+01  7.8809E+01  7.5856E+00 -6.2219E+00  8.6217E+00  7.8022E+00  1.2627E+01  2.0649E+01
             1.6660E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1568.88953623024        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:      188
 NPARAMETR:  1.0279E+00  7.5453E-01  3.2063E-01  1.1876E+00  4.2379E-01  1.0098E+00  4.5770E-01  1.4467E-01  1.2122E+00  4.1842E-01
             3.4277E+00
 PARAMETER:  1.2747E-01 -1.8167E-01 -1.0375E+00  2.7194E-01 -7.5853E-01  1.0978E-01 -6.8155E-01 -1.8333E+00  2.9246E-01 -7.7126E-01
             1.3319E+00
 GRADIENT:   2.0835E+01  6.7998E+01 -2.0274E+01  1.2184E+02 -1.7222E+01  5.3810E+00 -4.0399E+00  1.3241E-01 -9.7067E+00  1.2013E+00
             1.1220E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1591.48578633429        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      363
 NPARAMETR:  9.8632E-01  5.8511E-01  3.4091E-01  1.2082E+00  3.9265E-01  9.8995E-01  7.6883E-01  2.8219E-02  1.1638E+00  3.7614E-01
             2.8485E+00
 PARAMETER:  8.6229E-02 -4.3595E-01 -9.7613E-01  2.8914E-01 -8.3483E-01  8.9902E-02 -1.6288E-01 -3.4677E+00  2.5169E-01 -8.7780E-01
             1.1468E+00
 GRADIENT:  -3.6199E+01  1.4290E+01 -2.6343E+01  8.5865E+01  6.3277E+01 -1.2563E+00 -8.3630E+00 -1.4600E-02 -1.5756E+01 -1.1079E+01
            -1.4094E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1605.60172528710        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  1.0001E+00  4.7639E-01  2.6651E-01  1.1608E+00  3.1421E-01  9.9225E-01  8.2357E-01  1.0000E-02  1.2514E+00  5.8315E-01
             2.6268E+00
 PARAMETER:  1.0007E-01 -6.4151E-01 -1.2223E+00  2.4914E-01 -1.0577E+00  9.2217E-02 -9.4108E-02 -6.2712E+00  3.2428E-01 -4.3931E-01
             1.0658E+00
 GRADIENT:  -2.4618E-01  6.2799E+00 -2.0324E+01  4.6062E+01  4.0764E+01 -1.1148E+00 -9.3269E-01  0.0000E+00 -4.0045E+00  2.4346E+00
            -2.0217E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1609.64404419334        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      713
 NPARAMETR:  1.0007E+00  3.8476E-01  2.0952E-01  1.0778E+00  2.5454E-01  9.9473E-01  5.7704E-01  1.0000E-02  1.3194E+00  6.2991E-01
             2.6328E+00
 PARAMETER:  1.0069E-01 -8.5513E-01 -1.4629E+00  1.7488E-01 -1.2683E+00  9.4712E-02 -4.4985E-01 -8.1682E+00  3.7717E-01 -3.6218E-01
             1.0680E+00
 GRADIENT:   3.3095E-01 -5.8153E-01 -8.5345E-02 -1.1963E+00  6.8969E-01  1.4915E-02  7.1213E-01  0.0000E+00  1.0495E+00  1.0037E+00
             2.0609E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1609.87801655157        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      889
 NPARAMETR:  1.0006E+00  3.8137E-01  2.0107E-01  1.0675E+00  2.4835E-01  9.9541E-01  1.9523E-01  1.0000E-02  1.3363E+00  6.4969E-01
             2.6238E+00
 PARAMETER:  1.0055E-01 -8.6399E-01 -1.5041E+00  1.6532E-01 -1.2929E+00  9.5396E-02 -1.5336E+00 -7.1336E+00  3.8991E-01 -3.3125E-01
             1.0646E+00
 GRADIENT:  -2.0705E-01  2.6632E-01  8.3351E-01  5.2277E-01 -1.1101E+00  1.0176E-02  3.6812E-02  0.0000E+00  1.5577E-01  4.5358E-02
             1.6108E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1609.90073904293        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  1.0006E+00  3.8006E-01  1.9872E-01  1.0640E+00  2.4669E-01  9.9552E-01  2.1742E-02  1.0000E-02  1.3401E+00  6.5339E-01
             2.6217E+00
 PARAMETER:  1.0063E-01 -8.6743E-01 -1.5158E+00  1.6199E-01 -1.2996E+00  9.5514E-02 -3.7285E+00 -4.6285E+00  3.9273E-01 -3.2559E-01
             1.0638E+00
 GRADIENT:  -1.4880E-02 -6.1558E-02 -1.3992E-01  1.1323E-01  1.5427E-01  1.6357E-02  5.0480E-04  0.0000E+00  1.0911E-02 -4.9551E-02
            -2.5065E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1609.90093861175        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     1193
 NPARAMETR:  1.0007E+00  3.8017E-01  1.9877E-01  1.0638E+00  2.4667E-01  9.9549E-01  1.0000E-02  1.0000E-02  1.3404E+00  6.5360E-01
             2.6216E+00
 PARAMETER:  1.0070E-01 -8.6713E-01 -1.5156E+00  1.6185E-01 -1.2997E+00  9.5479E-02 -4.5839E+00 -4.9142E+00  3.9299E-01 -3.2527E-01
             1.0638E+00
 GRADIENT:   1.1040E-01  2.2080E-01  4.1586E-01 -4.5084E-02 -7.5348E-01  6.6775E-03  0.0000E+00  0.0000E+00  8.9655E-02 -4.6991E-03
            -5.3341E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1193
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2163E-03 -1.0470E-04  1.6060E-04 -8.3694E-03  1.0803E-03
 SE:             2.9159E-02  1.3403E-04  2.1240E-04  2.7393E-02  2.4294E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6673E-01  4.3468E-01  4.4956E-01  7.5996E-01  9.6453E-01

 ETASHRINKSD(%)  2.3130E+00  9.9551E+01  9.9288E+01  8.2299E+00  1.8614E+01
 ETASHRINKVR(%)  4.5726E+00  9.9998E+01  9.9995E+01  1.5782E+01  3.3763E+01
 EBVSHRINKSD(%)  2.2037E+00  9.9553E+01  9.9292E+01  6.4216E+00  1.9175E+01
 EBVSHRINKVR(%)  4.3587E+00  9.9998E+01  9.9995E+01  1.2431E+01  3.4673E+01
 RELATIVEINF(%)  9.5537E+01  3.2343E-04  3.4620E-04  4.8481E+01  2.4086E+00
 EPSSHRINKSD(%)  2.7275E+01
 EPSSHRINKVR(%)  4.7110E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1609.9009386117505     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -690.96240540707777     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1609.901       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  3.80E-01  1.99E-01  1.06E+00  2.47E-01  9.95E-01  1.00E-02  1.00E-02  1.34E+00  6.54E-01  2.62E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       98.632
Stop Time:
Sat Oct 23 22:29:42 CDT 2021
