Sat Oct 23 23:19:16 CDT 2021
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
$DATA ../../../../data/SD3/SL1/dat6.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m6.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2143.23695654727        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2934E+02 -4.2493E+01 -8.8611E+01  6.2856E+01  1.7096E+02  5.7801E+01  6.3654E+00  1.5112E+01  8.9230E+00  4.2649E+00
             2.5602E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2149.15507974868        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  1.0387E+00  1.0508E+00  1.1535E+00  1.0303E+00  8.7814E-01  9.7112E-01  9.6920E-01  9.2486E-01  1.0008E+00  9.2866E-01
             9.3574E-01
 PARAMETER:  1.3796E-01  1.4953E-01  2.4284E-01  1.2989E-01 -2.9944E-02  7.0693E-02  6.8714E-02  2.1892E-02  1.0079E-01  2.5984E-02
             3.3579E-02
 GRADIENT:   8.4793E+01  8.1999E+01  6.1673E+01  1.7582E+01 -1.2989E+02 -4.3178E+00  4.4149E+00  3.6514E+00 -2.6674E+00  3.9680E+00
            -3.7095E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2154.14346701800        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  1.0346E+00  9.3474E-01  9.9468E-01  1.0860E+00  8.0226E-01  9.8052E-01  7.0109E-01  3.8330E-01  9.9060E-01  1.0870E+00
             9.5852E-01
 PARAMETER:  1.3404E-01  3.2509E-02  9.4663E-02  1.8251E-01 -1.2032E-01  8.0333E-02 -2.5512E-01 -8.5894E-01  9.0552E-02  1.8338E-01
             5.7631E-02
 GRADIENT:   7.4307E+01  6.7894E+01  3.5916E+01  3.7787E+01 -1.0820E+02  2.6997E-01 -1.9489E+00  1.8384E+00 -1.1072E+01  2.8596E+01
            -5.8657E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2163.03294103958        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.9997E-01  8.3567E-01  1.0143E+00  1.1269E+00  8.3026E-01  9.7011E-01  9.2181E-01  3.8283E-01  9.4781E-01  9.2350E-01
             9.6441E-01
 PARAMETER:  9.9970E-02 -7.9523E-02  1.1420E-01  2.1947E-01 -8.6020E-02  6.9655E-02  1.8582E-02 -8.6016E-01  4.6402E-02  2.0412E-02
             6.3758E-02
 GRADIENT:  -3.2513E+00  6.4589E+00 -8.7879E-01  5.3501E+00  1.2119E+00 -1.3976E+00 -1.0517E+00  1.1414E+00 -2.2707E+00  3.4935E+00
            -1.1240E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2163.77386214527        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.9987E-01  6.9564E-01  9.9584E-01  1.2039E+00  7.7685E-01  9.7194E-01  1.1997E+00  2.2421E-01  8.8468E-01  8.8476E-01
             9.6507E-01
 PARAMETER:  9.9875E-02 -2.6292E-01  9.5828E-02  2.8553E-01 -1.5250E-01  7.1540E-02  2.8205E-01 -1.3952E+00 -2.2531E-02 -2.2441E-02
             6.4444E-02
 GRADIENT:  -8.5374E-02 -7.7205E-01 -8.1661E-01 -1.3151E+00  4.9196E-01 -8.1972E-02  1.5672E-01  1.8270E-01  1.1494E-01  2.5056E-01
             1.6153E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2163.77586923876        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      877
 NPARAMETR:  9.9981E-01  6.8780E-01  1.0009E+00  1.2094E+00  7.7685E-01  9.7207E-01  1.1970E+00  2.0321E-01  8.8327E-01  8.9093E-01
             9.6507E-01
 PARAMETER:  9.9807E-02 -2.7426E-01  1.0092E-01  2.9012E-01 -1.5251E-01  7.1672E-02  2.7980E-01 -1.4935E+00 -2.4120E-02 -1.5487E-02
             6.4441E-02
 GRADIENT:   4.2806E-02 -6.5696E-01 -5.7470E-01 -2.7540E-01  6.8078E-01  1.4680E-02  4.7703E-02  1.1586E-01  7.2103E-02  3.4562E-01
             6.3866E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2163.89553290381        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  1.0008E+00  7.4847E-01  9.5713E-01  1.1709E+00  7.7711E-01  9.7301E-01  1.1695E+00  1.3319E-02  8.9873E-01  8.8087E-01
             9.6457E-01
 PARAMETER:  1.0080E-01 -1.8972E-01  5.6182E-02  2.5779E-01 -1.5218E-01  7.2641E-02  2.5660E-01 -4.2186E+00 -6.7714E-03 -2.6844E-02
             6.3932E-02
 GRADIENT:   3.0762E-01  1.1560E+00 -1.1340E+00  1.5282E+00 -1.3734E+00  5.0258E-02  3.8571E-01  9.8227E-04  8.7241E-02  1.2391E+00
             3.0946E-01

0ITERATION NO.:   34    OBJECTIVE VALUE:  -2163.90012925980        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     1186
 NPARAMETR:  1.0014E+00  7.3802E-01  9.7041E-01  1.1768E+00  7.8003E-01  9.7307E-01  1.1600E+00  1.0000E-02  8.9765E-01  8.8497E-01
             9.6473E-01
 PARAMETER:  1.0137E-01 -2.0379E-01  6.9960E-02  2.6277E-01 -1.4842E-01  7.2698E-02  2.4845E-01 -4.6730E+00 -7.9712E-03 -2.2205E-02
             6.4095E-02
 GRADIENT:   2.0914E+00 -1.7080E-01  4.8843E-01 -1.0983E+00 -3.4102E-01  1.4647E-01  5.3770E-03  0.0000E+00 -2.5235E-02  3.9829E-02
             1.1986E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1186
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0906E-04 -4.7764E-03 -3.3801E-04 -3.0970E-04 -1.2388E-02
 SE:             2.9863E-02  1.5960E-02  1.9800E-04  2.6902E-02  2.5302E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9709E-01  7.6473E-01  8.7808E-02  9.9081E-01  6.2442E-01

 ETASHRINKSD(%)  1.0000E-10  4.6533E+01  9.9337E+01  9.8750E+00  1.5237E+01
 ETASHRINKVR(%)  1.0000E-10  7.1413E+01  9.9996E+01  1.8775E+01  2.8152E+01
 EBVSHRINKSD(%)  3.3659E-01  4.7175E+01  9.9373E+01  9.8451E+00  1.3905E+01
 EBVSHRINKVR(%)  6.7204E-01  7.2095E+01  9.9996E+01  1.8721E+01  2.5877E+01
 RELATIVEINF(%)  9.8568E+01  1.7817E+00  5.3840E-04  6.9133E+00  1.1420E+01
 EPSSHRINKSD(%)  3.3315E+01
 EPSSHRINKVR(%)  5.5531E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2163.9001292597995     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1244.9615960551268     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2163.900       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  7.38E-01  9.70E-01  1.18E+00  7.80E-01  9.73E-01  1.16E+00  1.00E-02  8.98E-01  8.85E-01  9.65E-01
 


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
 #CPUT: Total CPU Time in Seconds,       94.505
Stop Time:
Sat Oct 23 23:19:31 CDT 2021
