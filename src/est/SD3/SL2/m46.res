Sat Oct 23 23:57:35 CDT 2021
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
$DATA ../../../../data/SD3/SL2/dat46.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2101.30933546756        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6924E+02 -5.3347E+01 -4.4044E+01 -3.9319E+00  5.7647E+01  5.2899E+01  5.4222E+00  1.8730E+01  2.3680E+01  2.1816E+01
            -5.2406E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2112.72372447708        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.7131E-01  1.0759E+00  1.0927E+00  1.0348E+00  1.0459E+00  9.3572E-01  9.6997E-01  9.0790E-01  9.1278E-01  8.8869E-01
             1.0836E+00
 PARAMETER:  7.0890E-02  1.7318E-01  1.8861E-01  1.3423E-01  1.4484E-01  3.3562E-02  6.9513E-02  3.3840E-03  8.7342E-03 -1.8011E-02
             1.8032E-01
 GRADIENT:   9.1363E+00  4.6064E+00 -2.0743E+01  3.3219E+01  4.5462E+01 -6.6907E+00 -9.9504E-01  9.4143E+00 -1.2089E+00 -5.4256E+00
             6.2875E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2116.61082367495        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.7053E-01  8.9827E-01  1.0041E+00  1.1441E+00  9.2675E-01  8.9558E-01  1.0898E+00  3.5195E-01  8.8443E-01  9.2804E-01
             1.0550E+00
 PARAMETER:  7.0088E-02 -7.2871E-03  1.0412E-01  2.3462E-01  2.3932E-02 -1.0279E-02  1.8599E-01 -9.4428E-01 -2.2816E-02  2.5320E-02
             1.5353E-01
 GRADIENT:   1.1161E+01  6.8961E+00 -2.9477E+01  5.0884E+01  3.2867E+01 -2.4686E+01 -1.3572E+00  2.2000E+00  5.7209E+00  8.0160E+00
            -8.2948E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2119.22742094375        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.6601E-01  8.4274E-01  9.1537E-01  1.1466E+00  8.3954E-01  9.5597E-01  1.2327E+00  2.8336E-01  8.2950E-01  7.8735E-01
             1.0675E+00
 PARAMETER:  6.5419E-02 -7.1098E-02  1.1573E-02  2.3682E-01 -7.4898E-02  5.4974E-02  3.0918E-01 -1.1610E+00 -8.6934E-02 -1.3908E-01
             1.6532E-01
 GRADIENT:  -1.7388E+00  2.1414E+00 -5.9470E-01  4.3041E+00 -1.7338E+00  2.4482E+00 -6.3070E-01  1.1618E+00 -2.9427E-01  1.2215E+00
             8.1558E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2119.62767153528        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.6542E-01  7.5779E-01  9.5712E-01  1.1973E+00  8.3404E-01  9.4948E-01  1.3486E+00  9.5704E-02  8.0577E-01  8.1549E-01
             1.0669E+00
 PARAMETER:  6.4809E-02 -1.7735E-01  5.6177E-02  2.8011E-01 -8.1471E-02  4.8154E-02  3.9906E-01 -2.2465E+00 -1.1596E-01 -1.0397E-01
             1.6480E-01
 GRADIENT:  -2.7955E-01 -1.3205E+00 -1.7861E+00  1.8759E+00  2.7354E+00  4.0842E-01  9.8895E-02  1.1502E-01 -7.0905E-03  2.7576E-01
            -9.0465E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2119.81273046975        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  9.6693E-01  8.5275E-01  8.8689E-01  1.1356E+00  8.3090E-01  9.5016E-01  1.2336E+00  2.4431E-02  8.3559E-01  7.7476E-01
             1.0669E+00
 PARAMETER:  6.6372E-02 -5.9293E-02 -2.0029E-02  2.2715E-01 -8.5245E-02  4.8877E-02  3.0990E-01 -3.6119E+00 -7.9619E-02 -1.5520E-01
             1.6473E-01
 GRADIENT:  -2.8287E-02  1.7678E-01 -4.7008E-01  9.7336E-01  8.4178E-01 -3.6596E-02 -8.5433E-02  8.5661E-03  3.8146E-02 -4.5048E-02
            -3.3272E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2119.81385723377        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.6688E-01  8.4553E-01  8.8746E-01  1.1392E+00  8.2816E-01  9.5016E-01  1.2437E+00  1.7049E-02  8.3295E-01  7.7382E-01
             1.0674E+00
 PARAMETER:  6.6316E-02 -6.7787E-02 -1.9389E-02  2.3034E-01 -8.8554E-02  4.8870E-02  3.1805E-01 -3.9716E+00 -8.2783E-02 -1.5642E-01
             1.6526E-01
 GRADIENT:   5.3667E-04 -8.9189E-02  2.4500E-03 -2.6315E-01 -7.4666E-02  3.4748E-03  1.7302E-02  4.1756E-03 -9.7065E-03  4.4456E-02
             8.8494E-02

0ITERATION NO.:   34    OBJECTIVE VALUE:  -2119.81476599365        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1186
 NPARAMETR:  9.6686E-01  8.4467E-01  8.8846E-01  1.1399E+00  8.2835E-01  9.5013E-01  1.2445E+00  1.0000E-02  8.3271E-01  7.7440E-01
             1.0674E+00
 PARAMETER:  6.6301E-02 -6.8807E-02 -1.8263E-02  2.3094E-01 -8.8314E-02  4.8848E-02  3.1875E-01 -4.9458E+00 -8.3068E-02 -1.5567E-01
             1.6523E-01
 GRADIENT:   2.0153E-03 -1.2378E-02  3.9412E-02 -8.2408E-02 -9.5615E-02  2.5139E-03  1.2724E-02  0.0000E+00 -4.1257E-03  1.6723E-02
             4.1789E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1186
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.7074E-04  3.8085E-03 -3.8394E-04 -5.6398E-03 -1.0029E-02
 SE:             2.9844E-02  1.9940E-02  1.8961E-04  2.5331E-02  2.2910E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7405E-01  8.4853E-01  4.2882E-02  8.2381E-01  6.6157E-01

 ETASHRINKSD(%)  1.9424E-02  3.3198E+01  9.9365E+01  1.5138E+01  2.3249E+01
 ETASHRINKVR(%)  3.8844E-02  5.5375E+01  9.9996E+01  2.7985E+01  4.1093E+01
 EBVSHRINKSD(%)  4.2826E-01  3.3692E+01  9.9355E+01  1.5002E+01  2.2408E+01
 EBVSHRINKVR(%)  8.5469E-01  5.6033E+01  9.9996E+01  2.7754E+01  3.9794E+01
 RELATIVEINF(%)  9.8548E+01  3.0445E+00  4.1287E-04  6.2870E+00  6.1699E+00
 EPSSHRINKSD(%)  3.2339E+01
 EPSSHRINKVR(%)  5.4219E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2119.8147659936531     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1200.8762327889804     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2119.815       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.67E-01  8.45E-01  8.88E-01  1.14E+00  8.28E-01  9.50E-01  1.24E+00  1.00E-02  8.33E-01  7.74E-01  1.07E+00
 


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
 #CPUT: Total CPU Time in Seconds,       91.457
Stop Time:
Sat Oct 23 23:57:50 CDT 2021
