Sun Oct 24 02:35:41 CDT 2021
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
$DATA ../../../../data/SD4/A3/dat80.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   135.783969406214        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3026E+02  1.4566E+02  1.3993E+02  6.5491E+01  1.7579E+02  7.7510E+01 -4.5848E+01 -6.0090E+01 -9.2168E+01 -1.5964E+02
            -3.2299E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1022.22476432048        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0225E+00  7.9984E-01  6.4490E-01  1.0873E+00  6.3947E-01  8.4928E-01  9.5052E-01  9.0420E-01  1.1637E+00  1.0724E+00
             1.8016E+00
 PARAMETER:  1.2221E-01 -1.2335E-01 -3.3865E-01  1.8373E-01 -3.4712E-01 -6.3368E-02  4.9258E-02 -7.0429E-04  2.5165E-01  1.6991E-01
             6.8868E-01
 GRADIENT:   1.9048E+02  7.7023E+01  7.3088E+01  4.7177E+01 -1.6900E+01 -2.1359E-01 -9.8676E+00  3.5375E+00 -7.2198E+00  7.1205E+00
            -8.2311E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1026.04838009708        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:      205
 NPARAMETR:  1.0248E+00  6.8521E-01  3.9835E-01  1.1150E+00  4.4337E-01  8.6622E-01  8.3285E-01  7.7341E-01  9.4852E-01  9.6592E-01
             1.7665E+00
 PARAMETER:  1.2453E-01 -2.7803E-01 -8.2042E-01  2.0882E-01 -7.1336E-01 -4.3621E-02 -8.2905E-02 -1.5694E-01  4.7145E-02  6.5331E-02
             6.6898E-01
 GRADIENT:   2.5227E+01  1.1316E+02  8.4493E+01  6.0738E+01 -6.7352E+01 -4.8020E+00 -1.9104E+01 -2.3812E+00 -8.8133E+01  3.4026E+01
            -8.1442E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1203.50463392975        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      384
 NPARAMETR:  1.0658E+00  2.9826E-01  5.4720E-01  1.2909E+00  4.0381E-01  7.5072E-01  1.0000E-02  2.0181E-01  1.4574E+00  1.0409E+00
             2.6072E+00
 PARAMETER:  1.6374E-01 -1.1098E+00 -5.0294E-01  3.5536E-01 -8.0682E-01 -1.8672E-01 -5.3049E+00 -1.5005E+00  4.7666E-01  1.4011E-01
             1.0583E+00
 GRADIENT:   1.7315E+02  1.8426E+01  1.8826E+02 -2.4225E+01 -2.0528E+02 -5.5128E+01  0.0000E+00 -3.1408E-01  7.7482E+01  9.9948E+00
            -2.6245E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1294.53736415585        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      559
 NPARAMETR:  1.0370E+00  3.2435E-01  2.4196E-01  1.1351E+00  2.7206E-01  8.6246E-01  1.0000E-02  1.7441E-01  1.1137E+00  5.1108E-01
             3.5230E+00
 PARAMETER:  1.3631E-01 -1.0259E+00 -1.3190E+00  2.2670E-01 -1.2017E+00 -4.7964E-02 -6.7003E+00 -1.6464E+00  2.0772E-01 -5.7122E-01
             1.3593E+00
 GRADIENT:  -9.5487E+00  6.0756E+00  8.7950E+00  7.7538E+00 -2.8482E+01 -2.1550E+00  0.0000E+00  3.6445E-01 -4.2895E+00  8.1688E-01
             9.4563E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1295.31552320170        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  1.0338E+00  2.7962E-01  3.0867E-01  1.2202E+00  3.1212E-01  8.5781E-01  1.0000E-02  2.0836E-01  1.0496E+00  4.2169E-01
             3.6166E+00
 PARAMETER:  1.3320E-01 -1.1743E+00 -1.0755E+00  2.9899E-01 -1.0644E+00 -5.3374E-02 -7.4702E+00 -1.4685E+00  1.4843E-01 -7.6349E-01
             1.3855E+00
 GRADIENT:  -7.2648E+00  2.7348E+00  4.4445E+00  2.6564E+00 -9.9720E+00  1.0730E+00  0.0000E+00  5.1754E-01 -1.5015E+00  5.9887E-02
             6.3126E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1295.76462372552        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      910
 NPARAMETR:  1.0311E+00  2.0523E-01  3.2049E-01  1.2559E+00  3.1132E-01  8.5219E-01  1.0000E-02  7.0907E-02  1.0252E+00  4.6014E-01
             3.5823E+00
 PARAMETER:  1.3066E-01 -1.4836E+00 -1.0379E+00  3.2788E-01 -1.0669E+00 -5.9946E-02 -9.0134E+00 -2.5464E+00  1.2489E-01 -6.7622E-01
             1.3760E+00
 GRADIENT:  -1.6642E+00  8.1595E-02 -3.4536E+00  4.5177E-01  4.0092E+00  8.4536E-01  0.0000E+00  5.9467E-02 -2.9950E-01  1.0384E+00
             5.9677E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1295.86810320554        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1085
 NPARAMETR:  1.0269E+00  1.4861E-01  3.3838E-01  1.2913E+00  3.1778E-01  8.4412E-01  1.0000E-02  2.8188E-02  9.9493E-01  4.1917E-01
             3.6256E+00
 PARAMETER:  1.2651E-01 -1.8065E+00 -9.8358E-01  3.5568E-01 -1.0464E+00 -6.9456E-02 -1.1042E+01 -3.4689E+00  9.4914E-02 -7.6949E-01
             1.3880E+00
 GRADIENT:  -4.0175E-01  3.9275E-02  3.0359E-02  3.7473E-01 -1.2485E-01  7.8701E-02  0.0000E+00  6.6244E-03 -4.3623E-02 -8.2822E-02
             1.4884E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1295.87113755550        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1248
 NPARAMETR:  1.0268E+00  1.4409E-01  3.3794E-01  1.2919E+00  3.1686E-01  8.4380E-01  1.0000E-02  1.0000E-02  9.9405E-01  4.2382E-01
             3.6213E+00
 PARAMETER:  1.2642E-01 -1.8373E+00 -9.8490E-01  3.5615E-01 -1.0493E+00 -6.9838E-02 -1.1244E+01 -4.6238E+00  9.4033E-02 -7.5844E-01
             1.3868E+00
 GRADIENT:   1.5032E-01  1.5459E-02  7.2327E-01 -2.6744E-01 -9.0910E-01 -4.9282E-03  0.0000E+00  0.0000E+00 -3.3393E-02 -4.8566E-02
            -2.7661E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1248
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.4317E-04 -3.6412E-05  9.6407E-05 -1.3327E-02  1.1956E-03
 SE:             2.8203E-02  2.4982E-05  2.2558E-04  2.5666E-02  1.3733E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7332E-01  1.4497E-01  6.6911E-01  6.0358E-01  9.3063E-01

 ETASHRINKSD(%)  5.5171E+00  9.9916E+01  9.9244E+01  1.4017E+01  5.3991E+01
 ETASHRINKVR(%)  1.0730E+01  1.0000E+02  9.9994E+01  2.6069E+01  7.8832E+01
 EBVSHRINKSD(%)  5.3962E+00  9.9932E+01  9.9214E+01  1.2959E+01  5.4070E+01
 EBVSHRINKVR(%)  1.0501E+01  1.0000E+02  9.9994E+01  2.4238E+01  7.8904E+01
 RELATIVEINF(%)  6.8163E+01  3.2714E-06  1.6239E-04  1.6890E+01  3.5428E-01
 EPSSHRINKSD(%)  2.4326E+01
 EPSSHRINKVR(%)  4.2734E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1295.8711375555022     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -560.72031099176399     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1295.871       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.44E-01  3.38E-01  1.29E+00  3.17E-01  8.44E-01  1.00E-02  1.00E-02  9.94E-01  4.24E-01  3.62E+00
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       38.191
Stop Time:
Sun Oct 24 02:35:50 CDT 2021
