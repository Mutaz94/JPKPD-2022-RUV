Sat Oct 23 15:38:43 CDT 2021
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
$DATA ../../../../data/SD1/SL3/dat92.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      979
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      879
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -617.524635784727        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5251E+02  2.6273E+01  1.5670E+02  1.0991E+02  1.8513E+02  6.0327E+01 -1.2585E+02 -2.3632E+02 -1.2933E+02 -1.7902E+01
            -5.8234E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2680.29833292498        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0454E+00  1.3325E+00  8.8353E-01  9.4829E-01  1.0423E+00  8.4030E-01  1.3805E+00  1.0142E+00  1.2203E+00  9.4235E-01
             2.5990E+00
 PARAMETER:  1.4437E-01  3.8703E-01 -2.3830E-02  4.6900E-02  1.4143E-01 -7.3994E-02  4.2247E-01  1.1406E-01  2.9908E-01  4.0620E-02
             1.0551E+00
 GRADIENT:   2.1127E+02  1.2707E+02 -1.2936E+01  8.2003E+01 -1.3640E+01 -4.0368E+01  3.2415E+01  2.3868E+00  1.1692E+01  3.1472E+00
            -1.9754E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2684.05778735124        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      222
 NPARAMETR:  1.0546E+00  1.5002E+00  9.9273E-01  8.6847E-01  1.1966E+00  9.4264E-01  1.2401E+00  9.5318E-01  1.3129E+00  9.7213E-01
             2.6158E+00
 PARAMETER:  1.5317E-01  5.0561E-01  9.2699E-02 -4.1019E-02  2.7945E-01  4.0933E-02  3.1522E-01  5.2048E-02  3.7223E-01  7.1735E-02
             1.0616E+00
 GRADIENT:   1.1633E+02  8.7415E+01  7.7541E-01  7.7522E+01 -9.2495E+00 -7.8039E-01  1.5413E+01 -2.0515E+00  7.0398E+00 -1.2994E+01
            -2.6797E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2697.62051084309        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      399
 NPARAMETR:  1.0075E+00  1.6173E+00  1.2752E+00  7.1244E-01  1.4538E+00  9.3109E-01  9.4964E-01  1.7213E+00  1.2901E+00  1.2701E+00
             2.6251E+00
 PARAMETER:  1.0749E-01  5.8073E-01  3.4309E-01 -2.3906E-01  4.7421E-01  2.8603E-02  4.8333E-02  6.4310E-01  3.5472E-01  3.3908E-01
             1.0651E+00
 GRADIENT:   6.4607E+00 -5.3945E-01 -3.3480E-01  2.3325E+01  1.3553E+01  7.8412E-01 -7.2413E+00 -1.7842E+00 -5.1137E+00 -6.4372E-01
             2.3458E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2702.74587808165        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      574
 NPARAMETR:  1.0065E+00  2.1159E+00  1.1277E+00  3.8655E-01  1.7775E+00  9.2208E-01  8.4675E-01  2.1409E+00  1.8501E+00  1.5531E+00
             2.6201E+00
 PARAMETER:  1.0646E-01  8.4946E-01  2.2015E-01 -8.5050E-01  6.7520E-01  1.8874E-02 -6.6346E-02  8.6124E-01  7.1525E-01  5.4023E-01
             1.0632E+00
 GRADIENT:   4.7424E+00  1.9007E+01  1.9010E+00  8.2139E+00 -5.8583E+00 -2.7941E+00 -2.6568E+00  2.2927E-01 -2.8252E-01  2.5327E+00
             3.8273E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2706.12935975795        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      750
 NPARAMETR:  1.0049E+00  2.4882E+00  4.6487E-01  1.3749E-01  1.9538E+00  9.4279E-01  7.8069E-01  8.1499E-01  3.3857E+00  1.6673E+00
             2.5989E+00
 PARAMETER:  1.0490E-01  1.0116E+00 -6.6600E-01 -1.8842E+00  7.6975E-01  4.1092E-02 -1.4758E-01 -1.0458E-01  1.3196E+00  6.1121E-01
             1.0551E+00
 GRADIENT:   1.6423E+00  2.9980E+01  1.9470E+00  5.1141E+00 -7.2908E+00  5.0297E+00 -6.7581E+00  3.7692E-01  5.6315E+00 -8.9649E-01
            -8.9815E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2707.77153738944        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      927
 NPARAMETR:  1.0037E+00  2.5566E+00  2.7995E-01  8.1065E-02  2.0140E+00  9.2737E-01  7.9158E-01  4.8382E-01  3.9269E+00  1.7085E+00
             2.6052E+00
 PARAMETER:  1.0374E-01  1.0387E+00 -1.1731E+00 -2.4125E+00  8.0012E-01  2.4602E-02 -1.3373E-01 -6.2603E-01  1.4678E+00  6.3562E-01
             1.0575E+00
 GRADIENT:  -1.5432E+00  1.3462E+01 -1.2630E+00  1.3781E+00 -3.2866E-01 -8.3082E-01 -9.2436E-02  2.5614E-01 -5.1369E-01 -1.5600E+00
             1.9677E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2707.87112023708        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  1.0029E+00  2.5438E+00  2.9191E-01  7.9280E-02  2.0321E+00  9.2938E-01  7.9456E-01  1.2667E-01  3.9227E+00  1.7245E+00
             2.6011E+00
 PARAMETER:  1.0290E-01  1.0337E+00 -1.1313E+00 -2.4348E+00  8.0908E-01  2.6758E-02 -1.2997E-01 -1.9661E+00  1.4668E+00  6.4493E-01
             1.0559E+00
 GRADIENT:   6.3947E+01  4.7736E+02  1.0963E+00  2.0278E+00  3.6669E+01  6.2405E+00  6.9960E+00  1.7863E-02  9.8187E-02  9.3139E+00
             2.0480E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2707.97150273811        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1223
 NPARAMETR:  1.0044E+00  2.5463E+00  2.8782E-01  7.9942E-02  2.0276E+00  9.2913E-01  7.9066E-01  1.0727E-01  3.9494E+00  1.7143E+00
             2.6014E+00
 PARAMETER:  1.0438E-01  1.0346E+00 -1.1454E+00 -2.4265E+00  8.0687E-01  2.6494E-02 -1.3489E-01 -2.1324E+00  1.4736E+00  6.3902E-01
             1.0561E+00
 GRADIENT:   2.0992E-01  3.2033E+00 -5.1825E-02 -1.5094E+00  4.9289E+00 -1.3820E-01  4.1050E-01  1.2757E-02 -4.2086E+00  2.9296E-01
             1.0167E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2708.02353100020        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1415
 NPARAMETR:  1.0042E+00  2.5452E+00  3.0139E-01  8.0480E-02  2.0181E+00  9.2981E-01  7.9580E-01  5.4644E-02  3.9875E+00  1.7167E+00
             2.5998E+00
 PARAMETER:  1.0420E-01  1.0341E+00 -1.0994E+00 -2.4194E+00  8.0204E-01  2.7241E-02 -1.2842E-01 -2.7791E+00  1.4833E+00  6.4038E-01
             1.0554E+00
 GRADIENT:  -1.3633E-01 -2.8163E+01 -1.2777E-01  1.3677E+01 -4.7392E+01  1.8003E-01 -4.4248E-01  3.5103E-03  2.2683E+01 -4.0462E-01
            -3.2257E-01
 NUMSIGDIG:         3.3         2.6         2.8         2.4         2.3         2.3         2.5         0.5         2.4         2.8
                    3.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1415
 NO. OF SIG. DIGITS IN FINAL EST.:  0.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4159E-03 -2.3054E-02 -6.5198E-05  3.1804E-02 -2.3770E-02
 SE:             2.9340E-02  2.7164E-02  6.8274E-05  1.4312E-02  2.5863E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6151E-01  3.9604E-01  3.3960E-01  2.6271E-02  3.5806E-01

 ETASHRINKSD(%)  1.7072E+00  8.9986E+00  9.9771E+01  5.2052E+01  1.3355E+01
 ETASHRINKVR(%)  3.3853E+00  1.7187E+01  9.9999E+01  7.7010E+01  2.4926E+01
 EBVSHRINKSD(%)  1.8344E+00  6.8730E+00  9.9578E+01  6.6587E+01  9.1957E+00
 EBVSHRINKVR(%)  3.6351E+00  1.3274E+01  9.9998E+01  8.8836E+01  1.7546E+01
 RELATIVEINF(%)  9.6315E+01  3.8683E+01  1.5246E-03  4.6211E+00  6.7142E+01
 EPSSHRINKSD(%)  1.6577E+01
 EPSSHRINKVR(%)  3.0406E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          879
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1615.4939413738146     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2708.0235310002045     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1092.5295896263899     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.85
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2708.024       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.54E+00  3.01E-01  8.05E-02  2.02E+00  9.30E-01  7.96E-01  5.62E-02  3.99E+00  1.72E+00  2.60E+00
 


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
 #CPUT: Total CPU Time in Seconds,      102.272
Stop Time:
Sat Oct 23 15:39:00 CDT 2021
