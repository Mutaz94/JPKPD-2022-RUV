Sat Oct 23 15:07:29 CDT 2021
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
$DATA ../../../../data/SD1/SL2/dat42.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      996
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

 TOT. NO. OF OBS RECS:      896
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2615.86305735496        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0024E+02  1.6715E+02  7.2841E-01  2.0547E+02  1.6699E+02  5.4366E+01 -7.7968E+01 -1.1406E+02 -8.9283E+01 -2.4663E+01
            -2.2558E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3198.84793167486        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  1.0481E+00  9.7151E-01  1.0655E+00  9.5521E-01  9.4761E-01  1.0152E+00  1.0733E+00  1.0356E+00  1.1257E+00  1.0191E+00
             1.8017E+00
 PARAMETER:  1.4701E-01  7.1094E-02  1.6342E-01  5.4177E-02  4.6189E-02  1.1511E-01  1.7070E-01  1.3501E-01  2.1845E-01  1.1890E-01
             6.8871E-01
 GRADIENT:   2.1022E+02  2.8411E+00  6.4698E+00  7.8866E-01 -1.5309E+01  2.9685E+01 -1.8222E+00 -1.2085E+00  5.9712E+00  2.0346E+00
            -1.7576E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3201.91829774711        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0335E+00  1.0592E+00  1.2196E+00  9.3471E-01  1.0677E+00  9.4419E-01  1.0980E+00  1.3968E+00  1.0937E+00  1.1014E+00
             1.8249E+00
 PARAMETER:  1.3297E-01  1.5755E-01  2.9852E-01  3.2482E-02  1.6549E-01  4.2572E-02  1.9345E-01  4.3415E-01  1.8958E-01  1.9655E-01
             7.0152E-01
 GRADIENT:   1.4140E+02  2.9745E+01 -5.3145E+00  2.4919E+01  2.2938E+01 -2.0706E+00  1.2613E+01  6.2571E-01  2.2238E+00  4.1679E+00
             1.8781E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3209.34176057371        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      342
 NPARAMETR:  1.0721E+00  1.2322E+00  1.5978E+00  8.7131E-01  1.2767E+00  9.6486E-01  8.9030E-01  2.5555E+00  1.1040E+00  1.2147E+00
             1.8472E+00
 PARAMETER:  1.6962E-01  3.0877E-01  5.6865E-01 -3.7755E-02  3.4431E-01  6.4225E-02 -1.6192E-02  1.0382E+00  1.9897E-01  2.9451E-01
             7.1369E-01
 GRADIENT:   1.7451E+01  1.8187E+01 -2.2273E+01  8.3583E+00 -1.4536E+00 -7.1075E+00  2.7132E+00  5.7427E+00  1.7758E+00 -9.8094E+00
             5.2838E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3211.92636429140        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      525             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0698E+00  1.2470E+00  1.7649E+00  8.5738E-01  1.3275E+00  9.8091E-01  8.5488E-01  2.6661E+00  1.1103E+00  1.2918E+00
             1.8113E+00
 PARAMETER:  1.6751E-01  3.2074E-01  6.6809E-01 -5.3877E-02  3.8328E-01  8.0721E-02 -5.6792E-02  1.0806E+00  2.0461E-01  3.5602E-01
             6.9402E-01
 GRADIENT:   2.8689E+02  1.2036E+02 -5.1054E+00  2.4902E+01  5.5415E+01  1.3860E+01  2.8468E+00  8.5068E+00  7.1368E+00  5.3455E+00
             2.8254E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3212.17161266885        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      622
 NPARAMETR:  1.0505E+00  1.2058E+00  1.7649E+00  8.7247E-01  1.2896E+00  9.7143E-01  8.6252E-01  2.5925E+00  1.0907E+00  1.2657E+00
             1.8036E+00
 PARAMETER:  1.4924E-01  2.8716E-01  6.6807E-01 -3.6428E-02  3.5430E-01  7.1013E-02 -4.7892E-02  1.0526E+00  1.8680E-01  3.3560E-01
             6.8976E-01
 GRADIENT:  -2.8812E+01 -2.8492E-01 -1.2048E+01 -1.5521E+01 -1.0164E+00 -4.8667E+00 -4.1426E-01 -8.8374E-02 -5.4396E-01 -2.6516E+00
             1.0110E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3212.96272028970        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  1.0694E+00  1.1362E+00  1.7687E+00  9.2999E-01  1.2469E+00  9.8823E-01  9.3364E-01  2.6345E+00  1.0410E+00  1.2409E+00
             1.7857E+00
 PARAMETER:  1.6709E-01  2.2770E-01  6.7026E-01  2.7422E-02  3.2065E-01  8.8155E-02  3.1332E-02  1.0687E+00  1.4018E-01  3.1583E-01
             6.7981E-01
 GRADIENT:   1.3113E+01  6.9198E+00 -1.8401E+01  5.6239E+00 -1.4408E-01  2.1482E+00  5.0851E-01  3.4246E-02  1.0529E+00  1.1646E+00
            -1.8228E+00

0ITERATION NO.:   32    OBJECTIVE VALUE:  -3212.97449970003        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:      859
 NPARAMETR:  1.0691E+00  1.1345E+00  1.7688E+00  9.3080E-01  1.2460E+00  9.8779E-01  9.3498E-01  2.6345E+00  1.0400E+00  1.2393E+00
             1.7857E+00
 PARAMETER:  1.6682E-01  2.2623E-01  6.7032E-01  2.8288E-02  3.1994E-01  8.7712E-02  3.2773E-02  1.0687E+00  1.3922E-01  3.1458E-01
             6.7984E-01
 GRADIENT:  -5.4782E+02  2.1149E+02 -8.6340E+01  4.7124E+02  2.9337E+02  4.6768E+02  4.6621E+02  8.7795E+01 -6.6829E+02 -2.9589E+02
            -1.5595E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      859
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.8105E-03 -3.5376E-02 -1.0082E-02  1.4742E-02 -3.8181E-02
 SE:             2.9591E-02  1.9315E-02  2.0389E-02  2.5213E-02  2.3141E-02
 N:                     100         100         100         100         100

 P VAL.:         8.9754E-01  6.7021E-02  6.2094E-01  5.5874E-01  9.8961E-02

 ETASHRINKSD(%)  8.6485E-01  3.5293E+01  3.1695E+01  1.5533E+01  2.2474E+01
 ETASHRINKVR(%)  1.7222E+00  5.8130E+01  5.3345E+01  2.8654E+01  3.9897E+01
 EBVSHRINKSD(%)  7.3691E-01  3.5531E+01  3.4184E+01  1.6861E+01  1.9149E+01
 EBVSHRINKVR(%)  1.4684E+00  5.8438E+01  5.6682E+01  3.0879E+01  3.4632E+01
 RELATIVEINF(%)  9.8511E+01  1.2816E+01  2.8856E+01  2.4183E+01  3.4202E+01
 EPSSHRINKSD(%)  1.9430E+01
 EPSSHRINKVR(%)  3.5085E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          896
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1646.7378515027735     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3212.9744997000321     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1566.2366481972585     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.37
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3212.974       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.07E+00  1.13E+00  1.77E+00  9.31E-01  1.25E+00  9.88E-01  9.35E-01  2.63E+00  1.04E+00  1.24E+00  1.79E+00
 


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
 #CPUT: Total CPU Time in Seconds,       61.647
Stop Time:
Sat Oct 23 15:07:40 CDT 2021
