Sun Oct 24 03:21:42 CDT 2021
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
$DATA ../../../../data/SD4/SL2/dat75.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1643.09770514388        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2109E+02 -4.9899E+01  4.8932E+00 -7.2994E+01  6.7304E+00  4.1898E+01 -1.5303E+00  2.6002E-01 -1.0178E+01  7.3322E+00
            -3.0154E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1650.93934142283        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  9.9697E-01  1.0415E+00  1.0072E+00  1.0698E+00  1.0039E+00  1.0182E+00  1.0041E+00  1.0013E+00  1.0150E+00  9.8227E-01
             1.0727E+00
 PARAMETER:  9.6969E-02  1.4062E-01  1.0714E-01  1.6748E-01  1.0390E-01  1.1805E-01  1.0408E-01  1.0127E-01  1.1489E-01  8.2112E-02
             1.7022E-01
 GRADIENT:  -6.5735E-02  7.3511E-01  2.5152E-01 -1.9083E+00 -4.9680E+00 -2.9090E-01 -9.1492E-02 -7.7048E-01 -5.2464E-01  3.9258E+00
            -2.2234E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1651.02996956202        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  9.9759E-01  1.0217E+00  1.0444E+00  1.0852E+00  1.0111E+00  1.0135E+00  1.0133E+00  1.0744E+00  9.9582E-01  9.5164E-01
             1.0732E+00
 PARAMETER:  9.7585E-02  1.2152E-01  1.4345E-01  1.8174E-01  1.1106E-01  1.1342E-01  1.1325E-01  1.7174E-01  9.5809E-02  5.0434E-02
             1.7067E-01
 GRADIENT:   1.8884E+00  2.2424E+00  8.6062E-01 -5.9521E-01 -2.0957E+00 -1.9565E+00 -6.1414E-01 -5.8118E-02 -3.0222E+00 -7.5927E-01
            -9.2074E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1651.71798359267        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      532
 NPARAMETR:  9.9335E-01  9.2685E-01  1.2484E+00  1.1557E+00  1.0702E+00  1.0287E+00  7.2756E-01  1.2562E+00  1.0810E+00  1.0591E+00
             1.0819E+00
 PARAMETER:  9.3324E-02  2.4034E-02  3.2187E-01  2.4469E-01  1.6785E-01  1.2828E-01 -2.1806E-01  3.2806E-01  1.7789E-01  1.5738E-01
             1.7872E-01
 GRADIENT:  -3.4963E+00  4.9547E+00  4.2230E-01  1.0784E+01 -2.2218E+00  4.4524E+00  3.3008E+00 -1.1791E+00  8.9002E+00  1.4368E+00
             3.9423E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1653.97763816289        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  9.9421E-01  9.4144E-01  1.3505E+00  1.1299E+00  1.1287E+00  1.0182E+00  1.0849E-01  1.4547E+00  1.1231E+00  1.1068E+00
             1.0743E+00
 PARAMETER:  9.4198E-02  3.9657E-02  4.0047E-01  2.2211E-01  2.2104E-01  1.1799E-01 -2.1211E+00  4.7483E-01  2.1613E-01  2.0144E-01
             1.7171E-01
 GRADIENT:  -9.4259E-01 -3.6521E+00 -2.7626E+00 -3.5754E+00  5.1317E+00  5.9110E-01  1.6102E-01  5.0895E-01 -7.2648E-01 -8.2925E-02
             9.6911E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1654.07403920601        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.9487E-01  9.6003E-01  1.3453E+00  1.1183E+00  1.1275E+00  1.0167E+00  2.4464E-02  1.4589E+00  1.1415E+00  1.1067E+00
             1.0721E+00
 PARAMETER:  9.4857E-02  5.9211E-02  3.9664E-01  2.1185E-01  2.1998E-01  1.1660E-01 -3.6106E+00  4.7771E-01  2.3238E-01  2.0135E-01
             1.6962E-01
 GRADIENT:   1.4646E-01 -4.2043E-01  4.0059E-01 -1.0806E+00 -7.0943E-01 -6.3487E-03  1.0100E-02 -1.7020E-01  3.6306E-01  1.4632E-01
             7.1598E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1654.07637424094        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1063             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9556E-01  9.5153E-01  1.3500E+00  1.1246E+00  1.1264E+00  1.0172E+00  1.0000E-02  1.4589E+00  1.1347E+00  1.1052E+00
             1.0720E+00
 PARAMETER:  9.5553E-02  5.0321E-02  4.0011E-01  2.1743E-01  2.1902E-01  1.1706E-01 -4.5433E+00  4.7767E-01  2.2636E-01  2.0003E-01
             1.6954E-01
 GRADIENT:   3.6398E+02  2.9483E+01  5.0587E+00  1.3209E+02  1.0395E+01  5.1924E+01  5.2653E-04  1.4392E+00  1.6685E+01  1.4021E+00
             1.3523E+00

0ITERATION NO.:   32    OBJECTIVE VALUE:  -1654.07637424094        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:     1118
 NPARAMETR:  9.9556E-01  9.5153E-01  1.3500E+00  1.1246E+00  1.1264E+00  1.0172E+00  1.0000E-02  1.4589E+00  1.1347E+00  1.1052E+00
             1.0720E+00
 PARAMETER:  9.5553E-02  5.0321E-02  4.0011E-01  2.1743E-01  2.1902E-01  1.1706E-01 -4.5433E+00  4.7767E-01  2.2636E-01  2.0003E-01
             1.6954E-01
 GRADIENT:   1.7700E+00 -2.9703E-01  6.3033E-02 -4.5881E-01 -7.0353E-02  1.9341E-01  2.6327E-04 -1.7746E-04  5.5136E-02  1.8363E-02
            -4.0802E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1118
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.8262E-04 -6.5998E-04 -3.7580E-02 -2.7490E-03 -4.0522E-02
 SE:             2.9860E-02  1.8983E-04  1.7042E-02  2.9229E-02  2.1244E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7642E-01  5.0761E-04  2.7445E-02  9.2507E-01  5.6456E-02

 ETASHRINKSD(%)  1.0000E-10  9.9364E+01  4.2907E+01  2.0783E+00  2.8831E+01
 ETASHRINKVR(%)  1.0000E-10  9.9996E+01  6.7404E+01  4.1135E+00  4.9350E+01
 EBVSHRINKSD(%)  4.6620E-01  9.9421E+01  4.5939E+01  2.7088E+00  2.6303E+01
 EBVSHRINKVR(%)  9.3023E-01  9.9997E+01  7.0774E+01  5.3442E+00  4.5687E+01
 RELATIVEINF(%)  9.8768E+01  3.2958E-04  9.2112E+00  1.1377E+01  1.7185E+01
 EPSSHRINKSD(%)  4.4689E+01
 EPSSHRINKVR(%)  6.9407E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1654.0763742409356     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -918.92554767719741     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1654.076       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  9.52E-01  1.35E+00  1.12E+00  1.13E+00  1.02E+00  1.00E-02  1.46E+00  1.13E+00  1.11E+00  1.07E+00
 


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
 #CPUT: Total CPU Time in Seconds,       32.571
Stop Time:
Sun Oct 24 03:21:50 CDT 2021
