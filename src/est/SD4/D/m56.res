Sun Oct 24 04:21:00 CDT 2021
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
$DATA ../../../../data/SD4/D/dat56.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m56.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1637.90074320467        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7106E+02 -4.2937E+01 -3.3494E+00 -5.5184E+01 -2.7654E+01  2.4964E+00 -2.8593E+01  3.7945E+00 -2.7869E+01  2.6492E+01
             2.2568E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1648.33392456467        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0497E+00  1.1612E+00  1.1394E+00  1.0823E+00  1.1321E+00  1.3714E+00  1.3345E+00  1.0027E+00  1.2420E+00  8.3022E-01
             9.1684E-01
 PARAMETER:  1.4853E-01  2.4943E-01  2.3054E-01  1.7907E-01  2.2407E-01  4.1581E-01  3.8857E-01  1.0267E-01  3.1676E-01 -8.6067E-02
             1.3181E-02
 GRADIENT:   3.5858E+01  5.2174E+01  7.9001E+00  9.0507E+01  1.1234E+01  5.8714E+01 -2.3622E+00 -1.1295E+01  2.4793E+01 -5.1331E+00
            -1.8651E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1650.37445644424        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0368E+00  9.2951E-01  1.2054E+00  1.2350E+00  1.0274E+00  1.3415E+00  1.6276E+00  1.0874E+00  1.0970E+00  7.0304E-01
             9.0993E-01
 PARAMETER:  1.3610E-01  2.6904E-02  2.8683E-01  3.1104E-01  1.2706E-01  3.9379E-01  5.8713E-01  1.8377E-01  1.9259E-01 -2.5234E-01
             5.6081E-03
 GRADIENT:   2.3640E+01  5.3442E+01  1.7300E+01  1.0504E+02 -1.8986E+01  5.2885E+01 -1.1664E+00 -7.8675E+00  2.3738E+01 -1.0051E+01
            -2.2569E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1663.70196828705        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0187E+00  8.5660E-01  9.9601E-01  1.1585E+00  9.4419E-01  1.1537E+00  1.7545E+00  8.9644E-01  8.6030E-01  6.9140E-01
             9.5736E-01
 PARAMETER:  1.1856E-01 -5.4786E-02  9.5998E-02  2.4715E-01  4.2576E-02  2.4298E-01  6.6221E-01 -9.3261E-03 -5.0472E-02 -2.6903E-01
             5.6421E-02
 GRADIENT:  -1.1162E+00  1.5091E+00 -4.4741E+00  2.6304E+00  7.7769E+00 -1.7742E+00 -5.8976E-01  1.8282E-01 -8.7232E-01 -2.3651E-01
             2.0729E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.18768225623        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0165E+00  6.2889E-01  1.1534E+00  1.3088E+00  9.2525E-01  1.1553E+00  2.1713E+00  1.0115E+00  8.1755E-01  7.1972E-01
             9.4804E-01
 PARAMETER:  1.1639E-01 -3.6380E-01  2.4270E-01  3.6914E-01  2.2304E-02  2.4436E-01  8.7532E-01  1.1145E-01 -1.0145E-01 -2.2890E-01
             4.6641E-02
 GRADIENT:   1.1058E+00  5.6551E+00  4.6816E+00  6.2883E+00 -1.0993E+01  2.4727E-01  2.9604E-01 -5.7770E-01  2.3524E-01  1.1419E+00
            -1.2104E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1664.33908007448        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  1.0147E+00  5.4376E-01  1.2794E+00  1.3608E+00  9.5842E-01  1.1537E+00  2.3448E+00  1.1331E+00  8.0973E-01  7.4580E-01
             9.5059E-01
 PARAMETER:  1.1458E-01 -5.0925E-01  3.4640E-01  4.0805E-01  5.7527E-02  2.4296E-01  9.5222E-01  2.2496E-01 -1.1105E-01 -1.9330E-01
             4.9332E-02
 GRADIENT:   8.5094E-01  5.7390E-01  5.6820E-02 -3.3765E+00 -2.1384E-01  5.8984E-01  2.7655E-02 -2.6211E-01 -3.2158E-02  9.9183E-02
            -1.6948E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1664.34490329995        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1087             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0155E+00  5.4273E-01  1.2795E+00  1.3595E+00  9.5847E-01  1.1561E+00  2.3588E+00  1.1370E+00  8.0828E-01  7.4388E-01
             9.5094E-01
 PARAMETER:  1.1537E-01 -5.1114E-01  3.4646E-01  4.0710E-01  5.7588E-02  2.4509E-01  9.5814E-01  2.2838E-01 -1.1285E-01 -1.9587E-01
             4.9697E-02
 GRADIENT:   5.2521E+02  6.0372E+01  8.6014E+00  5.6289E+02  6.9270E+00  1.8045E+02  5.9319E+01  9.4845E-01  9.3043E+00  7.4475E-01
             8.0060E-01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -1664.34490329995        NO. OF FUNC. EVALS.:  61
 CUMULATIVE NO. OF FUNC. EVALS.:     1148
 NPARAMETR:  1.0155E+00  5.4273E-01  1.2795E+00  1.3595E+00  9.5847E-01  1.1561E+00  2.3588E+00  1.1370E+00  8.0828E-01  7.4388E-01
             9.5094E-01
 PARAMETER:  1.1537E-01 -5.1114E-01  3.4646E-01  4.0710E-01  5.7588E-02  2.4509E-01  9.5814E-01  2.2838E-01 -1.1285E-01 -1.9587E-01
             4.9697E-02
 GRADIENT:   5.8524E-04  5.2366E-02 -1.4965E-01 -3.8525E-02  2.0369E-01 -2.7866E-04 -2.5616E-02  4.4355E-02  1.8029E-02 -2.9992E-03
             9.8942E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1148
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.8143E-04  2.8411E-02 -3.9793E-02 -2.8153E-02 -2.5526E-02
 SE:             2.9948E-02  2.0606E-02  1.7806E-02  2.3024E-02  1.7787E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7386E-01  1.6796E-01  2.5433E-02  2.2141E-01  1.5127E-01

 ETASHRINKSD(%)  1.0000E-10  3.0967E+01  4.0347E+01  2.2868E+01  4.0410E+01
 ETASHRINKVR(%)  1.0000E-10  5.2344E+01  6.4415E+01  4.0507E+01  6.4490E+01
 EBVSHRINKSD(%)  3.0284E-01  3.4342E+01  4.2081E+01  2.0128E+01  3.8065E+01
 EBVSHRINKVR(%)  6.0477E-01  5.6891E+01  6.6454E+01  3.6205E+01  6.1640E+01
 RELATIVEINF(%)  9.8875E+01  6.4288E+00  6.4606E+00  1.0556E+01  7.0781E+00
 EPSSHRINKSD(%)  4.4965E+01
 EPSSHRINKVR(%)  6.9711E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1664.3449032999458     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -929.19407673620765     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.64
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1664.345       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  5.43E-01  1.28E+00  1.36E+00  9.58E-01  1.16E+00  2.36E+00  1.14E+00  8.08E-01  7.44E-01  9.51E-01
 


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
 #CPUT: Total CPU Time in Seconds,       35.989
Stop Time:
Sun Oct 24 04:21:08 CDT 2021
