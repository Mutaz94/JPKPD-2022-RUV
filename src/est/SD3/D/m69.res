Sun Oct 24 01:03:32 CDT 2021
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
$DATA ../../../../data/SD3/D/dat69.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2069.52123145763        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0407E+02  3.7879E+01 -5.1589E+01  1.4333E+02  1.0168E+02  1.2351E+01 -3.1950E+01 -1.5205E+00 -3.6188E+01 -2.2686E+00
            -1.3064E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2081.93112655833        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0250E+00  1.0419E+00  1.0888E+00  9.7193E-01  9.9439E-01  1.2301E+00  1.4332E+00  1.0427E+00  1.3506E+00  9.1775E-01
             9.8621E-01
 PARAMETER:  1.2466E-01  1.4101E-01  1.8512E-01  7.1528E-02  9.4374E-02  3.0708E-01  4.5988E-01  1.4182E-01  4.0051E-01  1.4171E-02
             8.6114E-02
 GRADIENT:  -2.5332E+01  3.6978E+00 -1.3331E+00  2.4726E+01  1.3951E+01  3.3033E+01  6.4683E+00 -2.9093E+00  3.8789E+01 -2.1376E+00
            -2.3613E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2084.10938601580        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0301E+00  9.1264E-01  1.1360E+00  1.0574E+00  9.5166E-01  1.2165E+00  1.7505E+00  1.1105E+00  1.1940E+00  8.1338E-01
             1.0209E+00
 PARAMETER:  1.2969E-01  8.5829E-03  2.2748E-01  1.5584E-01  5.0455E-02  2.9598E-01  6.5987E-01  2.0478E-01  2.7732E-01 -1.0656E-01
             1.2070E-01
 GRADIENT:  -1.7007E+01  1.3745E+01 -1.4541E+00  2.6759E+01  1.8491E+01  2.9475E+01  1.5411E+01 -8.9376E-01  3.1366E+01 -9.5023E+00
             2.6336E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2090.08881397042        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0414E+00  7.6664E-01  1.1347E+00  1.1206E+00  8.9876E-01  1.1218E+00  1.8058E+00  1.0817E+00  9.4848E-01  8.5668E-01
             1.0158E+00
 PARAMETER:  1.4057E-01 -1.6573E-01  2.2638E-01  2.1383E-01 -6.7338E-03  2.1491E-01  6.9100E-01  1.7853E-01  4.7103E-02 -5.4686E-02
             1.1569E-01
 GRADIENT:   1.3876E+00  2.4879E+00  1.4518E-01  1.7617E+00 -2.2329E+00 -4.2726E-01 -5.7724E-01  1.4096E-01 -5.0388E-01  5.6795E-01
             9.4358E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2090.22542235369        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      731             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0379E+00  6.5337E-01  1.2346E+00  1.1871E+00  8.9961E-01  1.1245E+00  2.0138E+00  1.1550E+00  9.1807E-01  8.6248E-01
             1.0161E+00
 PARAMETER:  1.3720E-01 -3.2561E-01  3.1076E-01  2.7154E-01 -5.7898E-03  2.1737E-01  8.0002E-01  2.4411E-01  1.4518E-02 -4.7948E-02
             1.1592E-01
 GRADIENT:   6.5866E+02  6.9201E+01  6.0314E+00  3.8295E+02  1.0579E+01  1.3026E+02  6.3961E+01  8.8092E-01  1.2290E+01  9.3306E-02
             1.3970E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2090.24013043033        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      910             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0386E+00  6.5030E-01  1.2366E+00  1.1917E+00  8.9678E-01  1.1233E+00  2.0274E+00  1.1538E+00  9.1656E-01  8.6726E-01
             1.0158E+00
 PARAMETER:  1.3789E-01 -3.3032E-01  3.1233E-01  2.7541E-01 -8.9425E-03  2.1627E-01  8.0677E-01  2.4306E-01  1.2872E-02 -4.2413E-02
             1.1572E-01
 GRADIENT:   6.6346E+02  7.1222E+01  7.3816E+00  3.9368E+02  6.3107E+00  1.2895E+02  6.4541E+01  8.4071E-01  1.2342E+01  1.0186E+00
             1.3407E+00

0ITERATION NO.:   28    OBJECTIVE VALUE:  -2090.24111753165        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     1005
 NPARAMETR:  1.0386E+00  6.5030E-01  1.2366E+00  1.1917E+00  8.9707E-01  1.1233E+00  2.0229E+00  1.1536E+00  9.1616E-01  8.6587E-01
             1.0158E+00
 PARAMETER:  1.3789E-01 -3.3032E-01  3.1233E-01  2.7541E-01 -8.6245E-03  2.1627E-01  8.0451E-01  2.4291E-01  1.2431E-02 -4.4020E-02
             1.1566E-01
 GRADIENT:  -5.8085E-01  1.6868E-01  6.1091E-01 -4.0719E+00 -1.1347E+00  5.8163E-01  3.7849E-01 -1.5509E-02  1.5820E-01  1.1471E-01
             7.3886E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1005
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8007E-03  1.5101E-02 -3.7774E-02 -1.5686E-02 -2.3896E-02
 SE:             2.9927E-02  2.0282E-02  1.6206E-02  2.3894E-02  2.0393E-02
 N:                     100         100         100         100         100

 P VAL.:         9.2544E-01  4.5652E-01  1.9755E-02  5.1153E-01  2.4129E-01

 ETASHRINKSD(%)  1.0000E-10  3.2054E+01  4.5709E+01  1.9951E+01  3.1681E+01
 ETASHRINKVR(%)  1.0000E-10  5.3833E+01  7.0525E+01  3.5921E+01  5.3325E+01
 EBVSHRINKSD(%)  2.8204E-01  3.3998E+01  4.8774E+01  1.8955E+01  2.8965E+01
 EBVSHRINKVR(%)  5.6328E-01  5.6438E+01  7.3759E+01  3.4316E+01  4.9541E+01
 RELATIVEINF(%)  9.8835E+01  4.7425E+00  5.0399E+00  7.6090E+00  1.1987E+01
 EPSSHRINKSD(%)  3.4431E+01
 EPSSHRINKVR(%)  5.7007E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2090.2411175316456     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1171.3025843269729     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2090.241       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  6.50E-01  1.24E+00  1.19E+00  8.97E-01  1.12E+00  2.02E+00  1.15E+00  9.16E-01  8.66E-01  1.02E+00
 


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
 #CPUT: Total CPU Time in Seconds,       36.516
Stop Time:
Sun Oct 24 01:03:40 CDT 2021
